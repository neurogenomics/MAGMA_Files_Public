pbs_logs <- function(dir="/Volumes/bms20/projects/neurogenomics-lab/live/Data/MAGMA_Files_Public/code/pbs_logs",
                     pattern="batch[0-9]+\\.pbs\\.",
                     show_plot = TRUE,
                     save_path = tempfile("pbs_report",fileext = ".csv")){ 
  
    requireNamespace("ggplot2")
  #### Make report out of pbs logs ####
  pbs_files <- list.files(dir, pattern, full.names = TRUE) 
  pbs_res <- parallel::mclapply(pbs_files, 
                    function(x){
      message("Processing: ",x)
      pbs <- readLines(x)
      pbs <- pbs[pbs!=""]
      #### Parse pbspro report ####
      if(grepl("Processing dataset :",pbs[1])){  
          i <- grep("Requested  :",pbs)
          data.table::data.table(
              id = utils::tail(strsplit(grep("Processing dataset : ",
                                      pbs,value = TRUE)," ")[[1]],2)[1],
              memory_gb_requested = strsplit(gsub("[ ]+"," ",pbs[i])," ")[[1]][4],
              memory_gb_peak = strsplit(gsub("[ ]+"," ",pbs[i+1])," ")[[1]][4],
              ncpus_requested = strsplit(gsub("[ ]+"," ",pbs[i])," ")[[1]][5],
              ncpus_ave = strsplit(gsub("[ ]+"," ",pbs[i+1])," ")[[1]][6],
              pbs_o = basename(x)
          ) 
      #### Parse MungeSumstats report ####
      } else if(any(grepl("Formatted summary statistics will be saved to",pbs))){
          munged_path <- tail(strsplit(grep(
              "^Formatted summary statistics will be saved to",
              pbs,value = TRUE)," ")[[1]],1)
          d <- grep(": Done in *.* minutes\\.", pbs,value = TRUE)
          if(length(d)>0) {
            id <- strsplit(d," ")[[1]][[1]]
            time <- strsplit(d," ")[[1]][[5]]
          } else {
            id <- time <- NA;
          }
          data.table::data.table(
              id = id,
              time = time,
              munged_path = munged_path,
              pbs_e = basename(x)
          ) 
      } else {
          pbs
      }
  }, mc.cores = 10) 
  #### Concat pbpro results ####
  probpro_report <- data.table::rbindlist(
      pbs_res[sapply(pbs_res,function(x){
          !is.na(x) && names(x)[2]=="memory_gb_requested"
          })]
  )
  cols <- c("memory_gb_requested","memory_gb_peak","ncpus_requested","ncpus_ave")
  probpro_report[, (cols) := lapply(.SD, as.numeric), .SDcols = cols]
  #### Concat MSS results ####
  mungesumstats_report <-  data.table::rbindlist(
      pbs_res[sapply(pbs_res,function(x){
          !is.na(x) && names(x)[2]=="time"
      })]
  )
  ## Get munged file sizes
  mungesumstats_report[,
      munged_size_mb:=file.size(
          remote_to_local(path = mungesumstats_report$munged_path)
          )/1e+6
  ]
  cols <- c("time")
  mungesumstats_report[, (cols) := lapply(.SD, as.numeric), .SDcols = cols]
  #### Merge reports ####
  report <- merge(probpro_report,mungesumstats_report,by = "id")
  #### Plot ####
  if(isTRUE(show_plot)){
      plot_dat <- data.table::melt.data.table(
          report,
          id.vars = c("id","munged_path","pbs_o","pbs_e"))
      ggp1 <- ggplot2::ggplot(plot_dat, ggplot2::aes(x = value)) + 
          ggplot2::geom_histogram() + 
          ggplot2::facet_wrap(~ variable, 
                              scales = "free", ncol = 2) +
          ggplot2::labs(title=paste(
              "MungeSumstats report:",
              formatC(length(unique(plot_dat$id)),big.mark = ","),
              "unique GWAS"))
      ggp1
  }
  #### Save ####
  if(!is.null(save_path)){
    message("Writing report --> ",save_path)
    dir.create(dirname(save_path),showWarnings = FALSE, recursive = TRUE)
    data.table::fwrite(report,save_path)
  }
  return(report)
} 
