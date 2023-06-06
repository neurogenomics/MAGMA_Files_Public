prepare_Vuckovic2020 <- function(save_dir = here::here("data","GWAS_munged"),
                                 nThread = 8){
  
  
  ## Get abbrev-to-trait mapping file
  trait_map <- readxl::read_xlsx(here::here("code/Vuckovic2020_mmc1.xlsx"), 
                                 skip = 2) |>
    data.table::data.table(check.names = TRUE) |> 
    tidyr::fill(`Cell.Type`) |>
    dplyr::mutate(file_abbreviation=tolower(
      gsub("RDW","rdw_cv",
           gsub("[%]","_p",
                gsub("#","",Standard.Abbreviation))
      )
    )
    )
  ## List all files on FTP server ####
  result <- (
    RCurl::getURL("ftp://ftp.sanger.ac.uk/pub/project/humgen/summary_statistics/UKBB_blood_cell_traits/",
                  verbose=TRUE,
                  ftp.use.epsv=TRUE,
                  dirlistonly = TRUE) |> 
      strsplit(split = "\r*\n")
  )[[1]] 
  ## Remove readme file
  result <- result[result!="nohup.out"]
  ## Merge metadata
  ftp_files <- data.table::data.table(
    file_abbreviation=gsub("\\.assoc$","",result),
    link=paste0("ftp://ftp.sanger.ac.uk/pub/project/humgen/summary_statistics/UKBB_blood_cell_traits/",result)
  ) |> merge(trait_map,by="file_abbreviation", all.x=TRUE)
  # ftp_files[is.na(Long.Name)] 
  
  #### Download raw files ####
  files <- downloadR::downloader(input_url = ftp_files$link, 
                                 output_path = file.path(here::here("data/GWAS_raw"),
                                                         paste0("Vuckovic2020.",
                                                                gsub("\\.assoc","",
                                                                     basename(ftp_files$link))
                                                         )
                                 )
  ) 
  names(files) <- basename(files) 
  
  #### Munge ####
  files_munged <- lapply(names(files),  
                         function(nm){
                           message("Processing: ",basename(nm)) 
                           save_dir_i <- here::here("data/GWAS_munged",nm)
                           MungeSumstats::format_sumstats(
                             path = files[[nm]], 
                             save_path = file.path(save_dir_i,paste0(nm,".tsv.gz")),
                             dbSNP = 144,
                             bi_allelic_filter = TRUE,
                             tabix_index = TRUE,
                             log_folder = file.path(save_dir_i,"logs"),
                             log_folder_ind = TRUE,
                             log_mungesumstats_msgs = TRUE,
                             nThread = nThread
                           )
                         }) 
  
  #### Map #### 
  # Unfortunately, per-SNP sample size is not included in any of the summary stats provided. 
  # Instead, we must assume a sample size of 408,112 for all blood trait GWAS, as reported in the original paper. 
  sumstats <- MungeSumstats::list_sumstats(save_dir = save_dir,
                                           pattern = "Vuckovic") 
  t1 <- Sys.time()
  magma_files <- parallel::mclapply(sumstats,
                                    function(x){ 
                                      MAGMA.Celltyping:::message_parallel("----- Processing: ",x," -----") 
                                      tryCatch(expr = {
                                        MAGMA.Celltyping::map_snps_to_genes(
                                          path_formatted = x,
                                          genome_build = "GRCH37",
                                          ## Assuming equal sample size across all GWAS 
                                          ## (only one N was reported in the paper)
                                          N = 408112,
                                          ## These UKB GWAS only use British European ancestry individuals 
                                          population = "EUR",
                                          upstream_kb = 35,  
                                          downstream_kb = 10, 
                                          force_new = FALSE
                                        )
                                      }, error = function(e) {message(e);NULL}) 
                                    }, mc.cores = min(length(sumstats),nThread) )
  t2 <- Sys.time()
  print(t2-t1) 
  #### Return ####
  return(
    list(gwas_paths=files_munged,
         magma_files=magma_files)
  )
}