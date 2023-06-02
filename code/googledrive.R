#' Google Drive upload
#' 
#' @source \href{https://github.com/tidyverse/googledrive/issues/345}{
#' drive_upload() hangs indefinitely #345}
googledrive_upload <- function(repo="/rds/general/project/neurogenomics-lab/ephemeral/Data/MAGMA_Files_Public",
                               gwas_dir=file.path(repo,"data","GWAS_sumstats"),
                               pattern="*\\.tsv",
                               drive_dir="~/GWAS_sumstats",
                               force_new=FALSE){ 
  
  # templateR:::args2vars(googledrive_upload) 
  files <- list.files(path = gwas_dir,
                      pattern = pattern,
                      full.names = TRUE, 
                      recursive = TRUE)
  snp_files <- c("dup_snp_id",
                 "np_not_found_from_chr_bp",
                 "snp_bi_allelic",
                 "snp_missing_rs",
                 "dup_base_pair_position",
                 "alleles_dont_match_ref_gen",
                 "snp_multi_colon")
  log_files <- list.files(path = gwas_dir,
                          pattern = ".txt$", 
                          full.names = TRUE,
                          recursive = TRUE)
  gwas_files <- grep(pattern = paste(snp_files,collapse = "|"),
                     x = files,
                     invert = TRUE,
                     value = TRUE)
  names(gwas_files) <- basename(dirname(gwas_files))
  message(formatC(length(gwas_files),big.mark = ",")," files found.")
  #### Check files in GoogleDrive ####
  drive_files <- googledrive::drive_ls(path = drive_dir, 
                                       type = c(".tsv",".tsv.bgz",".txt"),
                                       n_max = Inf,
                                       recursive = TRUE)
  #### Iterate over IDs ####
  all_ids <- stats::setNames(unique(names(gwas_files)),
                             unique(names(gwas_files)))
  drive_paths <- lapply(all_ids, 
                        function(id){
    message("Processing: ",id)
    id_files <- gwas_files[names(gwas_files)==id]
    tbx_files <- list("tsv"=grep(".tsv$",id_files, value = TRUE),
                      "tbx"=grep(".tsv.bgz$",id_files, value = TRUE),
                      "tbi"=grep(".tsv.bgz.tbi$",id_files, value = TRUE)) 
    #### Make sure there's a decompressed version ####
    if(length(tbx_files$tsv)==0 &&
       length(tbx_files$tbx)==1){
      tbx_files$tsv <- R.utils::gunzip(filename=tbx_files$tbx,
                                       destname=gsub("\\.bgz$","",tbx_files$tbx), 
                                       remove=FALSE,
                                       overwrite=TRUE) |> as.character()
    }
    #### Re-convert to tabix-indexed format ####
    # dat <- data.table::fread(tbx_files$tsv)
    if(length(tbx_files$tbi)==0){
      tbx_files2 <- echotabix::convert(target_path = tbx_files$tsv,
                                       chrom_col = "CHR",
                                       start_col = "BP",
                                       force_new = TRUE)
      tbx_files$tbx <- tbx_files2$path
      tbx_files$tbi <- tbx_files2$index
    } 
    #### Upload to googledrive ####
    upload_files <- unlist(tbx_files) 
    #### Make id directory ####
    message(" - Preparing subdirectory.")
    if(isTRUE(id %in% drive_files$name)){
      id_dir <- drive_files[drive_files$name==id,]
    } else {
      id_dir <- googledrive::drive_mkdir(name =  id,
                                         path = drive_dir)
    }
    #### Iterate over each file type ####
    tbx_files_drive <- lapply(upload_files,
                              function(x){                               
      #### Handle whether force new uploads or not ####
      if(isTRUE(basename(x) %in% drive_files$name) &&
         isFALSE(force_new)){
        message(" - Skipping: ",basename(x))
        return(subset(drive_files,name==basename(x)))
      } else {
        message(" - Uploading: ",basename(x))
        googledrive::drive_upload(media = x,
                                  path  = id_dir,
                                  name = basename(x))
      }
                                })   
    return(tbx_files_drive)
  })
  return(drive_paths)
}