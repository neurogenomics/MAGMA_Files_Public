prepare_Yengo2018 <- function(save_dir = here::here("data","GWAS_munged"),
                              id = "Yengo2018.Height",
                              nThread = 8){ 
  
  #### Download ####
  path <- downloadR::downloader("https://portals.broadinstitute.org/collaboration/giant/images/6/63/Meta-analysis_Wood_et_al%2BUKBiobank_2018.txt.gz",
                                output_dir = here::here("data","GWAS_raw"))
  #### Munge ####
  gwas_paths <- MungeSumstats::format_sumstats(
    path = path, 
    save_path = file.path(save_dir,id,paste0(id,".tsv.gz")),
    dbSNP = 144,
    ref_genome = "GRCh37",
    bi_allelic_filter = TRUE,
    tabix_index = TRUE,
    log_folder = file.path(save_dir,id,"logs"),
    log_folder_ind = TRUE,
    log_mungesumstats_msgs = TRUE,
    nThread = nThread
  ) 
  #### Map #### 
  t1 <- Sys.time()
  magma_files <-  MAGMA.Celltyping::map_snps_to_genes(
    path_formatted = gwas_paths$sumstats, 
    genome_build = "GRCH37",  
    population = "EUR",
    upstream_kb = 35,  
    downstream_kb = 10, 
    force_new = FALSE
  )
  t2 <- Sys.time()
  print(t2-t1)
  #### Return ####
  return(
    list(gwas_paths=gwas_paths,
         magma_files=magma_files)
  )
}
