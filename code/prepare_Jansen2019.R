prepare_Jansen2019 <- function(save_dir = here::here("data","GWAS_munged"),
                               id = "Jansen2019.AD",
                               nThread = 8){ 
  
  #### Download ####
  path <- downloadR::downloader("https://ctg.cncr.nl/documents/p1651/AD_sumstats_Jansenetal_2019sept.txt.gz",
                                output_dir = here::here("data","GWAS_raw"))
  #### Munge ####
  gwas_paths <- MungeSumstats::format_sumstats(
    path = path, 
    save_path = file.path(save_dir,id,paste0(id,".tsv.gz")),
    dbSNP = 144,
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
