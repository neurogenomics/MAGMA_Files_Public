prepare_Wightman2021 <- function(save_dir = here::here("data","GWAS_munged"),
                                 nThread = 10){
  
  # header: "chr" "PosGRCh37" "testedAllele" "otherAllele" "z" "p" "N"  
  mapping_file <- MungeSumstats:::sumstatsColHeaders
  mapping_file <- rbind(mapping_file,
                        data.frame("Uncorrected"=c("POSGRCH37","TESTEDALLELE"),
                                   "Corrected"=c("BP","A1"))) 
  path <- "https://ctg.cncr.nl/documents/p1651/PGCALZ2sumstatsExcluding23andMe.txt.gz"
  dat <- data.table::fread(path) 
  tmp <- file.path(tempdir(),basename(path))
  data.table::fwrite(dat, tmp, sep="\t", nThread = nThread)
  
  id <- "Wightman2021"
  #### Munge ####
  gwas_paths <- MungeSumstats::format_sumstats(
    path = tmp,
    save_path = here::here(save_dir,id,paste0(id,".tsv.gz")),
    ref_genome = "GRCh37",
    dbSNP = 144,
    bi_allelic_filter = TRUE,
    tabix_index = TRUE,
    nThread = nThread,
    force_new = FALSE,
    #### Record logs  
    log_folder = here::here(save_dir,id,"logs"),
    log_folder_ind = TRUE,
    log_mungesumstats_msgs = TRUE, 
    mapping_file = mapping_file) 
  # gwas_paths <- MungeSumstats::list_sumstats(save_dir = save_dir,
  #                                          pattern = "Wightman2021")
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