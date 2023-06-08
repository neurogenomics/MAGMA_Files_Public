prepare_Nalls2019_23andMe <- function(save_dir = here::here("data","GWAS_munged"),
                                      id = "Nalls2019_23andMe.PD",
                                      nThread = 8){ 
   
  #### Munge ####
  ## Previously munged by Alan Murphy.
  # ss <- data.table::fread("~/ParkinsonMeta2018.tbl --- Original From Julien")
  # reformatted <-
  #   MungeSumstats::format_sumstats(path=ss,
  #                                  ref_genome = "GRCh38",
  #                                  dbSNP="144",
  #                                  bi_allelic_filter = TRUE,
  #                                  force_new=TRUE,
  #                                  save_path = "~/PD_23andME_dbSNP144_no_non_bi.tsv.gz")
  path <- "/home/bms20/RDS/project/gwas_23andme/live/Munged/Munged_PD/PD_23andME_dbSNP144_no_non_bi.tsv.gz"
  tmp <- file.path(save_dir,id,paste0(id,".tsv.gz"))
  dir.create(dirname(tmp))
  file.copy(path,tmp)
  #### Map #### 
  ## Per-SNP sample size not provided in summary stats.
  ## Using reported sample size from paper instead: 
  # "37 688 cases, 18 618 UK Biobank proxy-cases 
  # (ie, individuals who do not have Parkinson's disease but have a first degree 
  # relative that does), and 1·4 million controls. 
  t1 <- Sys.time()
  magma_files <-  MAGMA.Celltyping::map_snps_to_genes(
    path_formatted = tmp, 
    genome_build = "GRCH38",  
    population = "EUR", 
    N = 37688 + 18618 + 1.4e6,
    upstream_kb = 35,  
    downstream_kb = 10, 
    force_new = FALSE
  )
  t2 <- Sys.time()
  print(t2-t1)
  #### Move MAGMA files to private repo ####
  magma_files_final <- file.path("/home/bms20/projects/MAGMA_Files")
  dir.create(magma_files_final)
  file.copy(dirname(magma_files),magma_files_final,recursive = TRUE)
  #### Update metadata file ####
  key_path <- "/home/bms20/projects/MAGMA_Files/study_key.xlsx"
  key <- xlsx::read.xlsx(key_path, sheetIndex = 1)
  key2 <- rbind(key,
        list(trait="Parkinson's disease",
             trait_source="Parkinson's disease (Nalls et al., 2019)",
             reference=NA,
             prefix=basename(dirname(magma_files)),
             restricted=TRUE
             )
        )
  xlsx::write.xlsx(key2,key_path)
  #### Delete temp files 
  unlink(file.path(save_dir,id),recursive = TRUE)  
  #### Return ####
  return(
    list(gwas_paths=NULL,
         magma_files=file.path(magma_files_final,basename(dirname(magma_files))))
  )
}
