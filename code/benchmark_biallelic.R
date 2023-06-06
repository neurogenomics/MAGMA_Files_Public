
root <- "/rds/general/project/neurogenomics-lab/ephemeral/MAGMA_benchmark/"
dir.create(root,showWarnings = FALSE, recursive = TRUE)

ids <- c("ieu-b-7"," ieu-b-5067","ieu-a-298")

#### Munge sumstats ####
gwas_paths <- MungeSumstats::import_sumstats(
  ids = ids, 
  save_dir = file.path(root,"multiallelic"), 
  nThread = 8,
  parallel_across_ids = FALSE, 
  force_new_vcf = TRUE,
  force_new = TRUE,
  vcf_download = TRUE,
  vcf_dir = file.path(root,"data","VCFs"),
  #### Filters ####
  bi_allelic_filter = FALSE,
  compute_z = TRUE,
  force_new_z = TRUE,
  tabix_index = TRUE,
  #### Record logs ####
  log_folder_ind = FALSE, # =TRUE takes up too much storage space 
  log_mungesumstats_msgs = TRUE
)  
gwas_paths <- MungeSumstats::import_sumstats(
  ids = ids, 
  save_dir = file.path(root,"biallelic"), 
  nThread = 8,
  parallel_across_ids = FALSE, 
  force_new_vcf = TRUE,
  force_new = TRUE,
  vcf_download = TRUE,
  vcf_dir = file.path(root,"data","VCFs"),
  #### Filters ####
  bi_allelic_filter = TRUE, # <---- Key difference
  compute_z = TRUE,
  force_new_z = TRUE,
  tabix_index = TRUE,
  #### Record logs ####
  log_folder_ind = FALSE, # =TRUE takes up too much storage space 
  log_mungesumstats_msgs = TRUE
)  


#### Gather metadata ####
source("code/utils.R")
meta <- gather_metadata(save_dir = root,
                        N_dict=c("Wightman2021"=1126563, 
                                 "Vuckovic2020"=408112)) 

#### Map SNPS to genes ####
t1 <- Sys.time()
magma_files <- parallel::mclapply(seq_len(nrow(meta)),
                                  function(i){ 
                                    EWCE:::message_parallel("----- ",i," : ",
                                                            meta$id[i]," -----") 
                                    tryCatch(expr = { 
                                      MAGMA.Celltyping::map_snps_to_genes(
                                        # version = "1.08",
                                        path_formatted = meta$munged_path[i],
                                        genome_build = meta$build_final[i],
                                        N = if(is.na(meta$N[i])) NULL else meta$N[i],
                                        population = meta$population_1KG[i],
                                        upstream_kb = 35,
                                        downstream_kb = 10,
                                        genes_only = FALSE, 
                                        force_new = FALSE
                                      )
                                    }, error = function(e) {EWCE:::message_parallel(e);NULL}) 
                                  }, mc.cores = min(nrow(meta),20) ) |> `names<-`(meta$id) 
t2 <- Sys.time()
print(t2-t1)
