#!/usr/bin/env Rscript
library("optparse")
 
option_list <-list(
    optparse::make_option(c("-i", "--id"), type="character", default=NULL, 
              help="OpenGWAS ID", metavar="character"),
    optparse::make_option( c("-n", "--ncpus"), type="integer", default=4,
              help="Number of CPUs to use.", metavar="character"),
    optparse::make_option( c("-p", "--population"), type="character", default="EUR",
                           help="Reference population", metavar="character")
);
opt_parser <- optparse::OptionParser(option_list=option_list)
opt <- optparse::parse_args(opt_parser)
root <- "/rds/general/project/neurogenomics-lab/ephemeral/MAGMA_Files_Public"

#### Prepare metadata: only process IDs not already done on OpenStack ####
# meta_done <- data.table::fread("https://github.com/neurogenomics/MAGMA_Files_Public/raw/master/metadata.csv")
# meta_done <- meta_done[meta_done$munged_path=="",]
# meta  <- MungeSumstats::find_sumstats(consortia = c("Neale Lab","UK Biobank"))
# meta_filt <- meta[!id %in% meta_done$id,]
# data.table::fwrite(meta_filt[,-("query")],
#                    "/Volumes/bms20/projects/neurogenomics-lab/live/Projects/phenome_decomposition/code/gwas_meta.tsv",
#                    sep="\t")

#### Munge ####
gwas_paths <- MungeSumstats::import_sumstats(
    ids = opt$id, 
    save_dir = file.path(root,"data","GWAS_munged"), 
    nThread = opt$ncpus,
    parallel_across_ids = FALSE, 
    force_new_vcf = FALSE,
    force_new = FALSE,
    vcf_download = TRUE,
    vcf_dir = file.path(root,"data","VCFs"),
    #### Filters ####
    dbSNP = 144,
    bi_allelic_filter = TRUE,
    tabix_index = TRUE,
    compute_z = TRUE,
    force_new_z = TRUE,
    #### Record logs ####
    log_folder_ind = TRUE,
    log_mungesumstats_msgs = TRUE
)  
#### Map ####
magma_paths <- lapply(gwas_paths, function(id){
  MAGMA.Celltyping::map_snps_to_genes(
    path_formatted = id$sumstats,
    genome_build = "GRCH37",  
    population = opt$population,
    upstream_kb = 35,  
    downstream_kb = 10, 
    force_new = FALSE
  )
})
