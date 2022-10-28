#!/usr/bin/env Rscript
library("optparse")
 
option_list <-list(
    optparse::make_option(c("-i", "--id"), type="character", default=NULL, 
              help="OpenGWAS ID", metavar="character"),
    optparse::make_option( c("-n", "--ncpus"), type="integer", default=4,
              help="OpenGWAS ID", metavar="character")
);
opt_parser <- optparse::OptionParser(option_list=option_list)
opt <- optparse::parse_args(opt_parser)
root <- "/rds/general/project/neurogenomics-lab/ephemeral/Data/MAGMA_Files_Public"

#### Prepare metadata: only process IDs not already done on OpenStack ####
# meta_done <- data.table::fread("https://github.com/neurogenomics/MAGMA_Files_Public/raw/master/metadata.csv")
# meta_done <- meta_done[meta_done$munged_path=="",]
# meta  <- MungeSumstats::find_sumstats(consortia = c("Neale Lab","UK Biobank"))
# meta_filt <- meta[!id %in% meta_done$id,]
# data.table::fwrite(meta_filt[,-("query")],
#                    "/Volumes/bms20/projects/neurogenomics-lab/live/Projects/phenome_decomposition/code/gwas_meta.tsv",
#                    sep="\t")

#### Munge sumstats ####
gwas_paths <- MungeSumstats::import_sumstats(
    ids = opt$id, 
    save_dir = file.path(root,"data/GWAS_sumstats"), 
    nThread = opt$ncpus,
    parallel_across_ids = FALSE, 
    force_new_vcf = TRUE,
    force_new = TRUE,
    vcf_download = TRUE,
    vcf_dir = file.path(root,"data/VCFs"),
    #### Filters ####
    bi_allelic_filter = TRUE,
    compute_z = TRUE,
    force_new_z = TRUE,
    tabix_index = TRUE,
    #### Record logs ####
    log_folder_ind = TRUE,
    log_mungesumstats_msgs = TRUE
)  

#### Remove temporary files ####
tmp_file <- file.path(root,"data","GWAS_sumstats",opt$id,paste0(opt$id,".tsv"))
if(file.exists(tmp_file)){
    #### Remove tsv if the compressed version also exists ####
    if(file.exists(gsub("\\.tsv",".bgz",tmp_file))){
        file.remove(tmp_file)
    }
    #### Remove vcf/index #####
    # vcf_file <- file.path(root,"data","VCFs",opt$id,paste0(opt$id,".vcf.gz"))
    # tbi_file <- paste0(vcf_file,".tbi")
    # if(file.exists(vcf_file)){
    #     file.remove(vcf_file)
    # }
    # if(file.exists(tbi_file)){
    #     file.remove(tbi_file)
    # }
}


