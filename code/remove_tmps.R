#' Remove temporary files
#' 
#' Remove temporary files created by MungeSumstats.
#' \strong{WARNING}:: Only do this after all jobs are done running. 
#' Otherwise, could accidentally delete tmp/vcf files as they're being used.
#' @export
#' @importFrom MungeSumstats list_sumstats
#' @export
#' @examples 
#' ### Find most recent vcfs
#' system("find MAGMA_Files_Public/data/VCFs/ -daystart -maxdepth 1 -mmin +1500 -type f")
#' remove_tmps(remove_vcfs=TRUE)
remove_tmps <- function(repo = "/rds/general/project/neurogenomics-lab/ephemeral/Data/MAGMA_Files_Public",
                        target_pattern = ".bgz",
                        residual_pattern = "\\.tsv$",
                        gwas_folder = file.path(repo,"data","GWAS_sumstats"),
                        vcf_folder = file.path(repo,"data","VCFs"),
                        remove_vcfs = FALSE){ 
    
    ### Remove residual temp files ####
    tmp_files <- MungeSumstats::list_sumstats(save_dir = gwas_folder,
                                              pattern = residual_pattern)
    tmp_files <- tmp_files[file.exists(paste0(tmp_files,target_pattern))]
    indexed_files <- file.exists(paste0(tmp_files,target_pattern,".tbi"))
    file.remove(tmp_files)
    if(isTRUE(remove_vcfs)){
        unlink(vcf_folder,force = TRUE, recursive = TRUE)   
    }
    return(list(tmp_files=tmp_files,
                indexed_files=indexed_files,
                vcf_folder=vcf_folder))
}
