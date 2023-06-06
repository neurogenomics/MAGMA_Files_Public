#!/usr/bin/env Rscript
library("optparse")
#### Example ####
# Rscript write_pbs.R -s 50 -r /rds/general/project/neurogenomics-lab/live/Projects/phenome_decomposition
option_list <-list(
  optparse::make_option(c("-s", "--batch_size"), type="character", default=50, 
                        help="Number of jobs per batch.", metavar="character"),
  optparse::make_option(c("-t", "--time"), type="character", default=36, 
                        help="Max number of hours requested.", metavar="character"),
  optparse::make_option(c("-r", "--repo"), type="character", default=NULL, 
                        help="Project repository path.", metavar="character")
);
opt_parser <- optparse::OptionParser(option_list=option_list)
opt <- optparse::parse_args(opt_parser) 

#' Write pbs
#' 
#' Write a series of pbs(pro) array job submission scripts to process all samples
#' listed in the \code{metadata}. Each batch is chained together, such that 
#' the previous job automatically launches the next one. This is useful
#' as our HPC limits the number of samples per array job to 50. This function
#' avoids having to manually launch each individual job script 
#' (useful for processing 100s-1000s of samples in parallel).
#' @param batch_size Number of samples per batch.
#' @param hours Max number of hours per job.
#' @param ncpus Number of CPUs per job.
#' @param mem_gb Memory size (gigabytes) per job.
#' @param conda_env Conda environment to run job script with.
#' @param repo Path to repository where files can be found.
#' @param metadata Path to tab-delimited metadata file path containing 
#' unique IDs for each sample in the first row. 
#' Used as the iterator for the array jobs.
#' @param r_script Path to Rscript to run.
#' @param remote_local_prefix Translate path prefixes from local mounted version
#' to actual path names on the server 
#' (e.g. \code{c("/path/on/server"="/Volumes/rds/path/on/local")}).
#' @param start_next On which iteration of a given job script do 
#' you want to submit the next job script? Set \code{start_next=1} if you want
#' to launch them as soon as the prior one is launched. Set to larger number 
#' if you want runs to be run (almost) sequentially. 
#' Note, you can still have more concurrently running jobs that 
#' \code{batch_size} even when \code{start_next=batch_size} 
#' (due to some previously launched jobs not finishing 
#' before the last job is launched),
#' @param submit Submit the first job to begin the job submission chain.
#' @export
#' @examples 
#' meta <- MungeSumstats::find_sumstats()
#' ## Sort by sample size
#' meta <- dplyr::arrange(meta,dplyr::desc(N))
#' metadata <- file.path("/rds/general/project/neurogenomics-lab/live/",
#'               "Data/MAGMA_Files_Public/code/gwas_meta.tsv")
#' data.table::fwrite(meta[seq(10001,20000),-c("query")],
#'                    metadata, sep="\t")
#' batch_files <- write_pbs(metadata = metadata)
write_pbs <- function(batch_size = 50,
                      ## Set hours <8 to run in 'short' queue.
                      hours = 36, 
                      ncpus = 4,
                      mem_gb = 50,
                      queue = "med-bio",
                      conda_env = "bioc",
                      repo = "/rds/general/project/neurogenomics-lab/live/Data/MAGMA_Files_Public",
                      metadata = file.path(repo,"code/gwas_meta.tsv"), 
                      r_script = file.path(repo,"code/mungesumstats.R"),
                      save_dir = file.path(repo,"code/pbs_scripts"),
                      remote_local_prefix = c("/rds/general/project"="/Volumes/bms20/projects"),
                      start_next = batch_size,
                      on_hpc = TRUE,
                      rsync_to = NULL,#"bms20@146.169.8.44:/shared/bms20/projects/MAGMA_Files_Public/data/GWAS_munged",
                      submit = FALSE){ 
  
  #### Find max rows ####
  max_rows <- nrow(
    data.table::fread(
      remote_to_local(path = metadata, 
                      remote_local_prefix=remote_local_prefix, 
                      invert = on_hpc))
  )
  batches <- seq(0,(ceiling(max_rows/batch_size)-1))
  batch_files <- lapply(stats::setNames(batches,
                                        paste0("batch",batches)), 
                        function(batch){ 
                          #### Determine whether njobs should be capped below 50 ####
                          end <- (batch_size*batch)+batch_size
                          if(end>max_rows){
                            max_J <- end%%max_rows
                          }else {
                            max_J <- batch_size
                          };
                          lines <- list(
                            args = c(
                              paste0("#PBS -l walltime=",hours,":00:00"),
                              paste0("#PBS -l select=1:ncpus=",ncpus,":mem=",mem_gb,"gb",
                                     if(!is.null(queue)) paste(" -q",queue)
                              ),
                              ## Must be at least 2
                              paste0("#PBS -J 1-",max(max_J,2)) 
                            ),
                            conda = c(
                              "module load anaconda3/personal",
                              paste("source activate",conda_env)
                            ),
                            wd = "cd $PBS_O_WORKDIR",
                            #read line from tab-delimited metadata file
                            comment1 = "\n\n#read line from tab-delimited metadata file",
                            repo = paste("repo",repo,sep="="),
                            file = paste("file",
                                         remote_to_local(
                                           path = metadata, 
                                           remote_local_prefix = remote_local_prefix,
                                           invert = on_hpc),
                                         sep="="),
                            #identify which row to use in metadata
                            comment2 = "\n#identify which row to use in metadata",
                            batch = paste("batch",batch,sep="="),
                            row = paste("row",
                                        paste0("$((batch*",batch_size,"+PBS_ARRAY_INDEX))"),
                                        sep="="),
                            #extract relevant columns from a specific row
                            comment3 = "\n#extract relevant columns from a specific row",
                            id = paste(
                              "id",
                              "`tail -n+2 $file | awk -F'\\t' -v row=$row '{if(NR==row)print $1}'`",
                              sep="="
                            ),
                            trait = paste(
                              "trait",
                              "`tail -n+2 $file | awk -F'\\t' -v row=$row '{if(NR==row)print $2}'`",
                              sep="="
                            ),
                            r_script = paste("Rscript",r_script,
                                             "-i $id",
                                             "-n",ncpus),
                            rsync_to = if(!is.null(rsync_to)){
                              paste0("rsync -r ",
                                     "/rds/general/project/neurogenomics-lab/ephemeral/Data/MAGMA_Files_Public/data/GWAS_munged ",
                                     rsync_to)
                            },
                            #### Submit the next batch when the prior one is done ####
                            next_batch = if(batch!=max(batches)){
                              paste0(
                                "\n#Submit the next batch when the prior one is done\n",
                                paste(paste0("if (( $PBS_ARRAY_INDEX == ",min(start_next,max_J)," ))"),
                                      "\tthen",
                                      paste(
                                        "\tqsub",
                                        file.path(save_dir,paste0("batch",batch+1,".pbs"))
                                      ),
                                      "fi",
                                      sep="\n")
                              )
                            }
                          ) 
                          #### Save #### 
                          save_path <- remote_to_local(
                            path = file.path(save_dir,paste0("batch",batch,".pbs")), 
                            remote_local_prefix = remote_local_prefix, 
                            invert = on_hpc)
                          dir.create(dirname(save_path),showWarnings = FALSE, recursive = TRUE)
                          writeLines(unlist(lines), save_path, useBytes = TRUE)
                          return(list(lines=lines, 
                                      save_path=save_path))
                        }) # end lapply
  #### Submit ####
  if(isTRUE(submit)){
    #### Only need to submit first one since they're chained together ####
    system(paste("qsub",batch_files[[1]]))
  } 
  return(batch_files)   
}

# write_pbs(batch_size = opt$batch_size,
#           hours = opt$time,
#           repo = opt$repo)
