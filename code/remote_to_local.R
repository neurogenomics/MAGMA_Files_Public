remote_to_local <- function(
    path, 
    remote_local_prefix = c("/rds/general/project"="/Volumes/bms20/projects"),
    invert=FALSE){
    if(isTRUE(invert)){
        gsub(paste0("^",unname(remote_local_prefix)),
             names(remote_local_prefix),
             path)
    }else {
        gsub(paste0("^",names(remote_local_prefix)),
             unname(remote_local_prefix),
             path)
    }
} 