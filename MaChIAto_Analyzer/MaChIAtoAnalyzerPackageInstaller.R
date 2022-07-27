list.of.packages <- c("BiocManager", "effsize", "feather", "hash", "seqinr", "ape", "tcltk2", "rDNAse", "ggplot2", "dplyr", "reshape2", "maptools", "scales", "plyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){
    install.packages(new.packages, repos = "https://cloud.r-project.org/")
}

list.of.bioc.packages <- c("Biostrings")
new.bioc.packages <- list.of.bioc.packages[!(list.of.bioc.packages %in% installed.packages()[,"Package"])]
for(new.bioc.package in new.bioc.packages){
    BiocManager::install(new.bioc.package)
}