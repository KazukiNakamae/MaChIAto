list.of.packages <- c("feather", "hash", "seqinr", "ape", "tcltk2", "rDNAse", "maptools", "ggpubr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){
    install.packages(new.packages, repos = "https://cran.ism.ac.jp/")
}

#list.of.bioc.packages <- c("Biostrings")
#new.bioc.packages <- list.of.bioc.packages[!(list.of.bioc.packages %in% installed.packages()[,"Package"])]
#for(new.bioc.package in new.bioc.packages){
#    source("http://bioconductor.org/biocLite.R")
#    biocLite(new.bioc.package)
#    #BiocManager::install(new.bioc.package)
#}