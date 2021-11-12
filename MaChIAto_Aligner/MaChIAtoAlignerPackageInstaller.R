

list.of.packages <- c("BiocManager", "grDevices", "network", "dplyr", "reshape2", "stringr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos='https://cloud.r-project.org/')

list.of.packages <- c("Biostrings", "Rsamtools", "GenomicRanges", "CrispRVariants")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)

# library("devtools")

# list.of.packages <- c("CrispRVariants")
# new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
# if(length(new.packages)) devtools::install_github("markrobinsonuzh/CrispRVariants", dependencies=FALSE)

# list.of.packages <- c("CrispRVariants")
# new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
# if(length(new.packages)) devtools::install_url("markrobinsonuzh/CrispRVariants")


# install.packages(new.packages, repos='https://cran.ism.ac.jp/')