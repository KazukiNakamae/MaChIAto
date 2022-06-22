list.of.packages <- c("effsize", "ggalluvial", "dplyr", "ggplot2", "reshape2", "stringr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos='https://cloud.r-project.org/')