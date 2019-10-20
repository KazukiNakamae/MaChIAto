# cd /Users/kazuki/Documents/Lab/MachinelearningForPITChprediction/20190624_develop_MaChIAto_Reviewer;
# source activate MaChIAto_Aligner;
# r
# source("./MaChIAto_Reviewer/MaChIAtoReviewerCore.R")
# quit(save = "no") # for debug
###################################################################################################################################################
# Version β2.0
#-- Setting ----------------------------------------------------#

message("----------------------------")
message("---Start MaChIAto_Reviewer---")
message("----------------------------")

debug.flag <- FALSE

# Get directory path
if(debug.flag == TRUE){
    script.dir <- "./MaChIAto_Reviewer"
}else{
    initial.options <- commandArgs(trailingOnly = FALSE)
    file.arg.name <- "--file="
    script.fn <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
    script.dir <- dirname(script.fn)
}

#-- Run Script ----------------------------------------------------#
message("---Load MaChIAto_Reviewer programs---")

source(file.path(script.dir, "MaChIAtoReviewerPackageInstaller.R"))
source(file.path(script.dir, "MaChIAtoReviewerFunctions.R"))
#source(file.path(script.dir, "RunLocalMutationAnalysis.R"))
#source(file.path(script.dir, "RunLocalKIAnalysis.R"))

#-- Load library ----------------------------------------------------#
message("---Load libraries---")

library(CrispRVariants)
#library(GenomicFeatures)
library(Biostrings)
#library(GenomicRanges)
#library(Rsamtools)
library(ggalluvial)
#library(grDevices)
#library(network)
library(dplyr)
library(ggplot2)
library(reshape2)
library(stringr)

if(debug.flag == TRUE){
    classifier.res.dir <- "/Users/kazuki/Documents/Lab/MachinelearningForPITChprediction/20190621_run_MaChIAto_Aligner/MaChIAto_v1.1.0_output_190710_2019y50data_th0995"
    # aligner.res.dir <- "/Volumes/databank1/MaChIAtoAligner_output_from_MaChIAto_output_190508"
    aligner.res.dir <- "/Volumes/databank1/MaChIAtoAlignerv12_output_190726_from_MaChIAto_v1.1.0_output_190710_2019y50data_th0995"
    output.prefix <- "./MaChIAtoReviewer_test_output"
    untreated.label <- "A"
    negative.ctrl.label <- "B"
    condition1.label <- "C"
    condition2.label <- "D"
    ignore.list <- "./ignore.list.csv"
}else{
    # Get commandline arguments
    classifier.res.dir <- GetAbsolutePath(commandArgs(trailingOnly=TRUE)[1])
    aligner.res.dir <- GetAbsolutePath(commandArgs(trailingOnly=TRUE)[2])
    output.prefix <- GetAbsolutePath(commandArgs(trailingOnly=TRUE)[3])
    untreated.label <- commandArgs(trailingOnly=TRUE)[4]
    negative.ctrl.label <- commandArgs(trailingOnly=TRUE)[5]
    condition1.label <- commandArgs(trailingOnly=TRUE)[6]
    condition2.label <- commandArgs(trailingOnly=TRUE)[7]
    ignore.list <- GetAbsolutePath(commandArgs(trailingOnly=TRUE)[8])
}

################################################################################################################################
### result file information data frame
################################################################################################################################
classifier.res.path.vec <- list.dirs(path = file.path(classifier.res.dir), full.names = TRUE, recursive = FALSE)
classifier.res.sample.name.vec <- sapply(classifier.res.path.vec, function(x){tail(strsplit(x, "_")[[1]], n=1)})
classifier.res.sample.target.vec <- sapply(classifier.res.sample.name.vec, function(x){head(strsplit(x, "-")[[1]], n=1)})
classifier.res.sample.label.vec <- sapply(classifier.res.sample.name.vec, function(x){tail(strsplit(x, "-")[[1]], n=1)})
classifier.res.file.df <- data.frame(
  classifier.res.dir.path = classifier.res.path.vec
  , sample.name = classifier.res.sample.name.vec
  , sample.target = classifier.res.sample.target.vec
  , sample.label = classifier.res.sample.label.vec
  , stringsAsFactors = FALSE
)
aligner.res.path.vec <- list.dirs(path = file.path(aligner.res.dir), full.names = TRUE, recursive = FALSE)
aligner.res.sample.name.vec <- sapply(aligner.res.path.vec, function(x){tail(strsplit(x, "_")[[1]], n=1)})
aligner.res.file.df <- data.frame(
  aligner.res.dir.path = aligner.res.path.vec
  , sample.name = aligner.res.sample.name.vec
  , stringsAsFactors = FALSE
)
res.file.df <- merge(classifier.res.file.df, aligner.res.file.df, by = "sample.name", all = TRUE)
################################################################################################################################
### Filtering
################################################################################################################################
message("Filtering data")
if(ignore.list != ""){
  remove.target.vec <- unlist(read.csv(ignore.list, header = FALSE, stringsAsFactors = FALSE)[1,])
}else{
  remove.target.vec <- character(0)
}

# Remove gene contains low number of reads
message("Check amount of reads...")
has.sufficint.reads.vec <- apply(res.file.df, MARGIN = 1, function(row){
  all.df <- read.csv(file.path(row["classifier.res.dir.path"], "ALL_dataframe.csv"), header = TRUE)
  #crispresso.unmodified.df <- filter(all.df, UNMODIFIED == "True")
  #crispresso.nhej.df <- filter(all.df, NHEJ == "True")
  #crispresso.mixed.df <- filter(all.df, UNMODIFIED == "False" & NHEJ == "False" & HDR == "False")
  #crispresso.hdr.df <- filter(all.df, HDR == "True")
  #machiato.unclassified.df <- filter(all.df, CRISPResso_reclassification_labels == "Unclassified")
  #machiato.unmodified.df <- filter(all.df, CRISPResso_reclassification_labels == "Unmodified")
  #machiato.nhej.df <- filter(all.df, CRISPResso_reclassification_labels == "NHEJ")
  #machiato.mixed.df <- filter(all.df, CRISPResso_reclassification_labels == "Mixed HDR-NHEJ")
  #machiato.hdr.df <- filter(all.df, CRISPResso_reclassification_labels == "HDR")

  if(sum(all.df$X.Reads) < 1000){
    message(paste0(row["sample.name"], " dosen't have sufficint reads (<1,000 reads)."))
    return(FALSE)
  }else{
    return(TRUE)
  }
})
remove.target.vec <- c(remove.target.vec, unique(res.file.df$sample.target[!has.sufficint.reads.vec]))

# Remove gene contains unusual Unmodified rate in untreated
message("Check conservative rate in untreated samples...")
is.conservative.sample.vec <- rep(TRUE, nrow(res.file.df))
if(untreated.label != ""){ # untreated label is entered.
  is.conservative.sample.vec <- apply(res.file.df, MARGIN = 1, function(row){
    if(row["sample.label"] != untreated.label){
      return(TRUE) # this is no target of filtering
    }

    all.df <- read.csv(file.path(row["classifier.res.dir.path"], "ALL_dataframe.csv"), header = TRUE)
    crispresso.unmodified.df <- filter(all.df, UNMODIFIED == "True")
    if(sum(crispresso.unmodified.df$X.Reads) / sum(all.df$X.Reads) < 0.8){
      message(paste0(row["sample.name"], " is not conservative (the rate of unmodified reads < 80%)."))
      return(FALSE)
    }else{
      return(TRUE)
    }
  })
}
remove.target.vec <- c(remove.target.vec, unique(res.file.df$sample.target[!is.conservative.sample.vec]))

# Remove gene contains no difference between wt and negative.ctrl on untreated
message("Check whether efficary of editing exists in negative ctrl...")
if(untreated.label != "" & negative.ctrl.label != ""){ # untreated label is entered.
  for(focus.sample.target in unique(res.file.df$sample.target)){

    # load untreated
    temp.untreated.table.path <- filter(res.file.df, sample.target == focus.sample.target & sample.label == untreated.label)$classifier.res.dir.path
    if(length(temp.untreated.table.path) > 1){
      stop(paste0(focus.sample.target, "-", untreated.label, " is not unique in ", classifier.res.dir, ". Check input."))
    }else if(length(temp.untreated.table.path) < 1){
      stop(paste0(focus.sample.target, "-", untreated.label, " dosen't exist in ", classifier.res.dir, ". Check input."))
    }
    temp.untreated.all.table <- read.csv(file.path(temp.untreated.table.path, "ALL_dataframe.csv"), header = TRUE)
    temp.untreated.crispresso.unmodified.df <- filter(temp.untreated.all.table, UNMODIFIED == "True")

    # load negative control
    temp.negative.ctrl.table.path <- filter(res.file.df, sample.target == focus.sample.target & sample.label == negative.ctrl.label)$classifier.res.dir.path
    if(length(temp.negative.ctrl.table.path) > 1){
      stop(paste0(focus.sample.target, "-", negative.ctrl.label, " is not unique in ", classifier.res.dir, ". Check input."))
    }else if(length(temp.negative.ctrl.table.path) < 1){
      stop(paste0(focus.sample.target, "-", negative.ctrl.label, " dosen't exist in ", classifier.res.dir, ". Check input."))
    }
    temp.negative.ctrl.all.table <- read.csv(file.path(temp.negative.ctrl.table.path, "ALL_dataframe.csv"), header = TRUE)
    temp.negative.ctrl.crispresso.unmodified.df <- filter(temp.negative.ctrl.all.table, UNMODIFIED == "True")

    temp.untreated.vs.negativectrl.c2.table <- matrix(c(
      sum(temp.untreated.crispresso.unmodified.df$X.Reads),
      sum(temp.untreated.all.table$X.Reads) - sum(temp.untreated.crispresso.unmodified.df$X.Reads),
      sum(temp.negative.ctrl.crispresso.unmodified.df$X.Reads),
      sum(temp.negative.ctrl.all.table$X.Reads) - sum(temp.negative.ctrl.crispresso.unmodified.df$X.Reads)),
      ncol=2,
      byrow=T
    )
    options(warn=-1)
    temp.untreated.vs.negativectrl.chisq.test.res <- chisq.test(temp.untreated.vs.negativectrl.c2.table)
    options(warn=0)
    if(is.nan(temp.untreated.vs.negativectrl.chisq.test.res$p.value) |
      temp.untreated.vs.negativectrl.chisq.test.res$p.value >= 0.001){
        message(paste0(focus.sample.target, " dosen't have significant change between untreated and negative.ctrl. (chisq.test, p >= 0.1%)"))
        remove.target.vec <- c(remove.target.vec, focus.sample.target)
    }
    
  }
}

# Remove gene contains no precise knock-in in condition1 & condition2 sample
message("Check whether precise knock-in exists in condition1 and condition2...")
if(condition1.label != "" & condition2.label != ""){ # untreated label is entered.
  for(focus.sample.target in unique(res.file.df$sample.target)){

    # load condition1
    temp.condition1.table.path <- filter(res.file.df, sample.target == focus.sample.target & sample.label == condition1.label)$classifier.res.dir.path
    if(length(temp.condition1.table.path) > 1){
      stop(paste0(focus.sample.target, "-", condition1.label, " is not unique in ", classifier.res.dir, ". Check input."))
    }else if(length(temp.condition1.table.path) < 1){
      stop(paste0(focus.sample.target, "-", condition1.label, " dosen't exist in ", classifier.res.dir, ". Check input."))
    }
    temp.condition1.all.table <- read.csv(file.path(temp.condition1.table.path, "ALL_dataframe.csv"), header = TRUE)
    temp.condition1.crispresso.hdr.df <- filter(temp.condition1.all.table, HDR == "True")

    # load condition2
    temp.condition2.table.path <- filter(res.file.df, sample.target == focus.sample.target & sample.label == condition2.label)$classifier.res.dir.path
    if(length(temp.condition2.table.path) > 1){
      stop(paste0(focus.sample.target, "-", condition2.label, " is not unique in ", classifier.res.dir, ". Check input."))
    }else if(length(temp.condition2.table.path) < 1){
      stop(paste0(focus.sample.target, "-", condition2.label, " dosen't exist in ", classifier.res.dir, ". Check input."))
    }
    temp.condition2.all.table <- read.csv(file.path(temp.condition2.table.path, "ALL_dataframe.csv"), header = TRUE)
    temp.condition2.crispresso.hdr.df <- filter(temp.condition2.all.table, HDR == "True")


    if(sum(temp.condition1.crispresso.hdr.df$X.Reads) < 1 | sum(temp.condition2.crispresso.hdr.df$X.Reads) < 1){
      message(paste0(focus.sample.target, " has no precise knock-in in condition1 or condition2."))
      remove.target.vec <- c(remove.target.vec, focus.sample.target)
    }

  }
}

remove.target.vec <- unique(remove.target.vec)
message("Remove...")
print(remove.target.vec)
fltr.res.file.df <- filter(res.file.df, !(sample.target %in% remove.target.vec))
message("Filtering is done.")

################################################################################################################################
### setting for saving
################################################################################################################################

#-- Make output directory ----------------------------------------------------#
message("---Make directory for the analysis---")

#time.id <- "20190404150926"
time.id <- format(Sys.time(), "%Y%m%d%H%M%S")
output.name <- paste("MaChIAtoReviewer_at_", time.id, sep = "")
dir.create(file.path(output.prefix), showWarnings=FALSE)
output.dir <- file.path(output.prefix, output.name)

message(paste("Creating Folder", output.dir, sep=" "))

### result data dir
comparison.analysis.dir <- file.path(output.dir, "[A]Validation_of_CRISPResso-MaChIAto_classification")
comparison.analysis.original.dir <- file.path(comparison.analysis.dir, "[a1]Piechart_of_CRISPResso-MaChIAto_classification")
comparison.analysis.machiato.dir <- file.path(comparison.analysis.dir, "[b1]Piechart_of_original_CRISPResso_classification_filtered_by_crispRvariants")
mutation.analysis.dir <- file.path(output.dir, "[B]Validation_of_mutation_among_endgenous_sequence")
imprecise.knockin.analysis.dir <- file.path(output.dir, "[C]Validation_of_junction_sequence_in_imprecise_knock-in")
###


for(dirname in c(output.dir, comparison.analysis.dir, comparison.analysis.original.dir, comparison.analysis.machiato.dir,
  mutation.analysis.dir, imprecise.knockin.analysis.dir)){
  dir.create(file.path(dirname), showWarnings=FALSE)
}

################################################################################################################################
### Make summary of classification
################################################################################################################################

########[a] Sum of CRISPResso Classification
message("[a]Make a summary of CRISPResso classification")

temp.crispresso.rate.table.mat <- apply(fltr.res.file.df, MARGIN = 1, function(row){
  temp.all.df <- read.csv(file.path(row["classifier.res.dir.path"], "ALL_dataframe.csv"), header = TRUE)
  temp.crispresso.unmodified.df <- filter(temp.all.df, UNMODIFIED == "True")
  temp.crispresso.nhej.df <- filter(temp.all.df, NHEJ == "True")
  temp.crispresso.mixed.df <- filter(temp.all.df, UNMODIFIED == "False" & NHEJ == "False" & HDR == "False")
  temp.crispresso.hdr.df <- filter(temp.all.df, HDR == "True")

  crispresso.class.count.sum.vec <- c(sum(temp.crispresso.unmodified.df$X.Reads)
    , sum(temp.crispresso.nhej.df$X.Reads)
    , sum(temp.crispresso.mixed.df$X.Reads)
    , sum(temp.crispresso.hdr.df$X.Reads)
  )
  names(crispresso.class.count.sum.vec) <- c("Unmodified", "InDel", "Imprecise knock-in", "Precise knock-in")
  crispresso.class.count.table <- data.frame(
    Class=names(crispresso.class.count.sum.vec)
    , Frequency=crispresso.class.count.sum.vec
    , Percentage=crispresso.class.count.sum.vec / sum(crispresso.class.count.sum.vec) * 100
  )

  SaveTable(crispresso.class.count.table, file.path(comparison.analysis.original.dir, paste0(row["sample.name"], "_CRISPResso_classification_table")))
  
  SavePieChart(crispresso.class.count.table
    , c("#D98F4E",  "#A65437", "#F2CA52", "#ADD4D9", "#F2E2CE")
    , "Reads"
    , file.path(comparison.analysis.original.dir, paste0(row["sample.name"], "_CRISPResso_classification_piechart.png"))
  )

  return(crispresso.class.count.table$Percentage)
})
rownames(temp.crispresso.rate.table.mat) <- c("Unmodified", "InDel", "Imprecise knock-in", "Precise knock-in")
colnames(temp.crispresso.rate.table.mat) <- fltr.res.file.df$sample.name

MakeSummaryClassBarplot(
  temp.crispresso.rate.table.mat,
  condition2.label,
  c("#BFBFBF", "#548C1C", "#A60321", "#F20C36"),
  paste0(untreated.label, ":Unmodified   ", negative.ctrl.label,":Negative Ctrl   ", condition1.label, ":Condition1   ", condition2.label, ":Condition2"),
  file.path(comparison.analysis.dir, "[a2]Summary_of_original_CRISPResso_classification.png")
)

########[b] Sum of MaChIAto Classification
message("[b]Make a summary of MaChIAto classification")

temp.machiato.rate.table.mat <- apply(fltr.res.file.df, MARGIN = 1, function(row){
  if(row["sample.label"] %in% c(untreated.label, negative.ctrl.label, condition1.label, condition2.label)){

    temp.all.df <- read.csv(file.path(row["classifier.res.dir.path"], "ALL_dataframe.csv"), header = TRUE)
    # temp.machiato.unclassified.df <- filter(temp.all.df, CRISPResso_reclassification_labels == "Unclassified")
    temp.machiato.unmodified.df <- filter(temp.all.df, CRISPResso_reclassification_labels == "Unmodified")
    temp.machiato.nhej.df <- filter(temp.all.df, CRISPResso_reclassification_labels == "NHEJ")
    temp.machiato.mixed.df <- filter(temp.all.df, CRISPResso_reclassification_labels == "Mixed HDR-NHEJ")
    temp.machiato.hdr.df <- filter(temp.all.df, CRISPResso_reclassification_labels == "HDR")

    machiato.class.count.sum.vec <- c(sum(temp.machiato.unmodified.df$X.Reads)
      , sum(temp.machiato.nhej.df$X.Reads)
      , sum(temp.machiato.mixed.df$X.Reads)
      , sum(temp.machiato.hdr.df$X.Reads)
    )
    names(machiato.class.count.sum.vec) <- c("Unmodified", "InDel", "Imprecise knock-in", "Precise knock-in")
    machiato.class.count.table <- data.frame(
      Class=names(machiato.class.count.sum.vec)
      , Frequency=machiato.class.count.sum.vec
      , Percentage=machiato.class.count.sum.vec / sum(machiato.class.count.sum.vec) * 100
    )

    SaveTable(machiato.class.count.table, file.path(comparison.analysis.machiato.dir, paste0(row["sample.name"], "_MaChIAto_classification_table")))
    
    SavePieChart(machiato.class.count.table
      , c("#D98F4E",  "#A65437", "#F2CA52", "#ADD4D9", "#F2E2CE")
      , "Reads"
      , file.path(comparison.analysis.machiato.dir, paste0(row["sample.name"], "_MaChIAto_classification_piechart.png"))
    )

    return(machiato.class.count.table$Percentage)

  }else{
    return(rep(NA, 4))
  }
})
temp.condition12.machiato.rate.table.mat <- temp.machiato.rate.table.mat[,complete.cases(t(temp.machiato.rate.table.mat))]
rownames(temp.condition12.machiato.rate.table.mat) <- c("Unmodified", "InDel", "Imprecise knock-in", "Precise knock-in")
colnames(temp.condition12.machiato.rate.table.mat) <- fltr.res.file.df$sample.name[complete.cases(t(temp.machiato.rate.table.mat))]

MakeSummaryClassBarplot(
  temp.condition12.machiato.rate.table.mat,
  condition2.label,
  c("#BFBFBF", "#548C1C", "#A60321", "#F20C36"),
  paste0(untreated.label, ":Unmodified   ", negative.ctrl.label,":Negative Ctrl   ", condition1.label, ":Condition1   ", condition2.label, ":Condition2"),
  file.path(comparison.analysis.dir, "[b2]Summary_of_CRISPResso_MaChIAto_classification.png")
)

########[c] Sum of CRISPResso - BWAMEM-CrispRvariants Classification
message("[c]Make a summary of CRISPResso - BWAMEM-CrispRvariants classification")

res.1IBi.path <- "result/[1]comparison_with_original_crispresso_data/[Ⅰ]crispresso_alignment_data/[B]crispresso_rate_data/[ⅰ]CRISPResso_Unmodified_NHEJ_HDR_table.rds"
temp.crispresso.bwa.rate.table.mat <- apply(fltr.res.file.df, MARGIN = 1, function(row){
  if(row["sample.label"] %in% c(untreated.label, negative.ctrl.label, condition1.label, condition2.label)){
    return(readRDS(file.path(row["aligner.res.dir.path"], res.1IBi.path))$Percentage)
  }else{
    return(rep(NA, 4))
  }
})
temp.condition12.crispresso.bwa.rate.table.mat <- temp.crispresso.bwa.rate.table.mat[,complete.cases(t(temp.crispresso.bwa.rate.table.mat))]
rownames(temp.condition12.crispresso.bwa.rate.table.mat) <- c("Unmodified", "InDel", "Imprecise knock-in", "Precise knock-in")
colnames(temp.condition12.crispresso.bwa.rate.table.mat) <- fltr.res.file.df$sample.name[complete.cases(t(temp.crispresso.bwa.rate.table.mat))]

MakeSummaryClassBarplot(
  temp.condition12.crispresso.bwa.rate.table.mat,
  condition2.label,
  c("#BFBFBF", "#548C1C", "#A60321", "#F20C36"),
  paste0(untreated.label, ":Unmodified   ", negative.ctrl.label,":Negative Ctrl   ", condition1.label, ":Condition1   ", condition2.label, ":Condition2"),
  file.path(comparison.analysis.dir, "[c]Summary_of_CRISPResso_BWA-MEM_CrispRVariants_classification.png")
)

########[d] Sum of MaChIAto - BWAMEM-CrispRvariants Classification
message("[d]Make a summary of MaChIAto - BWAMEM-CrispRvariants classification")

res.1IIBii.path <- "result/[1]comparison_with_original_crispresso_data/[Ⅱ]machiato_alignment_data/[B]machiato_rate_data/[ⅱ]CRISPResso-MaChIAto_Unmodified_NHEJ_HDR_table.rds"
temp.machiato.bwa.rate.table.mat <- apply(fltr.res.file.df, MARGIN = 1, function(row){
  if(row["sample.label"] %in% c(untreated.label, negative.ctrl.label, condition1.label, condition2.label)){
    return(readRDS(file.path(row["aligner.res.dir.path"], res.1IIBii.path))$Percentage)
  }else{
    return(rep(NA, 4))
  }
})
temp.condition12.machiato.bwa.rate.table.mat <- temp.machiato.bwa.rate.table.mat[,complete.cases(t(temp.machiato.bwa.rate.table.mat))]
rownames(temp.condition12.machiato.bwa.rate.table.mat) <- c("Unmodified", "InDel", "Imprecise knock-in", "Precise knock-in")
colnames(temp.condition12.machiato.bwa.rate.table.mat) <- fltr.res.file.df$sample.name[complete.cases(t(temp.machiato.bwa.rate.table.mat))]

MakeSummaryClassBarplot(
  temp.condition12.machiato.bwa.rate.table.mat,
  condition2.label,
  c("#BFBFBF", "#548C1C", "#A60321", "#F20C36"),
  paste0(untreated.label, ":Unmodified   ", negative.ctrl.label,":Negative Ctrl   ", condition1.label, ":Condition1   ", condition2.label, ":Condition2"),
  file.path(comparison.analysis.dir, "[d]Summary_of_MaChIAto_BWA-MEM_CrispRVariants_classification.png")
)

########[e] check defference between all and unclassified
message("[e]Check KL divergence with MaChIAto unclassified reads")

res.1IIAa.path <- "result/[1]comparison_with_original_crispresso_data/[Ⅱ]machiato_alignment_data/[Aa]machiato_all/crispr.set.rds"
machiato.all.crisprset.list <- MakeCrisprsetList(fltr.res.file.df, res.1IIAa.path)
res.1IIAb.path <- "result/[1]comparison_with_original_crispresso_data/[Ⅱ]machiato_alignment_data/[Ab]machiato_unmodified/crispr.set.rds"
machiato.unmodified.crisprset.list <- MakeCrisprsetList(fltr.res.file.df, res.1IIAb.path)
res.1IIAc.path <- "result/[1]comparison_with_original_crispresso_data/[Ⅱ]machiato_alignment_data/[Ac]machiato_nhej/crispr.set.rds"
machiato.nhej.crisprset.list <- MakeCrisprsetList(fltr.res.file.df, res.1IIAc.path)
res.1IIAd.path <- "result/[1]comparison_with_original_crispresso_data/[Ⅱ]machiato_alignment_data/[Ad]machiato_mixed_hdr_nhej/crispr.set.rds"
machiato.mixed_hdr_nhej.crisprset.list <- MakeCrisprsetList(fltr.res.file.df, res.1IIAd.path)
res.1IIAe.path <- "result/[1]comparison_with_original_crispresso_data/[Ⅱ]machiato_alignment_data/[Ae]machiato_hdr/crispr.set.rds"
machiato.hdr.crisprset.list <- MakeCrisprsetList(fltr.res.file.df, res.1IIAe.path)
res.1IIAf.path <- "result/[1]comparison_with_original_crispresso_data/[Ⅱ]machiato_alignment_data/[Af]machiato_unclassified/crispr.set.rds"
machiato.unclassified.crisprset.list <- MakeCrisprsetList(fltr.res.file.df, res.1IIAf.path)

unclassified.comp.kl.mat <- matrix(, nrow = 5)
for(ind in 1:length(machiato.unclassified.crisprset.list)){
  temp.sample.name <- names(machiato.unclassified.crisprset.list[[ind]])
  # search each index
  machiato.all.ind <- which(unlist(lapply(machiato.all.crisprset.list, names)) %in% temp.sample.name)
  machiato.unmodified.ind <- which(unlist(lapply(machiato.unmodified.crisprset.list, names)) %in% temp.sample.name)
  machiato.nhej.ind <- which(unlist(lapply(machiato.nhej.crisprset.list, names)) %in% temp.sample.name)
  machiato.mixed_hdr_nhej.ind <- which(unlist(lapply(machiato.mixed_hdr_nhej.crisprset.list, names)) %in% temp.sample.name)
  machiato.hdr.ind <- which(unlist(lapply(machiato.hdr.crisprset.list, names)) %in% temp.sample.name)
  # Calc divergence
  temp.kl.vec <- c(
    ALL = CalcKLdivergence(
      SafeGetElement(machiato.unclassified.crisprset.list[[ind]], temp.sample.name),
      SafeGetElement(machiato.all.crisprset.list[[machiato.all.ind]], temp.sample.name)
    ),
    Unmodified = CalcKLdivergence(
      SafeGetElement(machiato.unclassified.crisprset.list[[ind]], temp.sample.name),
      SafeGetElement(machiato.unmodified.crisprset.list[[machiato.unmodified.ind]], temp.sample.name)
    ),
    InDel = CalcKLdivergence(
      SafeGetElement(machiato.unclassified.crisprset.list[[ind]], temp.sample.name),
      SafeGetElement(machiato.nhej.crisprset.list[[machiato.nhej.ind]], temp.sample.name)
    ),
    Imprecise_knockin = CalcKLdivergence(
      SafeGetElement(machiato.unclassified.crisprset.list[[ind]], temp.sample.name),
      SafeGetElement(machiato.mixed_hdr_nhej.crisprset.list[[machiato.mixed_hdr_nhej.ind]], temp.sample.name)
    ),
    Precise_knockin = CalcKLdivergence(
      SafeGetElement(machiato.unclassified.crisprset.list[[ind]], temp.sample.name),
      SafeGetElement(machiato.hdr.crisprset.list[[machiato.hdr.ind]], temp.sample.name)
    )
  )
  unclassified.comp.kl.mat <- cbind(unclassified.comp.kl.mat, temp.kl.vec)
  colnames(unclassified.comp.kl.mat)[ind + 1] <- temp.sample.name
}
unclassified.comp.kl.mat <- unclassified.comp.kl.mat[,-1]
# Adjuct table for heatmap
melt.unclassified.comp.kl.mat <- melt(unclassified.comp.kl.mat)
melt.unclassified.comp.kl.mat$Var2 <- factor(melt.unclassified.comp.kl.mat$Var2, level = fltr.res.file.df$sample.name[complete.cases(t(temp.machiato.rate.table.mat))])
melt.unclassified.comp.kl.mat$Var1 <- factor(melt.unclassified.comp.kl.mat$Var1, level = c("ALL", "Unmodified", "InDel", "Imprecise_knockin", "Precise_knockin"))
# Make table for heatmap
options(warn=-1)
unclassified.comp.kl.p <- ggplot(melt.unclassified.comp.kl.mat
  , aes(x=Var1, y=Var2, fill=value)) +
geom_raster(aes(fill = value)) +
labs(title ="Symmetrized KL divergence with Unclassified Class"
  , x = "Class"
  , y = "Sample") +
theme(axis.ticks = element_blank()
  , axis.text.x = element_text(angle = 330, hjust = 0, size=20)
  , axis.text.y = element_text(size=7)
  , plot.title = element_text(size = 10)  
  , axis.title.x = element_text(size = 20)  
  , axis.title.y = element_text(size = 20)
  , legend.title = element_text(size = 20)
  , legend.text = element_text(size = 15)
) +
scale_fill_distiller(palette = "RdYlBu", name = "Symmetrized KL Divergence")

SaveTable(unclassified.comp.kl.mat, file.path(comparison.analysis.dir, "[e]Comparison_between_Unclassified_and_other_classes"))
ggsave(file = file.path(comparison.analysis.dir, "[e]Comparison_between_Unclassified_and_other_classes.png")
  , plot = unclassified.comp.kl.p, dpi = 300, width = 10, height = 15)
options(warn=0)

########[f] check consistency barplot by each class
message("[f]Check consistency between CRISPResso and MaChIAto")

consistency.mat <- apply(fltr.res.file.df, MARGIN = 1, function(row){
  temp.all.df <- read.csv(file.path(row["classifier.res.dir.path"], "ALL_dataframe.csv"), header = TRUE)
  # temp.machiato.unclassified.df <- filter(temp.all.df, CRISPResso_reclassification_labels == "Unclassified")
  temp.machiato.unmodified.df <- filter(temp.all.df, CRISPResso_reclassification_labels == "Unmodified")
  temp.machiato.nhej.df <- filter(temp.all.df, CRISPResso_reclassification_labels == "NHEJ")
  temp.machiato.mixed.df <- filter(temp.all.df, CRISPResso_reclassification_labels == "Mixed HDR-NHEJ")
  temp.machiato.hdr.df <- filter(temp.all.df, CRISPResso_reclassification_labels == "HDR")
  # temp.crispresso.unmodified.df <- filter(temp.all.df, UNMODIFIED == "True")
  # temp.crispresso.nhej.df <- filter(temp.all.df, NHEJ == "True")
  # temp.crispresso.mixed.df <- filter(temp.all.df, UNMODIFIED == "False" & NHEJ == "False" & HDR == "False")
  # temp.crispresso.hdr.df <- filter(temp.all.df, HDR == "True")

  temp.consistency <- c(sum(filter(temp.machiato.unmodified.df, UNMODIFIED == "True")$X.Reads) /
    sum(temp.machiato.unmodified.df$X.Reads) * 100,
    sum(filter(temp.machiato.nhej.df, NHEJ == "True")$X.Reads) /
    sum(temp.machiato.nhej.df$X.Reads) * 100,
    sum(filter(temp.machiato.mixed.df, UNMODIFIED == "False" & NHEJ == "False" & HDR == "False")$X.Reads) /
    sum(temp.machiato.mixed.df$X.Reads) * 100,
    sum(filter(temp.machiato.hdr.df, HDR == "True")$X.Reads) /
    sum(temp.machiato.hdr.df$X.Reads) * 100)

  return(temp.consistency)
})
fltred.consistency.mat <- consistency.mat[,complete.cases(t(consistency.mat))]
rownames(fltred.consistency.mat) <- c("Unmodified", "InDel", "Imprecise knock-in", "Precise knock-in")
colnames(fltred.consistency.mat) <- fltr.res.file.df$sample.name[complete.cases(t(consistency.mat))]
SaveTable(fltred.consistency.mat, file.path(comparison.analysis.dir, "[f]Boxplot_of_Rate_of_Consistency"))

melt.consistency.df <- melt(fltred.consistency.mat)
options(warn=-1)
consistency.box.plot <- ggplot(melt.consistency.df, aes(x=Var1, y=value, fill = factor(Var1))) + 
  geom_boxplot(width = 0.5, alpha = 0.6, show.legend = FALSE, outlier.shape = NA) + 
  geom_jitter(width = 0.25, alpha=0.9) +
  coord_cartesian(ylim = c(0, 100)) +
  scale_fill_manual(values = c("#BFBFBF", "#548C1C", "#A60321", "#F20C36")) +
  theme(legend.position = "none") +
  labs(x = "", y = "Percentage[%]") +
  scale_y_continuous(breaks=seq(0, c(0, 100)[2], 5))
ggsave(file = file.path(comparison.analysis.dir, "[f]Boxplot_of_Rate_of_Consistency.png"), plot = consistency.box.plot, dpi = 350, width = 10, height = 10)
options(warn=0)

############################################################

########[extra] Alluvial diagrams between CRISPResso and MaChIAto
message("[extra]Make an alluvial diagrams between CRISPResso and MaChIAto")

MakeAlluvialDiagramsCRISPRessoMaChIAto()


########[extra] CRISPResso and MaChIAto All Variants Barplot
message("[extra] CRISPResso and MaChIAto All Variants Barplot")

res.1IAa.path <- "result/[1]comparison_with_original_crispresso_data/[Ⅰ]crispresso_alignment_data/[Aa]crispresso_all/crispr.set.rds"
crispresso.all.variants.cnt.simple.list <- MakeCrisprsetListSimple(fltr.res.file.df, res.1IAa.path)
SaveVariantsBarPlot(crispresso.all.variants.cnt.simple.list, file.path(comparison.analysis.dir, "[extra]CRISPResso_ALL_Variants_Barplot"))

machiato.unmodified.variants.cnt.simple.list <- MakeCrisprsetListSimple(fltr.res.file.df, res.1IIAb.path)
machiato.nhej.variants.cnt.simple.list <- MakeCrisprsetListSimple(fltr.res.file.df, res.1IIAc.path)
machiato.mixed.variants.cnt.simple.list <- MakeCrisprsetListSimple(fltr.res.file.df, res.1IIAd.path)
machiato.hdr.variants.cnt.simple.list <- MakeCrisprsetListSimple(fltr.res.file.df, res.1IIAe.path)
SaveVariantsBarPlot(c(
  machiato.unmodified.variants.cnt.simple.list,
  machiato.nhej.variants.cnt.simple.list,
  machiato.mixed.variants.cnt.simple.list,
  machiato.hdr.variants.cnt.simple.list
  ),
  file.path(comparison.analysis.dir, "[extra]MaChIAto_ALL(Excluding_Unclassified)_Variants_Barplot")
)



############################################################

########[g] Alluvial diagrams between total crispresso mixed-HDR or InDels and crisprVariants class
message("[g]Make an alluvial diagrams between total crispresso mixed-HDR or InDels and crisprVariants class")

res.1IAb.path <- "result/[1]comparison_with_original_crispresso_data/[Ⅰ]crispresso_alignment_data/[Ab]crispresso_unmodified/crispr.set.rds"
crispresso.unmodified.including.unmodified.vec <- SummaryVariantsByUnmodified(fltr.res.file.df, res.1IAb.path)

res.1IAc.path <- "result/[1]comparison_with_original_crispresso_data/[Ⅰ]crispresso_alignment_data/[Ac]crispresso_nhej/crispr.set.rds"
crispresso.indel.including.unmodified.vec <- SummaryVariantsByUnmodified(fltr.res.file.df, res.1IAc.path)

res.1IAd.path <- "result/[1]comparison_with_original_crispresso_data/[Ⅰ]crispresso_alignment_data/[Ad]crispresso_mixed_hdr_nhej/crispr.set.rds"
crispresso.imprecise.knockin.including.unmodified.vec <- SummaryVariantsByUnmodified(fltr.res.file.df, res.1IAd.path)

res.1IAe.path <- "result/[1]comparison_with_original_crispresso_data/[Ⅰ]crispresso_alignment_data/[Ae]crispresso_hdr/crispr.set.rds"
crispresso.precise.knockin.including.unmodified.vec <- SummaryVariantsByUnmodified(fltr.res.file.df, res.1IAe.path)

editied.crispresso.crisprvariants.on.unmodified.melt.df <- melt(rbind(
  InDel = crispresso.indel.including.unmodified.vec ,
  Imprecise.knockin = crispresso.imprecise.knockin.including.unmodified.vec,
  Preicise.knockin = crispresso.precise.knockin.including.unmodified.vec))
SaveTable(editied.crispresso.crisprvariants.on.unmodified.melt.df, file.path(comparison.analysis.dir, "[g1]Alluvial_Diagrams_of_CRISPResso_Edited_Class"))

MakeAlluvialDiagrams(editied.crispresso.crisprvariants.on.unmodified.melt.df,
  c("#F20C36", "#A60321", "#548C1C", "#F28B66", "#BFA29B"),
  file.path(comparison.analysis.dir, "[g1]Alluvial_Diagrams_of_CRISPResso_Edited_Class.png"),
  "CRISPResso Classification", "BWAMEM-CrispRvariants Classification", "unmodified allele")

uneditied.crispresso.crisprvariants.on.unmodified.melt.df <- melt(rbind(
  Unmodified = crispresso.unmodified.including.unmodified.vec))
SaveTable(uneditied.crispresso.crisprvariants.on.unmodified.melt.df, file.path(comparison.analysis.dir, "[g2]Alluvial_Diagrams_of_CRISPResso_Unedited_Class"))

MakeAlluvialDiagrams(uneditied.crispresso.crisprvariants.on.unmodified.melt.df,
  c("#BFBFBF", "#F28B66", "#BFA29B"),
  file.path(comparison.analysis.dir, "[g2]Alluvial_Diagrams_of_CRISPResso_Unedited_Class.png"),
  "CRISPResso Classification", "BWAMEM-CrispRvariants Classification", "unmodified allele")

########[h] Alluvial diagrams between total machiato mixed-HDR or InDels and crisprVariants class
message("[h]Make an alluvial diagrams between total machiato mixed-HDR or InDels and crisprVariants class")

machiato.unmodified.including.unmodified.vec <- SummaryVariantsByUnmodified(fltr.res.file.df, res.1IIAb.path)

machiato.indel.including.unmodified.vec <- SummaryVariantsByUnmodified(fltr.res.file.df, res.1IIAc.path)

machiato.imprecise.knockin.including.unmodified.vec <- SummaryVariantsByUnmodified(fltr.res.file.df, res.1IIAd.path)

machiato.precise.knockin.including.unmodified.vec <- SummaryVariantsByUnmodified(fltr.res.file.df, res.1IIAe.path)

editied.machiato.crisprvariants.on.unmodified.melt.df <- melt(rbind(
  InDel = machiato.indel.including.unmodified.vec ,
  Imprecise.knockin = machiato.imprecise.knockin.including.unmodified.vec,
  Preicise.knockin = machiato.precise.knockin.including.unmodified.vec))
SaveTable(editied.machiato.crisprvariants.on.unmodified.melt.df, file.path(comparison.analysis.dir, "[h1]Alluvial_Diagrams_of_MaChIAto_Edited_Class.txt"))

MakeAlluvialDiagrams(editied.machiato.crisprvariants.on.unmodified.melt.df,
  c("#F20C36", "#A60321", "#548C1C", "#F28B66", "#BFA29B"),
  file.path(comparison.analysis.dir, "[h1]Alluvial_Diagrams_of_MaChIAto_Edited_Class.png"),
  "MaChIAto Classification", "BWAMEM-CrispRvariants Classification", "unmodified allele")

uneditied.machiato.crisprvariants.on.unmodified.melt.df <- melt(rbind(
  Unmodified = machiato.unmodified.including.unmodified.vec))
SaveTable(uneditied.machiato.crisprvariants.on.unmodified.melt.df, file.path(comparison.analysis.dir, "[h2]Alluvial_Diagrams_of_MaChIAto_Unedited_Class.txt"))

MakeAlluvialDiagrams(uneditied.machiato.crisprvariants.on.unmodified.melt.df,
  c("#BFBFBF", "#F28B66", "#BFA29B"),
  file.path(comparison.analysis.dir, "[h2]Alluvial_Diagrams_of_MaChIAto_Unedited_Class.png"),
  "MaChIAto Classification", "BWAMEM-CrispRvariants Classification", "unmodified allele")

########[i] Alluvial diagrams between 66I of crisprVariants class in total crispresso HDR
message("[i]Make an alluvial diagrams between 66I of crisprVariants class in total crispresso HDR class")

index.wt.path <- "index/wt.fa"
index.precise.knockin.path <- "index/mmej_ki.fa"

crispresso.precise.knockin.including.precise.knockin.length.vec <- SummaryVariantsByPreciseKnockin(fltr.res.file.df, res.1IAe.path, index.wt.path, index.precise.knockin.path)

knockin.crispresso.crisprvariants.on.precise.knockin.melt.df <- melt(rbind(
  Precise.knockin = crispresso.precise.knockin.including.precise.knockin.length.vec))
SaveTable(knockin.crispresso.crisprvariants.on.precise.knockin.melt.df, file.path(comparison.analysis.dir, "[i]Alluvial_Diagrams_of_CRISPResso_Knock-in_Class.txt"))

MakeAlluvialDiagrams(knockin.crispresso.crisprvariants.on.precise.knockin.melt.df,
  c("#F20C36", "#E5F38D", "#F27999"),
  file.path(comparison.analysis.dir, "[i]Alluvial_Diagrams_of_CRISPResso_Knock-in_Class.png"),
  "CRISPResso Classification", "BWAMEM-CrispRvariants Classification", "precise knock-in allele")

########[j] Alluvial diagrams between 66I of crisprVariants class in total machiato HDR
message("[j]Make an alluvial diagrams between 66I of crisprVariants class in total machiato HDR class")

machiato.precise.knockin.including.precise.knockin.length.vec <- SummaryVariantsByPreciseKnockin(fltr.res.file.df, res.1IIAe.path, index.wt.path, index.precise.knockin.path)

knockin.machiato.crisprvariants.on.precise.knockin.melt.df <- melt(rbind(
  Precise.knockin = machiato.precise.knockin.including.precise.knockin.length.vec))
SaveTable(knockin.machiato.crisprvariants.on.precise.knockin.melt.df, file.path(comparison.analysis.dir, "[j]Alluvial_Diagrams_of_MaChIAto_Knock-in_Class.txt"))

MakeAlluvialDiagrams(knockin.machiato.crisprvariants.on.precise.knockin.melt.df,
  c("#F20C36", "#E5F38D", "#F27999"),
  file.path(comparison.analysis.dir, "[j]Alluvial_Diagrams_of_MaChIAto_Knock-in_Class.png"),
  "MaChIAto Classification", "BWAMEM-CrispRvariants Classification", "precise knock-in allele")

########[k] Grouped and ranked by efficiency and enhancement
message("[k]Make groups and ranks by efficiency and enhancement")

machiato.rate.table <- data.frame(
  condition1.precise.knockin = MakeClassRateVec(temp.condition12.machiato.rate.table.mat, "Precise knock-in", condition1.label),
  condition1.imprecise.knockin = MakeClassRateVec(temp.condition12.machiato.rate.table.mat, "Imprecise knock-in", condition1.label),
  condition1.indel = MakeClassRateVec(temp.condition12.machiato.rate.table.mat, "InDel", condition1.label),
  condition1.unmodified = MakeClassRateVec(temp.condition12.machiato.rate.table.mat, "Unmodified", condition1.label),
  condition2.precise.knockin = MakeClassRateVec(temp.condition12.machiato.rate.table.mat, "Precise knock-in", condition2.label),
  condition2.imprecise.knockin = MakeClassRateVec(temp.condition12.machiato.rate.table.mat, "Imprecise knock-in", condition2.label),
  condition2.indel = MakeClassRateVec(temp.condition12.machiato.rate.table.mat, "InDel", condition2.label),
  condition2.unmodified = MakeClassRateVec(temp.condition12.machiato.rate.table.mat, "Unmodified", condition2.label),
  enhancement.1.2.precise.knockin = MakeClassRateVec(temp.condition12.machiato.rate.table.mat, "Precise knock-in", condition2.label) /
    MakeClassRateVec(temp.condition12.machiato.rate.table.mat, "Precise knock-in", condition1.label),
  reduction.1.2.imprecise.knockin = MakeClassRateVec(temp.condition12.machiato.rate.table.mat, "Imprecise knock-in", condition1.label) /
    MakeClassRateVec(temp.condition12.machiato.rate.table.mat, "Imprecise knock-in", condition2.label),
  reduction.1.2.indel = MakeClassRateVec(temp.condition12.machiato.rate.table.mat, "InDel", condition1.label) /
    MakeClassRateVec(temp.condition12.machiato.rate.table.mat, "InDel", condition2.label),
  reduction.1.2.unmodified = MakeClassRateVec(temp.condition12.machiato.rate.table.mat, "Unmodified", condition1.label) /
    MakeClassRateVec(temp.condition12.machiato.rate.table.mat, "Unmodified", condition2.label)
)
SaveTable(machiato.rate.table, file.path(comparison.analysis.dir, "[k1]Summary_of_Rate_of_Class_on_MaChIAto.txt"))

machiato.rank.table <- apply(-machiato.rate.table, MARGIN = 2, rank)
SaveTable(machiato.rank.table, file.path(comparison.analysis.dir, "[k2]Summary_of_Rank_of_Class_on_MaChIAto.txt"))

machiato.group.table <- apply(machiato.rank.table, MARGIN = 2, function(col){
  group.n <- round(length(col) / 3)
  high.vec <- col <= group.n
  low.vec <- col >= length(col) - group.n
  group.col <- col
  group.col[high.vec] <- "High"
  group.col[low.vec] <- "Low"
  group.col[!(high.vec | low.vec)] <- "Medium"
  return(group.col)
})
SaveTable(machiato.group.table, file.path(comparison.analysis.dir, "[k3]Summary_of_Group_of_Class_on_MaChIAto"))

# condition1

MakeBumpChart(machiato.rank.table,
  machiato.group.table,
  "condition1.precise.knockin",
  c("condition1.imprecise.knockin", "condition1.precise.knockin", "condition1.indel"),
  file.path(comparison.analysis.dir, "[k4]Relation_of_Rank_of_Class_on_MaChIAto_in_condition1.png")
)

# condition2

MakeBumpChart(machiato.rank.table,
  machiato.group.table,
  "condition2.precise.knockin",
  c("condition2.imprecise.knockin", "condition2.precise.knockin", "condition2.indel"),
  file.path(comparison.analysis.dir, "[k5]Relation_of_Rank_of_Class_on_MaChIAto_in_condition2.png")
)

# enhancement

MakeBumpChart(machiato.rank.table,
  machiato.group.table,
  "enhancement.1.2.precise.knockin",
  c("reduction.1.2.imprecise.knockin", "enhancement.1.2.precise.knockin", "reduction.1.2.indel"),
  file.path(comparison.analysis.dir, "[k6]Relation_of_Rank_of_Enhancement_on_MaChIAto.png")
)

MakePercentageBoxPlot(
  machiato.rate.table[,c("condition1.precise.knockin", "condition2.precise.knockin")],
  c("#F20C36", "#F20C36"),
  c(0, round(max(machiato.rate.table[,c("condition1.precise.knockin", "condition2.precise.knockin")]) + 1, digits = 0)),
  file.path(comparison.analysis.dir, "[k7]Barplot_of_Rate_of_Precise_knock-in.png")
)

max.imprecise.indel.unmodified <- max(cbind(machiato.rate.table[,c("condition1.imprecise.knockin", "condition2.imprecise.knockin")],
  machiato.rate.table[,c("condition1.indel", "condition2.indel")],
  machiato.rate.table[,c("condition1.unmodified", "condition2.unmodified")])
)

MakePercentageBoxPlot(
  machiato.rate.table[,c("condition1.imprecise.knockin", "condition2.imprecise.knockin")],
  c("#A60321", "#A60321"),
  c(0, round(max.imprecise.indel.unmodified + 1, digits = 0)),
  file.path(comparison.analysis.dir, "[k8]Barplot_of_Rate_of_Imprecise_knock-in.png")
)

MakePercentageBoxPlot(
  machiato.rate.table[,c("condition1.indel", "condition2.indel")],
  c("#548C1C", "#548C1C"),
  c(0, round(max.imprecise.indel.unmodified + 1, digits = 0)),
  file.path(comparison.analysis.dir, "[k9]Barplot_of_Rate_of_InDel.png")
)

MakePercentageBoxPlot(
  machiato.rate.table[,c("condition1.unmodified", "condition2.unmodified")],
  c("#BFBFBF", "#BFBFBF"),
  c(0, round(max.imprecise.indel.unmodified + 1, digits = 0)),
  file.path(comparison.analysis.dir, "[k10]Barplot_of_Rate_of_Unmodified.png")
)

MakeFoldBoxPlot(
  machiato.rate.table[,c("enhancement.1.2.precise.knockin", "reduction.1.2.imprecise.knockin", "reduction.1.2.indel", "reduction.1.2.unmodified")],
  c("#FF540D", "#FF0DFF", "#FF0DFF", "#FF0DFF"),
  c(FALSE, TRUE, TRUE, TRUE),
  c(round(min(-log2(machiato.rate.table[,c("enhancement.1.2.precise.knockin", "reduction.1.2.imprecise.knockin", "reduction.1.2.indel", "reduction.1.2.unmodified")])) - 0.5, digits = 0),
    round(max(log2(machiato.rate.table[,c("enhancement.1.2.precise.knockin", "reduction.1.2.imprecise.knockin", "reduction.1.2.indel", "reduction.1.2.unmodified")])) + 0.5, digits = 0)),
  file.path(comparison.analysis.dir, "[k11]Barplot_of_Enhancement_Reduction.png")
)

# T-test on precise knock-in
MakePairedStaticalProfile(
  machiato.rate.table$condition1.precise.knockin,
  machiato.rate.table$condition2.precise.knockin,
  file.path(comparison.analysis.dir, "[k12]Paired_t-test_on_precise_knock-in.txt")
)

# T-test on imprecise knock-in
MakePairedStaticalProfile(
  machiato.rate.table$condition1.imprecise.knockin,
  machiato.rate.table$condition2.imprecise.knockin,
  file.path(comparison.analysis.dir, "[k13]Paired_t-test_on_imprecise_knock-in.txt")
)

# T-test on indel
MakePairedStaticalProfile(
  machiato.rate.table$condition1.indel,
  machiato.rate.table$condition2.indel,
  file.path(comparison.analysis.dir, "[k14]Paired_t-test_on_indel.txt")
)

# T-test on unmodified
MakePairedStaticalProfile(
  machiato.rate.table$condition1.unmodified,
  machiato.rate.table$condition2.unmodified,
  file.path(comparison.analysis.dir, "[k15]Paired_t-test_on_unmodified.txt")
)


################################################################################################################################
message("---Run Mutation analysis---")
########[a] Position analysis in indel class
message("[a]Position analysis in indel class")

############################ Insert Position ############################

res.2nhejIIaia.path <- "result/[2]machiato_local_alignment/NHEJ/[Ⅱa]analysis_of_mutation_from_endogenous_locus/[ⅰa]edited_mutation.insertion.position.table.rds"

MakePositionRateSumTablePlots(
  res.2nhejIIaia.path,
  "Insert",
  c(
    file.path(mutation.analysis.dir, "[a1-1]Barplot_of_insert_position_of_condition1_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[a1-2]Barplot_of_insert_position_of_condition1_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[a1-3]Barplot_of_insert_position_of_condition1_reads_by_enhancement_group"),
    file.path(mutation.analysis.dir, "[a1-4]Barplot_of_insert_position_of_condition2_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[a1-5]Barplot_of_insert_position_of_condition2_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[a1-6]Barplot_of_insert_position_of_condition2_reads_by_enhancement_group")
  )
)

############################ Deletion Position ############################

res.2nhejIIaib.path <- "result/[2]machiato_local_alignment/NHEJ/[Ⅱa]analysis_of_mutation_from_endogenous_locus/[ⅰb]edited_mutation.deletion.position.table.rds"

MakePositionRateSumTablePlots(
  res.2nhejIIaib.path,
  "Deletion",
  c(
    file.path(mutation.analysis.dir, "[a2-1]Barplot_of_deletion_position_of_condition1_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[a2-2]Barplot_of_deletion_position_of_condition1_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[a2-3]Barplot_of_deletion_position_of_condition1_reads_by_enhancement_group"),
    file.path(mutation.analysis.dir, "[a2-4]Barplot_of_deletion_position_of_condition2_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[a2-5]Barplot_of_deletion_position_of_condition2_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[a2-6]Barplot_of_deletion_position_of_condition2_reads_by_enhancement_group")
  )
)


############################ No Overlap Deletion Position ############################

res.2nhejIIaic.path <- "result/[2]machiato_local_alignment/NHEJ/[Ⅱa]analysis_of_mutation_from_endogenous_locus/[ⅰc]edited_mutation.deletion.no.overlap.position.table.rds"

MakePositionRateSumTablePlots(
  res.2nhejIIaic.path,
  "No Overlap Deletion",
  c(
    file.path(mutation.analysis.dir, "[a3-1]Barplot_of_no_overlap_deletion_position_of_condition1_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[a3-2]Barplot_of_no_overlap_deletion_position_of_condition1_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[a3-3]Barplot_of_no_overlap_deletion_position_of_condition1_reads_by_enhancement_group"),
    file.path(mutation.analysis.dir, "[a3-4]Barplot_of_no_overlap_deletion_position_of_condition2_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[a3-5]Barplot_of_no_overlap_deletion_position_of_condition2_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[a3-6]Barplot_of_no_overlap_deletion_position_of_condition2_reads_by_enhancement_group")
  )
)

############################ Mutation Substitution Position ############################

res.2nhejIIaie.path <- "result/[2]machiato_local_alignment/NHEJ/[Ⅱa]analysis_of_mutation_from_endogenous_locus/[ⅰe]mutation.substitution.position.table.rds"

MakePositionRateSumTablePlots(
  res.2nhejIIaie.path ,
  "Mutation Substitution",
  c(
    file.path(mutation.analysis.dir, "[a4-1]Barplot_of_mutation_substitution_position_of_condition1_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[a4-2]Barplot_of_mutation_substitution_position_of_condition1_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[a4-3]Barplot_of_mutation_substitution_position_of_condition1_reads_by_enhancement_group"),
    file.path(mutation.analysis.dir, "[a4-4]Barplot_of_mutation_substitution_position_of_condition2_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[a4-5]Barplot_of_mutation_substitution_position_of_condition2_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[a4-6]Barplot_of_mutation_substitution_position_of_condition2_reads_by_enhancement_group")
  )
)

########[b] Indel distribution analysis in indel class
message("[b]Indel distribution analysis in indel class")

res.2nhejIIaiia.path <- "result/[2]machiato_local_alignment/NHEJ/[Ⅱa]analysis_of_mutation_from_endogenous_locus/[ⅱa]Distribution_of_indel_size_on_target_site.table.rds"

MakeRateSumTablePlots(
  res.2nhejIIaiia.path ,
  "InDel size in Indel Class",
  c(
    file.path(mutation.analysis.dir, "[b-1]Barplot_of_mutation_indel_size_of_condition1_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[b-2]Barplot_of_mutation_indel_size_of_condition1_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[b-3]Barplot_of_mutation_indel_size_of_condition1_reads_by_enhancement_group"),
    file.path(mutation.analysis.dir, "[b-4]Barplot_of_mutation_indel_size_of_condition2_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[b-5]Barplot_of_mutation_indel_size_of_condition2_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[b-6]Barplot_of_mutation_indel_size_of_condition2_reads_by_enhancement_group")
  ),
  TRUE
)

########[c] Substitution Number analysis in indel class
message("[c]Substitution Number analysis in indel class")

res.2nhejIIaiib.path <- "result/[2]machiato_local_alignment/NHEJ/[Ⅱa]analysis_of_mutation_from_endogenous_locus/[ⅱb]Distribution_of_substitution_number_on_target_site.table.rds"

MakeRateSumTablePlots(
  res.2nhejIIaiib.path ,
  "Substitution number in Indel Class",
  c(
    file.path(mutation.analysis.dir, "[c-1]Barplot_of_mutation_substitution_number_of_condition1_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[c-2]Barplot_of_mutation_substitution_number_of_condition1_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[c-3]Barplot_of_mutation_substitution_number_of_condition1_reads_by_enhancement_group"),
    file.path(mutation.analysis.dir, "[c-4]Barplot_of_mutation_substitution_number_of_condition2_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[c-5]Barplot_of_mutation_substitution_number_of_condition2_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[c-6]Barplot_of_mutation_substitution_number_of_condition2_reads_by_enhancement_group")
  ),
  FALSE
)

########[d] Major Indel analysis in indel class
message("[d]Major Indel analysis in indel class")

# major size class

res.2nhejIIaiiib.path <- "result/[2]machiato_local_alignment/NHEJ/[Ⅱa]analysis_of_mutation_from_endogenous_locus/[ⅲb]Rate_of_mutation_indel_size.rds"

MakeMaxSumInDelTablePlots(
  res.2nhejIIaiiib.path,
  c(
    file.path(mutation.analysis.dir, "[d-1]Piechart_of_max_indel_size_class_in_condition1_reads_among_condition1"),
    file.path(mutation.analysis.dir, "[d-2]Piechart_of_max_indel_size_class_in_condition1_reads_among_condition2"),
    file.path(mutation.analysis.dir, "[d-3]Piechart_of_max_indel_size_class_in_condition1_reads_among_enhancement"),
    file.path(mutation.analysis.dir, "[d-4]Piechart_of_max_indel_size_class_in_condition2_reads_among_condition1"),
    file.path(mutation.analysis.dir, "[d-5]Piechart_of_max_indel_size_class_in_condition2_reads_among_condition2"),
    file.path(mutation.analysis.dir, "[d-6]Piechart_of_max_indel_size_class_in_condition2_reads_among_enhancement")
  ), type = "indel.size.indel"
)

########[e] Major Indel analysis in MMEJ indel class
message("[e]Major Indel analysis in MMEJ indel class")

res.2nhejIIaivmmej.path <- "result/[2]machiato_local_alignment/NHEJ/[Ⅱa]analysis_of_mutation_from_endogenous_locus/[ⅳ]Rate_of_MMEJ_mut_indel_size.rds"

MakeMaxSumInDelTablePlots(
  res.2nhejIIaivmmej.path,
  c(
    file.path(mutation.analysis.dir, "[e-1]Piechart_of_max_mmej_indel_size_class_in_condition1_reads_among_condition1"),
    file.path(mutation.analysis.dir, "[e-2]Piechart_of_max_mmej_indel_size_class_in_condition1_reads_among_condition2"),
    file.path(mutation.analysis.dir, "[e-3]Piechart_of_max_mmej_indel_size_class_in_condition1_reads_among_enhancement"),
    file.path(mutation.analysis.dir, "[e-4]Piechart_of_max_mmej_indel_size_class_in_condition2_reads_among_condition1"),
    file.path(mutation.analysis.dir, "[e-5]Piechart_of_max_mmej_indel_size_class_in_condition2_reads_among_condition2"),
    file.path(mutation.analysis.dir, "[e-6]Piechart_of_max_mmej_indel_size_class_in_condition2_reads_among_enhancement")
  ), type = "ej.indel.size.indel"
)

########[f] Major Indel analysis in NHEJ indel class
message("[f]Major Indel analysis in NHEJ indel class")

res.2nhejIIaivnhej.path <- "result/[2]machiato_local_alignment/NHEJ/[Ⅱa]analysis_of_mutation_from_endogenous_locus/[ⅳ]Rate_of_NHEJ_mut_indel_size.rds"

MakeMaxSumInDelTablePlots(
  res.2nhejIIaivnhej.path,
  c(
    file.path(mutation.analysis.dir, "[f-1]Piechart_of_max_nhej_indel_size_class_in_condition1_reads_among_condition1"),
    file.path(mutation.analysis.dir, "[f-2]Piechart_of_max_nhej_indel_size_class_in_condition1_reads_among_condition2"),
    file.path(mutation.analysis.dir, "[f-3]Piechart_of_max_nhej_indel_size_class_in_condition1_reads_among_enhancement"),
    file.path(mutation.analysis.dir, "[f-4]Piechart_of_max_nhej_indel_size_class_in_condition2_reads_among_condition1"),
    file.path(mutation.analysis.dir, "[f-5]Piechart_of_max_nhej_indel_size_class_in_condition2_reads_among_condition2"),
    file.path(mutation.analysis.dir, "[f-6]Piechart_of_max_nhej_indel_size_class_in_condition2_reads_among_enhancement")
  ), type = "ej.indel.size.indel"
)

########[g] Total Reads Indel analysis in indel class
message("[g]Total Reads Indel analysis in indel class")

MakeReadsSumInDelTablePlots(
  res.2nhejIIaiiib.path,
  c(
    file.path(mutation.analysis.dir, "[g-1]Piechart_of_total_reads_indel_size_class_in_condition1_reads_among_condition1"),
    file.path(mutation.analysis.dir, "[g-2]Piechart_of_total_reads_indel_size_class_in_condition1_reads_among_condition2"),
    file.path(mutation.analysis.dir, "[g-3]Piechart_of_total_reads_indel_size_class_in_condition1_reads_among_enhancement"),
    file.path(mutation.analysis.dir, "[g-4]Piechart_of_total_reads_indel_size_class_in_condition2_reads_among_condition1"),
    file.path(mutation.analysis.dir, "[g-5]Piechart_of_total_reads_indel_size_class_in_condition2_reads_among_condition2"),
    file.path(mutation.analysis.dir, "[g-6]Piechart_of_total_reads_indel_size_class_in_condition2_reads_among_enhancement")
  ), type = "indel.size.indel"
)

########[h] Total Reads Indel analysis in MMEJ indel class
message("[h]Total Reads Indel analysis in MMEJ indel class")

MakeReadsSumInDelTablePlots(
  res.2nhejIIaivmmej.path,
  c(
    file.path(mutation.analysis.dir, "[h-1]Piechart_of_total_reads_mmej_indel_size_class_in_condition1_reads_among_condition1"),
    file.path(mutation.analysis.dir, "[h-2]Piechart_of_total_reads_mmej_indel_size_class_in_condition1_reads_among_condition2"),
    file.path(mutation.analysis.dir, "[h-3]Piechart_of_total_reads_mmej_indel_size_class_in_condition1_reads_among_enhancement"),
    file.path(mutation.analysis.dir, "[h-4]Piechart_of_total_reads_mmej_indel_size_class_in_condition2_reads_among_condition1"),
    file.path(mutation.analysis.dir, "[h-5]Piechart_of_total_reads_mmej_indel_size_class_in_condition2_reads_among_condition2"),
    file.path(mutation.analysis.dir, "[h-6]Piechart_of_total_reads_mmej_indel_size_class_in_condition2_reads_among_enhancement")
  ), type = "ej.indel.size.indel"
)

########[i] Total Reads Indel analysis in NHEJ indel class
message("[i]Total Reads Indel analysis in NHEJ indel class")

MakeReadsSumInDelTablePlots(
  res.2nhejIIaivnhej.path,
  c(
    file.path(mutation.analysis.dir, "[i-1]Piechart_of_total_reads_nhej_indel_size_class_in_condition1_reads_among_condition1"),
    file.path(mutation.analysis.dir, "[i-2]Piechart_of_total_reads_nhej_indel_size_class_in_condition1_reads_among_condition2"),
    file.path(mutation.analysis.dir, "[i-3]Piechart_of_total_reads_nhej_indel_size_class_in_condition1_reads_among_enhancement"),
    file.path(mutation.analysis.dir, "[i-4]Piechart_of_total_reads_nhej_indel_size_class_in_condition2_reads_among_condition1"),
    file.path(mutation.analysis.dir, "[i-5]Piechart_of_total_reads_nhej_indel_size_class_in_condition2_reads_among_condition2"),
    file.path(mutation.analysis.dir, "[i-6]Piechart_of_total_reads_nhej_indel_size_class_in_condition2_reads_among_enhancement")
  ), type = "ej.indel.size.indel"
)

########[j] MMEJ indel size analysis in indel class
message("[j]MMEJ indel size analysis in indel class")

res.2nhejIIavia.path <- "result/[2]machiato_local_alignment/NHEJ/[Ⅱa]analysis_of_mutation_from_endogenous_locus/[ⅵa]Distribution_of_MMEJ_indel_size_on_target_site.table.rds"

MakeRateSumTablePlots(
  res.2nhejIIavia.path,
  "InDel size in MMEJ Indel Class",
  c(
    file.path(mutation.analysis.dir, "[j-1]Barplot_of_mutation_MMEJ-mediated_indel_size_of_condition1_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[j-2]Barplot_of_mutation_MMEJ-mediated_indel_size_of_condition1_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[j-3]Barplot_of_mutation_MMEJ-mediated_indel_size_of_condition1_reads_by_enhancement_group"),
    file.path(mutation.analysis.dir, "[j-4]Barplot_of_mutation_MMEJ-mediated_indel_size_of_condition2_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[j-5]Barplot_of_mutation_MMEJ-mediated_indel_size_of_condition2_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[j-6]Barplot_of_mutation_MMEJ-mediated_indel_size_of_condition2_reads_by_enhancement_group")
  ),
  TRUE
)

########[k] NHEJ indel size analysis in indel class
message("[k]NHEJ indel size analysis in indel class")

res.2nhejIIavib.path <- "result/[2]machiato_local_alignment/NHEJ/[Ⅱa]analysis_of_mutation_from_endogenous_locus/[ⅵb]Distribution_of_NHEJ_indel_size_on_target_site.table.rds"

MakeRateSumTablePlots(
  res.2nhejIIavib.path,
  "InDel size in NHEJ Indel Class",
  c(
    file.path(mutation.analysis.dir, "[k-1]Barplot_of_mutation_NHEJ-mediated_indel_size_of_condition1_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[k-2]Barplot_of_mutation_NHEJ-mediated_indel_size_of_condition1_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[k-3]Barplot_of_mutation_NHEJ-mediated_indel_size_of_condition1_reads_by_enhancement_group"),
    file.path(mutation.analysis.dir, "[k-4]Barplot_of_mutation_NHEJ-mediated_indel_size_of_condition2_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[k-5]Barplot_of_mutation_NHEJ-mediated_indel_size_of_condition2_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[k-6]Barplot_of_mutation_NHEJ-mediated_indel_size_of_condition2_reads_by_enhancement_group")
  ),
  TRUE
)

########[l] Microhomology size analysis in MMEJ indel class
message("[l]Microhomology size analysis in MMEJ indel class")

res.2nhejIIaviia.path <- "result/[2]machiato_local_alignment/NHEJ/[Ⅱa]analysis_of_mutation_from_endogenous_locus/[ⅶa]Distribution_of_MMEJ_microhomology_length_on_target_site.table.rds"


MakeRateSumTablePlots(
  res.2nhejIIaviia.path,
  "Microhomology size in MMEJ Indel Class",
  c(
    file.path(mutation.analysis.dir, "[l-1]Barplot_of_mutation_microhomology_size_of_condition1_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[l-2]Barplot_of_mutation_microhomology_size_of_condition1_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[l-3]Barplot_of_mutation_microhomology_size_of_condition1_reads_by_enhancement_group"),
    file.path(mutation.analysis.dir, "[l-4]Barplot_of_mutation_microhomology_size_of_condition2_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[l-5]Barplot_of_mutation_microhomology_size_of_condition2_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[l-6]Barplot_of_mutation_microhomology_size_of_condition2_reads_by_enhancement_group")
  ),
  FALSE
)

########[m] intervening size analysis in MMEJ indel class
message("[m]Intervening size analysis in MMEJ indel class")

res.2nhejIIaviiia.path <- "result/[2]machiato_local_alignment/NHEJ/[Ⅱa]analysis_of_mutation_from_endogenous_locus/[ⅷa]Distribution_of_MMEJ_Trimmed_Seq_length_on_target_site.table.rds"


MakeRateSumTablePlots(
  res.2nhejIIaviiia.path,
  "Intervening size in MMEJ Indel Class",
  c(
    file.path(mutation.analysis.dir, "[m-1]Barplot_of_mutation_MMEJ-mediated_intervening_size_of_condition1_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[m-2]Barplot_of_mutation_MMEJ-mediated_intervening_size_of_condition1_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[m-3]Barplot_of_mutation_MMEJ-mediated_intervening_size_of_condition1_reads_by_enhancement_group"),
    file.path(mutation.analysis.dir, "[m-4]Barplot_of_mutation_MMEJ-mediated_intervening_size_of_condition2_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[m-5]Barplot_of_mutation_MMEJ-mediated_intervening_size_of_condition2_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[m-6]Barplot_of_mutation_MMEJ-mediated_intervening_size_of_condition2_reads_by_enhancement_group")
  ),
  FALSE
)

########[n] intervening size analysis in NHEJ indel class
message("[n]Intervening size analysis in NHEJ indel class")

res.2nhejIIaviiib.path <- "result/[2]machiato_local_alignment/NHEJ/[Ⅱa]analysis_of_mutation_from_endogenous_locus/[ⅷb]Distribution_of_NHEJ_Trimmed_Seq_length_on_target_site.table.rds"

MakeRateSumTablePlots(
  res.2nhejIIaviiib.path,
  "Intervening size in NHEJ Indel Class",
  c(
    file.path(mutation.analysis.dir, "[n-1]Barplot_of_mutation_NHEJ-mediated_intervening_size_of_condition1_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[n-2]Barplot_of_mutation_NHEJ-mediated_intervening_size_of_condition1_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[n-3]Barplot_of_mutation_NHEJ-mediated_intervening_size_of_condition1_reads_by_enhancement_group"),
    file.path(mutation.analysis.dir, "[n-4]Barplot_of_mutation_NHEJ-mediated_intervening_size_of_condition2_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[n-5]Barplot_of_mutation_NHEJ-mediated_intervening_size_of_condition2_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[n-6]Barplot_of_mutation_NHEJ-mediated_intervening_size_of_condition2_reads_by_enhancement_group")
  ),
  FALSE
)

########[o] microhomology size vs intervening size analysis in MMEJ indel class
message("[o]Microhomology size vs intervening size analysis in MMEJ indel class")

res.2nhejIIavix.path <- "result/[2]machiato_local_alignment/NHEJ/[Ⅱa]analysis_of_mutation_from_endogenous_locus/[ⅸ]Distribution_of_mut_MMEJ_Microhomology_length_Trimmed_Seq_length.table.rds"

MakeScatterPlots(
  res.2nhejIIavix.path,
  "Correlation between microhomology and intervening size in MMEJ Indel Class",
  c(
    file.path(mutation.analysis.dir, "[o-1]Countplot_of_mutation_microhomology_intervening_size_of_condition1_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[o-2]Countplot_of_mutation_microhomology_intervening_size_of_condition1_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[o-3]Countplot_of_mutation_microhomology_intervening_size_of_condition1_reads_by_enhancement_group"),
    file.path(mutation.analysis.dir, "[o-4]Countplot_of_mutation_microhomology_intervening_size_of_condition2_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[o-5]Countplot_of_mutation_microhomology_intervening_size_of_condition2_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[o-6]Countplot_of_mutation_microhomology_intervening_size_of_condition2_reads_by_enhancement_group")
  )
)

########[p] Distribution of total microhomology sequence vs frequency analysis in indel class
message("[p]Distribution of total microhomology sequence vs frequency analysis in indel class")

mono.nuc.seq.vec <- colnames(oligonucleotideFrequency(DNAStringSet(""), width = 1))
di.nuc.seq.vec <- colnames(oligonucleotideFrequency(DNAStringSet(""), width = 2))

res.2nhejIIaxa.path <- "result/[2]machiato_local_alignment/NHEJ/[Ⅱa]analysis_of_mutation_from_endogenous_locus/[ⅹa]Distribution_of_mut_MMEJ_Total_Microhomology_Sequence_Frequency.table.rds"

MakeSeqRateSumTablePlots(
  res.2nhejIIaxa.path ,
  "Sequence Frequency in Indel Class",
  "sequences",
  c(
    file.path(mutation.analysis.dir, "[p-1]Barplot_of_microhomology_sequence_frequency_of_condition1_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[p-2]Barplot_of_microhomology_sequence_frequency_of_condition1_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[p-3]Barplot_of_microhomology_sequence_frequency_of_condition1_reads_by_enhancement_group"),
    file.path(mutation.analysis.dir, "[p-4]Barplot_of_microhomology_sequence_frequency_of_condition2_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[p-5]Barplot_of_microhomology_sequence_frequency_of_condition2_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[p-6]Barplot_of_microhomology_sequence_frequency_of_condition2_reads_by_enhancement_group")
  ),
  0.2, 0.3
)

########[q] Distribution of mono-microhomology sequence vs frequency analysis in indel class
message("[q]Distribution of mono-microhomology sequence vs frequency analysis in indel class")

res.2nhejIIaxb.path <- "result/[2]machiato_local_alignment/NHEJ/[Ⅱa]analysis_of_mutation_from_endogenous_locus/[ⅹb]Distribution_of_mut_MMEJ_Mono_Microhomology_Sequence_Frequency.table.rds"

MakeSeqRateSumTablePlots(
  res.2nhejIIaxb.path ,
  "Sequence Frequency in Indel Class",
  "reads",
  c(
    file.path(mutation.analysis.dir, "[q-1]Barplot_of_mono-microhomology_sequence_frequency_of_condition1_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[q-2]Barplot_of_mono-microhomology_sequence_frequency_of_condition1_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[q-3]Barplot_of_mono-microhomology_sequence_frequency_of_condition1_reads_by_enhancement_group"),
    file.path(mutation.analysis.dir, "[q-4]Barplot_of_mono-microhomology_sequence_frequency_of_condition2_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[q-5]Barplot_of_mono-microhomology_sequence_frequency_of_condition2_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[q-6]Barplot_of_mono-microhomology_sequence_frequency_of_condition2_reads_by_enhancement_group")
  ),
  0.5, 0.8
)

########[r] Distribution of di-microhomology sequence vs frequency analysis in indel class
message("[r]Distribution of di-microhomology sequence vs frequency analysis in indel class")

res.2nhejIIaxc.path <- "result/[2]machiato_local_alignment/NHEJ/[Ⅱa]analysis_of_mutation_from_endogenous_locus/[ⅹc]Distribution_of_mut_MMEJ_Di_Microhomology_Sequence_Frequency.table.rds"

MakeSeqRateSumTablePlots(
  res.2nhejIIaxc.path ,
  "Sequence Frequency in Indel Class",
  "reads",
  c(
    file.path(mutation.analysis.dir, "[r-1]Barplot_of_di-microhomology_sequence_frequency_of_condition1_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[r-2]Barplot_of_di-microhomology_sequence_frequency_of_condition1_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[r-3]Barplot_of_di-microhomology_sequence_frequency_of_condition1_reads_by_enhancement_group"),
    file.path(mutation.analysis.dir, "[r-4]Barplot_of_di-microhomology_sequence_frequency_of_condition2_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[r-5]Barplot_of_di-microhomology_sequence_frequency_of_condition2_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[r-6]Barplot_of_di-microhomology_sequence_frequency_of_condition2_reads_by_enhancement_group")
  ),
  0.2, 0.3
)

########[s] Distribution of intervening length vs microhomology sequence frequency analysis in indel class
message("[s]Distribution of intervening length vs microhomology sequence frequency analysis in indel class")

res.2nhejIIaxd.path <- "result/[2]machiato_local_alignment/NHEJ/[Ⅱa]analysis_of_mutation_from_endogenous_locus/[ⅹd]Distribution_of_mut_MMEJ_Intervening_Length_Microhomology_Sequence_Mono_Frequency.table.rds"

MakeSeqLengthRateSumTablePlots(
  res.2nhejIIaxd.path ,
  "Sequence Frequency in Indel Class",
  "sequences",
  c(
    file.path(mutation.analysis.dir, "[s-1]Barplot_of_intervening_length_microhomology_sequence_frequency_of_condition1_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[s-2]Barplot_of_intervening_length_microhomology_sequence_frequency_of_condition1_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[s-3]Barplot_of_intervening_length_microhomology_sequence_frequency_of_condition1_reads_by_enhancement_group"),
    file.path(mutation.analysis.dir, "[s-4]Barplot_of_intervening_length_microhomology_sequence_frequency_of_condition2_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[s-5]Barplot_of_intervening_length_microhomology_sequence_frequency_of_condition2_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[s-6]Barplot_of_intervening_length_microhomology_sequence_frequency_of_condition2_reads_by_enhancement_group")
  ),
  0.05, 0.1
)

########[t] Distribution of intervening length vs microhomology di-sequence frequency analysis in indel class
message("[t]Distribution of intervening length vs microhomology sequence frequency (per 2nt) analysis in indel class")

res.2nhejIIaxe.path <- "result/[2]machiato_local_alignment/NHEJ/[Ⅱa]analysis_of_mutation_from_endogenous_locus/[ⅹe]Distribution_of_mut_MMEJ_Intervening_Length_Microhomology_Sequence_Di_Frequency.table.rds"

MakeSeqLengthRateSumTablePlots(
  res.2nhejIIaxe.path ,
  "Sequence Frequency in Indel Class",
  "sequences",
  c(
    file.path(mutation.analysis.dir, "[t-1]Barplot_of_intervening_length_microhomology_sequence_frequency(per 2nt)_of_condition1_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[t-2]Barplot_of_intervening_length_microhomology_sequence_frequency(per 2nt)_of_condition1_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[t-3]Barplot_of_intervening_length_microhomology_sequence_frequency(per 2nt)_of_condition1_reads_by_enhancement_group"),
    file.path(mutation.analysis.dir, "[t-4]Barplot_of_intervening_length_microhomology_sequence_frequency(per 2nt)_of_condition2_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[t-5]Barplot_of_intervening_length_microhomology_sequence_frequency(per 2nt)_of_condition2_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[t-6]Barplot_of_intervening_length_microhomology_sequence_frequency(per 2nt)_of_condition2_reads_by_enhancement_group")
  ),
  0.025, 0.05
)

########[u] Distribution of total intervening sequence vs frequency analysis in indel class
message("[u]Distribution of total intervening sequence vs frequency analysis in indel class")

res.2nhejIIaxia.path <- "result/[2]machiato_local_alignment/NHEJ/[Ⅱa]analysis_of_mutation_from_endogenous_locus/[xia]Distribution_of_mut_MMEJ_Total_Intervening_Sequence_Frequency.table.rds"

MakeSeqRateSumTablePlots(
  res.2nhejIIaxia.path ,
  "Sequence Frequency in Indel Class",
  "sequences",
  c(
    file.path(mutation.analysis.dir, "[u-1]Barplot_of_intervening_sequence_frequency_of_condition1_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[u-2]Barplot_of_intervening_sequence_frequency_of_condition1_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[u-3]Barplot_of_intervening_sequence_frequency_of_condition1_reads_by_enhancement_group"),
    file.path(mutation.analysis.dir, "[u-4]Barplot_of_intervening_sequence_frequency_of_condition2_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[u-5]Barplot_of_intervening_sequence_frequency_of_condition2_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[u-6]Barplot_of_intervening_sequence_frequency_of_condition2_reads_by_enhancement_group")
  ),
  0.2, 0.3
)

########[v] Distribution of mono-intervening sequence vs frequency analysis in indel class
message("[v]Distribution of mono-intervening sequence vs frequency analysis in indel class")

res.2nhejIIaxib.path <- "result/[2]machiato_local_alignment/NHEJ/[Ⅱa]analysis_of_mutation_from_endogenous_locus/[xib]Distribution_of_mut_MMEJ_Mono_Intervening_Sequence_Frequency.table.rds"

MakeSeqRateSumTablePlots(
  res.2nhejIIaxib.path ,
  "Sequence Frequency in Indel Class",
  "reads",
  c(
    file.path(mutation.analysis.dir, "[v-1]Barplot_of_mono-intervening_sequence_frequency_of_condition1_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[v-2]Barplot_of_mono-intervening_sequence_frequency_of_condition1_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[v-3]Barplot_of_mono-intervening_sequence_frequency_of_condition1_reads_by_enhancement_group"),
    file.path(mutation.analysis.dir, "[v-4]Barplot_of_mono-intervening_sequence_frequency_of_condition2_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[v-5]Barplot_of_mono-intervening_sequence_frequency_of_condition2_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[v-6]Barplot_of_mono-intervening_sequence_frequency_of_condition2_reads_by_enhancement_group")
  ),
  0.5, 0.8
)

########[w] Distribution of di-intervening sequence vs frequency analysis in indel class
message("[w]Distribution of di-intervening sequence vs frequency analysis in indel class")

res.2nhejIIaxic.path <- "result/[2]machiato_local_alignment/NHEJ/[Ⅱa]analysis_of_mutation_from_endogenous_locus/[xic]Distribution_of_mut_MMEJ_Di_Intervening_Sequence_Frequency.table.rds"

MakeSeqRateSumTablePlots(
  res.2nhejIIaxic.path ,
  "Sequence Frequency in Indel Class",
  "reads",
  c(
    file.path(mutation.analysis.dir, "[w-1]Barplot_of_di-intervening_sequence_frequency_of_condition1_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[w-2]Barplot_of_di-intervening_sequence_frequency_of_condition1_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[w-3]Barplot_of_di-intervening_sequence_frequency_of_condition1_reads_by_enhancement_group"),
    file.path(mutation.analysis.dir, "[w-4]Barplot_of_di-intervening_sequence_frequency_of_condition2_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[w-5]Barplot_of_di-intervening_sequence_frequency_of_condition2_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[w-6]Barplot_of_di-intervening_sequence_frequency_of_condition2_reads_by_enhancement_group")
  ),
  0.2, 0.3
)

########[x] Distribution of intervening length vs intervening sequence frequency analysis in indel class
message("[x]Distribution of intervening length vs intervening sequence frequency analysis in indel class")

res.2nhejIIaxid.path <- "result/[2]machiato_local_alignment/NHEJ/[Ⅱa]analysis_of_mutation_from_endogenous_locus/[xid]Distribution_of_mut_MMEJ_Intervening_Length_Intervening_Sequence_Mono_Frequency.table.rds"

MakeSeqLengthRateSumTablePlots(
  res.2nhejIIaxid.path ,
  "Sequence Frequency in Indel Class",
  "sequences",
  c(
    file.path(mutation.analysis.dir, "[x-1]Barplot_of_intervening_length_intervening_sequence_frequency_of_condition1_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[x-2]Barplot_of_intervening_length_intervening_sequence_frequency_of_condition1_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[x-3]Barplot_of_intervening_length_intervening_sequence_frequency_of_condition1_reads_by_enhancement_group"),
    file.path(mutation.analysis.dir, "[x-4]Barplot_of_intervening_length_intervening_sequence_frequency_of_condition2_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[x-5]Barplot_of_intervening_length_intervening_sequence_frequency_of_condition2_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[x-6]Barplot_of_intervening_length_intervening_sequence_frequency_of_condition2_reads_by_enhancement_group")
  ),
  0.05, 0.1
)

########[y] Distribution of intervening length vs intervening di-sequence frequency analysis in indel class
message("[y]Distribution of intervening length vs intervening sequence frequency (per 2nt) analysis in indel class")

res.2nhejIIaxie.path <- "result/[2]machiato_local_alignment/NHEJ/[Ⅱa]analysis_of_mutation_from_endogenous_locus/[xie]Distribution_of_mut_MMEJ_Intervening_Length_Intervening_Sequence_Di_Frequency.table.rds"

MakeSeqLengthRateSumTablePlots(
  res.2nhejIIaxie.path ,
  "Sequence Frequency in Indel Class",
  "sequences",
  c(
    file.path(mutation.analysis.dir, "[y-1]Barplot_of_intervening_length_intervening_sequence_frequency(per 2nt)_of_condition1_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[y-2]Barplot_of_intervening_length_intervening_sequence_frequency(per 2nt)_of_condition1_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[y-3]Barplot_of_intervening_length_intervening_sequence_frequency(per 2nt)_of_condition1_reads_by_enhancement_group"),
    file.path(mutation.analysis.dir, "[y-4]Barplot_of_intervening_length_intervening_sequence_frequency(per 2nt)_of_condition2_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[y-5]Barplot_of_intervening_length_intervening_sequence_frequency(per 2nt)_of_condition2_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[y-6]Barplot_of_intervening_length_intervening_sequence_frequency(per 2nt)_of_condition2_reads_by_enhancement_group")
  ),
  0.025, 0.05
)






################################################################################################################################
message("---Run Knock-in analysis---")
########[a] Indel distribution analysis in imprecise knock-in class
message("[a]Indel distribution analysis in imprecise knock-in class")

res.2mixedIIbialeft.path <- "result/[2]machiato_local_alignment/Mixed/[Ⅱb]analysis_of_knock-in_junction/[ⅰa]Distribution_of_indel_size_on_correct_left_ki_junction.table.rds"
res.2mixedIIbiaright.path <- "result/[2]machiato_local_alignment/Mixed/[Ⅱb]analysis_of_knock-in_junction/[ⅰa]Distribution_of_indel_size_on_correct_right_ki_junction.table.rds"


MakeRateSumTablePlots(
  res.2mixedIIbialeft.path,
  "InDel size in Imprecise Knock-In Class",
  c(
    file.path(imprecise.knockin.analysis.dir, "[a-1]Barplot_of_left_junction_of_indel_size_of_condition1_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[a-2]Barplot_of_left_junction_of_indel_size_of_condition1_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[a-3]Barplot_of_left_junction_of_indel_size_of_condition1_reads_by_enhancement_group"),
    file.path(imprecise.knockin.analysis.dir, "[a-4]Barplot_of_left_junction_of_indel_size_of_condition2_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[a-5]Barplot_of_left_junction_of_indel_size_of_condition2_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[a-6]Barplot_of_left_junction_of_indel_size_of_condition2_reads_by_enhancement_group")
  ),
  FALSE
)

MakeRateSumTablePlots(
  res.2mixedIIbiaright.path,
  "InDel size in Imprecise Knock-In Class",
  c(
    file.path(imprecise.knockin.analysis.dir, "[a-1]Barplot_of_right_junction_of_indel_size_of_condition1_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[a-2]Barplot_of_right_junction_of_indel_size_of_condition1_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[a-3]Barplot_of_right_junction_of_indel_size_of_condition1_reads_by_enhancement_group"),
    file.path(imprecise.knockin.analysis.dir, "[a-4]Barplot_of_right_junction_of_indel_size_of_condition2_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[a-5]Barplot_of_right_junction_of_indel_size_of_condition2_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[a-6]Barplot_of_right_junction_of_indel_size_of_condition2_reads_by_enhancement_group")
  ),
  FALSE
)


########[b] Substitution Number analysis in indel class
message("[b]Substitution Number analysis in indel class")

res.2mixedIIbibleft.path <- "result/[2]machiato_local_alignment/Mixed/[Ⅱb]analysis_of_knock-in_junction/[ⅰb]Distribution_of_substitution_number_on_correct_left_ki_junction.table.rds"
res.2mixedIIbibright.path <- "result/[2]machiato_local_alignment/Mixed/[Ⅱb]analysis_of_knock-in_junction/[ⅰb]Distribution_of_substitution_number_on_correct_right_ki_junction.table.rds"

MakeRateSumTablePlots(
  res.2mixedIIbibleft.path ,
  "Substitution number in Imprecise Knock-In Class",
  c(
    file.path(imprecise.knockin.analysis.dir, "[b-1]Barplot_of_left_junction_substitution_number_of_condition1_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[b-2]Barplot_of_left_junction_substitution_number_of_condition1_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[b-3]Barplot_of_left_junction_substitution_number_of_condition1_reads_by_enhancement_group"),
    file.path(imprecise.knockin.analysis.dir, "[b-4]Barplot_of_left_junction_substitution_number_of_condition2_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[b-5]Barplot_of_left_junction_substitution_number_of_condition2_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[b-6]Barplot_of_left_junction_substitution_number_of_condition2_reads_by_enhancement_group")
  ),
  FALSE
)

MakeRateSumTablePlots(
  res.2mixedIIbibleft.path ,
  "Substitution number in Imprecise Knock-In Class",
  c(
    file.path(imprecise.knockin.analysis.dir, "[b-1]Barplot_of_right_junction_substitution_number_of_condition1_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[b-2]Barplot_of_right_junction_substitution_number_of_condition1_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[b-3]Barplot_of_right_junction_substitution_number_of_condition1_reads_by_enhancement_group"),
    file.path(imprecise.knockin.analysis.dir, "[b-4]Barplot_of_right_junction_substitution_number_of_condition2_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[b-5]Barplot_of_right_junction_substitution_number_of_condition2_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[b-6]Barplot_of_right_junction_substitution_number_of_condition2_reads_by_enhancement_group")
  ),
  FALSE
)


########[c] Major Indel analysis in imprecise knock-in class
message("[c] Major Indel analysis in imprecise knock-in class")

# major size class

res.2mixedIIbiibleft.path <- "result/[2]machiato_local_alignment/Mixed/[Ⅱb]analysis_of_knock-in_junction/[ⅱb]left_ki_Rate_of_knock-in_indel_size.rds"
res.2mixedIIbiibright.path <- "result/[2]machiato_local_alignment/Mixed/[Ⅱb]analysis_of_knock-in_junction/[ⅱb]right_ki_Rate_of_knock-in_indel_size.rds"


MakeMaxSumInDelTablePlots(
  res.2mixedIIbiibleft.path,
  c(
    file.path(imprecise.knockin.analysis.dir, "[c-1]Piechart_of_max_indel_size_class_on_left_junction_in_condition1_reads_among_condition1"),
    file.path(imprecise.knockin.analysis.dir, "[c-2]Piechart_of_max_indel_size_class_on_left_junction_in_condition1_reads_among_condition2"),
    file.path(imprecise.knockin.analysis.dir, "[c-3]Piechart_of_max_indel_size_class_on_left_junction_in_condition1_reads_among_enhancement"),
    file.path(imprecise.knockin.analysis.dir, "[c-4]Piechart_of_max_indel_size_class_on_left_junction_in_condition2_reads_among_condition1"),
    file.path(imprecise.knockin.analysis.dir, "[c-5]Piechart_of_max_indel_size_class_on_left_junction_in_condition2_reads_among_condition2"),
    file.path(imprecise.knockin.analysis.dir, "[c-6]Piechart_of_max_indel_size_class_on_left_junction_in_condition2_reads_among_enhancement")
  ), type = "indel.size.impreciseki"
)

MakeMaxSumInDelTablePlots(
  res.2mixedIIbiibright.path,
  c(
    file.path(imprecise.knockin.analysis.dir, "[c-1]Piechart_of_max_indel_size_class_on_right_junction_in_condition1_reads_among_condition1"),
    file.path(imprecise.knockin.analysis.dir, "[c-2]Piechart_of_max_indel_size_class_on_right_junction_in_condition1_reads_among_condition2"),
    file.path(imprecise.knockin.analysis.dir, "[c-3]Piechart_of_max_indel_size_class_on_right_junction_in_condition1_reads_among_enhancement"),
    file.path(imprecise.knockin.analysis.dir, "[c-4]Piechart_of_max_indel_size_class_on_right_junction_in_condition2_reads_among_condition1"),
    file.path(imprecise.knockin.analysis.dir, "[c-5]Piechart_of_max_indel_size_class_on_right_junction_in_condition2_reads_among_condition2"),
    file.path(imprecise.knockin.analysis.dir, "[c-6]Piechart_of_max_indel_size_class_on_right_junction_in_condition2_reads_among_enhancement")
  ), type = "indel.size.impreciseki"
)


########[d] Major Indel analysis in MMEJ imprecise knock-in class
message("[d] Major Indel analysis in MMEJ imprecise knock-in class")

# major size class

res.2mixedIIbiiimmejleft.path <- "result/[2]machiato_local_alignment/Mixed/[Ⅱb]analysis_of_knock-in_junction/[ⅲ]Rate_of_MMEJ_left_ki_indel_size.rds"
res.2mixedIIbiiimmejright.path <- "result/[2]machiato_local_alignment/Mixed/[Ⅱb]analysis_of_knock-in_junction/[ⅲ]Rate_of_MMEJ_right_ki_indel_size.rds"


MakeMaxSumInDelTablePlots(
  res.2mixedIIbiiimmejleft.path,
  c(
    file.path(imprecise.knockin.analysis.dir, "[d-1]Piechart_of_max_mmej_indel_size_class_on_left_junction_in_condition1_reads_among_condition1"),
    file.path(imprecise.knockin.analysis.dir, "[d-2]Piechart_of_max_mmej_indel_size_class_on_left_junction_in_condition1_reads_among_condition2"),
    file.path(imprecise.knockin.analysis.dir, "[d-3]Piechart_of_max_mmej_indel_size_class_on_left_junction_in_condition1_reads_among_enhancement"),
    file.path(imprecise.knockin.analysis.dir, "[d-4]Piechart_of_max_mmej_indel_size_class_on_left_junction_in_condition2_reads_among_condition1"),
    file.path(imprecise.knockin.analysis.dir, "[d-5]Piechart_of_max_mmej_indel_size_class_on_left_junction_in_condition2_reads_among_condition2"),
    file.path(imprecise.knockin.analysis.dir, "[d-6]Piechart_of_max_mmej_indel_size_class_on_left_junction_in_condition2_reads_among_enhancement")
  ), type = "ej.indel.size.impreciseki"
)

MakeMaxSumInDelTablePlots(
  res.2mixedIIbiiimmejright.path,
  c(
    file.path(imprecise.knockin.analysis.dir, "[d-1]Piechart_of_max_mmej_indel_size_class_on_right_junction_in_condition1_reads_among_condition1"),
    file.path(imprecise.knockin.analysis.dir, "[d-2]Piechart_of_max_mmej_indel_size_class_on_right_junction_in_condition1_reads_among_condition2"),
    file.path(imprecise.knockin.analysis.dir, "[d-3]Piechart_of_max_mmej_indel_size_class_on_right_junction_in_condition1_reads_among_enhancement"),
    file.path(imprecise.knockin.analysis.dir, "[d-4]Piechart_of_max_mmej_indel_size_class_on_right_junction_in_condition2_reads_among_condition1"),
    file.path(imprecise.knockin.analysis.dir, "[d-5]Piechart_of_max_mmej_indel_size_class_on_right_junction_in_condition2_reads_among_condition2"),
    file.path(imprecise.knockin.analysis.dir, "[d-6]Piechart_of_max_mmej_indel_size_class_on_right_junction_in_condition2_reads_among_enhancement")
  ), type = "ej.indel.size.impreciseki"
)

########[e] Major Indel analysis in NHEJ imprecise knock-in class
message("[e] Major Indel analysis in NHEJ imprecise knock-in class")

# major size class

res.2mixedIIbiiinhejleft.path <- "result/[2]machiato_local_alignment/Mixed/[Ⅱb]analysis_of_knock-in_junction/[ⅲ]Rate_of_NHEJ_left_ki_indel_size.rds"
res.2mixedIIbiiinhejright.path <- "result/[2]machiato_local_alignment/Mixed/[Ⅱb]analysis_of_knock-in_junction/[ⅲ]Rate_of_NHEJ_right_ki_indel_size.rds"


MakeMaxSumInDelTablePlots(
  res.2mixedIIbiiinhejleft.path,
  c(
    file.path(imprecise.knockin.analysis.dir, "[e-1]Piechart_of_max_nhej_indel_size_class_on_left_junction_in_condition1_reads_among_condition1"),
    file.path(imprecise.knockin.analysis.dir, "[e-2]Piechart_of_max_nhej_indel_size_class_on_left_junction_in_condition1_reads_among_condition2"),
    file.path(imprecise.knockin.analysis.dir, "[e-3]Piechart_of_max_nhej_indel_size_class_on_left_junction_in_condition1_reads_among_enhancement"),
    file.path(imprecise.knockin.analysis.dir, "[e-4]Piechart_of_max_nhej_indel_size_class_on_left_junction_in_condition2_reads_among_condition1"),
    file.path(imprecise.knockin.analysis.dir, "[e-5]Piechart_of_max_nhej_indel_size_class_on_left_junction_in_condition2_reads_among_condition2"),
    file.path(imprecise.knockin.analysis.dir, "[e-6]Piechart_of_max_nhej_indel_size_class_on_left_junction_in_condition2_reads_among_enhancement")
  ), type = "ej.indel.size.impreciseki"
)

MakeMaxSumInDelTablePlots(
  res.2mixedIIbiiinhejright.path,
  c(
    file.path(imprecise.knockin.analysis.dir, "[e-1]Piechart_of_max_nhej_indel_size_class_on_right_junction_in_condition1_reads_among_condition1"),
    file.path(imprecise.knockin.analysis.dir, "[e-2]Piechart_of_max_nhej_indel_size_class_on_right_junction_in_condition1_reads_among_condition2"),
    file.path(imprecise.knockin.analysis.dir, "[e-3]Piechart_of_max_nhej_indel_size_class_on_right_junction_in_condition1_reads_among_enhancement"),
    file.path(imprecise.knockin.analysis.dir, "[e-4]Piechart_of_max_nhej_indel_size_class_on_right_junction_in_condition2_reads_among_condition1"),
    file.path(imprecise.knockin.analysis.dir, "[e-5]Piechart_of_max_nhej_indel_size_class_on_right_junction_in_condition2_reads_among_condition2"),
    file.path(imprecise.knockin.analysis.dir, "[e-6]Piechart_of_max_nhej_indel_size_class_on_right_junction_in_condition2_reads_among_enhancement")
  ), type = "ej.indel.size.impreciseki"
)

########[f] Total Reads Indel analysis in imprecise knock-in class
message("[f]Total Reads Indel analysis in imprecise knock-in class")

MakeReadsSumInDelTablePlots(
  res.2mixedIIbiibleft.path,
  c(
    file.path(imprecise.knockin.analysis.dir, "[f-1]Piechart_of_total_reads_indel_size_class_on_left_junction_in_condition1_reads_among_condition1"),
    file.path(imprecise.knockin.analysis.dir, "[f-2]Piechart_of_total_reads_indel_size_class_on_left_junction_in_condition1_reads_among_condition2"),
    file.path(imprecise.knockin.analysis.dir, "[f-3]Piechart_of_total_reads_indel_size_class_on_left_junction_in_condition1_reads_among_enhancement"),
    file.path(imprecise.knockin.analysis.dir, "[f-4]Piechart_of_total_reads_indel_size_class_on_left_junction_in_condition2_reads_among_condition1"),
    file.path(imprecise.knockin.analysis.dir, "[f-5]Piechart_of_total_reads_indel_size_class_on_left_junction_in_condition2_reads_among_condition2"),
    file.path(imprecise.knockin.analysis.dir, "[f-6]Piechart_of_total_reads_indel_size_class_on_left_junction_in_condition2_reads_among_enhancement")
  ), type = "indel.size.impreciseki"
)

MakeReadsSumInDelTablePlots(
  res.2mixedIIbiibright.path,
  c(
    file.path(imprecise.knockin.analysis.dir, "[f-1]Piechart_of_total_reads_indel_size_class_on_right_junction_in_condition1_reads_among_condition1"),
    file.path(imprecise.knockin.analysis.dir, "[f-2]Piechart_of_total_reads_indel_size_class_on_right_junction_in_condition1_reads_among_condition2"),
    file.path(imprecise.knockin.analysis.dir, "[f-3]Piechart_of_total_reads_indel_size_class_on_right_junction_in_condition1_reads_among_enhancement"),
    file.path(imprecise.knockin.analysis.dir, "[f-4]Piechart_of_total_reads_indel_size_class_on_right_junction_in_condition2_reads_among_condition1"),
    file.path(imprecise.knockin.analysis.dir, "[f-5]Piechart_of_total_reads_indel_size_class_on_right_junction_in_condition2_reads_among_condition2"),
    file.path(imprecise.knockin.analysis.dir, "[f-6]Piechart_of_total_reads_indel_size_class_on_right_junction_in_condition2_reads_among_enhancement")
  ), type = "indel.size.impreciseki"
)

########[g] Total Reads Indel analysis in MMEJ imprecise knock-in class
message("[g]Total Reads Indel analysis in MMEJ imprecise knock-in class")

MakeReadsSumInDelTablePlots(
  res.2mixedIIbiiimmejleft.path,
  c(
    file.path(imprecise.knockin.analysis.dir, "[g-1]Piechart_of_total_reads_mmej_indel_size_class_on_left_junction_in_condition1_reads_among_condition1"),
    file.path(imprecise.knockin.analysis.dir, "[g-2]Piechart_of_total_reads_mmej_indel_size_class_on_left_junction_in_condition1_reads_among_condition2"),
    file.path(imprecise.knockin.analysis.dir, "[g-3]Piechart_of_total_reads_mmej_indel_size_class_on_left_junction_in_condition1_reads_among_enhancement"),
    file.path(imprecise.knockin.analysis.dir, "[g-4]Piechart_of_total_reads_mmej_indel_size_class_on_left_junction_in_condition2_reads_among_condition1"),
    file.path(imprecise.knockin.analysis.dir, "[g-5]Piechart_of_total_reads_mmej_indel_size_class_on_left_junction_in_condition2_reads_among_condition2"),
    file.path(imprecise.knockin.analysis.dir, "[g-6]Piechart_of_total_reads_mmej_indel_size_class_on_left_junction_in_condition2_reads_among_enhancement")
  ), type = "ej.indel.size.impreciseki"
)

MakeReadsSumInDelTablePlots(
  res.2mixedIIbiiimmejright.path,
  c(
    file.path(imprecise.knockin.analysis.dir, "[g-1]Piechart_of_total_reads_mmej_indel_size_class_on_right_junction_in_condition1_reads_among_condition1"),
    file.path(imprecise.knockin.analysis.dir, "[g-2]Piechart_of_total_reads_mmej_indel_size_class_on_right_junction_in_condition1_reads_among_condition2"),
    file.path(imprecise.knockin.analysis.dir, "[g-3]Piechart_of_total_reads_mmej_indel_size_class_on_right_junction_in_condition1_reads_among_enhancement"),
    file.path(imprecise.knockin.analysis.dir, "[g-4]Piechart_of_total_reads_mmej_indel_size_class_on_right_junction_in_condition2_reads_among_condition1"),
    file.path(imprecise.knockin.analysis.dir, "[g-5]Piechart_of_total_reads_mmej_indel_size_class_on_right_junction_in_condition2_reads_among_condition2"),
    file.path(imprecise.knockin.analysis.dir, "[g-6]Piechart_of_total_reads_mmej_indel_size_class_on_right_junction_in_condition2_reads_among_enhancement")
  ), type = "ej.indel.size.impreciseki"
)

########[h] Total Reads Indel analysis in NHEJ imprecise knock-in class
message("[h]Total Reads Indel analysis in NHEJ imprecise knock-in class")

MakeReadsSumInDelTablePlots(
  res.2mixedIIbiiinhejleft.path,
  c(
    file.path(imprecise.knockin.analysis.dir, "[h-1]Piechart_of_total_reads_nhej_indel_size_class_on_left_junction_in_condition1_reads_among_condition1"),
    file.path(imprecise.knockin.analysis.dir, "[h-2]Piechart_of_total_reads_nhej_indel_size_class_on_left_junction_in_condition1_reads_among_condition2"),
    file.path(imprecise.knockin.analysis.dir, "[h-3]Piechart_of_total_reads_nhej_indel_size_class_on_left_junction_in_condition1_reads_among_enhancement"),
    file.path(imprecise.knockin.analysis.dir, "[h-4]Piechart_of_total_reads_nhej_indel_size_class_on_left_junction_in_condition2_reads_among_condition1"),
    file.path(imprecise.knockin.analysis.dir, "[h-5]Piechart_of_total_reads_nhej_indel_size_class_on_left_junction_in_condition2_reads_among_condition2"),
    file.path(imprecise.knockin.analysis.dir, "[h-6]Piechart_of_total_reads_nhej_indel_size_class_on_left_junction_in_condition2_reads_among_enhancement")
  ), type = "ej.indel.size.impreciseki"
)

MakeReadsSumInDelTablePlots(
  res.2mixedIIbiiinhejright.path,
  c(
    file.path(imprecise.knockin.analysis.dir, "[h-1]Piechart_of_total_reads_nhej_indel_size_class_on_right_junction_in_condition1_reads_among_condition1"),
    file.path(imprecise.knockin.analysis.dir, "[h-2]Piechart_of_total_reads_nhej_indel_size_class_on_right_junction_in_condition1_reads_among_condition2"),
    file.path(imprecise.knockin.analysis.dir, "[h-3]Piechart_of_total_reads_nhej_indel_size_class_on_right_junction_in_condition1_reads_among_enhancement"),
    file.path(imprecise.knockin.analysis.dir, "[h-4]Piechart_of_total_reads_nhej_indel_size_class_on_right_junction_in_condition2_reads_among_condition1"),
    file.path(imprecise.knockin.analysis.dir, "[h-5]Piechart_of_total_reads_nhej_indel_size_class_on_right_junction_in_condition2_reads_among_condition2"),
    file.path(imprecise.knockin.analysis.dir, "[h-6]Piechart_of_total_reads_nhej_indel_size_class_on_right_junction_in_condition2_reads_among_enhancement")
  ), type = "ej.indel.size.impreciseki"
)



########[i]Major Id-Ip indel rate analysis in imprecise knock-in class
message("[i]Major Id-Ip indel rate analysis in imprecise knock-in class")

res.2mixedIIbivleft.path <- "result/[2]machiato_local_alignment/Mixed/[Ⅱb]analysis_of_knock-in_junction/[ⅳ]Rate_of_Id-Ip_left_ki_endjoining_type.rds"
res.2mixedIIbivright.path <- "result/[2]machiato_local_alignment/Mixed/[Ⅱb]analysis_of_knock-in_junction/[ⅳ]Rate_of_Id-Ip_right_ki_endjoining_type.rds"

MakeMaxSumInDelTablePlots(
  res.2mixedIIbivleft.path,
  c(
    file.path(imprecise.knockin.analysis.dir, "[i-1]Piechart_of_max_Id-Ip_indel_size_class_on_left_junction_in_condition1_reads_among_condition1"),
    file.path(imprecise.knockin.analysis.dir, "[i-2]Piechart_of_max_Id-Ip_indel_size_class_on_left_junction_in_condition1_reads_among_condition2"),
    file.path(imprecise.knockin.analysis.dir, "[i-3]Piechart_of_max_Id-Ip_indel_size_class_on_left_junction_in_condition1_reads_among_enhancement"),
    file.path(imprecise.knockin.analysis.dir, "[i-4]Piechart_of_max_Id-Ip_indel_size_class_on_left_junction_in_condition2_reads_among_condition1"),
    file.path(imprecise.knockin.analysis.dir, "[i-5]Piechart_of_max_Id-Ip_indel_size_class_on_left_junction_in_condition2_reads_among_condition2"),
    file.path(imprecise.knockin.analysis.dir, "[i-6]Piechart_of_max_Id-Ip_indel_size_class_on_left_junction_in_condition2_reads_among_enhancement")
  ), type = "id-ip.ej.indel.size.impreciseki"
)

MakeMaxSumInDelTablePlots(
  res.2mixedIIbivright.path,
  c(
    file.path(imprecise.knockin.analysis.dir, "[i-1]Piechart_of_max_Id-Ip_indel_size_class_on_right_junction_in_condition1_reads_among_condition1"),
    file.path(imprecise.knockin.analysis.dir, "[i-2]Piechart_of_max_Id-Ip_indel_size_class_on_right_junction_in_condition1_reads_among_condition2"),
    file.path(imprecise.knockin.analysis.dir, "[i-3]Piechart_of_max_Id-Ip_indel_size_class_on_right_junction_in_condition1_reads_among_enhancement"),
    file.path(imprecise.knockin.analysis.dir, "[i-4]Piechart_of_max_Id-Ip_indel_size_class_on_right_junction_in_condition2_reads_among_condition1"),
    file.path(imprecise.knockin.analysis.dir, "[i-5]Piechart_of_max_Id-Ip_indel_size_class_on_right_junction_in_condition2_reads_among_condition2"),
    file.path(imprecise.knockin.analysis.dir, "[i-6]Piechart_of_max_Id-Ip_indel_size_class_on_right_junction_in_condition2_reads_among_enhancement")
  ), type = "id-ip.ej.indel.size.impreciseki"
)


########[j]Total Reads  Id-Ip indel rate analysis in imprecise knock-in class
message("[j]Total Reads Id-Ip indel rate analysis in imprecise knock-in class")

MakeReadsSumInDelTablePlots(
  res.2mixedIIbivleft.path,
  c(
    file.path(imprecise.knockin.analysis.dir, "[j-1]Piechart_of_total_reads_Id-Ip_indel_size_class_on_left_junction_in_condition1_reads_among_condition1"),
    file.path(imprecise.knockin.analysis.dir, "[j-2]Piechart_of_total_reads_Id-Ip_indel_size_class_on_left_junction_in_condition1_reads_among_condition2"),
    file.path(imprecise.knockin.analysis.dir, "[j-3]Piechart_of_total_reads_Id-Ip_indel_size_class_on_left_junction_in_condition1_reads_among_enhancement"),
    file.path(imprecise.knockin.analysis.dir, "[j-4]Piechart_of_total_reads_Id-Ip_indel_size_class_on_left_junction_in_condition2_reads_among_condition1"),
    file.path(imprecise.knockin.analysis.dir, "[j-5]Piechart_of_total_reads_Id-Ip_indel_size_class_on_left_junction_in_condition2_reads_among_condition2"),
    file.path(imprecise.knockin.analysis.dir, "[j-6]Piechart_of_total_reads_Id-Ip_indel_size_class_on_left_junction_in_condition2_reads_among_enhancement")
  ), type = "id-ip.ej.indel.size.impreciseki"
)

MakeReadsSumInDelTablePlots(
  res.2mixedIIbivright.path,
  c(
    file.path(imprecise.knockin.analysis.dir, "[j-1]Piechart_of_total_reads_Id-Ip_indel_size_class_on_right_junction_in_condition1_reads_among_condition1"),
    file.path(imprecise.knockin.analysis.dir, "[j-2]Piechart_of_total_reads_Id-Ip_indel_size_class_on_right_junction_in_condition1_reads_among_condition2"),
    file.path(imprecise.knockin.analysis.dir, "[j-3]Piechart_of_total_reads_Id-Ip_indel_size_class_on_right_junction_in_condition1_reads_among_enhancement"),
    file.path(imprecise.knockin.analysis.dir, "[j-4]Piechart_of_total_reads_Id-Ip_indel_size_class_on_right_junction_in_condition2_reads_among_condition1"),
    file.path(imprecise.knockin.analysis.dir, "[j-5]Piechart_of_total_reads_Id-Ip_indel_size_class_on_right_junction_in_condition2_reads_among_condition2"),
    file.path(imprecise.knockin.analysis.dir, "[j-6]Piechart_of_total_reads_Id-Ip_indel_size_class_on_right_junction_in_condition2_reads_among_enhancement")
  ), type = "id-ip.ej.indel.size.impreciseki"
)

########[k] MMEJ Correct Direction Indel distribution analysis in imprecise knock-in class
message("[k]MMEJ Correct Direction Indel distribution analysis in imprecise knock-in class")

res.2mixedIIbiiibleft.path <- "result/[2]machiato_local_alignment/Mixed/[Ⅱb]analysis_of_knock-in_junction/[iiib]Distribution_of_MMEJ_deletion_size_on_correct_left_junction.table.rds"
res.2mixedIIbiiibright.path <- "result/[2]machiato_local_alignment/Mixed/[Ⅱb]analysis_of_knock-in_junction/[iiib]Distribution_of_MMEJ_deletion_size_on_correct_right_junction.table.rds"

MakeRateSumTablePlots(
  res.2mixedIIbiiibleft.path,
  "MMEJ Correct Direction InDel size in Imprecise Knock-In Class",
  c(
    file.path(imprecise.knockin.analysis.dir, "[k-1]Barplot_of_left_junction_of_mmej_correct_direction_indel_size_of_condition1_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[k-2]Barplot_of_left_junction_of_mmej_correct_direction_indel_size_of_condition1_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[k-3]Barplot_of_left_junction_of_mmej_correct_direction_indel_size_of_condition1_reads_by_enhancement_group"),
    file.path(imprecise.knockin.analysis.dir, "[k-4]Barplot_of_left_junction_of_mmej_correct_direction_indel_size_of_condition2_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[k-5]Barplot_of_left_junction_of_mmej_correct_direction_indel_size_of_condition2_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[k-6]Barplot_of_left_junction_of_mmej_correct_direction_indel_size_of_condition2_reads_by_enhancement_group")
  ),
  FALSE
)

MakeRateSumTablePlots(
  res.2mixedIIbiiibright.path,
  "MMEJ Correct Direction InDel size in Imprecise Knock-In Class",
  c(
    file.path(imprecise.knockin.analysis.dir, "[k-1]Barplot_of_right_junction_of_mmej_correct_direction_indel_size_of_condition1_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[k-2]Barplot_of_right_junction_of_mmej_correct_direction_indel_size_of_condition1_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[k-3]Barplot_of_right_junction_of_mmej_correct_direction_indel_size_of_condition1_reads_by_enhancement_group"),
    file.path(imprecise.knockin.analysis.dir, "[k-4]Barplot_of_right_junction_of_mmej_correct_direction_indel_size_of_condition2_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[k-5]Barplot_of_right_junction_of_mmej_correct_direction_indel_size_of_condition2_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[k-6]Barplot_of_right_junction_of_mmej_correct_direction_indel_size_of_condition2_reads_by_enhancement_group")
  ),
  FALSE
)


########[l] NHEJ Correct Direction Indel distribution analysis in imprecise knock-in class
message("[l]NHEJ Correct Direction Indel distribution analysis in imprecise knock-in class")

res.2mixedIIbiiicleft.path <- "result/[2]machiato_local_alignment/Mixed/[Ⅱb]analysis_of_knock-in_junction/[iiic]Distribution_of_NHEJ_deletion_size_on_correct_left_junction.table.rds"
res.2mixedIIbiiicright.path <- "result/[2]machiato_local_alignment/Mixed/[Ⅱb]analysis_of_knock-in_junction/[iiic]Distribution_of_NHEJ_deletion_size_on_correct_right_junction.table.rds"

MakeRateSumTablePlots(
  res.2mixedIIbiiicleft.path,
  "NHEJ Correct Direction InDel size in Imprecise Knock-In Class",
  c(
    file.path(imprecise.knockin.analysis.dir, "[l-1]Barplot_of_left_junction_of_nhej_correct_direction_indel_size_of_condition1_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[l-2]Barplot_of_left_junction_of_nhej_correct_direction_indel_size_of_condition1_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[l-3]Barplot_of_left_junction_of_nhej_correct_direction_indel_size_of_condition1_reads_by_enhancement_group"),
    file.path(imprecise.knockin.analysis.dir, "[l-4]Barplot_of_left_junction_of_nhej_correct_direction_indel_size_of_condition2_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[l-5]Barplot_of_left_junction_of_nhej_correct_direction_indel_size_of_condition2_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[l-6]Barplot_of_left_junction_of_nhej_correct_direction_indel_size_of_condition2_reads_by_enhancement_group")
  ),
  FALSE
)

MakeRateSumTablePlots(
  res.2mixedIIbiiicright.path,
  "NHEJ Correct Direction InDel size in Imprecise Knock-In Class",
  c(
    file.path(imprecise.knockin.analysis.dir, "[l-1]Barplot_of_right_junction_of_nhej_correct_direction_indel_size_of_condition1_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[l-2]Barplot_of_right_junction_of_nhej_correct_direction_indel_size_of_condition1_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[l-3]Barplot_of_right_junction_of_nhej_correct_direction_indel_size_of_condition1_reads_by_enhancement_group"),
    file.path(imprecise.knockin.analysis.dir, "[l-4]Barplot_of_right_junction_of_nhej_correct_direction_indel_size_of_condition2_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[l-5]Barplot_of_right_junction_of_nhej_correct_direction_indel_size_of_condition2_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[l-6]Barplot_of_right_junction_of_nhej_correct_direction_indel_size_of_condition2_reads_by_enhancement_group")
  ),
  FALSE
)

########[m] MMEJ Reverse Direction Indel distribution analysis in imprecise knock-in class
message("[m]MMEJ Reverse Direction Indel distribution analysis in imprecise knock-in class")

res.2mixedIIbiiidleft.path <- "result/[2]machiato_local_alignment/Mixed/[Ⅱb]analysis_of_knock-in_junction/[iiid]Distribution_of_MMEJ_deletion_size_on_reverse_left_junction.table.rds"
res.2mixedIIbiiidright.path <- "result/[2]machiato_local_alignment/Mixed/[Ⅱb]analysis_of_knock-in_junction/[iiid]Distribution_of_MMEJ_deletion_size_on_reverse_right_junction.table.rds"

MakeRateSumTablePlots(
  res.2mixedIIbiiidleft.path,
  "MMEJ Reverse Direction InDel size in Imprecise Knock-In Class",
  c(
    file.path(imprecise.knockin.analysis.dir, "[m-1]Barplot_of_left_junction_of_mmej_reverse_direction_indel_size_of_condition1_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[m-2]Barplot_of_left_junction_of_mmej_reverse_direction_indel_size_of_condition1_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[m-3]Barplot_of_left_junction_of_mmej_reverse_direction_indel_size_of_condition1_reads_by_enhancement_group"),
    file.path(imprecise.knockin.analysis.dir, "[m-4]Barplot_of_left_junction_of_mmej_reverse_direction_indel_size_of_condition2_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[m-5]Barplot_of_left_junction_of_mmej_reverse_direction_indel_size_of_condition2_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[m-6]Barplot_of_left_junction_of_mmej_reverse_direction_indel_size_of_condition2_reads_by_enhancement_group")
  ),
  FALSE
)

MakeRateSumTablePlots(
  res.2mixedIIbiiidright.path,
  "MMEJ Reverse Direction InDel size in Imprecise Knock-In Class",
  c(
    file.path(imprecise.knockin.analysis.dir, "[m-1]Barplot_of_right_junction_of_mmej_reverse_direction_indel_size_of_condition1_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[m-2]Barplot_of_right_junction_of_mmej_reverse_direction_indel_size_of_condition1_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[m-3]Barplot_of_right_junction_of_mmej_reverse_direction_indel_size_of_condition1_reads_by_enhancement_group"),
    file.path(imprecise.knockin.analysis.dir, "[m-4]Barplot_of_right_junction_of_mmej_reverse_direction_indel_size_of_condition2_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[m-5]Barplot_of_right_junction_of_mmej_reverse_direction_indel_size_of_condition2_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[m-6]Barplot_of_right_junction_of_mmej_reverse_direction_indel_size_of_condition2_reads_by_enhancement_group")
  ),
  FALSE
)

########[n] NHEJ Reverse Direction Indel distribution analysis in imprecise knock-in class
message("[n]NHEJ Reverse Direction Indel distribution analysis in imprecise knock-in class")

res.2mixedIIbiiieleft.path <- "result/[2]machiato_local_alignment/Mixed/[Ⅱb]analysis_of_knock-in_junction/[iiie]Distribution_of_NHEJ_deletion_size_on_reverse_left_junction.table.rds"
res.2mixedIIbiiieright.path <- "result/[2]machiato_local_alignment/Mixed/[Ⅱb]analysis_of_knock-in_junction/[iiie]Distribution_of_NHEJ_deletion_size_on_reverse_left_junction.table.rds"

MakeRateSumTablePlots(
  res.2mixedIIbiiieleft.path,
  "NHEJ Reverse Direction InDel size in Imprecise Knock-In Class",
  c(
    file.path(imprecise.knockin.analysis.dir, "[n-1]Barplot_of_left_junction_of_nhej_reverse_direction_indel_size_of_condition1_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[n-2]Barplot_of_left_junction_of_nhej_reverse_direction_indel_size_of_condition1_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[n-3]Barplot_of_left_junction_of_nhej_reverse_direction_indel_size_of_condition1_reads_by_enhancement_group"),
    file.path(imprecise.knockin.analysis.dir, "[n-4]Barplot_of_left_junction_of_nhej_reverse_direction_indel_size_of_condition2_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[n-5]Barplot_of_left_junction_of_nhej_reverse_direction_indel_size_of_condition2_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[n-6]Barplot_of_left_junction_of_nhej_reverse_direction_indel_size_of_condition2_reads_by_enhancement_group")
  ),
  FALSE
)

MakeRateSumTablePlots(
  res.2mixedIIbiiieright.path,
  "NHEJ Reverse Direction InDel size in Imprecise Knock-In Class",
  c(
    file.path(imprecise.knockin.analysis.dir, "[n-1]Barplot_of_right_junction_of_nhej_reverse_direction_indel_size_of_condition1_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[n-2]Barplot_of_right_junction_of_nhej_reverse_direction_indel_size_of_condition1_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[n-3]Barplot_of_right_junction_of_nhej_reverse_direction_indel_size_of_condition1_reads_by_enhancement_group"),
    file.path(imprecise.knockin.analysis.dir, "[n-4]Barplot_of_right_junction_of_nhej_reverse_direction_indel_size_of_condition2_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[n-5]Barplot_of_right_junction_of_nhej_reverse_direction_indel_size_of_condition2_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[n-6]Barplot_of_right_junction_of_nhej_reverse_direction_indel_size_of_condition2_reads_by_enhancement_group")
  ),
  FALSE
)

########[o] Microhomology size analysis in MMEJ Imprecise Knock-In Class
message("[o]Microhomology size analysis in MMEJ Imprecise Knock-In Class")

res.2mixedIIbvaleft.path <- "result/[2]machiato_local_alignment/Mixed/[Ⅱb]analysis_of_knock-in_junction/[ⅴa]Distribution_of_MMEJ_microhomology_length_on_left_junction.table.rds"
res.2mixedIIbvaright.path <- "result/[2]machiato_local_alignment/Mixed/[Ⅱb]analysis_of_knock-in_junction/[ⅴa]Distribution_of_MMEJ_microhomology_length_on_right_junction.table.rds"


MakeRateSumTablePlots(
  res.2mixedIIbvaleft.path,
  "Microhomology size in MMEJ Imprecise Knock-In Class",
  c(
    file.path(imprecise.knockin.analysis.dir, "[o-1]Barplot_of_left_junction_microhomology_size_of_condition1_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[o-2]Barplot_of_left_junction_microhomology_size_of_condition1_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[o-3]Barplot_of_left_junction_microhomology_size_of_condition1_reads_by_enhancement_group"),
    file.path(imprecise.knockin.analysis.dir, "[o-4]Barplot_of_left_junction_microhomology_size_of_condition2_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[o-5]Barplot_of_left_junction_microhomology_size_of_condition2_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[o-6]Barplot_of_left_junction_microhomology_size_of_condition2_reads_by_enhancement_group")
  ),
  FALSE
)

MakeRateSumTablePlots(
  res.2mixedIIbvaright.path,
  "Microhomology size in MMEJ Imprecise Knock-In Class",
  c(
    file.path(imprecise.knockin.analysis.dir, "[o-1]Barplot_of_right_junction_microhomology_size_of_condition1_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[o-2]Barplot_of_right_junction_microhomology_size_of_condition1_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[o-3]Barplot_of_right_junction_microhomology_size_of_condition1_reads_by_enhancement_group"),
    file.path(imprecise.knockin.analysis.dir, "[o-4]Barplot_of_right_junction_microhomology_size_of_condition2_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[o-5]Barplot_of_right_junction_microhomology_size_of_condition2_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[o-6]Barplot_of_right_junction_microhomology_size_of_condition2_reads_by_enhancement_group")
  ),
  FALSE
)

########[p] intervening size analysis in MMEJ Imprecise Knock-In Class
message("[p]Intervening size analysis in MMEJ Imprecise Knock-In Class")

res.2mixedIIbvialeft.path <- "result/[2]machiato_local_alignment/Mixed/[Ⅱb]analysis_of_knock-in_junction/[ⅵa]Distribution_of_MMEJ_Trimmed_Seq_length_on_left_junction.table.rds"
res.2mixedIIbviaright.path <- "result/[2]machiato_local_alignment/Mixed/[Ⅱb]analysis_of_knock-in_junction/[ⅵa]Distribution_of_MMEJ_Trimmed_Seq_length_on_right_junction.table.rds"


MakeRateSumTablePlots(
  res.2mixedIIbvialeft.path,
  "Intervening size in MMEJ Imprecise Knock-In Class",
  c(
    file.path(imprecise.knockin.analysis.dir, "[p-1]Barplot_of_left_junction_MMEJ-mediated_intervening_size_of_condition1_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[p-2]Barplot_of_left_junction_MMEJ-mediated_intervening_size_of_condition1_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[p-3]Barplot_of_left_junction_MMEJ-mediated_intervening_size_of_condition1_reads_by_enhancement_group"),
    file.path(imprecise.knockin.analysis.dir, "[p-4]Barplot_of_left_junction_MMEJ-mediated_intervening_size_of_condition2_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[p-5]Barplot_of_left_junction_MMEJ-mediated_intervening_size_of_condition2_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[p-6]Barplot_of_left_junction_MMEJ-mediated_intervening_size_of_condition2_reads_by_enhancement_group")
  ),
  FALSE
)

MakeRateSumTablePlots(
  res.2mixedIIbviaright.path,
  "Intervening size in MMEJ Imprecise Knock-In Class",
  c(
    file.path(imprecise.knockin.analysis.dir, "[p-1]Barplot_of_right_junction_MMEJ-mediated_intervening_size_of_condition1_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[p-2]Barplot_of_right_junction_MMEJ-mediated_intervening_size_of_condition1_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[p-3]Barplot_of_right_junction_MMEJ-mediated_intervening_size_of_condition1_reads_by_enhancement_group"),
    file.path(imprecise.knockin.analysis.dir, "[p-4]Barplot_of_right_junction_MMEJ-mediated_intervening_size_of_condition2_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[p-5]Barplot_of_right_junction_MMEJ-mediated_intervening_size_of_condition2_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[p-6]Barplot_of_right_junction_MMEJ-mediated_intervening_size_of_condition2_reads_by_enhancement_group")
  ),
  FALSE
)

########[q] intervening size analysis in MMEJ Imprecise Knock-In Class
#message("[q]Intervening size analysis in MMEJ Imprecise Knock-In Class")

#res.2mixedIIbvialeft.path <- "result/[2]machiato_local_alignment/Mixed/[Ⅱb]analysis_of_knock-in_junction/[ⅵa]Distribution_of_MMEJ_Trimmed_Seq_length_on_left_junction.table.rds"
#res.2mixedIIbviaright.path <- "result/[2]machiato_local_alignment/Mixed/[Ⅱb]analysis_of_knock-in_junction/[ⅵa]Distribution_of_MMEJ_Trimmed_Seq_length_on_right_junction.table.rds"


#MakeRateSumTablePlots(
#  res.2mixedIIbvialeft.path,
#  "Intervening size in MMEJ Imprecise Knock-In Class",
#  c(
#    file.path(imprecise.knockin.analysis.dir, "[q-1]Barplot_of_left_junction_MMEJ-mediated_intervening_size_of_condition1_reads_by_condition1_group"),
#    file.path(imprecise.knockin.analysis.dir, "[q-2]Barplot_of_left_junction_MMEJ-mediated_intervening_size_of_condition1_reads_by_condition2_group"),
#    file.path(imprecise.knockin.analysis.dir, "[q-3]Barplot_of_left_junction_MMEJ-mediated_intervening_size_of_condition1_reads_by_enhancement_group"),
#    file.path(imprecise.knockin.analysis.dir, "[q-4]Barplot_of_left_junction_MMEJ-mediated_intervening_size_of_condition2_reads_by_condition1_group"),
#    file.path(imprecise.knockin.analysis.dir, "[q-5]Barplot_of_left_junction_MMEJ-mediated_intervening_size_of_condition2_reads_by_condition2_group"),
#    file.path(imprecise.knockin.analysis.dir, "[q-6]Barplot_of_left_junction_MMEJ-mediated_intervening_size_of_condition2_reads_by_enhancement_group")
#  ),
#  FALSE
#)

#MakeRateSumTablePlots(
#  res.2mixedIIbviaright.path,
#  "Intervening size in MMEJ Imprecise Knock-In Class",
#  c(
#    file.path(imprecise.knockin.analysis.dir, "[q-1]Barplot_of_right_junction_MMEJ-mediated_intervening_size_of_condition1_reads_by_condition1_group"),
#    file.path(imprecise.knockin.analysis.dir, "[q-2]Barplot_of_right_junction_MMEJ-mediated_intervening_size_of_condition1_reads_by_condition2_group"),
#    file.path(imprecise.knockin.analysis.dir, "[q-3]Barplot_of_right_junction_MMEJ-mediated_intervening_size_of_condition1_reads_by_enhancement_group"),
#    file.path(imprecise.knockin.analysis.dir, "[q-4]Barplot_of_right_junction_MMEJ-mediated_intervening_size_of_condition2_reads_by_condition1_group"),
#    file.path(imprecise.knockin.analysis.dir, "[q-5]Barplot_of_right_junction_MMEJ-mediated_intervening_size_of_condition2_reads_by_condition2_group"),
#    file.path(imprecise.knockin.analysis.dir, "[q-6]Barplot_of_right_junction_MMEJ-mediated_intervening_size_of_condition2_reads_by_enhancement_group")
#  ),
#  FALSE
#)

########[r] intervening size analysis in NHEJ Imprecise Knock-In Class
message("[r]Intervening size analysis in NHEJ Imprecise Knock-In Class")

res.2mixedIIbvibleft.path <- "result/[2]machiato_local_alignment/Mixed/[Ⅱb]analysis_of_knock-in_junction/[ⅵb]Distribution_of_NHEJ_Trimmed_Seq_length_on_left_junction.table.rds"
res.2mixedIIbvibright.path <- "result/[2]machiato_local_alignment/Mixed/[Ⅱb]analysis_of_knock-in_junction/[ⅵb]Distribution_of_NHEJ_Trimmed_Seq_length_on_right_junction.table.rds"

MakeRateSumTablePlots(
  res.2mixedIIbvibleft.path,
  "Intervening size in NHEJ Imprecise Knock-In Class",
  c(
    file.path(imprecise.knockin.analysis.dir, "[r-1]Barplot_of_left_junction_NHEJ-mediated_intervening_size_of_condition1_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[r-2]Barplot_of_left_junction_NHEJ-mediated_intervening_size_of_condition1_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[r-3]Barplot_of_left_junction_NHEJ-mediated_intervening_size_of_condition1_reads_by_enhancement_group"),
    file.path(imprecise.knockin.analysis.dir, "[r-4]Barplot_of_left_junction_NHEJ-mediated_intervening_size_of_condition2_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[r-5]Barplot_of_left_junction_NHEJ-mediated_intervening_size_of_condition2_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[r-6]Barplot_of_left_junction_NHEJ-mediated_intervening_size_of_condition2_reads_by_enhancement_group")
  ),
  FALSE
)

MakeRateSumTablePlots(
  res.2mixedIIbvibright.path,
  "Intervening size in NHEJ Imprecise Knock-In Class",
  c(
    file.path(imprecise.knockin.analysis.dir, "[r-1]Barplot_of_right_junction_NHEJ-mediated_intervening_size_of_condition1_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[r-2]Barplot_of_right_junction_NHEJ-mediated_intervening_size_of_condition1_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[r-3]Barplot_of_right_junction_NHEJ-mediated_intervening_size_of_condition1_reads_by_enhancement_group"),
    file.path(imprecise.knockin.analysis.dir, "[r-4]Barplot_of_right_junction_NHEJ-mediated_intervening_size_of_condition2_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[r-5]Barplot_of_right_junction_NHEJ-mediated_intervening_size_of_condition2_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[r-6]Barplot_of_right_junction_NHEJ-mediated_intervening_size_of_condition2_reads_by_enhancement_group")
  ),
  FALSE
)

########[s] microhomology size vs intervening size analysis in MMEJ Imprecise Knock-In Class
message("[s]Microhomology size vs intervening size analysis in MMEJ Imprecise Knock-In Class")

res.2mixedIIbviileft.path <- "result/[2]machiato_local_alignment/Mixed/[Ⅱb]analysis_of_knock-in_junction/[ⅶ]Distribution_of_left_ki_MMEJ_Microhomology_length_Trimmed_Seq_length.table.rds"
res.2mixedIIbviiright.path <- "result/[2]machiato_local_alignment/Mixed/[Ⅱb]analysis_of_knock-in_junction/[ⅶ]Distribution_of_right_ki_MMEJ_Microhomology_length_Trimmed_Seq_length.table.rds"

MakeScatterPlots(
  res.2mixedIIbviileft.path,
  "Correlation between microhomology and intervening size in MMEJ Imprecise Knock-In Class",
  c(
    file.path(imprecise.knockin.analysis.dir, "[s-1]Countplot_of_left_junction_microhomology_intervening_size_of_condition1_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[s-2]Countplot_of_left_junction_microhomology_intervening_size_of_condition1_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[s-3]Countplot_of_left_junction_microhomology_intervening_size_of_condition1_reads_by_enhancement_group"),
    file.path(imprecise.knockin.analysis.dir, "[s-4]Countplot_of_left_junction_microhomology_intervening_size_of_condition2_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[s-5]Countplot_of_left_junction_microhomology_intervening_size_of_condition2_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[s-6]Countplot_of_left_junction_microhomology_intervening_size_of_condition2_reads_by_enhancement_group")
  )
)

MakeScatterPlots(
  res.2mixedIIbviiright.path,
  "Correlation between microhomology and intervening size in MMEJ Imprecise Knock-In Class",
  c(
    file.path(imprecise.knockin.analysis.dir, "[s-1]Countplot_of_right_junction_microhomology_intervening_size_of_condition1_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[s-2]Countplot_of_right_junction_microhomology_intervening_size_of_condition1_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[s-3]Countplot_of_right_junction_microhomology_intervening_size_of_condition1_reads_by_enhancement_group"),
    file.path(imprecise.knockin.analysis.dir, "[s-4]Countplot_of_right_junction_microhomology_intervening_size_of_condition2_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[s-5]Countplot_of_right_junction_microhomology_intervening_size_of_condition2_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[s-6]Countplot_of_right_junction_microhomology_intervening_size_of_condition2_reads_by_enhancement_group")
  )
)

########[t] Distribution of total microhomology sequence vs frequency analysis in Imprecise Knock-In Class
message("[t]Distribution of total microhomology sequence vs frequency analysis in Imprecise Knock-In Class")

res.2mixedIIbviiialeft.path <- "result/[2]machiato_local_alignment/Mixed/[Ⅱb]analysis_of_knock-in_junction/[ⅷa]Distribution_of_left_ki_MMEJ_Total_Microhomology_Sequence_Frequency.table.rds"
res.2mixedIIbviiiaright.path <- "result/[2]machiato_local_alignment/Mixed/[Ⅱb]analysis_of_knock-in_junction/[ⅷa]Distribution_of_right_ki_MMEJ_Total_Microhomology_Sequence_Frequency.table.rds"

MakeSeqRateSumTablePlots(
  res.2mixedIIbviiialeft.path ,
  "Sequence Frequency in Imprecise Knock-In Class",
  "sequences",
  c(
    file.path(imprecise.knockin.analysis.dir, "[t-1]Barplot_of_left_junction_microhomology_sequence_frequency_of_condition1_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[t-2]Barplot_of_left_junction_microhomology_sequence_frequency_of_condition1_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[t-3]Barplot_of_left_junction_microhomology_sequence_frequency_of_condition1_reads_by_enhancement_group"),
    file.path(imprecise.knockin.analysis.dir, "[t-4]Barplot_of_left_junction_microhomology_sequence_frequency_of_condition2_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[t-5]Barplot_of_left_junction_microhomology_sequence_frequency_of_condition2_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[t-6]Barplot_of_left_junction_microhomology_sequence_frequency_of_condition2_reads_by_enhancement_group")
  ),
  0.2, 0.3
)

MakeSeqRateSumTablePlots(
  res.2mixedIIbviiiaright.path ,
  "Sequence Frequency in Imprecise Knock-In Class",
  "sequences",
  c(
    file.path(imprecise.knockin.analysis.dir, "[t-1]Barplot_of_right_junction_microhomology_sequence_frequency_of_condition1_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[t-2]Barplot_of_right_junction_microhomology_sequence_frequency_of_condition1_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[t-3]Barplot_of_right_junction_microhomology_sequence_frequency_of_condition1_reads_by_enhancement_group"),
    file.path(imprecise.knockin.analysis.dir, "[t-4]Barplot_of_right_junction_microhomology_sequence_frequency_of_condition2_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[t-5]Barplot_of_right_junction_microhomology_sequence_frequency_of_condition2_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[t-6]Barplot_of_right_junction_microhomology_sequence_frequency_of_condition2_reads_by_enhancement_group")
  ),
  0.2, 0.3
)

########[u] Distribution of mono-microhomology sequence vs frequency analysis in Imprecise Knock-In Class
message("[u]Distribution of mono-microhomology sequence vs frequency analysis in Imprecise Knock-In Class")

res.2mixedIIbviiibleft.path <- "result/[2]machiato_local_alignment/Mixed/[Ⅱb]analysis_of_knock-in_junction/[ⅷb]Distribution_of_left_ki_MMEJ_Mono_Microhomology_Sequence_Frequency.table.rds"
res.2mixedIIbviiibright.path <- "result/[2]machiato_local_alignment/Mixed/[Ⅱb]analysis_of_knock-in_junction/[ⅷb]Distribution_of_right_ki_MMEJ_Mono_Microhomology_Sequence_Frequency.table.rds"

MakeSeqRateSumTablePlots(
  res.2mixedIIbviiibleft.path ,
  "Sequence Frequency in Imprecise Knock-In Class",
  "reads",
  c(
    file.path(imprecise.knockin.analysis.dir, "[u-1]Barplot_of_left_junction_mono-microhomology_sequence_frequency_of_condition1_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[u-2]Barplot_of_left_junction_mono-microhomology_sequence_frequency_of_condition1_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[u-3]Barplot_of_left_junction_mono-microhomology_sequence_frequency_of_condition1_reads_by_enhancement_group"),
    file.path(imprecise.knockin.analysis.dir, "[u-4]Barplot_of_left_junction_mono-microhomology_sequence_frequency_of_condition2_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[u-5]Barplot_of_left_junction_mono-microhomology_sequence_frequency_of_condition2_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[u-6]Barplot_of_left_junction_mono-microhomology_sequence_frequency_of_condition2_reads_by_enhancement_group")
  ),
  0.5, 0.8
)

MakeSeqRateSumTablePlots(
  res.2mixedIIbviiibright.path ,
  "Sequence Frequency in Imprecise Knock-In Class",
  "reads",
  c(
    file.path(imprecise.knockin.analysis.dir, "[u-1]Barplot_of_right_junction_mono-microhomology_sequence_frequency_of_condition1_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[u-2]Barplot_of_right_junction_mono-microhomology_sequence_frequency_of_condition1_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[u-3]Barplot_of_right_junction_mono-microhomology_sequence_frequency_of_condition1_reads_by_enhancement_group"),
    file.path(imprecise.knockin.analysis.dir, "[u-4]Barplot_of_right_junction_mono-microhomology_sequence_frequency_of_condition2_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[u-5]Barplot_of_right_junction_mono-microhomology_sequence_frequency_of_condition2_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[u-6]Barplot_of_right_junction_mono-microhomology_sequence_frequency_of_condition2_reads_by_enhancement_group")
  ),
  0.5, 0.8
)

########[v] Distribution of di-microhomology sequence vs frequency analysis in Imprecise Knock-In Class
message("[v]Distribution of di-microhomology sequence vs frequency analysis in Imprecise Knock-In Class")

res.2mixedIIbviiicleft.path <- "result/[2]machiato_local_alignment/Mixed/[Ⅱb]analysis_of_knock-in_junction/[ⅷc]Distribution_of_left_ki_MMEJ_Di_Microhomology_Sequence_Frequency.table.rds"
res.2mixedIIbviiicright.path <- "result/[2]machiato_local_alignment/Mixed/[Ⅱb]analysis_of_knock-in_junction/[ⅷc]Distribution_of_right_ki_MMEJ_Di_Microhomology_Sequence_Frequency.table.rds"

MakeSeqRateSumTablePlots(
  res.2mixedIIbviiicleft.path ,
  "Sequence Frequency in Imprecise Knock-In Class",
  "reads",
  c(
    file.path(imprecise.knockin.analysis.dir, "[v-1]Barplot_of_left_junction_di-microhomology_sequence_frequency_of_condition1_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[v-2]Barplot_of_left_junction_di-microhomology_sequence_frequency_of_condition1_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[v-3]Barplot_of_left_junction_di-microhomology_sequence_frequency_of_condition1_reads_by_enhancement_group"),
    file.path(imprecise.knockin.analysis.dir, "[v-4]Barplot_of_left_junction_di-microhomology_sequence_frequency_of_condition2_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[v-5]Barplot_of_left_junction_di-microhomology_sequence_frequency_of_condition2_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[v-6]Barplot_of_left_junction_di-microhomology_sequence_frequency_of_condition2_reads_by_enhancement_group")
  ),
  0.2, 0.3
)

MakeSeqRateSumTablePlots(
  res.2mixedIIbviiicright.path ,
  "Sequence Frequency in Imprecise Knock-In Class",
  "reads",
  c(
    file.path(imprecise.knockin.analysis.dir, "[v-1]Barplot_of_right_junction_di-microhomology_sequence_frequency_of_condition1_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[v-2]Barplot_of_right_junction_di-microhomology_sequence_frequency_of_condition1_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[v-3]Barplot_of_right_junction_di-microhomology_sequence_frequency_of_condition1_reads_by_enhancement_group"),
    file.path(imprecise.knockin.analysis.dir, "[v-4]Barplot_of_right_junction_di-microhomology_sequence_frequency_of_condition2_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[v-5]Barplot_of_right_junction_di-microhomology_sequence_frequency_of_condition2_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[v-6]Barplot_of_right_junction_di-microhomology_sequence_frequency_of_condition2_reads_by_enhancement_group")
  ),
  0.2, 0.3
)

########[w] Distribution of intervening length vs microhomology sequence frequency analysis in Imprecise Knock-In Class
message("[w]Distribution of intervening length vs microhomology sequence frequency analysis in Imprecise Knock-In Class")

res.2mixedIIbviiidleft.path <- "result/[2]machiato_local_alignment/Mixed/[Ⅱb]analysis_of_knock-in_junction/[ⅷd]Distribution_of_left_ki_MMEJ_Intervening_Length_Microhomology_Sequence_Mono_Frequency.table.rds"
res.2mixedIIbviiidright.path <- "result/[2]machiato_local_alignment/Mixed/[Ⅱb]analysis_of_knock-in_junction/[ⅷd]Distribution_of_right_ki_MMEJ_Intervening_Length_Microhomology_Sequence_Mono_Frequency.table.rds"

MakeSeqLengthRateSumTablePlots(
  res.2mixedIIbviiidleft.path ,
  "Sequence Frequency in Imprecise Knock-In Class",
  "sequences",
  c(
    file.path(imprecise.knockin.analysis.dir, "[w-1]Barplot_of_left_junction_intervening_length_microhomology_sequence_frequency_of_condition1_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[w-2]Barplot_of_left_junction_intervening_length_microhomology_sequence_frequency_of_condition1_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[w-3]Barplot_of_left_junction_intervening_length_microhomology_sequence_frequency_of_condition1_reads_by_enhancement_group"),
    file.path(imprecise.knockin.analysis.dir, "[w-4]Barplot_of_left_junction_intervening_length_microhomology_sequence_frequency_of_condition2_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[w-5]Barplot_of_left_junction_intervening_length_microhomology_sequence_frequency_of_condition2_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[w-6]Barplot_of_left_junction_intervening_length_microhomology_sequence_frequency_of_condition2_reads_by_enhancement_group")
  ),
  0.1, 0.15
)

MakeSeqLengthRateSumTablePlots(
  res.2mixedIIbviiidright.path ,
  "Sequence Frequency in Imprecise Knock-In Class",
  "sequences",
  c(
    file.path(imprecise.knockin.analysis.dir, "[w-1]Barplot_of_right_junction_intervening_length_microhomology_sequence_frequency_of_condition1_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[w-2]Barplot_of_right_junction_intervening_length_microhomology_sequence_frequency_of_condition1_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[w-3]Barplot_of_right_junction_intervening_length_microhomology_sequence_frequency_of_condition1_reads_by_enhancement_group"),
    file.path(imprecise.knockin.analysis.dir, "[w-4]Barplot_of_right_junction_intervening_length_microhomology_sequence_frequency_of_condition2_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[w-5]Barplot_of_right_junction_intervening_length_microhomology_sequence_frequency_of_condition2_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[w-6]Barplot_of_right_junction_intervening_length_microhomology_sequence_frequency_of_condition2_reads_by_enhancement_group")
  ),
  0.1, 0.15
)

########[x] Distribution of intervening length vs microhomology di-sequence frequency analysis in Imprecise Knock-In Class
message("[x]Distribution of intervening length vs microhomology sequence frequency(per 2nt) analysis in Imprecise Knock-In Class")

res.2mixedIIbviiieleft.path <- "result/[2]machiato_local_alignment/Mixed/[Ⅱb]analysis_of_knock-in_junction/[ⅷe]Distribution_of_left_ki_MMEJ_Intervening_Length_Microhomology_Sequence_Di_Frequency.table.rds"
res.2mixedIIbviiieright.path <- "result/[2]machiato_local_alignment/Mixed/[Ⅱb]analysis_of_knock-in_junction/[ⅷe]Distribution_of_right_ki_MMEJ_Intervening_Length_Microhomology_Sequence_Di_Frequency.table.rds"

MakeSeqLengthRateSumTablePlots(
  res.2mixedIIbviiieleft.path ,
  "Sequence Frequency in Imprecise Knock-In Class",
  "sequences",
  c(
    file.path(imprecise.knockin.analysis.dir, "[x-1]Barplot_of_left_junction_intervening_length_microhomology_sequence_frequency(per 2nt)_of_condition1_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[x-2]Barplot_of_left_junction_intervening_length_microhomology_sequence_frequency(per 2nt)_of_condition1_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[x-3]Barplot_of_left_junction_intervening_length_microhomology_sequence_frequency(per 2nt)_of_condition1_reads_by_enhancement_group"),
    file.path(imprecise.knockin.analysis.dir, "[x-4]Barplot_of_left_junction_intervening_length_microhomology_sequence_frequency(per 2nt)_of_condition2_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[x-5]Barplot_of_left_junction_intervening_length_microhomology_sequence_frequency(per 2nt)_of_condition2_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[x-6]Barplot_of_left_junction_intervening_length_microhomology_sequence_frequency(per 2nt)_of_condition2_reads_by_enhancement_group")
  ),
  0.05, 0.05
)

MakeSeqLengthRateSumTablePlots(
  res.2mixedIIbviiieright.path ,
  "Sequence Frequency in Imprecise Knock-In Class",
  "sequences",
  c(
    file.path(imprecise.knockin.analysis.dir, "[x-1]Barplot_of_right_junction_intervening_length_microhomology_sequence_frequency(per 2nt)_of_condition1_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[x-2]Barplot_of_right_junction_intervening_length_microhomology_sequence_frequency(per 2nt)_of_condition1_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[x-3]Barplot_of_right_junction_intervening_length_microhomology_sequence_frequency(per 2nt)_of_condition1_reads_by_enhancement_group"),
    file.path(imprecise.knockin.analysis.dir, "[x-4]Barplot_of_right_junction_intervening_length_microhomology_sequence_frequency(per 2nt)_of_condition2_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[x-5]Barplot_of_right_junction_intervening_length_microhomology_sequence_frequency(per 2nt)_of_condition2_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[x-6]Barplot_of_right_junction_intervening_length_microhomology_sequence_frequency(per 2nt)_of_condition2_reads_by_enhancement_group")
  ),
  0.05, 0.05
)

########[y] Distribution of total intervening sequence vs frequency analysis in Imprecise Knock-In Class
message("[y]Distribution of total intervening sequence vs frequency analysis in Imprecise Knock-In Class")

res.2mixedIIbvixaleft.path <- "result/[2]machiato_local_alignment/Mixed/[Ⅱb]analysis_of_knock-in_junction/[ixa]Distribution_of_left_ki_MMEJ_Total_Intervening_Sequence_Frequency.table.rds"
res.2mixedIIbvixaright.path <- "result/[2]machiato_local_alignment/Mixed/[Ⅱb]analysis_of_knock-in_junction/[ixa]Distribution_of_right_ki_MMEJ_Total_Intervening_Sequence_Frequency.table.rds"

MakeSeqRateSumTablePlots(
  res.2mixedIIbvixaleft.path ,
  "Sequence Frequency in Imprecise Knock-In Class",
  "sequences",
  c(
    file.path(imprecise.knockin.analysis.dir, "[y-1]Barplot_of_left_junction_intervening_sequence_frequency_of_condition1_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[y-2]Barplot_of_left_junction_intervening_sequence_frequency_of_condition1_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[y-3]Barplot_of_left_junction_intervening_sequence_frequency_of_condition1_reads_by_enhancement_group"),
    file.path(imprecise.knockin.analysis.dir, "[y-4]Barplot_of_left_junction_intervening_sequence_frequency_of_condition2_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[y-5]Barplot_of_left_junction_intervening_sequence_frequency_of_condition2_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[y-6]Barplot_of_left_junction_intervening_sequence_frequency_of_condition2_reads_by_enhancement_group")
  ),
  0.2, 0.3
)

MakeSeqRateSumTablePlots(
  res.2mixedIIbvixaright.path ,
  "Sequence Frequency in Imprecise Knock-In Class",
  "sequences",
  c(
    file.path(imprecise.knockin.analysis.dir, "[y-1]Barplot_of_right_junction_intervening_sequence_frequency_of_condition1_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[y-2]Barplot_of_right_junction_intervening_sequence_frequency_of_condition1_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[y-3]Barplot_of_right_junction_intervening_sequence_frequency_of_condition1_reads_by_enhancement_group"),
    file.path(imprecise.knockin.analysis.dir, "[y-4]Barplot_of_right_junction_intervening_sequence_frequency_of_condition2_reads_by_condition1_group"),
    file.path(imprecise.knockin.analysis.dir, "[y-5]Barplot_of_right_junction_intervening_sequence_frequency_of_condition2_reads_by_condition2_group"),
    file.path(imprecise.knockin.analysis.dir, "[y-6]Barplot_of_right_junction_intervening_sequence_frequency_of_condition2_reads_by_enhancement_group")
  ),
  0.2, 0.3
)




################################################################################################################################
message("---Run optional analysis---")
########[extra1]Comparison between InDelphi InDel Variants and detected InDel Variants
message("[extra1]Comparison between InDelphi InDel Variants and detected InDel Variants")

indelphi.data.list <- apply(res.file.df, MARGIN = 1, function(row){
  temp.file.name <- paste0("inDelphi_", as.character(row["sample.target"]), ".csv")
  if((row["sample.label"] %in% c(condition1.label)) & file.exists(file.path("/Users/kazuki/Documents/Lab/MachinelearningForPITChprediction/20190624_develop_MaChIAto_Reviewer/InDelphi_output_190716", temp.file.name))){
    temp.data <- read.csv(file.path("/Users/kazuki/Documents/Lab/MachinelearningForPITChprediction/20190624_develop_MaChIAto_Reviewer/InDelphi_output_190716", temp.file.name))
    temp.data$Category <- gsub("del", "D", temp.data$Category)
    temp.data$Category <- gsub("ins", "I", temp.data$Category)
    aggregated.freq.table <- aggregate(reads ~ label, data = data.frame(label = paste0(temp.data$Length, temp.data$Category),  reads = temp.data$Predicted.frequency), sum)
    temp.df <- matrix(
      aggregated.freq.table$reads,
      ncol = 1
    )
    rownames(temp.df) <- aggregated.freq.table$label
    return(temp.df)
  }else if(!(row["sample.label"] %in% c(condition1.label))){
    return(NULL)
  }else{
    message(paste0("Ignore ", row["sample.name"], " because ", temp.file.name, " was not found."))
    return(NULL)
  }
})
indelphi.data.name.vec <- res.file.df$sample.name[!sapply(indelphi.data.list, is.null)]
indelphi.data.list <- Filter(Negate(is.null), indelphi.data.list)
names(indelphi.data.list) <- indelphi.data.name.vec

# untreated
indelphi.untreated.data.list <- indelphi.data.list
names(indelphi.untreated.data.list) <- gsub("-.*", paste0("-", untreated.label), names(indelphi.untreated.data.list))
machiato.untreated.nhej.crisprset.list <- MakeSpecificCrisprsetList(fltr.res.file.df, res.1IIAc.path, untreated.label)
# negative.ctrl
indelphi.negative.ctrl.data.list <- indelphi.data.list
names(indelphi.negative.ctrl.data.list) <- gsub("-.*", paste0("-", negative.ctrl.label), names(indelphi.negative.ctrl.data.list))
machiato.negative.ctrl.nhej.crisprset.list <- MakeSpecificCrisprsetList(fltr.res.file.df, res.1IIAc.path, negative.ctrl.label)
# condition1
machiato.condition1.nhej.crisprset.list <- MakeSpecificCrisprsetList(fltr.res.file.df, res.1IIAc.path, condition1.label)
# condition2
indelphi.condition2.data.list <- indelphi.data.list
names(indelphi.condition2.data.list) <- gsub("-.*", paste0("-", condition2.label), names(indelphi.condition2.data.list))
machiato.condition2.nhej.crisprset.list <- MakeSpecificCrisprsetList(fltr.res.file.df, res.1IIAc.path, condition2.label)

MakeHeatMap(
  indelphi.untreated.data.list,
  machiato.untreated.nhej.crisprset.list,
  "Predicted InDel Variants based on InDelphi",
  "Detected InDel Variants in untreated sample (MaChIAto Classification)",
  file.path(comparison.analysis.dir, "[extra1]Comparison_between_InDelphi_and_MaChIAto_InDel_class_in_untreated_samples"),
  untreated.label
)
MakeHeatMap(
  indelphi.negative.ctrl.data.list,
  machiato.negative.ctrl.nhej.crisprset.list,
  "Predicted InDel Variants based on InDelphi",
  "Detected InDel Variants in negative control sample (MaChIAto Classification)",
  file.path(comparison.analysis.dir, "[extra1]Comparison_between_InDelphi_and_MaChIAto_InDel_class_in_negative_control_samples"),
  negative.ctrl.label
)
MakeHeatMap(
  indelphi.data.list,
  machiato.condition1.nhej.crisprset.list,
  "Predicted InDel Variants based on InDelphi",
  "Detected InDel Variants in condition1 sample (MaChIAto Classification)",
  file.path(comparison.analysis.dir, "[extra1]Comparison_between_InDelphi_and_MaChIAto_InDel_class_in_condition1_samples"),
  condition1.label
)
MakeHeatMap(
  indelphi.condition2.data.list,
  machiato.condition2.nhej.crisprset.list,
  "Predicted InDel Variants based on InDelphi",
  "Detected InDel Variants in condition2 sample (MaChIAto Classification)",
  file.path(comparison.analysis.dir, "[extra1]Comparison_between_InDelphi_and_MaChIAto_InDel_class_in_condition2_samples"),
  condition2.label
)


########[extra2]Comparison between FORECasT InDel Variants and detected InDel Variants
message("[extra2]Comparison between FORECasT InDel Variants and detected InDel Variants")


forecast.table <- data.frame(
  label = gsub("_[-A-Za-z0-9]+$" , "", read.table("/Users/kazuki/Documents/Lab/MachinelearningForPITChprediction/20190624_develop_MaChIAto_Reviewer/FORECasT/output_predictedindelsummary.txt"
  , fill = TRUE , header = FALSE)$V1),
  reads = read.table("/Users/kazuki/Documents/Lab/MachinelearningForPITChprediction/20190624_develop_MaChIAto_Reviewer/FORECasT/output_predictedindelsummary.txt"
  , fill = TRUE , header = FALSE)$V3
)
forecast.table <- forecast.table[!forecast.table$label %in% c("-"), ]

forecast.list <- list()
temp.table <- data.frame(reads = numeric(0))
for(row.ind in 1:nrow(forecast.table)){
  if(length(grep("^@@@", forecast.table[row.ind, 1])) > 0){
    if(nrow(temp.table) > 0){
      temp.mat <- matrix(
        aggregate(reads ~ label, data = temp.table, sum)$reads,
        ncol = 1
      )
      temp.name.vec <- paste0(
        gsub("[[:alpha:]]", "", aggregate(reads ~ label, data = temp.table, sum)$label),
        gsub("[[:digit:]]", "", aggregate(reads ~ label, data = temp.table, sum)$label))
      rownames(temp.mat) <- temp.name.vec # add row name
      temp.list <- list(temp.mat)
      names(temp.list) <- temp.target.name # add list name
      forecast.list <- c(forecast.list, temp.list) # add list
      temp.target.name <- paste0(gsub("^@@@", "", forecast.table[row.ind, 1]), "-", condition1.label) # update list name
    }else if(length(forecast.list) < 1){ # initial process
      temp.target.name <- paste0(gsub("^@@@", "", forecast.table[row.ind, 1]), "-", condition1.label)
    }
    temp.table <- data.frame(label = character(0), reads = numeric(0))
    next
  }else{
    temp.table <- rbind(temp.table, forecast.table[row.ind, ])
  }
}

# untreated
forecast.untreated.list <- forecast.list
names(forecast.untreated.list) <- gsub("-.*", paste0("-", untreated.label), names(forecast.untreated.list))
# negative.ctrl
forecast.negative.ctrl.list <- forecast.list
names(forecast.negative.ctrl.list) <- gsub("-.*", paste0("-", negative.ctrl.label), names(forecast.negative.ctrl.list))
# condition2
forecast.condition2.list <- forecast.list
names(forecast.condition2.list) <- gsub("-.*", paste0("-", condition2.label), names(forecast.condition2.list))


MakeHeatMap(
  forecast.untreated.list,
  machiato.untreated.nhej.crisprset.list,
  "Predicted InDel Variants based on FORECasT",
  "Detected InDel Variants in untreated sample (MaChIAto Classification)",
  file.path(comparison.analysis.dir, "[extra2]Comparison_between_FORECasT_and_MaChIAto_InDel_class_in_untreated_samples"),
  untreated.label
)
MakeHeatMap(
  forecast.negative.ctrl.list,
  machiato.negative.ctrl.nhej.crisprset.list,
  "Predicted InDel Variants based on FORECasT",
  "Detected InDel Variants in negative control sample (MaChIAto Classification)",
  file.path(comparison.analysis.dir, "[extra2]Comparison_between_FORECasT_and_MaChIAto_InDel_class_in_negative_control_samples"),
  negative.ctrl.label
)
MakeHeatMap(
  forecast.list,
  machiato.condition1.nhej.crisprset.list,
  "Predicted InDel Variants based on FORECasT",
  "Detected InDel Variants in condition1 sample (MaChIAto Classification)",
  file.path(comparison.analysis.dir, "[extra2]Comparison_between_FORECasT_and_MaChIAto_InDel_class_in_condition1_samples"),
  condition1.label
)
MakeHeatMap(
  forecast.condition2.list,
  machiato.condition2.nhej.crisprset.list,
  "Predicted InDel Variants based on FORECasT",
  "Detected InDel Variants in condition2 sample (MaChIAto Classification)",
  file.path(comparison.analysis.dir, "[extra2]Comparison_between_FORECasT_and_MaChIAto_InDel_class_in_condition2_samples"),
  condition2.label
)



# TODO make randomie value in each target
# source(file.path(script.dir, "MaChIAtoReviewerFunctions.R"))



browser()



message("---All process is completed.---")
quit(save = "no")










