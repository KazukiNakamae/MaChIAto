debug.flag <- FALSE

message("
888b     d888         .d8888b. 888     8888888       d8888888            
8888b   d8888        d88P  Y88b888       888        d88888888            
88888b.d88888        888    888888       888       d88P888888            
888Y88888P888 8888b. 888       88888b.   888      d88P 888888888 .d88b.  
888 Y888P 888    \"88b888       888 \"88b  888     d88P  888888   d88\"\"88b 
888  Y8P  888.d888888888    888888  888  888    d88P   888888   888  888 
888   \"   888888  888Y88b  d88P888  888  888   d8888888888Y88b. Y88..88P 
888       888\"Y888888 \"Y8888P\" 888  8888888888d88P     888 \"Y888 \"Y88P\"  
                                                                         
                                                                         
                                                                         
8888888b.                 d8b                                     
888   Y88b                Y8P                                     
888    888                                                        
888   d88P .d88b. 888  888888 .d88b. 888  888  888 .d88b. 888d888 
8888888P\" d8P  Y8b888  888888d8P  Y8b888  888  888d8P  Y8b888P\"   
888 T88b  88888888Y88  88P88888888888888  888  88888888888888     
888  T88b Y8b.     Y8bd8P 888Y8b.    Y88b 888 d88PY8b.    888     
888   T88b \"Y8888   Y88P  888 \"Y8888  \"Y8888888P\"  \"Y8888 888     
                                                                  
                                                                      
version 1.0
")

help.message <- "

-------------------------------------------------------------------------------
This is not correct input. The description below would be helpful.
-------------------------------------------------------------------------------


-------------------------------------------------------------------------------
EXAMPLE INPUT
-------------------------------------------------------------------------------
Rscript MaChIAto_Reviewer/MaChIAtoReviewer.R <Summary directory> <Output prefix>

-------------------------------------------------------------------------------
DESPRIPTION
-------------------------------------------------------------------------------

Summary directory:
The summary directory generated with “collect_MaChIAto_data.py”.

Output prefix:
The directory into which the output directory is saved.

-------------------------------------------------------------------------------

"

#-- Setting ----------------------------------------------------#

message("----------------------------")
message("---Start MaChIAto_Reviewer---")
message("----------------------------")

# Get directory path
if(debug.flag){
    script.dir <- "/Volumes/databank2/temp/MaChIAto/MaChIAto_Reviewer"
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

#-- Load library ----------------------------------------------------#
message("---Load libraries---")

library(CrispRVariants)
library(Biostrings)
library(effsize)
library(ggalluvial)
library(dplyr)
library(ggplot2)
library(reshape2)
library(stringr)

if(debug.flag == TRUE){
    classifier.res.dir <- "./example/MaChIAto_Classifier_output"
    aligner.res.dir <- "./example/MaChIAto_Aligner_output"
    summary.res.dir <- "./example/collections/Double_knockin_analysis"
    output.prefix <- "./example/MaChIAto_Reviewer_output"
}else{
    if(length(commandArgs(trailingOnly=TRUE)) != 4){
      message(help.message)
      q()
    }
    # Get commandline arguments
    classifier.res.dir <- GetAbsolutePath(commandArgs(trailingOnly=TRUE)[1])
    aligner.res.dir <- GetAbsolutePath(commandArgs(trailingOnly=TRUE)[2])
    summary.res.dir <- GetAbsolutePath(commandArgs(trailingOnly=TRUE)[3])
    output.prefix <- GetAbsolutePath(commandArgs(trailingOnly=TRUE)[4])
}

################################################################################################################################
### read summary data
################################################################################################################################
all.count.table <- read.csv(file.path(summary.res.dir, "sum_ALL_read_countdata.csv"), stringsAsFactors=FALSE)
all.count.table <- cbind(all.count.table,
    sample.target = as.character(sapply(all.count.table$name, function(x){strsplit(x, "-")[[1]]})[1,]),
    sample.label = as.character(sapply(all.count.table$name, function(x){strsplit(x, "-")[[1]]})[2,]),
    stringsAsFactors=FALSE)
all.rate.table <- read.csv(file.path(summary.res.dir, "sum_ALL_read_rates.csv"), stringsAsFactors=FALSE)
label.table <- read.csv(file.path(summary.res.dir, "label_sample_type.csv"), stringsAsFactors=FALSE)

################################################################################################################################
### recognize analysis type
################################################################################################################################
### Check analysis mode
message("Checking analysis mode...")
# Count label number
n.ko.sample <- 0
n.ki.sample <- 0
for(label.vec in strsplit(label.table$sample.type, "_")){
    if(label.vec[1] == "knock-out"){
        n.ko.sample <- n.ko.sample + 1
    }else if(label.vec[1] == "knock-in"){
        n.ki.sample <- n.ki.sample + 1
    }
}
# recognize analysis mode
message("---------------------------------------------")
message("")
if((n.ko.sample == 1) & (n.ki.sample == 2)){
    message("ANALYSIS TYPE: Double knock-in analysis")
    analysis.mode = 1
    untreated.label <- label.table$label[which(label.table$sample.type == "untreated")]
    negative.ctrl.label <- label.table$label[which(label.table$sample.type == "knock-out_1")]
    condition1.label <- label.table$label[which(label.table$sample.type == "knock-in_1")]
    condition2.label <- label.table$label[which(label.table$sample.type == "knock-in_2")]
    x.label <- paste0("Unmodified: ", untreated.label, "Negative Ctrl: ", negative.ctrl.label, "Condition1: ", condition1.label, "Condition2: ", condition2.label)
}else if((n.ko.sample == 1) & (n.ki.sample == 1)){
    message("ANALYSIS TYPE: Single knock-in analysis")
    analysis.mode = 2
    untreated.label <- label.table$label[which(label.table$sample.type == "untreated")]
    negative.ctrl.label <- label.table$label[which(label.table$sample.type == "knock-out_1")]
    condition1.label <- label.table$label[which(label.table$sample.type == "knock-in_1")]
    condition2.label <- NULL
    x.label <- paste0("Unmodified: ", untreated.label, "Negative Ctrl: ", negative.ctrl.label, "Condition1: ", condition1.label)
}else if((n.ko.sample == 0) & (n.ki.sample == 1)){
    message("ANALYSIS TYPE: Simple knock-in analysis")
    analysis.mode = 3
    untreated.label <- label.table$label[which(label.table$sample.type == "untreated")]
    negative.ctrl.label <- NULL
    condition1.label <- label.table$label[which(label.table$sample.type == "knock-in_1")]
    condition2.label <- NULL
    x.label <- paste0("Unmodified: ", untreated.label, "Condition1: ", condition1.label)
}else if((n.ko.sample == 2) & (n.ki.sample == 0)){
    message("ANALYSIS TYPE: Double knock-out analysis")
    analysis.mode = 4
    untreated.label <- label.table$label[which(label.table$sample.type == "untreated")]
    negative.ctrl.label <- NULL
    condition1.label <- label.table$label[which(label.table$sample.type == "knock-out_1")]
    condition2.label <- label.table$label[which(label.table$sample.type == "knock-out_2")]
    x.label <- paste0("Unmodified: ", untreated.label, "Negative Ctrl: ", negative.ctrl.label, "Condition1: ", condition1.label)
}else if((n.ko.sample == 1) & (n.ki.sample == 0)){
    message("ANALYSIS TYPE: Single knock-out analysis")
    analysis.mode = 5
    untreated.label <- label.table$label[which(label.table$sample.type == "untreated")]
    negative.ctrl.label <- NULL
    condition1.label <- label.table$label[which(label.table$sample.type == "knock-out_1")]
    condition2.label <- NULL
    x.label <- paste0("Unmodified: ", untreated.label, "Condition1: ", condition1.label)
}
message("")
message("---------------------------------------------")
# Check whether machiato_dummy_sample label exists
is.dummy.untreated <- FALSE
if("machiato_dummy_sample" %in% label.table$label){
    is.dummy.untreated <- TRUE
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

remove.target.vec <- character(0)

### filtering sample included in the ignore_list
ignore.target <- setdiff(res.file.df$sample.target, all.rate.table$group.name)
if(length(ignore.target) > 0){
  message("Remove samples included in the ignore_list...")
  remove.target.vec <- c(remove.target.vec, ignore.target)
}

### filtering sample has "other label"
evaluated.label <- label.table$label
evaluated.count.table <- NULL
if(length(evaluated.label) > 0){
    # Remove gene has "other" label
    message("Remove sample has 'other' label...")
    # count table
    temp.evaluated.count.df <- filter(all.count.table, sample.label %in% evaluated.label)
    if(nrow(temp.evaluated.count.df) > 0){
        if(is.null(evaluated.count.table)){
            evaluated.count.table <- temp.evaluated.count.df
        }else{
            evaluated.count.table <- rbind(evaluated.count.table, temp.evaluated.count.df)
        }
        all.count.table <- evaluated.count.table
    }
    # rate table
    all.rate.table <- filter(all.rate.table, !(sample.type %in% "other"))
}

# Remove gene contains low number of reads
message("Check amount of reads...")
has.sufficint.reads.vec <- apply(all.count.table, MARGIN = 1, function(row){
  if(as.numeric(row["total"]) < 1000){
    message(paste0(row["name"], " dosen't have sufficint reads (< 1,000 reads)."))
    return(FALSE)
  }else{
    return(TRUE)
  }
})
remove.target.vec <- c(remove.target.vec, unique(as.character(all.count.table$sample.target)[!has.sufficint.reads.vec]))

# Remove gene contains low number of "classified" reads [Point!!!]
message("Check amount of reads...")
has.classified.reads.vec <- apply(all.count.table, MARGIN = 1, function(row){
  classified.read.count <- (as.numeric(row["UNTREATED"])
      + as.numeric(row["INSERT_EDITING"])
      + as.numeric(row["SUBSTITUTION_EDITING"])
      + as.numeric(row["DELETION_EDITING"])
      + as.numeric(row["IMPRECISE_KNOCK_IN"])
      + as.numeric(row["PRECISE_KNOCK_IN"]))
  if(classified.read.count < 1000){
    message(paste0(row["name"], " dosen't have sufficint classified reads (< 1,000 reads)."))
    return(FALSE)
  }else{
    return(TRUE)
  }
})
remove.target.vec <- c(remove.target.vec, unique(as.character(all.count.table$sample.target)[!has.classified.reads.vec]))

# Remove gene contains unusual Unmodified rate in untreated
message("Check conservative rate in untreated samples...")
is.conservative.sample.vec <- apply(all.count.table, MARGIN = 1, function(row){
    if(row["sample.label"] != untreated.label){
        return(TRUE) # this is no target of filtering
    }
    if(as.numeric(row["untreated"]) / as.numeric(row["total"]) < 0.8){
        message(paste0(row["name"], " is not conservative (the rate of unmodified reads < 80%)."))
        return(FALSE)
    }else{
        return(TRUE)
    }
    })
remove.target.vec <- c(remove.target.vec, unique(as.character(all.count.table$sample.target)[!is.conservative.sample.vec]))

if(analysis.mode != 3){
    # Remove gene contains no difference between wt and negative.ctrl on untreated
    message("Check whether editing efficary exists in knock-out sample...")
    knock.out.1.label <- label.table$label[which(label.table$sample.type == "knock-out_1")]
    for(focus.sample.target in unique(as.character(all.count.table$sample.target))){

        # load untreated
        temp.unmodified.df <- filter(all.count.table, sample.target == focus.sample.target & sample.label == untreated.label)

        # load negative control
        temp.negative.ctrl.df <- filter(all.count.table, sample.target == focus.sample.target & sample.label == knock.out.1.label)
        
        temp.untreated.vs.negativectrl.c2.table <- matrix(c(
            temp.unmodified.df$untreated,
            temp.unmodified.df$total - temp.unmodified.df$untreated,
            temp.negative.ctrl.df$untreated,
            temp.negative.ctrl.df$total - temp.negative.ctrl.df$untreated),
            ncol=2,
            byrow=T
        )
        options(warn=-1)
        temp.untreated.vs.negativectrl.chisq.test.res <- chisq.test(temp.untreated.vs.negativectrl.c2.table)
        options(warn=0)
        if(is.nan(temp.untreated.vs.negativectrl.chisq.test.res$p.value) |
            temp.untreated.vs.negativectrl.chisq.test.res$p.value >= 0.001){
            message(paste0(focus.sample.target, " dosen't have significant change between untreated and knock-out. (chisq.test, p >= 0.1%)"))
            remove.target.vec <- c(remove.target.vec, focus.sample.target)
        }
    }
}

if((analysis.mode == 1) | (analysis.mode == 4)){
    # Remove gene contains no expected editing read in both samples.
    message("Check whether expected editing read exists in condition1 and condition2...")

    # Read label
    if(analysis.mode == 1){
        condition.1.label <- label.table$label[which(label.table$sample.type == "knock-in_1")]
        condition.2.label <- label.table$label[which(label.table$sample.type == "knock-in_2")]
    }else if(analysis.mode == 4){
        condition.1.label <- label.table$label[which(label.table$sample.type == "knock-out_1")]
        condition.2.label <- label.table$label[which(label.table$sample.type == "knock-out_2")]
    }

    for(focus.sample.target in unique(as.character(all.count.table$sample.target))){

        # Load count table of condition 1
        temp.condition1.df <- filter(all.count.table, sample.target == focus.sample.target & sample.label == condition.1.label)
        # Load count table of condition 2
        temp.condition2.df <- filter(all.count.table, sample.target == focus.sample.target & sample.label == condition.2.label)

        if(analysis.mode == 1){
                cnt.condition.1 <- temp.condition1.df$PRECISE_KNOCK_IN
                cnt.condition.2 <- temp.condition2.df$PRECISE_KNOCK_IN
        }else if(analysis.mode == 4){
                cnt.condition.1 <- temp.condition1.df$DELETION_EDITING + temp.condition1.df$INSERT_EDITING
                cnt.condition.2 <- temp.condition2.df$DELETION_EDITING + temp.condition2.df$INSERT_EDITING
        }

        if(cnt.condition.1 < 1 | cnt.condition.2 < 1){
            message(paste0(focus.sample.target, " has no expected editing read in the both samples."))
            remove.target.vec <- c(remove.target.vec, focus.sample.target)
        }
    }
}

remove.target.vec <- unique(remove.target.vec)
message("Remove...")
print(remove.target.vec)
fltr.temp.res.file.df <- filter(res.file.df, (sample.label %in% label.table$label))
fltr.res.file.df <- filter(fltr.temp.res.file.df, !(sample.target %in% remove.target.vec))
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
  condition1.label,
  c("#BFBFBF", "#548C1C", "#A60321", "#F20C36"),
  x.label,
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
  condition1.label,
  c("#BFBFBF", "#548C1C", "#A60321", "#F20C36"),
  x.label,
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
  condition1.label,
  c("#BFBFBF", "#548C1C", "#A60321", "#F20C36"),
  x.label,
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
  condition1.label,
  c("#BFBFBF", "#548C1C", "#A60321", "#F20C36"),
  x.label,
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
fltred.consistency.mat <- tidyr::replace_na(consistency.mat, 100)
rownames(fltred.consistency.mat) <- c("Unmodified", "InDel", "Imprecise knock-in", "Precise knock-in")
# Is it not needed?
# colnames(fltred.consistency.mat) <- fltr.res.file.df$sample.name[complete.cases(t(fltred.consistency.mat))]
colnames(fltred.consistency.mat) <- fltr.res.file.df$sample.name
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
## TODO unclassifiedの出力がおかしい、データをちゃんと読み込めてるか確認
########[g] Alluvial diagrams between CRISPResso and MaChIAto
message("[g]Make an alluvial diagrams between CRISPResso and MaChIAto")
MakeAlluvialDiagramsCRISPRessoMaChIAto()

############################################################

if(analysis.mode < 4){
  ########[h] Alluvial diagrams between total crispresso mixed-HDR or InDels and crisprVariants class
  message("[h]Make an alluvial diagrams between total crispresso mixed-HDR or InDels and crisprVariants class")

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
  SaveTable(editied.crispresso.crisprvariants.on.unmodified.melt.df, file.path(comparison.analysis.dir, "[h1]Alluvial_Diagrams_of_CRISPResso_Edited_Class"))

  MakeAlluvialDiagrams(editied.crispresso.crisprvariants.on.unmodified.melt.df,
    c("#F20C36", "#A60321", "#548C1C", "#F28B66", "#BFA29B"),
    file.path(comparison.analysis.dir, "[h1]Alluvial_Diagrams_of_CRISPResso_Edited_Class.png"),
    "CRISPResso Classification", "BWAMEM-CrispRvariants Classification", "unmodified allele")

  uneditied.crispresso.crisprvariants.on.unmodified.melt.df <- melt(rbind(
    Unmodified = crispresso.unmodified.including.unmodified.vec))
  SaveTable(uneditied.crispresso.crisprvariants.on.unmodified.melt.df, file.path(comparison.analysis.dir, "[h2]Alluvial_Diagrams_of_CRISPResso_Unedited_Class"))

  MakeAlluvialDiagrams(uneditied.crispresso.crisprvariants.on.unmodified.melt.df,
    c("#BFBFBF", "#F28B66", "#BFA29B"),
    file.path(comparison.analysis.dir, "[h2]Alluvial_Diagrams_of_CRISPResso_Unedited_Class.png"),
    "CRISPResso Classification", "BWAMEM-CrispRvariants Classification", "unmodified allele")

  ########[i] Alluvial diagrams between total machiato mixed-HDR or InDels and crisprVariants class
  message("[h]Make an alluvial diagrams between total machiato mixed-HDR or InDels and crisprVariants class")

  machiato.unmodified.including.unmodified.vec <- SummaryVariantsByUnmodified(fltr.res.file.df, res.1IIAb.path)

  machiato.indel.including.unmodified.vec <- SummaryVariantsByUnmodified(fltr.res.file.df, res.1IIAc.path)

  machiato.imprecise.knockin.including.unmodified.vec <- SummaryVariantsByUnmodified(fltr.res.file.df, res.1IIAd.path)

  machiato.precise.knockin.including.unmodified.vec <- SummaryVariantsByUnmodified(fltr.res.file.df, res.1IIAe.path)

  editied.machiato.crisprvariants.on.unmodified.melt.df <- melt(rbind(
    InDel = machiato.indel.including.unmodified.vec ,
    Imprecise.knockin = machiato.imprecise.knockin.including.unmodified.vec,
    Preicise.knockin = machiato.precise.knockin.including.unmodified.vec))
  SaveTable(editied.machiato.crisprvariants.on.unmodified.melt.df, file.path(comparison.analysis.dir, "[i1]Alluvial_Diagrams_of_MaChIAto_Edited_Class.txt"))

  MakeAlluvialDiagrams(editied.machiato.crisprvariants.on.unmodified.melt.df,
    c("#F20C36", "#A60321", "#548C1C", "#F28B66", "#BFA29B"),
    file.path(comparison.analysis.dir, "[i1]Alluvial_Diagrams_of_MaChIAto_Edited_Class.png"),
    "MaChIAto Classification", "BWAMEM-CrispRvariants Classification", "unmodified allele")

  uneditied.machiato.crisprvariants.on.unmodified.melt.df <- melt(rbind(
    Unmodified = machiato.unmodified.including.unmodified.vec))
  SaveTable(uneditied.machiato.crisprvariants.on.unmodified.melt.df, file.path(comparison.analysis.dir, "[i2]Alluvial_Diagrams_of_MaChIAto_Unedited_Class.txt"))

  MakeAlluvialDiagrams(uneditied.machiato.crisprvariants.on.unmodified.melt.df,
    c("#BFBFBF", "#F28B66", "#BFA29B"),
    file.path(comparison.analysis.dir, "[i2]Alluvial_Diagrams_of_MaChIAto_Unedited_Class.png"),
    "MaChIAto Classification", "BWAMEM-CrispRvariants Classification", "unmodified allele")

  ########[j] Alluvial diagrams between CRISPResso HDR class and precise knock-in size of crisprVariants class
  message("[j]Make an alluvial diagrams between HDR CRISPResso HDR class and precise knock-in size of crisprVariants class")

  index.wt.path <- "index/wt.fa"
  index.precise.knockin.path <- "index/mmej_ki.fa"

  knock.in.fltr.res.file.df <- filter(fltr.res.file.df, (sample.label %in% c(condition1.label, condition2.label)))
  crispresso.precise.knockin.including.precise.knockin.length.vec <- SummaryVariantsByPreciseKnockin(knock.in.fltr.res.file.df, res.1IAe.path, index.wt.path, index.precise.knockin.path)

  knockin.crispresso.crisprvariants.on.precise.knockin.melt.df <- melt(rbind(
    Precise.knockin = crispresso.precise.knockin.including.precise.knockin.length.vec))
  SaveTable(knockin.crispresso.crisprvariants.on.precise.knockin.melt.df, file.path(comparison.analysis.dir, "[j]Alluvial_Diagrams_of_CRISPResso_Knock-in_Class.txt"))

  MakeAlluvialDiagrams(knockin.crispresso.crisprvariants.on.precise.knockin.melt.df,
    c("#F20C36", "#E5F38D", "#F27999"),
    file.path(comparison.analysis.dir, "[j]Alluvial_Diagrams_of_CRISPResso_Knock-in_Class.png"),
    "CRISPResso Classification", "BWAMEM-CrispRvariants Classification", "precise knock-in allele")

  ########[k] Alluvial diagrams between MaChIAto HDR class and precise knock-in size of crisprVariants class
  message("[k]Make an alluvial diagrams between MaChIAto HDR class and precise knock-in size of crisprVariants class")

  machiato.precise.knockin.including.precise.knockin.length.vec <- SummaryVariantsByPreciseKnockin(knock.in.fltr.res.file.df, res.1IIAe.path, index.wt.path, index.precise.knockin.path)

  knockin.machiato.crisprvariants.on.precise.knockin.melt.df <- melt(rbind(
    Precise.knockin = machiato.precise.knockin.including.precise.knockin.length.vec))
  SaveTable(knockin.machiato.crisprvariants.on.precise.knockin.melt.df, file.path(comparison.analysis.dir, "[k]Alluvial_Diagrams_of_MaChIAto_Knock-in_Class.txt"))

  MakeAlluvialDiagrams(knockin.machiato.crisprvariants.on.precise.knockin.melt.df,
    c("#F20C36", "#E5F38D", "#F27999"),
    file.path(comparison.analysis.dir, "[k]Alluvial_Diagrams_of_MaChIAto_Knock-in_Class.png"),
    "MaChIAto Classification", "BWAMEM-CrispRvariants Classification", "precise knock-in allele")
}


########[l] Grouped and ranked by efficiency and enhancement
message("[l]Make groups and ranks by efficiency and enhancement")
if(analysis.mode %in% c(1, 4)){
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
}else if(analysis.mode %in% c(2, 3, 5)){
  machiato.rate.table <- data.frame(
    condition1.precise.knockin = MakeClassRateVec(temp.condition12.machiato.rate.table.mat, "Precise knock-in", condition1.label),
    condition1.imprecise.knockin = MakeClassRateVec(temp.condition12.machiato.rate.table.mat, "Imprecise knock-in", condition1.label),
    condition1.indel = MakeClassRateVec(temp.condition12.machiato.rate.table.mat, "InDel", condition1.label),
    condition1.unmodified = MakeClassRateVec(temp.condition12.machiato.rate.table.mat, "Unmodified", condition1.label)
  )
}else{
  stop("Unexpected Error: an unknown analysis.mode ")
}

if(analysis.mode %in% c(1, 2, 3)){
  condition1.focus.variable <- "condition1.precise.knockin"
  condition2.focus.variable <- "condition2.precise.knockin"
  enhancement.1.2.focus.variable <- "enhancement.1.2.precise.knockin"
}else{
  condition1.focus.variable <- "condition1.indel"
  condition2.focus.variable <- "condition2.indel"
  enhancement.1.2.focus.variable <- "reduction.1.2.indel"
}

# note: The null function applies the table to modify mode setting of the table.
SaveTable(apply(machiato.rate.table, MARGIN = 2, function(col){return(col)}) , file.path(comparison.analysis.dir, "[l1]Summary_of_Rate_of_Class_on_MaChIAto"))

machiato.rank.table <- apply(-machiato.rate.table, MARGIN = 2, rank)
SaveTable(machiato.rank.table, file.path(comparison.analysis.dir, "[l2]Summary_of_Rank_of_Class_on_MaChIAto"))

machiato.group.table <- apply(machiato.rank.table, MARGIN = 2, function(col){
  group.n <- round(length(col) / 3)
  high.vec <- col <= group.n
  low.vec <- col > length(col) - group.n
  group.col <- col
  group.col[high.vec] <- "High"
  group.col[low.vec] <- "Low"
  group.col[!(high.vec | low.vec)] <- "Medium"
  return(group.col)
})
SaveTable(machiato.group.table, file.path(comparison.analysis.dir, "[l3]Summary_of_Group_of_Class_on_MaChIAto"))

# condition1
if(analysis.mode %in% c(1, 2, 3)){
  MakeBumpChart(machiato.rank.table,
    machiato.group.table,
    condition1.focus.variable,
    c("condition1.imprecise.knockin", "condition1.precise.knockin", "condition1.indel"),
    file.path(comparison.analysis.dir, "[l4]Relation_of_Rank_of_Class_on_MaChIAto_in_condition1.png")
  )
}else{
  message("[l4]Relation_of_Rank_of_Class_on_MaChIAto_in_condition1.png is skipped")
}

if(analysis.mode %in% c(1)){
  max.imprecise.indel.unmodified <- max(cbind(machiato.rate.table[,c("condition1.imprecise.knockin", "condition2.imprecise.knockin")],
    machiato.rate.table[,c("condition1.indel", "condition2.indel")],
    machiato.rate.table[,c("condition1.unmodified", "condition2.unmodified")])
  )
}else if(analysis.mode %in% c(4)){
  max.imprecise.indel.unmodified <- max(cbind(machiato.rate.table[,c("condition1.indel", "condition2.indel")],
    machiato.rate.table[,c("condition1.unmodified", "condition2.unmodified")])
  )
}

if(analysis.mode %in% c(1)){
  # condition2

  MakeBumpChart(machiato.rank.table,
    machiato.group.table,
    condition2.focus.variable,
    c("condition2.imprecise.knockin", "condition2.precise.knockin", "condition2.indel"),
    file.path(comparison.analysis.dir, "[l5]Relation_of_Rank_of_Class_on_MaChIAto_in_condition2.png")
  )

  # enhancement

  MakeBumpChart(machiato.rank.table,
    machiato.group.table,
    enhancement.1.2.focus.variable,
    c("reduction.1.2.imprecise.knockin", "enhancement.1.2.precise.knockin", "reduction.1.2.indel"),
    file.path(comparison.analysis.dir, "[l6]Relation_of_Rank_of_Enhancement_on_MaChIAto.png")
  )

  MakePercentageBoxPlot(
    machiato.rate.table[,c("condition1.precise.knockin", "condition2.precise.knockin")],
    c("#F20C36", "#F20C36"),
    c(0, round(max(machiato.rate.table[,c("condition1.precise.knockin", "condition2.precise.knockin")]) + 1, digits = 0)),
    file.path(comparison.analysis.dir, "[l7]Barplot_of_Rate_of_Precise_knock-in.png")
  )

  MakePercentageBoxPlot(
    machiato.rate.table[,c("condition1.imprecise.knockin", "condition2.imprecise.knockin")],
    c("#A60321", "#A60321"),
    c(0, round(max.imprecise.indel.unmodified + 1, digits = 0)),
    file.path(comparison.analysis.dir, "[l8]Barplot_of_Rate_of_Imprecise_knock-in.png")
  )
}else{
  message("[l5]Relation_of_Rank_of_Class_on_MaChIAto_in_condition2.png is skipped")
  message("[l6]Relation_of_Rank_of_Enhancement_on_MaChIAto.png is skipped")
  message("[l7]Barplot_of_Rate_of_Precise_knock-in.png is skipped")
  message("[l8]Barplot_of_Rate_of_Imprecise_knock-in.png is skipped")
}

if(analysis.mode %in% c(1, 4)){

  MakePercentageBoxPlot(
    machiato.rate.table[,c("condition1.indel", "condition2.indel")],
    c("#548C1C", "#548C1C"),
    c(0, round(max.imprecise.indel.unmodified + 1, digits = 0)),
    file.path(comparison.analysis.dir, "[l9]Barplot_of_Rate_of_InDel.png")
  )

  MakePercentageBoxPlot(
    machiato.rate.table[,c("condition1.unmodified", "condition2.unmodified")],
    c("#BFBFBF", "#BFBFBF"),
    c(0, round(max.imprecise.indel.unmodified + 1, digits = 0)),
    file.path(comparison.analysis.dir, "[l10]Barplot_of_Rate_of_Unmodified.png")
  )
}else{
  message("[l9]Barplot_of_Rate_of_InDel.png is skipped")
  message("[l10]Barplot_of_Rate_of_Unmodified.png is skipped")
}

if(analysis.mode %in% c(1)){
  MakeFoldBoxPlot(
    machiato.rate.table[,c("enhancement.1.2.precise.knockin", "reduction.1.2.imprecise.knockin", "reduction.1.2.indel", "reduction.1.2.unmodified")],
    c("#FF540D", "#FF0DFF", "#FF0DFF", "#FF0DFF"),
    c(FALSE, TRUE, TRUE, TRUE),
    c(round(min(-log2(machiato.rate.table[,c("enhancement.1.2.precise.knockin", "reduction.1.2.imprecise.knockin", "reduction.1.2.indel", "reduction.1.2.unmodified")])) - 0.5, digits = 0),
      round(max(log2(machiato.rate.table[,c("enhancement.1.2.precise.knockin", "reduction.1.2.imprecise.knockin", "reduction.1.2.indel", "reduction.1.2.unmodified")])) + 0.5, digits = 0)),
    file.path(comparison.analysis.dir, "[l11]Barplot_of_Enhancement_Reduction.png")
  )
}else if(analysis.mode %in% c(4)){
  MakeFoldBoxPlot(
    machiato.rate.table[,c("reduction.1.2.indel", "reduction.1.2.unmodified")],
    c("#FF0DFF", "#FF0DFF"),
    c(TRUE, TRUE),
    c(round(min(-log2(machiato.rate.table[,c("reduction.1.2.indel", "reduction.1.2.unmodified")])) - 0.5, digits = 0),
      round(max(log2(machiato.rate.table[,c("reduction.1.2.indel", "reduction.1.2.unmodified")])) + 0.5, digits = 0)),
    file.path(comparison.analysis.dir, "[l11]Barplot_of_Enhancement_Reduction.png")
  )
}else{
  message("[l11]Barplot_of_Enhancement_Reduction.png is skipped")
}

if(analysis.mode %in% c(1)){
  # T-test on precise knock-in
  MakePairedStaticalProfile(
    machiato.rate.table$condition1.precise.knockin,
    machiato.rate.table$condition2.precise.knockin,
    file.path(comparison.analysis.dir, "[l12]Paired_t-test_on_precise_knock-in.txt")
  )

  # T-test on imprecise knock-in
  MakePairedStaticalProfile(
    machiato.rate.table$condition1.imprecise.knockin,
    machiato.rate.table$condition2.imprecise.knockin,
    file.path(comparison.analysis.dir, "[l13]Paired_t-test_on_imprecise_knock-in.txt")
  )
}else{
  message("[l12]Paired_t-test_on_precise_knock-in.txt is skipped")
  message("[l13]Paired_t-test_on_imprecise_knock-in.txt is skipped")
}

if(analysis.mode %in% c(1, 4)){
  # T-test on indel
  MakePairedStaticalProfile(
    machiato.rate.table$condition1.indel,
    machiato.rate.table$condition2.indel,
    file.path(comparison.analysis.dir, "[l14]Paired_t-test_on_indel.txt")
  )

  # T-test on unmodified
  MakePairedStaticalProfile(
    machiato.rate.table$condition1.unmodified,
    machiato.rate.table$condition2.unmodified,
    file.path(comparison.analysis.dir, "[l15]Paired_t-test_on_unmodified.txt")
  )
}else{
  message("[l14]Paired_t-test_on_indel.txt is skipped")
  message("[l15]Paired_t-test_on_unmodified.txt is skipped")
}

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
    file.path(mutation.analysis.dir, "[a1-6]Barplot_of_insert_position_of_condition2_reads_by_enhancement_group"),
    file.path(mutation.analysis.dir, "[a1-7]Barplot_of_insert_position_of_negative.ctrl_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[a1-8]Barplot_of_insert_position_of_negative.ctrl_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[a1-9]Barplot_of_insert_position_of_negative.ctrl_reads_by_enhancement_group")
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
    file.path(mutation.analysis.dir, "[a2-6]Barplot_of_deletion_position_of_condition2_reads_by_enhancement_group"),
    file.path(mutation.analysis.dir, "[a2-7]Barplot_of_deletion_position_of_negative.ctrl_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[a2-8]Barplot_of_deletion_position_of_negative.ctrl_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[a2-9]Barplot_of_deletion_position_of_negative.ctrl_reads_by_enhancement_group")
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
    file.path(mutation.analysis.dir, "[a3-6]Barplot_of_no_overlap_deletion_position_of_condition2_reads_by_enhancement_group"),
    file.path(mutation.analysis.dir, "[a3-7]Barplot_of_no_overlap_deletion_position_of_negative.ctrl_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[a3-8]Barplot_of_no_overlap_deletion_position_of_negative.ctrl_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[a3-9]Barplot_of_no_overlap_deletion_position_of_negative.ctrl_reads_by_enhancement_group")
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
    file.path(mutation.analysis.dir, "[a4-6]Barplot_of_mutation_substitution_position_of_condition2_reads_by_enhancement_group"),
    file.path(mutation.analysis.dir, "[a4-7]Barplot_of_mutation_substitution_position_of_negative.ctrl_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[a4-8]Barplot_of_mutation_substitution_position_of_negative.ctrl_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[a4-9]Barplot_of_mutation_substitution_position_of_negative.ctrl_reads_by_enhancement_group")
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
    file.path(mutation.analysis.dir, "[b-6]Barplot_of_mutation_indel_size_of_condition2_reads_by_enhancement_group"),
    file.path(mutation.analysis.dir, "[b-7]Barplot_of_mutation_indel_size_of_negative.ctrl_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[b-8]Barplot_of_mutation_indel_size_of_negative.ctrl_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[b-9]Barplot_of_mutation_indel_size_of_negative.ctrl_reads_by_enhancement_group")
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
    file.path(mutation.analysis.dir, "[c-6]Barplot_of_mutation_substitution_number_of_condition2_reads_by_enhancement_group"),
    file.path(mutation.analysis.dir, "[c-7]Barplot_of_mutation_substitution_number_of_negative.ctrl_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[c-8]Barplot_of_mutation_substitution_number_of_negative.ctrl_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[c-9]Barplot_of_mutation_substitution_number_of_negative.ctrl_reads_by_enhancement_group")
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
    file.path(mutation.analysis.dir, "[d-6]Piechart_of_max_indel_size_class_in_condition2_reads_among_enhancement"),
    file.path(mutation.analysis.dir, "[d-7]Piechart_of_max_indel_size_class_in_negative.ctrl_reads_among_condition1"),
    file.path(mutation.analysis.dir, "[d-8]Piechart_of_max_indel_size_class_in_negative.ctrl_reads_among_condition2"),
    file.path(mutation.analysis.dir, "[d-9]Piechart_of_max_indel_size_class_in_negative.ctrl_reads_among_enhancement")
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
    file.path(mutation.analysis.dir, "[e-6]Piechart_of_max_mmej_indel_size_class_in_condition2_reads_among_enhancement"),
    file.path(mutation.analysis.dir, "[e-7]Piechart_of_max_mmej_indel_size_class_in_negative.ctrl_reads_among_condition1"),
    file.path(mutation.analysis.dir, "[e-8]Piechart_of_max_mmej_indel_size_class_in_negative.ctrl_reads_among_condition2"),
    file.path(mutation.analysis.dir, "[e-9]Piechart_of_max_mmej_indel_size_class_in_negative.ctrl_reads_among_enhancement")
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
    file.path(mutation.analysis.dir, "[f-6]Piechart_of_max_nhej_indel_size_class_in_condition2_reads_among_enhancement"),
    file.path(mutation.analysis.dir, "[f-7]Piechart_of_max_nhej_indel_size_class_in_negative.ctrl_reads_among_condition1"),
    file.path(mutation.analysis.dir, "[f-8]Piechart_of_max_nhej_indel_size_class_in_negative.ctrl_reads_among_condition2"),
    file.path(mutation.analysis.dir, "[f-9]Piechart_of_max_nhej_indel_size_class_in_negative.ctrl_reads_among_enhancement")
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
    file.path(mutation.analysis.dir, "[g-6]Piechart_of_total_reads_indel_size_class_in_condition2_reads_among_enhancement"),
    file.path(mutation.analysis.dir, "[g-7]Piechart_of_total_reads_indel_size_class_in_negative.ctrl_reads_among_condition1"),
    file.path(mutation.analysis.dir, "[g-8]Piechart_of_total_reads_indel_size_class_in_negative.ctrl_reads_among_condition2"),
    file.path(mutation.analysis.dir, "[g-9]Piechart_of_total_reads_indel_size_class_in_negative.ctrl_reads_among_enhancement")
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
    file.path(mutation.analysis.dir, "[h-6]Piechart_of_total_reads_mmej_indel_size_class_in_condition2_reads_among_enhancement"),
    file.path(mutation.analysis.dir, "[h-7]Piechart_of_total_reads_mmej_indel_size_class_in_negative.ctrl_reads_among_condition1"),
    file.path(mutation.analysis.dir, "[h-8]Piechart_of_total_reads_mmej_indel_size_class_in_negative.ctrl_reads_among_condition2"),
    file.path(mutation.analysis.dir, "[h-9]Piechart_of_total_reads_mmej_indel_size_class_in_negative.ctrl_reads_among_enhancement")
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
    file.path(mutation.analysis.dir, "[i-6]Piechart_of_total_reads_nhej_indel_size_class_in_condition2_reads_among_enhancement"),
    file.path(mutation.analysis.dir, "[i-7]Piechart_of_total_reads_nhej_indel_size_class_in_negative.ctrl_reads_among_condition1"),
    file.path(mutation.analysis.dir, "[i-8]Piechart_of_total_reads_nhej_indel_size_class_in_negative.ctrl_reads_among_condition2"),
    file.path(mutation.analysis.dir, "[i-9]Piechart_of_total_reads_nhej_indel_size_class_in_negative.ctrl_reads_among_enhancement")
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
    file.path(mutation.analysis.dir, "[j-6]Barplot_of_mutation_MMEJ-mediated_indel_size_of_condition2_reads_by_enhancement_group"),
    file.path(mutation.analysis.dir, "[j-7]Barplot_of_mutation_MMEJ-mediated_indel_size_of_negative.ctrl_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[j-8]Barplot_of_mutation_MMEJ-mediated_indel_size_of_negative.ctrl_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[j-9]Barplot_of_mutation_MMEJ-mediated_indel_size_of_negative.ctrl_reads_by_enhancement_group")
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
    file.path(mutation.analysis.dir, "[k-6]Barplot_of_mutation_NHEJ-mediated_indel_size_of_condition2_reads_by_enhancement_group"),
    file.path(mutation.analysis.dir, "[k-7]Barplot_of_mutation_NHEJ-mediated_indel_size_of_negative.ctrl_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[k-8]Barplot_of_mutation_NHEJ-mediated_indel_size_of_negative.ctrl_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[k-9]Barplot_of_mutation_NHEJ-mediated_indel_size_of_negative.ctrl_reads_by_enhancement_group")
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
    file.path(mutation.analysis.dir, "[l-6]Barplot_of_mutation_microhomology_size_of_condition2_reads_by_enhancement_group"),
    file.path(mutation.analysis.dir, "[l-7]Barplot_of_mutation_microhomology_size_of_negative.ctrl_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[l-8]Barplot_of_mutation_microhomology_size_of_negative.ctrl_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[l-9]Barplot_of_mutation_microhomology_size_of_negative.ctrl_reads_by_enhancement_group")
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
    file.path(mutation.analysis.dir, "[m-6]Barplot_of_mutation_MMEJ-mediated_intervening_size_of_condition2_reads_by_enhancement_group"),
    file.path(mutation.analysis.dir, "[m-7]Barplot_of_mutation_MMEJ-mediated_intervening_size_of_negative.ctrl_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[m-8]Barplot_of_mutation_MMEJ-mediated_intervening_size_of_negative.ctrl_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[m-9]Barplot_of_mutation_MMEJ-mediated_intervening_size_of_negative.ctrl_reads_by_enhancement_group")
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
    file.path(mutation.analysis.dir, "[n-6]Barplot_of_mutation_NHEJ-mediated_intervening_size_of_condition2_reads_by_enhancement_group"),
    file.path(mutation.analysis.dir, "[n-7]Barplot_of_mutation_NHEJ-mediated_intervening_size_of_negative.ctrl_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[n-8]Barplot_of_mutation_NHEJ-mediated_intervening_size_of_negative.ctrl_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[n-9]Barplot_of_mutation_NHEJ-mediated_intervening_size_of_negative.ctrl_reads_by_enhancement_group")
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
    file.path(mutation.analysis.dir, "[o-6]Countplot_of_mutation_microhomology_intervening_size_of_condition2_reads_by_enhancement_group"),
    file.path(mutation.analysis.dir, "[o-7]Countplot_of_mutation_microhomology_intervening_size_of_negative.ctrl_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[o-8]Countplot_of_mutation_microhomology_intervening_size_of_negative.ctrl_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[o-9]Countplot_of_mutation_microhomology_intervening_size_of_negative.ctrl_reads_by_enhancement_group")
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
    file.path(mutation.analysis.dir, "[p-6]Barplot_of_microhomology_sequence_frequency_of_condition2_reads_by_enhancement_group"),
    file.path(mutation.analysis.dir, "[p-7]Barplot_of_microhomology_sequence_frequency_of_negative.ctrl_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[p-8]Barplot_of_microhomology_sequence_frequency_of_negative.ctrl_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[p-9]Barplot_of_microhomology_sequence_frequency_of_negative.ctrl_reads_by_enhancement_group")
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
    file.path(mutation.analysis.dir, "[q-6]Barplot_of_mono-microhomology_sequence_frequency_of_condition2_reads_by_enhancement_group"),
    file.path(mutation.analysis.dir, "[q-7]Barplot_of_mono-microhomology_sequence_frequency_of_negative.ctrl_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[q-8]Barplot_of_mono-microhomology_sequence_frequency_of_negative.ctrl_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[q-9]Barplot_of_mono-microhomology_sequence_frequency_of_negative.ctrl_reads_by_enhancement_group")
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
    file.path(mutation.analysis.dir, "[r-6]Barplot_of_di-microhomology_sequence_frequency_of_condition2_reads_by_enhancement_group"),
    file.path(mutation.analysis.dir, "[r-7]Barplot_of_di-microhomology_sequence_frequency_of_negative.ctrl_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[r-8]Barplot_of_di-microhomology_sequence_frequency_of_negative.ctrl_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[r-9]Barplot_of_di-microhomology_sequence_frequency_of_negative.ctrl_reads_by_enhancement_group")
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
    file.path(mutation.analysis.dir, "[s-6]Barplot_of_intervening_length_microhomology_sequence_frequency_of_condition2_reads_by_enhancement_group"),
    file.path(mutation.analysis.dir, "[s-7]Barplot_of_intervening_length_microhomology_sequence_frequency_of_negative.ctrl_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[s-8]Barplot_of_intervening_length_microhomology_sequence_frequency_of_negative.ctrl_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[s-9]Barplot_of_intervening_length_microhomology_sequence_frequency_of_negative.ctrl_reads_by_enhancement_group")
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
    file.path(mutation.analysis.dir, "[t-6]Barplot_of_intervening_length_microhomology_sequence_frequency(per 2nt)_of_condition2_reads_by_enhancement_group"),
    file.path(mutation.analysis.dir, "[t-7]Barplot_of_intervening_length_microhomology_sequence_frequency(per 2nt)_of_negative.ctrl_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[t-8]Barplot_of_intervening_length_microhomology_sequence_frequency(per 2nt)_of_negative.ctrl_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[t-9]Barplot_of_intervening_length_microhomology_sequence_frequency(per 2nt)_of_negative.ctrl_reads_by_enhancement_group")
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
    file.path(mutation.analysis.dir, "[u-6]Barplot_of_intervening_sequence_frequency_of_condition2_reads_by_enhancement_group"),
    file.path(mutation.analysis.dir, "[u-7]Barplot_of_intervening_sequence_frequency_of_negative.ctrl_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[u-8]Barplot_of_intervening_sequence_frequency_of_negative.ctrl_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[u-9]Barplot_of_intervening_sequence_frequency_of_negative.ctrl_reads_by_enhancement_group")
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
    file.path(mutation.analysis.dir, "[v-6]Barplot_of_mono-intervening_sequence_frequency_of_condition2_reads_by_enhancement_group"),
    file.path(mutation.analysis.dir, "[v-7]Barplot_of_mono-intervening_sequence_frequency_of_negative.ctrl_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[v-8]Barplot_of_mono-intervening_sequence_frequency_of_negative.ctrl_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[v-9]Barplot_of_mono-intervening_sequence_frequency_of_negative.ctrl_reads_by_enhancement_group")
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
    file.path(mutation.analysis.dir, "[w-6]Barplot_of_di-intervening_sequence_frequency_of_condition2_reads_by_enhancement_group"),
    file.path(mutation.analysis.dir, "[w-7]Barplot_of_di-intervening_sequence_frequency_of_negative.ctrl_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[w-8]Barplot_of_di-intervening_sequence_frequency_of_negative.ctrl_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[w-9]Barplot_of_di-intervening_sequence_frequency_of_negative.ctrl_reads_by_enhancement_group")
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
    file.path(mutation.analysis.dir, "[x-6]Barplot_of_intervening_length_intervening_sequence_frequency_of_condition2_reads_by_enhancement_group"),
    file.path(mutation.analysis.dir, "[x-7]Barplot_of_intervening_length_intervening_sequence_frequency_of_negative.ctrl_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[x-8]Barplot_of_intervening_length_intervening_sequence_frequency_of_negative.ctrl_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[x-9]Barplot_of_intervening_length_intervening_sequence_frequency_of_negative.ctrl_reads_by_enhancement_group")
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
    file.path(mutation.analysis.dir, "[y-6]Barplot_of_intervening_length_intervening_sequence_frequency(per 2nt)_of_condition2_reads_by_enhancement_group"),
    file.path(mutation.analysis.dir, "[y-7]Barplot_of_intervening_length_intervening_sequence_frequency(per 2nt)_of_negative.ctrl_reads_by_condition1_group"),
    file.path(mutation.analysis.dir, "[y-8]Barplot_of_intervening_length_intervening_sequence_frequency(per 2nt)_of_negative.ctrl_reads_by_condition2_group"),
    file.path(mutation.analysis.dir, "[y-9]Barplot_of_intervening_length_intervening_sequence_frequency(per 2nt)_of_negative.ctrl_reads_by_enhancement_group")
  ),
  0.025, 0.05
)



if(analysis.mode %in% c(1, 2, 3)){

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
      file.path(imprecise.knockin.analysis.dir, "[a-6]Barplot_of_left_junction_of_indel_size_of_condition2_reads_by_enhancement_group"),
      file.path(imprecise.knockin.analysis.dir, "[a-7]Barplot_of_left_junction_of_indel_size_of_negative.ctrl_reads_by_condition1_group"),
      file.path(imprecise.knockin.analysis.dir, "[a-8]Barplot_of_left_junction_of_indel_size_of_negative.ctrl_reads_by_condition2_group"),
      file.path(imprecise.knockin.analysis.dir, "[a-9]Barplot_of_left_junction_of_indel_size_of_negative.ctrl_reads_by_enhancement_group")
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
      file.path(imprecise.knockin.analysis.dir, "[a-6]Barplot_of_right_junction_of_indel_size_of_condition2_reads_by_enhancement_group"),
      file.path(imprecise.knockin.analysis.dir, "[a-7]Barplot_of_right_junction_of_indel_size_of_negative.ctrl_reads_by_condition1_group"),
      file.path(imprecise.knockin.analysis.dir, "[a-8]Barplot_of_right_junction_of_indel_size_of_negative.ctrl_reads_by_condition2_group"),
      file.path(imprecise.knockin.analysis.dir, "[a-9]Barplot_of_right_junction_of_indel_size_of_negative.ctrl_reads_by_enhancement_group")
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
      file.path(imprecise.knockin.analysis.dir, "[b-6]Barplot_of_left_junction_substitution_number_of_condition2_reads_by_enhancement_group"),
      file.path(imprecise.knockin.analysis.dir, "[b-7]Barplot_of_left_junction_substitution_number_of_negative.ctrl_reads_by_condition1_group"),
      file.path(imprecise.knockin.analysis.dir, "[b-8]Barplot_of_left_junction_substitution_number_of_negative.ctrl_reads_by_condition2_group"),
      file.path(imprecise.knockin.analysis.dir, "[b-9]Barplot_of_left_junction_substitution_number_of_negative.ctrl_reads_by_enhancement_group")
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
      file.path(imprecise.knockin.analysis.dir, "[b-6]Barplot_of_right_junction_substitution_number_of_condition2_reads_by_enhancement_group"),
      file.path(imprecise.knockin.analysis.dir, "[b-7]Barplot_of_right_junction_substitution_number_of_negative.ctrl_reads_by_condition1_group"),
      file.path(imprecise.knockin.analysis.dir, "[b-8]Barplot_of_right_junction_substitution_number_of_negative.ctrl_reads_by_condition2_group"),
      file.path(imprecise.knockin.analysis.dir, "[b-9]Barplot_of_right_junction_substitution_number_of_negative.ctrl_reads_by_enhancement_group")
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
      file.path(imprecise.knockin.analysis.dir, "[c-6]Piechart_of_max_indel_size_class_on_left_junction_in_condition2_reads_among_enhancement"),
      file.path(imprecise.knockin.analysis.dir, "[c-7]Piechart_of_max_indel_size_class_on_left_junction_in_negative.ctrl_reads_among_condition1"),
      file.path(imprecise.knockin.analysis.dir, "[c-8]Piechart_of_max_indel_size_class_on_left_junction_in_negative.ctrl_reads_among_condition2"),
      file.path(imprecise.knockin.analysis.dir, "[c-9]Piechart_of_max_indel_size_class_on_left_junction_in_negative.ctrl_reads_among_enhancement")
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
      file.path(imprecise.knockin.analysis.dir, "[c-6]Piechart_of_max_indel_size_class_on_right_junction_in_condition2_reads_among_enhancement"),
      file.path(imprecise.knockin.analysis.dir, "[c-7]Piechart_of_max_indel_size_class_on_right_junction_in_negative.ctrl_reads_among_condition1"),
      file.path(imprecise.knockin.analysis.dir, "[c-8]Piechart_of_max_indel_size_class_on_right_junction_in_negative.ctrl_reads_among_condition2"),
      file.path(imprecise.knockin.analysis.dir, "[c-9]Piechart_of_max_indel_size_class_on_right_junction_in_negative.ctrl_reads_among_enhancement")
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
      file.path(imprecise.knockin.analysis.dir, "[d-6]Piechart_of_max_mmej_indel_size_class_on_left_junction_in_condition2_reads_among_enhancement"),
      file.path(imprecise.knockin.analysis.dir, "[d-7]Piechart_of_max_mmej_indel_size_class_on_left_junction_in_negative.ctrl_reads_among_condition1"),
      file.path(imprecise.knockin.analysis.dir, "[d-8]Piechart_of_max_mmej_indel_size_class_on_left_junction_in_negative.ctrl_reads_among_condition2"),
      file.path(imprecise.knockin.analysis.dir, "[d-9]Piechart_of_max_mmej_indel_size_class_on_left_junction_in_negative.ctrl_reads_among_enhancement")
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
      file.path(imprecise.knockin.analysis.dir, "[d-6]Piechart_of_max_mmej_indel_size_class_on_right_junction_in_condition2_reads_among_enhancement"),
      file.path(imprecise.knockin.analysis.dir, "[d-7]Piechart_of_max_mmej_indel_size_class_on_right_junction_in_negative.ctrl_reads_among_condition1"),
      file.path(imprecise.knockin.analysis.dir, "[d-8]Piechart_of_max_mmej_indel_size_class_on_right_junction_in_negative.ctrl_reads_among_condition2"),
      file.path(imprecise.knockin.analysis.dir, "[d-9]Piechart_of_max_mmej_indel_size_class_on_right_junction_in_negative.ctrl_reads_among_enhancement")
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
      file.path(imprecise.knockin.analysis.dir, "[e-6]Piechart_of_max_nhej_indel_size_class_on_left_junction_in_condition2_reads_among_enhancement"),
      file.path(imprecise.knockin.analysis.dir, "[e-7]Piechart_of_max_nhej_indel_size_class_on_left_junction_in_negative.ctrl_reads_among_condition1"),
      file.path(imprecise.knockin.analysis.dir, "[e-8]Piechart_of_max_nhej_indel_size_class_on_left_junction_in_negative.ctrl_reads_among_condition2"),
      file.path(imprecise.knockin.analysis.dir, "[e-9]Piechart_of_max_nhej_indel_size_class_on_left_junction_in_negative.ctrl_reads_among_enhancement")
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
      file.path(imprecise.knockin.analysis.dir, "[e-6]Piechart_of_max_nhej_indel_size_class_on_right_junction_in_condition2_reads_among_enhancement"),
      file.path(imprecise.knockin.analysis.dir, "[e-7]Piechart_of_max_nhej_indel_size_class_on_right_junction_in_negative.ctrl_reads_among_condition1"),
      file.path(imprecise.knockin.analysis.dir, "[e-8]Piechart_of_max_nhej_indel_size_class_on_right_junction_in_negative.ctrl_reads_among_condition2"),
      file.path(imprecise.knockin.analysis.dir, "[e-9]Piechart_of_max_nhej_indel_size_class_on_right_junction_in_negative.ctrl_reads_among_enhancement")
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
      file.path(imprecise.knockin.analysis.dir, "[f-6]Piechart_of_total_reads_indel_size_class_on_left_junction_in_condition2_reads_among_enhancement"),
      file.path(imprecise.knockin.analysis.dir, "[f-7]Piechart_of_total_reads_indel_size_class_on_left_junction_in_negative.ctrl_reads_among_condition1"),
      file.path(imprecise.knockin.analysis.dir, "[f-8]Piechart_of_total_reads_indel_size_class_on_left_junction_in_negative.ctrl_reads_among_condition2"),
      file.path(imprecise.knockin.analysis.dir, "[f-9]Piechart_of_total_reads_indel_size_class_on_left_junction_in_negative.ctrl_reads_among_enhancement")
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
      file.path(imprecise.knockin.analysis.dir, "[f-6]Piechart_of_total_reads_indel_size_class_on_right_junction_in_condition2_reads_among_enhancement"),
      file.path(imprecise.knockin.analysis.dir, "[f-7]Piechart_of_total_reads_indel_size_class_on_right_junction_in_negative.ctrl_reads_among_condition1"),
      file.path(imprecise.knockin.analysis.dir, "[f-8]Piechart_of_total_reads_indel_size_class_on_right_junction_in_negative.ctrl_reads_among_condition2"),
      file.path(imprecise.knockin.analysis.dir, "[f-9]Piechart_of_total_reads_indel_size_class_on_right_junction_in_negative.ctrl_reads_among_enhancement")
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
      file.path(imprecise.knockin.analysis.dir, "[g-6]Piechart_of_total_reads_mmej_indel_size_class_on_left_junction_in_condition2_reads_among_enhancement"),
      file.path(imprecise.knockin.analysis.dir, "[g-7]Piechart_of_total_reads_mmej_indel_size_class_on_left_junction_in_negative.ctrl_reads_among_condition1"),
      file.path(imprecise.knockin.analysis.dir, "[g-8]Piechart_of_total_reads_mmej_indel_size_class_on_left_junction_in_negative.ctrl_reads_among_condition2"),
      file.path(imprecise.knockin.analysis.dir, "[g-9]Piechart_of_total_reads_mmej_indel_size_class_on_left_junction_in_negative.ctrl_reads_among_enhancement")
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
      file.path(imprecise.knockin.analysis.dir, "[g-6]Piechart_of_total_reads_mmej_indel_size_class_on_right_junction_in_condition2_reads_among_enhancement"),
      file.path(imprecise.knockin.analysis.dir, "[g-7]Piechart_of_total_reads_mmej_indel_size_class_on_right_junction_in_negative.ctrl_reads_among_condition1"),
      file.path(imprecise.knockin.analysis.dir, "[g-8]Piechart_of_total_reads_mmej_indel_size_class_on_right_junction_in_negative.ctrl_reads_among_condition2"),
      file.path(imprecise.knockin.analysis.dir, "[g-9]Piechart_of_total_reads_mmej_indel_size_class_on_right_junction_in_negative.ctrl_reads_among_enhancement")
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
      file.path(imprecise.knockin.analysis.dir, "[h-6]Piechart_of_total_reads_nhej_indel_size_class_on_left_junction_in_condition2_reads_among_enhancement"),
      file.path(imprecise.knockin.analysis.dir, "[h-7]Piechart_of_total_reads_nhej_indel_size_class_on_left_junction_in_negative.ctrl_reads_among_condition1"),
      file.path(imprecise.knockin.analysis.dir, "[h-8]Piechart_of_total_reads_nhej_indel_size_class_on_left_junction_in_negative.ctrl_reads_among_condition2"),
      file.path(imprecise.knockin.analysis.dir, "[h-9]Piechart_of_total_reads_nhej_indel_size_class_on_left_junction_in_negative.ctrl_reads_among_enhancement")
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
      file.path(imprecise.knockin.analysis.dir, "[h-6]Piechart_of_total_reads_nhej_indel_size_class_on_right_junction_in_condition2_reads_among_enhancement"),
      file.path(imprecise.knockin.analysis.dir, "[h-7]Piechart_of_total_reads_nhej_indel_size_class_on_right_junction_in_negative.ctrl_reads_among_condition1"),
      file.path(imprecise.knockin.analysis.dir, "[h-8]Piechart_of_total_reads_nhej_indel_size_class_on_right_junction_in_negative.ctrl_reads_among_condition2"),
      file.path(imprecise.knockin.analysis.dir, "[h-9]Piechart_of_total_reads_nhej_indel_size_class_on_right_junction_in_negative.ctrl_reads_among_enhancement")
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
      file.path(imprecise.knockin.analysis.dir, "[i-6]Piechart_of_max_Id-Ip_indel_size_class_on_left_junction_in_condition2_reads_among_enhancement"),
      file.path(imprecise.knockin.analysis.dir, "[i-7]Piechart_of_max_Id-Ip_indel_size_class_on_left_junction_in_negative.ctrl_reads_among_condition1"),
      file.path(imprecise.knockin.analysis.dir, "[i-8]Piechart_of_max_Id-Ip_indel_size_class_on_left_junction_in_negative.ctrl_reads_among_condition2"),
      file.path(imprecise.knockin.analysis.dir, "[i-9]Piechart_of_max_Id-Ip_indel_size_class_on_left_junction_in_negative.ctrl_reads_among_enhancement")
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
      file.path(imprecise.knockin.analysis.dir, "[i-6]Piechart_of_max_Id-Ip_indel_size_class_on_right_junction_in_condition2_reads_among_enhancement"),
      file.path(imprecise.knockin.analysis.dir, "[i-7]Piechart_of_max_Id-Ip_indel_size_class_on_right_junction_in_negative.ctrl_reads_among_condition1"),
      file.path(imprecise.knockin.analysis.dir, "[i-8]Piechart_of_max_Id-Ip_indel_size_class_on_right_junction_in_negative.ctrl_reads_among_condition2"),
      file.path(imprecise.knockin.analysis.dir, "[i-9]Piechart_of_max_Id-Ip_indel_size_class_on_right_junction_in_negative.ctrl_reads_among_enhancement")
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
      file.path(imprecise.knockin.analysis.dir, "[j-6]Piechart_of_total_reads_Id-Ip_indel_size_class_on_left_junction_in_condition2_reads_among_enhancement"),
      file.path(imprecise.knockin.analysis.dir, "[j-7]Piechart_of_total_reads_Id-Ip_indel_size_class_on_left_junction_in_negative.ctrl_reads_among_condition1"),
      file.path(imprecise.knockin.analysis.dir, "[j-8]Piechart_of_total_reads_Id-Ip_indel_size_class_on_left_junction_in_negative.ctrl_reads_among_condition2"),
      file.path(imprecise.knockin.analysis.dir, "[j-9]Piechart_of_total_reads_Id-Ip_indel_size_class_on_left_junction_in_negative.ctrl_reads_among_enhancement")
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
      file.path(imprecise.knockin.analysis.dir, "[j-6]Piechart_of_total_reads_Id-Ip_indel_size_class_on_right_junction_in_condition2_reads_among_enhancement"),
      file.path(imprecise.knockin.analysis.dir, "[j-7]Piechart_of_total_reads_Id-Ip_indel_size_class_on_right_junction_in_negative.ctrl_reads_among_condition1"),
      file.path(imprecise.knockin.analysis.dir, "[j-8]Piechart_of_total_reads_Id-Ip_indel_size_class_on_right_junction_in_negative.ctrl_reads_among_condition2"),
      file.path(imprecise.knockin.analysis.dir, "[j-9]Piechart_of_total_reads_Id-Ip_indel_size_class_on_right_junction_in_negative.ctrl_reads_among_enhancement")
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
      file.path(imprecise.knockin.analysis.dir, "[k-6]Barplot_of_left_junction_of_mmej_correct_direction_indel_size_of_condition2_reads_by_enhancement_group"),
      file.path(imprecise.knockin.analysis.dir, "[k-7]Barplot_of_left_junction_of_mmej_correct_direction_indel_size_of_negative.ctrl_reads_by_condition1_group"),
      file.path(imprecise.knockin.analysis.dir, "[k-8]Barplot_of_left_junction_of_mmej_correct_direction_indel_size_of_negative.ctrl_reads_by_condition2_group"),
      file.path(imprecise.knockin.analysis.dir, "[k-9]Barplot_of_left_junction_of_mmej_correct_direction_indel_size_of_negative.ctrl_reads_by_enhancement_group")
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
      file.path(imprecise.knockin.analysis.dir, "[k-6]Barplot_of_right_junction_of_mmej_correct_direction_indel_size_of_condition2_reads_by_enhancement_group"),
      file.path(imprecise.knockin.analysis.dir, "[k-7]Barplot_of_right_junction_of_mmej_correct_direction_indel_size_of_negative.ctrl_reads_by_condition1_group"),
      file.path(imprecise.knockin.analysis.dir, "[k-8]Barplot_of_right_junction_of_mmej_correct_direction_indel_size_of_negative.ctrl_reads_by_condition2_group"),
      file.path(imprecise.knockin.analysis.dir, "[k-9]Barplot_of_right_junction_of_mmej_correct_direction_indel_size_of_negative.ctrl_reads_by_enhancement_group")
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
      file.path(imprecise.knockin.analysis.dir, "[l-6]Barplot_of_left_junction_of_nhej_correct_direction_indel_size_of_condition2_reads_by_enhancement_group"),
      file.path(imprecise.knockin.analysis.dir, "[l-7]Barplot_of_left_junction_of_nhej_correct_direction_indel_size_of_negative.ctrl_reads_by_condition1_group"),
      file.path(imprecise.knockin.analysis.dir, "[l-8]Barplot_of_left_junction_of_nhej_correct_direction_indel_size_of_negative.ctrl_reads_by_condition2_group"),
      file.path(imprecise.knockin.analysis.dir, "[l-9]Barplot_of_left_junction_of_nhej_correct_direction_indel_size_of_negative.ctrl_reads_by_enhancement_group")
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
      file.path(imprecise.knockin.analysis.dir, "[l-6]Barplot_of_right_junction_of_nhej_correct_direction_indel_size_of_condition2_reads_by_enhancement_group"),
      file.path(imprecise.knockin.analysis.dir, "[l-7]Barplot_of_right_junction_of_nhej_correct_direction_indel_size_of_negative.ctrl_reads_by_condition1_group"),
      file.path(imprecise.knockin.analysis.dir, "[l-8]Barplot_of_right_junction_of_nhej_correct_direction_indel_size_of_negative.ctrl_reads_by_condition2_group"),
      file.path(imprecise.knockin.analysis.dir, "[l-9]Barplot_of_right_junction_of_nhej_correct_direction_indel_size_of_negative.ctrl_reads_by_enhancement_group")
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
      file.path(imprecise.knockin.analysis.dir, "[m-6]Barplot_of_left_junction_of_mmej_reverse_direction_indel_size_of_condition2_reads_by_enhancement_group"),
      file.path(imprecise.knockin.analysis.dir, "[m-7]Barplot_of_left_junction_of_mmej_reverse_direction_indel_size_of_negative.ctrl_reads_by_condition1_group"),
      file.path(imprecise.knockin.analysis.dir, "[m-8]Barplot_of_left_junction_of_mmej_reverse_direction_indel_size_of_negative.ctrl_reads_by_condition2_group"),
      file.path(imprecise.knockin.analysis.dir, "[m-9]Barplot_of_left_junction_of_mmej_reverse_direction_indel_size_of_negative.ctrl_reads_by_enhancement_group")
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
      file.path(imprecise.knockin.analysis.dir, "[m-6]Barplot_of_right_junction_of_mmej_reverse_direction_indel_size_of_condition2_reads_by_enhancement_group"),
      file.path(imprecise.knockin.analysis.dir, "[m-7]Barplot_of_right_junction_of_mmej_reverse_direction_indel_size_of_negative.ctrl_reads_by_condition1_group"),
      file.path(imprecise.knockin.analysis.dir, "[m-8]Barplot_of_right_junction_of_mmej_reverse_direction_indel_size_of_negative.ctrl_reads_by_condition2_group"),
      file.path(imprecise.knockin.analysis.dir, "[m-9]Barplot_of_right_junction_of_mmej_reverse_direction_indel_size_of_negative.ctrl_reads_by_enhancement_group")
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
      file.path(imprecise.knockin.analysis.dir, "[n-6]Barplot_of_left_junction_of_nhej_reverse_direction_indel_size_of_condition2_reads_by_enhancement_group"),
      file.path(imprecise.knockin.analysis.dir, "[n-7]Barplot_of_left_junction_of_nhej_reverse_direction_indel_size_of_negative.ctrl_reads_by_condition1_group"),
      file.path(imprecise.knockin.analysis.dir, "[n-8]Barplot_of_left_junction_of_nhej_reverse_direction_indel_size_of_negative.ctrl_reads_by_condition2_group"),
      file.path(imprecise.knockin.analysis.dir, "[n-9]Barplot_of_left_junction_of_nhej_reverse_direction_indel_size_of_negative.ctrl_reads_by_enhancement_group")
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
      file.path(imprecise.knockin.analysis.dir, "[n-6]Barplot_of_right_junction_of_nhej_reverse_direction_indel_size_of_condition2_reads_by_enhancement_group"),
      file.path(imprecise.knockin.analysis.dir, "[n-7]Barplot_of_right_junction_of_nhej_reverse_direction_indel_size_of_negative.ctrl_reads_by_condition1_group"),
      file.path(imprecise.knockin.analysis.dir, "[n-8]Barplot_of_right_junction_of_nhej_reverse_direction_indel_size_of_negative.ctrl_reads_by_condition2_group"),
      file.path(imprecise.knockin.analysis.dir, "[n-9]Barplot_of_right_junction_of_nhej_reverse_direction_indel_size_of_negative.ctrl_reads_by_enhancement_group")
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
      file.path(imprecise.knockin.analysis.dir, "[o-6]Barplot_of_left_junction_microhomology_size_of_condition2_reads_by_enhancement_group"),
      file.path(imprecise.knockin.analysis.dir, "[o-7]Barplot_of_left_junction_microhomology_size_of_negative.ctrl_reads_by_condition1_group"),
      file.path(imprecise.knockin.analysis.dir, "[o-8]Barplot_of_left_junction_microhomology_size_of_negative.ctrl_reads_by_condition2_group"),
      file.path(imprecise.knockin.analysis.dir, "[o-9]Barplot_of_left_junction_microhomology_size_of_negative.ctrl_reads_by_enhancement_group")
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
      file.path(imprecise.knockin.analysis.dir, "[o-6]Barplot_of_right_junction_microhomology_size_of_condition2_reads_by_enhancement_group"),
      file.path(imprecise.knockin.analysis.dir, "[o-7]Barplot_of_right_junction_microhomology_size_of_negative.ctrl_reads_by_condition1_group"),
      file.path(imprecise.knockin.analysis.dir, "[o-8]Barplot_of_right_junction_microhomology_size_of_negative.ctrl_reads_by_condition2_group"),
      file.path(imprecise.knockin.analysis.dir, "[o-9]Barplot_of_right_junction_microhomology_size_of_negative.ctrl_reads_by_enhancement_group")
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
      file.path(imprecise.knockin.analysis.dir, "[p-6]Barplot_of_left_junction_MMEJ-mediated_intervening_size_of_condition2_reads_by_enhancement_group"),
      file.path(imprecise.knockin.analysis.dir, "[p-7]Barplot_of_left_junction_MMEJ-mediated_intervening_size_of_negative.ctrl_reads_by_condition1_group"),
      file.path(imprecise.knockin.analysis.dir, "[p-8]Barplot_of_left_junction_MMEJ-mediated_intervening_size_of_negative.ctrl_reads_by_condition2_group"),
      file.path(imprecise.knockin.analysis.dir, "[p-9]Barplot_of_left_junction_MMEJ-mediated_intervening_size_of_negative.ctrl_reads_by_enhancement_group")
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
      file.path(imprecise.knockin.analysis.dir, "[p-6]Barplot_of_right_junction_MMEJ-mediated_intervening_size_of_condition2_reads_by_enhancement_group"),
      file.path(imprecise.knockin.analysis.dir, "[p-7]Barplot_of_right_junction_MMEJ-mediated_intervening_size_of_negative.ctrl_reads_by_condition1_group"),
      file.path(imprecise.knockin.analysis.dir, "[p-8]Barplot_of_right_junction_MMEJ-mediated_intervening_size_of_negative.ctrl_reads_by_condition2_group"),
      file.path(imprecise.knockin.analysis.dir, "[p-9]Barplot_of_right_junction_MMEJ-mediated_intervening_size_of_negative.ctrl_reads_by_enhancement_group")
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
      file.path(imprecise.knockin.analysis.dir, "[r-6]Barplot_of_left_junction_NHEJ-mediated_intervening_size_of_condition2_reads_by_enhancement_group"),
      file.path(imprecise.knockin.analysis.dir, "[r-7]Barplot_of_left_junction_NHEJ-mediated_intervening_size_of_negative.ctrl_reads_by_condition1_group"),
      file.path(imprecise.knockin.analysis.dir, "[r-8]Barplot_of_left_junction_NHEJ-mediated_intervening_size_of_negative.ctrl_reads_by_condition2_group"),
      file.path(imprecise.knockin.analysis.dir, "[r-9]Barplot_of_left_junction_NHEJ-mediated_intervening_size_of_negative.ctrl_reads_by_enhancement_group")
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
      file.path(imprecise.knockin.analysis.dir, "[r-6]Barplot_of_right_junction_NHEJ-mediated_intervening_size_of_condition2_reads_by_enhancement_group"),
      file.path(imprecise.knockin.analysis.dir, "[r-7]Barplot_of_right_junction_NHEJ-mediated_intervening_size_of_negative.ctrl_reads_by_condition1_group"),
      file.path(imprecise.knockin.analysis.dir, "[r-8]Barplot_of_right_junction_NHEJ-mediated_intervening_size_of_negative.ctrl_reads_by_condition2_group"),
      file.path(imprecise.knockin.analysis.dir, "[r-9]Barplot_of_right_junction_NHEJ-mediated_intervening_size_of_negative.ctrl_reads_by_enhancement_group")
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
      file.path(imprecise.knockin.analysis.dir, "[s-6]Countplot_of_left_junction_microhomology_intervening_size_of_condition2_reads_by_enhancement_group"),
      file.path(imprecise.knockin.analysis.dir, "[s-7]Countplot_of_left_junction_microhomology_intervening_size_of_negative.ctrl_reads_by_condition1_group"),
      file.path(imprecise.knockin.analysis.dir, "[s-8]Countplot_of_left_junction_microhomology_intervening_size_of_negative.ctrl_reads_by_condition2_group"),
      file.path(imprecise.knockin.analysis.dir, "[s-9]Countplot_of_left_junction_microhomology_intervening_size_of_negative.ctrl_reads_by_enhancement_group")
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
      file.path(imprecise.knockin.analysis.dir, "[s-6]Countplot_of_right_junction_microhomology_intervening_size_of_condition2_reads_by_enhancement_group"),
      file.path(imprecise.knockin.analysis.dir, "[s-7]Countplot_of_right_junction_microhomology_intervening_size_of_negative.ctrl_reads_by_condition1_group"),
      file.path(imprecise.knockin.analysis.dir, "[s-8]Countplot_of_right_junction_microhomology_intervening_size_of_negative.ctrl_reads_by_condition2_group"),
      file.path(imprecise.knockin.analysis.dir, "[s-9]Countplot_of_right_junction_microhomology_intervening_size_of_negative.ctrl_reads_by_enhancement_group")
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
      file.path(imprecise.knockin.analysis.dir, "[t-6]Barplot_of_left_junction_microhomology_sequence_frequency_of_condition2_reads_by_enhancement_group"),
      file.path(imprecise.knockin.analysis.dir, "[t-7]Barplot_of_left_junction_microhomology_sequence_frequency_of_negative.ctrl_reads_by_condition1_group"),
      file.path(imprecise.knockin.analysis.dir, "[t-8]Barplot_of_left_junction_microhomology_sequence_frequency_of_negative.ctrl_reads_by_condition2_group"),
      file.path(imprecise.knockin.analysis.dir, "[t-9]Barplot_of_left_junction_microhomology_sequence_frequency_of_negative.ctrl_reads_by_enhancement_group")
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
      file.path(imprecise.knockin.analysis.dir, "[t-6]Barplot_of_right_junction_microhomology_sequence_frequency_of_condition2_reads_by_enhancement_group"),
      file.path(imprecise.knockin.analysis.dir, "[t-7]Barplot_of_right_junction_microhomology_sequence_frequency_of_negative.ctrl_reads_by_condition1_group"),
      file.path(imprecise.knockin.analysis.dir, "[t-8]Barplot_of_right_junction_microhomology_sequence_frequency_of_negative.ctrl_reads_by_condition2_group"),
      file.path(imprecise.knockin.analysis.dir, "[t-9]Barplot_of_right_junction_microhomology_sequence_frequency_of_negative.ctrl_reads_by_enhancement_group")
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
      file.path(imprecise.knockin.analysis.dir, "[u-6]Barplot_of_left_junction_mono-microhomology_sequence_frequency_of_condition2_reads_by_enhancement_group"),
      file.path(imprecise.knockin.analysis.dir, "[u-7]Barplot_of_left_junction_mono-microhomology_sequence_frequency_of_negative.ctrl_reads_by_condition1_group"),
      file.path(imprecise.knockin.analysis.dir, "[u-8]Barplot_of_left_junction_mono-microhomology_sequence_frequency_of_negative.ctrl_reads_by_condition2_group"),
      file.path(imprecise.knockin.analysis.dir, "[u-9]Barplot_of_left_junction_mono-microhomology_sequence_frequency_of_negative.ctrl_reads_by_enhancement_group")
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
      file.path(imprecise.knockin.analysis.dir, "[u-6]Barplot_of_right_junction_mono-microhomology_sequence_frequency_of_condition2_reads_by_enhancement_group"),
      file.path(imprecise.knockin.analysis.dir, "[u-7]Barplot_of_right_junction_mono-microhomology_sequence_frequency_of_negative.ctrl_reads_by_condition1_group"),
      file.path(imprecise.knockin.analysis.dir, "[u-8]Barplot_of_right_junction_mono-microhomology_sequence_frequency_of_negative.ctrl_reads_by_condition2_group"),
      file.path(imprecise.knockin.analysis.dir, "[u-9]Barplot_of_right_junction_mono-microhomology_sequence_frequency_of_negative.ctrl_reads_by_enhancement_group")
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
      file.path(imprecise.knockin.analysis.dir, "[v-6]Barplot_of_left_junction_di-microhomology_sequence_frequency_of_condition2_reads_by_enhancement_group"),
      file.path(imprecise.knockin.analysis.dir, "[v-7]Barplot_of_left_junction_di-microhomology_sequence_frequency_of_negative.ctrl_reads_by_condition1_group"),
      file.path(imprecise.knockin.analysis.dir, "[v-8]Barplot_of_left_junction_di-microhomology_sequence_frequency_of_negative.ctrl_reads_by_condition2_group"),
      file.path(imprecise.knockin.analysis.dir, "[v-9]Barplot_of_left_junction_di-microhomology_sequence_frequency_of_negative.ctrl_reads_by_enhancement_group")
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
      file.path(imprecise.knockin.analysis.dir, "[v-6]Barplot_of_right_junction_di-microhomology_sequence_frequency_of_condition2_reads_by_enhancement_group"),
      file.path(imprecise.knockin.analysis.dir, "[v-7]Barplot_of_right_junction_di-microhomology_sequence_frequency_of_negative.ctrl_reads_by_condition1_group"),
      file.path(imprecise.knockin.analysis.dir, "[v-8]Barplot_of_right_junction_di-microhomology_sequence_frequency_of_negative.ctrl_reads_by_condition2_group"),
      file.path(imprecise.knockin.analysis.dir, "[v-9]Barplot_of_right_junction_di-microhomology_sequence_frequency_of_negative.ctrl_reads_by_enhancement_group")
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
      file.path(imprecise.knockin.analysis.dir, "[w-6]Barplot_of_left_junction_intervening_length_microhomology_sequence_frequency_of_condition2_reads_by_enhancement_group"),
      file.path(imprecise.knockin.analysis.dir, "[w-7]Barplot_of_left_junction_intervening_length_microhomology_sequence_frequency_of_negative.ctrl_reads_by_condition1_group"),
      file.path(imprecise.knockin.analysis.dir, "[w-8]Barplot_of_left_junction_intervening_length_microhomology_sequence_frequency_of_negative.ctrl_reads_by_condition2_group"),
      file.path(imprecise.knockin.analysis.dir, "[w-9]Barplot_of_left_junction_intervening_length_microhomology_sequence_frequency_of_negative.ctrl_reads_by_enhancement_group")
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
      file.path(imprecise.knockin.analysis.dir, "[w-6]Barplot_of_right_junction_intervening_length_microhomology_sequence_frequency_of_condition2_reads_by_enhancement_group"),
      file.path(imprecise.knockin.analysis.dir, "[w-7]Barplot_of_right_junction_intervening_length_microhomology_sequence_frequency_of_negative.ctrl_reads_by_condition1_group"),
      file.path(imprecise.knockin.analysis.dir, "[w-8]Barplot_of_right_junction_intervening_length_microhomology_sequence_frequency_of_negative.ctrl_reads_by_condition2_group"),
      file.path(imprecise.knockin.analysis.dir, "[w-9]Barplot_of_right_junction_intervening_length_microhomology_sequence_frequency_of_negative.ctrl_reads_by_enhancement_group")
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
      file.path(imprecise.knockin.analysis.dir, "[x-6]Barplot_of_left_junction_intervening_length_microhomology_sequence_frequency(per 2nt)_of_condition2_reads_by_enhancement_group"),
      file.path(imprecise.knockin.analysis.dir, "[x-7]Barplot_of_left_junction_intervening_length_microhomology_sequence_frequency(per 2nt)_of_negative.ctrl_reads_by_condition1_group"),
      file.path(imprecise.knockin.analysis.dir, "[x-8]Barplot_of_left_junction_intervening_length_microhomology_sequence_frequency(per 2nt)_of_negative.ctrl_reads_by_condition2_group"),
      file.path(imprecise.knockin.analysis.dir, "[x-9]Barplot_of_left_junction_intervening_length_microhomology_sequence_frequency(per 2nt)_of_negative.ctrl_reads_by_enhancement_group")
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
      file.path(imprecise.knockin.analysis.dir, "[x-6]Barplot_of_right_junction_intervening_length_microhomology_sequence_frequency(per 2nt)_of_condition2_reads_by_enhancement_group"),
      file.path(imprecise.knockin.analysis.dir, "[x-7]Barplot_of_right_junction_intervening_length_microhomology_sequence_frequency(per 2nt)_of_negative.ctrl_reads_by_condition1_group"),
      file.path(imprecise.knockin.analysis.dir, "[x-8]Barplot_of_right_junction_intervening_length_microhomology_sequence_frequency(per 2nt)_of_negative.ctrl_reads_by_condition2_group"),
      file.path(imprecise.knockin.analysis.dir, "[x-9]Barplot_of_right_junction_intervening_length_microhomology_sequence_frequency(per 2nt)_of_negative.ctrl_reads_by_enhancement_group")
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
      file.path(imprecise.knockin.analysis.dir, "[y-6]Barplot_of_left_junction_intervening_sequence_frequency_of_condition2_reads_by_enhancement_group"),
      file.path(imprecise.knockin.analysis.dir, "[y-7]Barplot_of_left_junction_intervening_sequence_frequency_of_negative.ctrl_reads_by_condition1_group"),
      file.path(imprecise.knockin.analysis.dir, "[y-8]Barplot_of_left_junction_intervening_sequence_frequency_of_negative.ctrl_reads_by_condition2_group"),
      file.path(imprecise.knockin.analysis.dir, "[y-9]Barplot_of_left_junction_intervening_sequence_frequency_of_negative.ctrl_reads_by_enhancement_group")
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
      file.path(imprecise.knockin.analysis.dir, "[y-6]Barplot_of_right_junction_intervening_sequence_frequency_of_condition2_reads_by_enhancement_group"),
      file.path(imprecise.knockin.analysis.dir, "[y-7]Barplot_of_right_junction_intervening_sequence_frequency_of_negative.ctrl_reads_by_condition1_group"),
      file.path(imprecise.knockin.analysis.dir, "[y-8]Barplot_of_right_junction_intervening_sequence_frequency_of_negative.ctrl_reads_by_condition2_group"),
      file.path(imprecise.knockin.analysis.dir, "[y-9]Barplot_of_right_junction_intervening_sequence_frequency_of_negative.ctrl_reads_by_enhancement_group")
    ),
    0.2, 0.3
  )
}

message("---All process is completed.---")
quit(save = "no")










