

debug.flag <- FALSE

# version 1.2

message("
                                                                            
,-.-.     ,---.|    |,---.|         ,---.          |                        
| | |,---.|    |---.||---||--- ,---.|---|,---.,---.|    ,   .,---,,---.,---.
| | |,---||    |   |||   ||    |   ||   ||   |,---||    |   | .-' |---'|    
` ' '`---^`---'`   '``   '`---'`---'`   '`   '`---^`---'`---|'---'`---'`    
                                                        `---'               
version 1.2.beta
")

# Get directory path
if(debug.flag == TRUE){
    script.dir <- "./MaChIAtoAnalyzer"
}else{
    initial.options <- commandArgs(trailingOnly = FALSE)
    file.arg.name <- "--file="
    script.fn <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
    script.dir <- dirname(script.fn)
}
ref.dir <- file.path(script.dir, "ref")
source(file.path(script.dir, "MaChIAtoAnalyzerPackageInstaller.R"))
source(file.path(script.dir, "MaChIAtoAnalyzerFunctions.R"))
source(file.path(script.dir, "CalcMicroHomologyPropPK-class.R"))
source(file.path(script.dir, "CalcMHScore-class.R"))

# library(randomForest)
library(ggplot2)
library(dplyr)
library(pryr)
library(reshape2)
library(maptools)
library(scales)
library(ggpubr)

ignore.list.fn <- character(0)
if(debug.flag == TRUE){
    input.dir <- "./20190722_MaChIAto2_data_analysis"
    output.prefix <- "./MaChIAtoAnalyzer_output"
    extra.fn <- "./20190722_MaChIAto2_data_analysis/extra_data.csv"
    ignore.list.fn <- "./20190722_MaChIAto2_data_analysis/ignore_list.csv"
}else{
    # Get commandline arguments
    input.dir <- GetAbsolutePath(commandArgs(trailingOnly=TRUE)[1])
    # input.fn <- GetAbsolutePath(commandArgs(trailingOnly=TRUE)[1])
    output.prefix <- GetAbsolutePath(commandArgs(trailingOnly=TRUE)[2])
    extra.fn <- GetAbsolutePath(commandArgs(trailingOnly=TRUE)[3])
    ignore.list.fn <- GetAbsolutePath(commandArgs(trailingOnly=TRUE)[4])
}

# Make output directory
time.id <- format(Sys.time(), "%Y%m%d%H%M%S")
output.name <- paste("MaChIAtoAnalyzer_at_", time.id, sep = "")
output.dir <- file.path(output.prefix, output.name)
figure.dir <- file.path(output.dir, "fig")
system(paste("mkdir", output.dir, sep=" "))
system(paste("mkdir", figure.dir, sep=" "))
message(paste("Creating Folder", output.dir, sep=" "))
message(paste("done"))

all.count.table <- read.csv(file.path(input.dir, "sum_ALL_read_countdata.csv"), stringsAsFactors=FALSE)
all.count.table <- cbind(all.count.table,
    sample.target = as.character(sapply(all.count.table$name, function(x){strsplit(x, "-")[[1]]})[1,]),
    sample.label = as.character(sapply(all.count.table$name, function(x){strsplit(x, "-")[[1]]})[2,]),
    stringsAsFactors=FALSE)
all.rate.table <- read.csv(file.path(input.dir, "sum_ALL_read_rates.csv"), stringsAsFactors=FALSE)
extra.table <- read.csv(extra.fn, stringsAsFactors=FALSE)

################################################################################################################################
### Filtering
################################################################################################################################
message("Filtering data")
if(ignore.list.fn != ""){
    remove.target.vec <- gsub(" ", "",unlist(read.csv(ignore.list.fn, header = FALSE, stringsAsFactors = FALSE)[1,]))
}else{
    remove.target.vec <- character(0)
}

# Remove gene contains low number of reads
message("Check amount of reads...")
has.sufficint.reads.vec <- apply(all.count.table, MARGIN = 1, function(row){
  if(as.numeric(row["total"]) < 1000){
    message(paste0(row["name"], " dosen't have sufficint reads (<1,000 reads)."))
    return(FALSE)
  }else{
    return(TRUE)
  }
})
remove.target.vec <- c(remove.target.vec, unique(as.character(all.count.table$sample.target)[!has.sufficint.reads.vec]))

# Remove gene contains unusual Unmodified rate in untreated
message("Check conservative rate in untreated samples...")
is.conservative.sample.vec <- apply(all.count.table, MARGIN = 1, function(row){
    if(row["sample.label"] != "A"){
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

# Remove gene contains no difference between wt and negative.ctrl on untreated
message("Check whether efficary of editing exists in negative ctrl...")
for(focus.sample.target in unique(as.character(all.count.table$sample.target))){

    # load untreated
    temp.unmodified.df <- filter(all.count.table, sample.target == focus.sample.target & sample.label == "A")

    # load negative control
    temp.negative.ctrl.df <- filter(all.count.table, sample.target == focus.sample.target & sample.label == "B")
    
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
        message(paste0(focus.sample.target, " dosen't have significant change between untreated and negative.ctrl. (chisq.test, p >= 0.1%)"))
        remove.target.vec <- c(remove.target.vec, focus.sample.target)
    }
}

# Remove gene contains no precise knock-in in condition1 & condition2 sample
message("Check whether precise knock-in exists in condition1 and condition2...")
for(focus.sample.target in unique(as.character(all.count.table$sample.target))){

    # load untreated
    temp.condition1.df <- filter(all.count.table, sample.target == focus.sample.target & sample.label == "C")

    # load negative control
    temp.condition2.df <- filter(all.count.table, sample.target == focus.sample.target & sample.label == "D")
    
    if(temp.condition1.df$PRECISE_KNOCK_IN < 1 | temp.condition2.df$PRECISE_KNOCK_IN < 1){
      message(paste0(focus.sample.target, " has no precise knock-in in condition1 or condition2."))
      remove.target.vec <- c(remove.target.vec, focus.sample.target)
    }
}

message("Remove...")
remove.target.vec <- unique(remove.target.vec)
fltr.rate.table <- all.rate.table
fltr.count.table <- all.count.table
for(remove.target in remove.target.vec){
    message(remove.target)
    fltr.rate.table <- fltr.rate.table[as.character(fltr.rate.table[, "group.name"]) != as.character(remove.target),]
    fltr.count.table <- fltr.count.table[as.character(fltr.count.table[, "sample.target"]) != as.character(remove.target),]
}
message("Filtering is done.")
#############################################################################################################################################################################################################3

# split dataframe into list
fltr.rate.list <- split(fltr.rate.table, fltr.rate.table$group.name)

# Make Value set by sample.type
ctrl.value.set <- as.data.frame(MakeValueSet(fltr.rate.list, "ctrl"), stringsAsFactors = FALSE)
cas9.1.value.set <- as.data.frame(MakeValueSet(fltr.rate.list, "Cas9_1"), stringsAsFactors = FALSE)
ki.1.value.set <- as.data.frame(MakeValueSet(fltr.rate.list, "KI_1"), stringsAsFactors = FALSE)
ki.2.value.set <- as.data.frame(MakeValueSet(fltr.rate.list, "KI_2"), stringsAsFactors = FALSE)

value.rowname <- c("Accuracy" ,"Conservative" ,"Editing.Efficiency"
    , "Indels.Editing.Efficiency" ,"Integration.zygomorphy" ,"Knock.in.Efficiency"
    , "Left.MMEJ.tendency" ,"Left.MMEJ.NHEJ.joinability" ,"Left.NHEJ.tendency"
    , "Left.accuracy" ,"Left.efficiency" ,"Left.occupancy" ,"Left.reversibility"
    , "Left.unpredictability" ,"Left.unpredictable.knock.in.deletion"
    , "Left.unpredictable.knock.in.insert" ,"NHEJ.assisted.Knock.in.zygomorphy"
    , "Occupancy" ,"Precise.Knock.in.Efficiency" ,"Precise.Knock.in.zygomorphy"
    , "Right.MMEJ.tendency" ,"Right.MMEJ.NHEJ.joinability" ,"Right.NHEJ.tendency"
    , "Right.accuracy" ,"Right.efficiency" ,"Right.occupancy" ,"Right.reversibility"
    , "Right.unpredictability" ,"Right.unpredictable.knock.in.deletion"
    , "Right.unpredictable.knock.in.insert"
)

# second filterling
# non negaive
left.eff.ki.1.value.set <- filter(ki.1.value.set, Left.efficiency > 0)
right.eff.ki.1.value.set <- filter(ki.1.value.set, Right.efficiency > 0)
left.acc.ki.1.value.set <- filter(ki.1.value.set, Left.accuracy > 0)
right.acc.ki.1.value.set <- filter(ki.1.value.set, Right.accuracy > 0)
left.eff.ki.2.value.set <- filter(ki.2.value.set, Left.efficiency > 0)
right.eff.ki.2.value.set <- filter(ki.2.value.set, Right.efficiency > 0)
left.acc.ki.2.value.set <- filter(ki.2.value.set, Left.accuracy > 0)
right.acc.ki.2.value.set <- filter(ki.2.value.set, Right.accuracy > 0)

whole.eff.ki.1.value.set <- filter(ki.1.value.set, Precise.Knock.in.Efficiency > 0)
whole.acc.ki.1.value.set <- filter(ki.1.value.set, Accuracy > 0)
whole.eff.ki.2.value.set <- filter(ki.2.value.set, Precise.Knock.in.Efficiency > 0)
whole.acc.ki.2.value.set <- filter(ki.2.value.set, Accuracy > 0)

### Convert value for the following calculation
# raw
left.eff.ki.1.value.set[, "Left.efficiency"] <- left.eff.ki.1.value.set[, "Left.efficiency"]
right.eff.ki.1.value.set[, "Right.efficiency"] <- right.eff.ki.1.value.set[, "Right.efficiency"]
left.acc.ki.1.value.set[, "Left.accuracy"] <- left.acc.ki.1.value.set[, "Left.accuracy"]
right.acc.ki.1.value.set[, "Right.accuracy"] <- right.acc.ki.1.value.set[, "Right.accuracy"]
left.eff.ki.2.value.set[, "Left.efficiency"] <- left.eff.ki.2.value.set[, "Left.efficiency"]
right.eff.ki.2.value.set[, "Right.efficiency"] <- right.eff.ki.2.value.set[, "Right.efficiency"]
left.acc.ki.2.value.set[, "Left.accuracy"] <- left.acc.ki.2.value.set[, "Left.accuracy"]
right.acc.ki.2.value.set[, "Right.accuracy"] <- right.acc.ki.2.value.set[, "Right.accuracy"]
# log
whole.eff.ki.1.value.set[, "Precise.Knock.in.Efficiency"] <- whole.eff.ki.1.value.set[, "Precise.Knock.in.Efficiency"]
whole.acc.ki.1.value.set[, "Accuracy"] <- whole.acc.ki.1.value.set[, "Accuracy"]
whole.eff.ki.2.value.set[, "Precise.Knock.in.Efficiency"] <- whole.eff.ki.2.value.set[, "Precise.Knock.in.Efficiency"]
whole.acc.ki.2.value.set[, "Accuracy"] <- whole.acc.ki.2.value.set[, "Accuracy"]

enh.store.vec <- c("name", "lmh", "rmh", "bmh", "protospacer", "fourtybase.region.surrounding.targetseq", "sample.name", "group.name")
# enhancement including InDel and unmodified
whole.eff.enh.value.set <- data.frame(whole.eff.ki.2.value.set[, enh.store.vec]
    , Efficiency.Enhancement = log2(whole.eff.ki.2.value.set[, "Precise.Knock.in.Efficiency"]) - log2(whole.eff.ki.1.value.set[, "Precise.Knock.in.Efficiency"])
    , stringsAsFactors = FALSE
)
# enhancement excluding InDel and unmodified
whole.acc.enh.value.set <- data.frame(whole.acc.ki.2.value.set[, enh.store.vec]
    , Accuracy.Enhancement = log2(whole.acc.ki.2.value.set[, "Accuracy"]) - log2(whole.acc.ki.1.value.set[, "Accuracy"])
    , stringsAsFactors = FALSE
)

####################################################################################################################################

# Make annotation set

# left.eff.ki.1.value.set <- Removeoutliers(left.eff.ki.1.value.set, "Left.efficiency")
options(warn = -1)
left.eff.ki.1.annotation.set <- apply(left.eff.ki.1.value.set, MARGIN=2, MakeKmeansAnnotationSetDF2)
left.eff.ki.1.value.set <- left.eff.ki.1.value.set[c(left.eff.ki.1.annotation.set[, "Left.efficiency"] != "out"), ]
left.eff.ki.1.annotation.set <- left.eff.ki.1.annotation.set[c(left.eff.ki.1.annotation.set[, "Left.efficiency"] != "out"), ]
SaveLabeledGroupPlot(left.eff.ki.1.value.set, left.eff.ki.1.annotation.set, "Left.efficiency", figure.dir, "PITCh")

# right.eff.ki.1.value.set <- Removeoutliers(right.eff.ki.1.value.set, "Right.efficiency")
right.eff.ki.1.annotation.set <- apply(right.eff.ki.1.value.set, MARGIN=2, MakeKmeansAnnotationSetDF2)
right.eff.ki.1.value.set <- right.eff.ki.1.value.set[c(right.eff.ki.1.annotation.set[, "Right.efficiency"] != "out"), ]
right.eff.ki.1.annotation.set <- right.eff.ki.1.annotation.set[c(right.eff.ki.1.annotation.set[, "Right.efficiency"] != "out"), ]
SaveLabeledGroupPlot(right.eff.ki.1.value.set, right.eff.ki.1.annotation.set, "Right.efficiency", figure.dir, "PITCh")

# left.acc.ki.1.value.set <- Removeoutliers(left.acc.ki.1.value.set, "Left.accuracy")
left.acc.ki.1.annotation.set <- apply(left.acc.ki.1.value.set, MARGIN=2, MakeKmeansAnnotationSetDF2)
left.acc.ki.1.value.set <- left.acc.ki.1.value.set[c(left.acc.ki.1.annotation.set[, "Left.accuracy"] != "out"), ]
left.acc.ki.1.annotation.set <- left.acc.ki.1.annotation.set[c(left.acc.ki.1.annotation.set[, "Left.accuracy"] != "out"), ]
SaveLabeledGroupPlot(left.acc.ki.1.value.set, left.acc.ki.1.annotation.set, "Left.accuracy", figure.dir, "PITCh")

# right.acc.ki.1.value.set <- Removeoutliers(right.acc.ki.1.value.set, "Right.accuracy")
right.acc.ki.1.annotation.set <- apply(right.acc.ki.1.value.set, MARGIN=2, MakeKmeansAnnotationSetDF2)
right.acc.ki.1.value.set <- right.acc.ki.1.value.set[c(right.acc.ki.1.annotation.set[, "Right.accuracy"] != "out"), ]
right.acc.ki.1.annotation.set <- right.acc.ki.1.annotation.set[c(right.acc.ki.1.annotation.set[, "Right.accuracy"] != "out"), ]
SaveLabeledGroupPlot(right.acc.ki.1.value.set, right.acc.ki.1.annotation.set, "Right.accuracy", figure.dir, "PITCh")

# left.eff.ki.2.value.set <- Removeoutliers(left.eff.ki.2.value.set, "Left.efficiency")
left.eff.ki.2.annotation.set <- apply(left.eff.ki.2.value.set, MARGIN=2, MakeKmeansAnnotationSetDF2)
left.eff.ki.2.value.set <- left.eff.ki.2.value.set[c(left.eff.ki.2.annotation.set[, "Left.efficiency"] != "out"), ]
left.eff.ki.2.annotation.set <- left.eff.ki.2.annotation.set[c(left.eff.ki.2.annotation.set[, "Left.efficiency"] != "out"), ]
SaveLabeledGroupPlot(left.eff.ki.2.value.set , left.eff.ki.2.annotation.set, "Left.efficiency", figure.dir, "LoAD")

# right.eff.ki.2.value.set <- Removeoutliers(right.eff.ki.2.value.set, "Right.efficiency")
right.eff.ki.2.annotation.set <- apply(right.eff.ki.2.value.set, MARGIN=2, MakeKmeansAnnotationSetDF2)
right.eff.ki.2.value.set <- right.eff.ki.2.value.set[c(right.eff.ki.2.annotation.set[, "Right.efficiency"] != "out"), ]
right.eff.ki.2.annotation.set <- right.eff.ki.2.annotation.set[c(right.eff.ki.2.annotation.set[, "Right.efficiency"] != "out"), ]
SaveLabeledGroupPlot(right.eff.ki.2.value.set , right.eff.ki.2.annotation.set, "Right.efficiency", figure.dir, "LoAD")

# left.acc.ki.2.value.set <- Removeoutliers(left.acc.ki.2.value.set, "Left.accuracy")
left.acc.ki.2.annotation.set <- apply(left.acc.ki.2.value.set, MARGIN=2, MakeKmeansAnnotationSetDF2)
left.acc.ki.2.value.set <- left.acc.ki.2.value.set[c(left.acc.ki.2.annotation.set[, "Left.accuracy"] != "out"), ]
left.acc.ki.2.annotation.set <- left.acc.ki.2.annotation.set[c(left.acc.ki.2.annotation.set[, "Left.accuracy"] != "out"), ]
SaveLabeledGroupPlot(left.acc.ki.2.value.set , left.acc.ki.2.annotation.set, "Left.accuracy", figure.dir, "LoAD")

# right.acc.ki.2.value.set <- Removeoutliers(right.acc.ki.2.value.set, "Right.accuracy")
right.acc.ki.2.annotation.set <- apply(right.acc.ki.2.value.set, MARGIN=2, MakeKmeansAnnotationSetDF2)
right.acc.ki.2.value.set <- right.acc.ki.2.value.set[c(right.acc.ki.2.annotation.set[, "Right.accuracy"] != "out"), ]
right.acc.ki.2.annotation.set <- right.acc.ki.2.annotation.set[c(right.acc.ki.2.annotation.set[, "Right.accuracy"] != "out"), ]
SaveLabeledGroupPlot(right.acc.ki.2.value.set , right.acc.ki.2.annotation.set, "Right.accuracy", figure.dir, "LoAD")

# whole.eff.ki.1.value.set <- Removeoutliers(whole.eff.ki.1.value.set, "Precise.Knock.in.Efficiency")
whole.eff.ki.1.annotation.set <- apply(whole.eff.ki.1.value.set, MARGIN=2, MakeKmeansAnnotationSetDF2)
whole.eff.ki.1.value.set <- whole.eff.ki.1.value.set[c(whole.eff.ki.1.annotation.set[, "Precise.Knock.in.Efficiency"] != "out"), ]
whole.eff.ki.1.annotation.set <- whole.eff.ki.1.annotation.set[c(whole.eff.ki.1.annotation.set[, "Precise.Knock.in.Efficiency"] != "out"), ]
SaveLabeledGroupPlot(whole.eff.ki.1.value.set , whole.eff.ki.1.annotation.set, "Precise.Knock.in.Efficiency", figure.dir, "PITCh")

# whole.acc.ki.1.value.set <- Removeoutliers(whole.acc.ki.1.value.set, "Accuracy")
whole.acc.ki.1.annotation.set <- apply(whole.acc.ki.1.value.set, MARGIN=2, MakeKmeansAnnotationSetDF2)
whole.acc.ki.1.value.set <- whole.acc.ki.1.value.set[c(whole.acc.ki.1.annotation.set[, "Accuracy"] != "out"), ]
whole.acc.ki.1.annotation.set <- whole.acc.ki.1.annotation.set[c(whole.acc.ki.1.annotation.set[, "Accuracy"] != "out"), ]
SaveLabeledGroupPlot(whole.acc.ki.1.value.set , whole.acc.ki.1.annotation.set, "Accuracy", figure.dir, "PITCh")

# whole.eff.ki.2.value.set <- Removeoutliers(whole.eff.ki.2.value.set, "Precise.Knock.in.Efficiency")
whole.eff.ki.2.annotation.set <- apply(whole.eff.ki.2.value.set, MARGIN=2, MakeKmeansAnnotationSetDF2)
whole.eff.ki.2.value.set <- whole.eff.ki.2.value.set[c(whole.eff.ki.2.annotation.set[, "Precise.Knock.in.Efficiency"] != "out"), ]
whole.eff.ki.2.annotation.set <- whole.eff.ki.2.annotation.set[c(whole.eff.ki.2.annotation.set[, "Precise.Knock.in.Efficiency"] != "out"), ]
SaveLabeledGroupPlot(whole.eff.ki.2.value.set , whole.eff.ki.2.annotation.set, "Precise.Knock.in.Efficiency", figure.dir, "LoAD")

# whole.acc.ki.2.value.set <- Removeoutliers(whole.acc.ki.2.value.set, "Accuracy")
whole.acc.ki.2.annotation.set <- apply(whole.acc.ki.2.value.set, MARGIN=2, MakeKmeansAnnotationSetDF2)
whole.acc.ki.2.value.set <- whole.acc.ki.2.value.set[c(whole.acc.ki.2.annotation.set[, "Accuracy"] != "out"), ]
whole.acc.ki.2.annotation.set <- whole.acc.ki.2.annotation.set[c(whole.acc.ki.2.annotation.set[, "Accuracy"] != "out"), ]
SaveLabeledGroupPlot(whole.acc.ki.2.value.set , whole.acc.ki.2.annotation.set, "Accuracy", figure.dir, "LoAD")

# whole.eff.enh.value.set <- Removeoutliers(whole.eff.enh.value.set, "Efficiency.Enhancement")
whole.eff.enh.annotation.set <- apply(whole.eff.enh.value.set, MARGIN=2, MakeKmeansAnnotationSetDF2)
whole.eff.enh.value.set <- whole.eff.enh.value.set[c(whole.eff.enh.annotation.set[, "Efficiency.Enhancement"] != "out"), ]
whole.eff.enh.annotation.set <- whole.eff.enh.annotation.set[c(whole.eff.enh.annotation.set[, "Efficiency.Enhancement"] != "out"), ]
SaveLabeledGroupPlot(whole.eff.enh.value.set , whole.eff.enh.annotation.set, "Efficiency.Enhancement", figure.dir, "Enhancement")

# whole.acc.enh.value.set <- Removeoutliers(whole.acc.enh.value.set, "Accuracy.Enhancement")
whole.acc.enh.annotation.set <- apply(whole.acc.enh.value.set, MARGIN=2, MakeKmeansAnnotationSetDF2)
whole.acc.enh.value.set <- whole.acc.enh.value.set[c(whole.acc.enh.annotation.set[, "Accuracy.Enhancement"] != "out"), ]
whole.acc.enh.annotation.set <- whole.acc.enh.annotation.set[c(whole.acc.enh.annotation.set[, "Accuracy.Enhancement"] != "out"), ]
SaveLabeledGroupPlot(whole.acc.enh.value.set , whole.acc.enh.annotation.set, "Accuracy.Enhancement", figure.dir, "Enhancement")
options(warn = 0)

# make filtered extra table
# browser()
fltr.extra.table <- merge(cas9.1.value.set[,c("Editing.Efficiency", "lmh", "rmh", "bmh", "group.name")], extra.table, by = "group.name", all=TRUE)
fltr.extra.table <- fltr.extra.table[which(!(fltr.extra.table[, "bmh"] %in% NA)) ,] # remove row including NA because these column is filtered above.
# MHScore
mp.score.value.table <- MakeMHScoreTable(fltr.extra.table, script.dir)
fltr.extra.table <- merge(fltr.extra.table, mp.score.value.table, by = "bmh", all=TRUE)
# browser()

scaffold.seq <- left.eff.ki.1.value.set[, "scaffold"]
# scaffold.seq <- "gttttagagctaggccaacatgaggatcacccatgtctgcagggcctagcaagttaaaataaggctagtccgttatcaacttggccaacatgaggatcacccatgtctgcagggccaagtggcaccgagtcggtgc"
# Calc microhomology proparty
# on efficiency
# browser()
pitch.left.microhomology.40.40.CalcMicroHomologyPropPK <- CalcMHProp(
    microhomology.seqset = left.eff.ki.1.value.set[, "lmh"],
    focus.length = 40,
    scaffoldset = paste(left.eff.ki.1.value.set[, "protospacer"], left.eff.ki.1.value.set[, "scaffold"], sep=""),
    annotationset = left.eff.ki.1.annotation.set[, "Left.efficiency"],
    is.left.microhomology.logical = TRUE,
    out.align.file = "pitch.left.microhomology.40.40.sence.sence.align.rds",
    out.prop.file = "pitch.left.microhomology.40.40.CalcMicroHomologyPropPK.rds"
)

pitch.right.microhomology.40.40.CalcMicroHomologyPropPK <- CalcMHProp(
    microhomology.seqset = right.eff.ki.1.value.set[, "rmh"],
    focus.length = 40,
    scaffoldset = paste(right.eff.ki.1.value.set[, "protospacer"], right.eff.ki.1.value.set[, "scaffold"], sep=""),
    annotationset = right.eff.ki.1.annotation.set[, "Right.efficiency"],
    is.left.microhomology.logical = FALSE,
    out.align.file = "pitch.right.microhomology.40.40.sence.sence.align.rds",
    out.prop.file = "pitch.right.microhomology.40.40.CalcMicroHomologyPropPK.rds"
)

load.left.microhomology.40.40.CalcMicroHomologyPropPK <- CalcMHProp(
    microhomology.seqset = left.eff.ki.2.value.set[, "lmh"],
    focus.length = 40,
    scaffoldset = paste(left.eff.ki.2.value.set[, "protospacer"], left.eff.ki.2.value.set[, "scaffold"], sep=""),
    annotationset = left.eff.ki.2.annotation.set[, "Left.efficiency"],
    is.left.microhomology.logical = TRUE,
    out.align.file = "load.left.microhomology.40.40.sence.sence.align.rds",
    out.prop.file = "load.left.microhomology.40.40.CalcMicroHomologyPropPK.rds"
)

load.right.microhomology.40.40.CalcMicroHomologyPropPK <- CalcMHProp(
    microhomology.seqset = right.eff.ki.2.value.set[, "rmh"],
    focus.length = 40,
    scaffoldset = paste(right.eff.ki.2.value.set[, "protospacer"], right.eff.ki.2.value.set[, "scaffold"], sep=""),
    annotationset = right.eff.ki.2.annotation.set[, "Right.efficiency"],
    is.left.microhomology.logical = FALSE,
    out.align.file = "load.right.microhomology.40.40.sence.sence.align.rds",
    out.prop.file = "load.right.microhomology.40.40.CalcMicroHomologyPropPK.rds"
)

# on accuracy

pitch.acc.left.microhomology.40.40.CalcMicroHomologyPropPK <- CalcMHProp(
    microhomology.seqset = left.acc.ki.1.value.set[, "lmh"],
    focus.length = 40,
    scaffoldset = paste(left.acc.ki.1.value.set[, "protospacer"], left.acc.ki.1.value.set[, "scaffold"], sep=""),
    annotationset = left.acc.ki.1.annotation.set[, "Left.accuracy"],
    is.left.microhomology.logical = TRUE,
    out.align.file = "pitch.acc.left.microhomology.40.40.sence.sence.align.rds",
    out.prop.file = "pitch.acc.left.microhomology.40.40.CalcMicroHomologyPropPK.rds"
)

pitch.acc.right.microhomology.40.40.CalcMicroHomologyPropPK <- CalcMHProp(
    microhomology.seqset = right.acc.ki.1.value.set[, "rmh"],
    focus.length = 40,
    scaffoldset = paste(right.acc.ki.1.value.set[, "protospacer"], right.acc.ki.1.value.set[, "scaffold"], sep=""),
    annotationset = right.acc.ki.1.annotation.set[, "Right.accuracy"],
    is.left.microhomology.logical = FALSE,
    out.align.file = "pitch.acc.right.microhomology.40.40.sence.sence.align.rds",
    out.prop.file = "pitch.acc.right.microhomology.40.40.CalcMicroHomologyPropPK.rds"
)

load.acc.left.microhomology.40.40.CalcMicroHomologyPropPK <- CalcMHProp(
    microhomology.seqset = left.acc.ki.2.value.set[, "lmh"],
    focus.length = 40,
    scaffoldset = paste(left.acc.ki.2.value.set[, "protospacer"], left.acc.ki.2.value.set[, "scaffold"], sep=""),
    annotationset = left.acc.ki.2.annotation.set[, "Left.accuracy"],
    is.left.microhomology.logical = TRUE,
    out.align.file = "load.acc.left.microhomology.40.40.sence.sence.align.rds",
    out.prop.file = "load.acc.left.microhomology.40.40.CalcMicroHomologyPropPK.rds"
)

load.acc.right.microhomology.40.40.CalcMicroHomologyPropPK <- CalcMHProp(
    microhomology.seqset = right.acc.ki.2.value.set[, "rmh"],
    focus.length = 40,
    scaffoldset = paste(right.acc.ki.2.value.set[, "protospacer"], right.acc.ki.2.value.set[, "scaffold"], sep=""),
    annotationset = right.acc.ki.2.annotation.set[, "Right.accuracy"],
    is.left.microhomology.logical = FALSE,
    out.align.file = "load.acc.right.microhomology.40.40.sence.sence.align.rds",
    out.prop.file = "load.acc.right.microhomology.40.40.CalcMicroHomologyPropPK.rds"
)

pitch.eff.whole.microhomology.40.40.CalcMicroHomologyPropPK <- CalcMHProp(
    microhomology.seqset = whole.eff.ki.1.value.set[, "bmh"],
    focus.length = 80,
    scaffoldset = paste(whole.eff.ki.1.value.set[, "protospacer"], whole.eff.ki.1.value.set[, "scaffold"], sep=""),
    annotationset = whole.eff.ki.1.annotation.set[, "Precise.Knock.in.Efficiency"],
    is.left.microhomology.logical = TRUE, # No problem
    out.align.file = "pitch.eff.whole.microhomology.40.40.sence.sence.align.rds",
    out.prop.file = "pitch.eff.whole.microhomology.40.40.CalcMicroHomologyPropPK.rds"
)
pitch.acc.whole.microhomology.40.40.CalcMicroHomologyPropPK <- CalcMHProp(
    microhomology.seqset = whole.acc.ki.1.value.set[, "bmh"],
    focus.length = 80,
    scaffoldset = paste(whole.acc.ki.1.value.set[, "protospacer"], whole.acc.ki.1.value.set[, "scaffold"], sep=""),
    annotationset = whole.acc.ki.1.annotation.set[, "Accuracy"],
    is.left.microhomology.logical = TRUE, # No problem
    out.align.file = "pitch.acc.whole.microhomology.40.40.sence.sence.align.rds",
    out.prop.file = "pitch.acc.whole.microhomology.40.40.CalcMicroHomologyPropPK.rds"
)
load.eff.whole.microhomology.40.40.CalcMicroHomologyPropPK <- CalcMHProp(
    microhomology.seqset = whole.eff.ki.2.value.set[, "bmh"],
    focus.length = 80,
    scaffoldset = paste(whole.eff.ki.2.value.set[, "protospacer"], whole.eff.ki.2.value.set[, "scaffold"], sep=""),
    annotationset = whole.eff.ki.2.annotation.set[, "Precise.Knock.in.Efficiency"],
    is.left.microhomology.logical = TRUE, # No problem
    out.align.file = "load.eff.whole.microhomology.40.40.sence.sence.align.rds",
    out.prop.file = "load.eff.whole.microhomology.40.40.CalcMicroHomologyPropPK.rds"
)
load.acc.whole.microhomology.40.40.CalcMicroHomologyPropPK <- CalcMHProp(
    microhomology.seqset = whole.acc.ki.2.value.set[, "bmh"],
    focus.length = 80,
    scaffoldset = paste(whole.acc.ki.2.value.set[, "protospacer"], whole.acc.ki.2.value.set[, "scaffold"], sep=""),
    annotationset = whole.acc.ki.2.annotation.set[, "Accuracy"],
    is.left.microhomology.logical = TRUE, # No problem
    out.align.file = "load.acc.whole.microhomology.40.40.sence.sence.align.rds",
    out.prop.file = "load.acc.whole.microhomology.40.40.CalcMicroHomologyPropPK.rds"
)

whole.eff.enh.microhomology.40.40.CalcMicroHomologyPropPK <- CalcMHProp(
    microhomology.seqset = whole.eff.enh.value.set[, "bmh"],
    focus.length = 80,
    scaffoldset = paste(whole.eff.enh.value.set[, "protospacer"], whole.eff.enh.value.set[, "scaffold"], sep=""),
    annotationset = whole.eff.enh.annotation.set[, "Efficiency.Enhancement"],
    is.left.microhomology.logical = TRUE, # No problem
    out.align.file = "whole.eff.enh.microhomology.40.40.sence.sence.align.rds",
    out.prop.file = "whole.eff.enh.microhomology.40.40.CalcMicroHomologyPropPK.rds"
)
whole.acc.enh.microhomology.40.40.CalcMicroHomologyPropPK <- CalcMHProp(
    microhomology.seqset = whole.acc.enh.value.set[, "bmh"],
    focus.length = 80,
    scaffoldset = paste(whole.acc.enh.value.set[, "protospacer"], whole.acc.enh.value.set[, "scaffold"], sep=""),
    annotationset = whole.acc.enh.annotation.set[, "Accuracy.Enhancement"],
    is.left.microhomology.logical = TRUE, # No problem
    out.align.file = "whole.acc.enh.microhomology.40.40.sence.sence.align.rds",
    out.prop.file = "whole.acc.enh.microhomology.40.40.CalcMicroHomologyPropPK.rds"
)

####################################################################

# add extra data set to prop data
left.value.names.vec <- c("Editing.Efficiency", "Chromatin_size", "ExonNumber", "Position.aa.", "protein_length.aa.", "Position.relative.", "CNV_HMM_293T", "Expression"
    , "Left_accessibility_20bp", "Left_accessibility_40bp", "Left_accessibility_100bp", "Left_accessibility_500bp", "Left_accessibility_1kbp"
    , "Endogenous.MP.Score", "Endogenous.MPOF.Score", "LeftRevKI.MP.Score", "LeftRevKI.MPOF.Score", "LeftKI.MP.Score", "LeftKI.MPOF.Score"
    , "Endogenous.InDelphi.MS.Score", "Endogenous.InDelphi.MS.Percentile", "Endogenous.InDelphi.F.Frequency", "Endogenous.InDelphi.F.Percentile")
right.value.names.vec <- c("Editing.Efficiency", "Chromatin_size", "ExonNumber", "Position.aa.", "protein_length.aa.", "Position.relative.", "CNV_HMM_293T", "Expression"
    , "Right_accessibility_20bp", "Right_accessibility_40bp", "Right_accessibility_100bp", "Right_accessibility_500bp", "Right_accessibility_1kbp"
    , "Endogenous.MP.Score", "Endogenous.MPOF.Score", "RightRevKI.MP.Score", "RightRevKI.MPOF.Score", "RightKI.MP.Score", "RightKI.MPOF.Score"
    , "Endogenous.InDelphi.MS.Score", "Endogenous.InDelphi.MS.Percentile", "Endogenous.InDelphi.F.Frequency", "Endogenous.InDelphi.F.Percentile")
#####
whole.eff.ki.1.fltr.extra.table <- cbind(
    fltr.extra.table
    , Accessibility_20bp = rowSums(fltr.extra.table[, c("Left_accessibility_20bp", "Right_accessibility_20bp")])
    , Accessibility_40bp = rowSums(fltr.extra.table[, c("Left_accessibility_40bp", "Right_accessibility_40bp")])
    , Accessibility_100bp = rowSums(fltr.extra.table[, c("Left_accessibility_100bp", "Right_accessibility_100bp")])
    , Accessibility_500bp = rowSums(fltr.extra.table[, c("Left_accessibility_500bp", "Right_accessibility_500bp")])
    , Accessibility_1kbp = rowSums(fltr.extra.table[, c("Left_accessibility_1kbp", "Right_accessibility_1kbp")])
    , sep="")

whole.acc.ki.1.fltr.extra.table <- cbind(
    fltr.extra.table
    , Accessibility_20bp = rowSums(fltr.extra.table[, c("Left_accessibility_20bp", "Right_accessibility_20bp")])
    , Accessibility_40bp = rowSums(fltr.extra.table[, c("Left_accessibility_40bp", "Right_accessibility_40bp")])
    , Accessibility_100bp = rowSums(fltr.extra.table[, c("Left_accessibility_100bp", "Right_accessibility_100bp")])
    , Accessibility_500bp = rowSums(fltr.extra.table[, c("Left_accessibility_500bp", "Right_accessibility_500bp")])
    , Accessibility_1kbp = rowSums(fltr.extra.table[, c("Left_accessibility_1kbp", "Right_accessibility_1kbp")])
    , sep="")

whole.eff.ki.2.fltr.extra.table <- cbind(
    fltr.extra.table
    , Accessibility_20bp = rowSums(fltr.extra.table[, c("Left_accessibility_20bp", "Right_accessibility_20bp")])
    , Accessibility_40bp = rowSums(fltr.extra.table[, c("Left_accessibility_40bp", "Right_accessibility_40bp")])
    , Accessibility_100bp = rowSums(fltr.extra.table[, c("Left_accessibility_100bp", "Right_accessibility_100bp")])
    , Accessibility_500bp = rowSums(fltr.extra.table[, c("Left_accessibility_500bp", "Right_accessibility_500bp")])
    , Accessibility_1kbp = rowSums(fltr.extra.table[, c("Left_accessibility_1kbp", "Right_accessibility_1kbp")])
    , sep="")

whole.acc.ki.2.fltr.extra.table <- cbind(
    fltr.extra.table
    , Accessibility_20bp = rowSums(fltr.extra.table[, c("Left_accessibility_20bp", "Right_accessibility_20bp")])
    , Accessibility_40bp = rowSums(fltr.extra.table[, c("Left_accessibility_40bp", "Right_accessibility_40bp")])
    , Accessibility_100bp = rowSums(fltr.extra.table[, c("Left_accessibility_100bp", "Right_accessibility_100bp")])
    , Accessibility_500bp = rowSums(fltr.extra.table[, c("Left_accessibility_500bp", "Right_accessibility_500bp")])
    , Accessibility_1kbp = rowSums(fltr.extra.table[, c("Left_accessibility_1kbp", "Right_accessibility_1kbp")])
    , sep="")

whole.eff.enh.fltr.extra.table <- cbind(
    fltr.extra.table
    , Accessibility_20bp = rowSums(fltr.extra.table[, c("Left_accessibility_20bp", "Right_accessibility_20bp")])
    , Accessibility_40bp = rowSums(fltr.extra.table[, c("Left_accessibility_40bp", "Right_accessibility_40bp")])
    , Accessibility_100bp = rowSums(fltr.extra.table[, c("Left_accessibility_100bp", "Right_accessibility_100bp")])
    , Accessibility_500bp = rowSums(fltr.extra.table[, c("Left_accessibility_500bp", "Right_accessibility_500bp")])
    , Accessibility_1kbp = rowSums(fltr.extra.table[, c("Left_accessibility_1kbp", "Right_accessibility_1kbp")])
    , sep="")

whole.acc.enh.fltr.extra.table <- cbind(
    fltr.extra.table
    , Accessibility_20bp = rowSums(fltr.extra.table[, c("Left_accessibility_20bp", "Right_accessibility_20bp")])
    , Accessibility_40bp = rowSums(fltr.extra.table[, c("Left_accessibility_40bp", "Right_accessibility_40bp")])
    , Accessibility_100bp = rowSums(fltr.extra.table[, c("Left_accessibility_100bp", "Right_accessibility_100bp")])
    , Accessibility_500bp = rowSums(fltr.extra.table[, c("Left_accessibility_500bp", "Right_accessibility_500bp")])
    , Accessibility_1kbp = rowSums(fltr.extra.table[, c("Left_accessibility_1kbp", "Right_accessibility_1kbp")])
    , sep="")
#####
whole.value.names.vec <- c("Editing.Efficiency", "Chromatin_size", "ExonNumber", "Position.aa.", "protein_length.aa.", "Position.relative.", "CNV_HMM_293T", "Expression"
    , "Right_accessibility_20bp", "Right_accessibility_40bp", "Right_accessibility_100bp", "Right_accessibility_500bp", "Right_accessibility_1kbp"
    , "Left_accessibility_20bp", "Left_accessibility_40bp", "Left_accessibility_100bp", "Left_accessibility_500bp", "Left_accessibility_1kbp"
    , "Accessibility_20bp", "Accessibility_40bp", "Accessibility_100bp", "Accessibility_500bp", "Accessibility_1kbp"
    , "Endogenous.MP.Score", "Endogenous.MPOF.Score", "LeftRevKI.MP.Score", "LeftRevKI.MPOF.Score", "RightRevKI.MP.Score", "RightRevKI.MPOF.Score", "LeftKI.MP.Score"
    , "LeftKI.MPOF.Score", "RightKI.MP.Score", "RightKI.MPOF.Score"
    , "Endogenous.InDelphi.MS.Score", "Endogenous.InDelphi.MS.Percentile", "Endogenous.InDelphi.F.Frequency", "Endogenous.InDelphi.F.Percentile"
)
all.pitch.eff.left.future.data.df <- MergeExtraData(pitch.left.microhomology.40.40.CalcMicroHomologyPropPK$ShowDataValue("all"), fltr.extra.table, "lmh", left.value.names.vec)
all.pitch.eff.right.future.data.df <- MergeExtraData(pitch.right.microhomology.40.40.CalcMicroHomologyPropPK$ShowDataValue("all"), fltr.extra.table, "rmh", right.value.names.vec)
all.load.eff.left.future.data.df <- MergeExtraData(load.left.microhomology.40.40.CalcMicroHomologyPropPK$ShowDataValue("all"), fltr.extra.table, "lmh", left.value.names.vec)
all.load.eff.right.future.data.df <- MergeExtraData(load.right.microhomology.40.40.CalcMicroHomologyPropPK$ShowDataValue("all"), fltr.extra.table, "rmh", right.value.names.vec)
all.pitch.acc.left.future.data.df <- MergeExtraData(pitch.acc.left.microhomology.40.40.CalcMicroHomologyPropPK$ShowDataValue("all"), fltr.extra.table, "lmh", left.value.names.vec)
all.pitch.acc.right.future.data.df <- MergeExtraData(pitch.acc.right.microhomology.40.40.CalcMicroHomologyPropPK$ShowDataValue("all"), fltr.extra.table, "rmh", right.value.names.vec)
all.load.acc.left.future.data.df <- MergeExtraData(load.acc.left.microhomology.40.40.CalcMicroHomologyPropPK$ShowDataValue("all"), fltr.extra.table, "lmh", left.value.names.vec)
all.load.acc.right.future.data.df <- MergeExtraData(load.acc.right.microhomology.40.40.CalcMicroHomologyPropPK$ShowDataValue("all"), fltr.extra.table, "rmh", right.value.names.vec)

all.pitch.eff.whole.future.data.df <- MergeExtraData(pitch.eff.whole.microhomology.40.40.CalcMicroHomologyPropPK$ShowDataValue("all"), whole.eff.ki.1.fltr.extra.table, "bmh", whole.value.names.vec)
all.load.eff.whole.future.data.df <- MergeExtraData(pitch.acc.whole.microhomology.40.40.CalcMicroHomologyPropPK$ShowDataValue("all"), whole.acc.ki.1.fltr.extra.table, "bmh", whole.value.names.vec)
all.pitch.acc.whole.future.data.df <- MergeExtraData(load.eff.whole.microhomology.40.40.CalcMicroHomologyPropPK$ShowDataValue("all"), whole.eff.ki.2.fltr.extra.table, "bmh", whole.value.names.vec)
all.load.acc.whole.future.data.df <- MergeExtraData(load.acc.whole.microhomology.40.40.CalcMicroHomologyPropPK$ShowDataValue("all"), whole.acc.ki.2.fltr.extra.table, "bmh", whole.value.names.vec)

# continuous future only
continuous.pitch.eff.left.future.data.df <- MergeExtraData(pitch.left.microhomology.40.40.CalcMicroHomologyPropPK$ShowDataValue("continuous"), fltr.extra.table, "lmh", left.value.names.vec)
continuous.pitch.eff.right.future.data.df <- MergeExtraData(pitch.right.microhomology.40.40.CalcMicroHomologyPropPK$ShowDataValue("continuous"), fltr.extra.table, "rmh", right.value.names.vec)
continuous.load.eff.left.future.data.df <- MergeExtraData(load.left.microhomology.40.40.CalcMicroHomologyPropPK$ShowDataValue("continuous"), fltr.extra.table, "lmh", left.value.names.vec)
continuous.load.eff.right.future.data.df <- MergeExtraData(load.right.microhomology.40.40.CalcMicroHomologyPropPK$ShowDataValue("continuous"), fltr.extra.table, "rmh", right.value.names.vec)
continuous.pitch.acc.left.future.data.df <- MergeExtraData(pitch.acc.left.microhomology.40.40.CalcMicroHomologyPropPK$ShowDataValue("continuous"), fltr.extra.table, "lmh", left.value.names.vec)
continuous.pitch.acc.right.future.data.df <- MergeExtraData(pitch.acc.right.microhomology.40.40.CalcMicroHomologyPropPK$ShowDataValue("continuous"), fltr.extra.table, "rmh", right.value.names.vec)
continuous.load.acc.left.future.data.df <- MergeExtraData(load.acc.left.microhomology.40.40.CalcMicroHomologyPropPK$ShowDataValue("continuous"), fltr.extra.table, "lmh", left.value.names.vec)
continuous.load.acc.right.future.data.df <- MergeExtraData(load.acc.right.microhomology.40.40.CalcMicroHomologyPropPK$ShowDataValue("continuous"), fltr.extra.table, "rmh", right.value.names.vec)

continuous.pitch.eff.whole.future.data.df <- MergeExtraData(pitch.eff.whole.microhomology.40.40.CalcMicroHomologyPropPK$ShowDataValue("continuous"), whole.eff.ki.1.fltr.extra.table, "bmh", whole.value.names.vec)
continuous.load.eff.whole.future.data.df <- MergeExtraData(pitch.acc.whole.microhomology.40.40.CalcMicroHomologyPropPK$ShowDataValue("continuous"), whole.acc.ki.1.fltr.extra.table, "bmh", whole.value.names.vec)
continuous.pitch.acc.whole.future.data.df <- MergeExtraData(load.eff.whole.microhomology.40.40.CalcMicroHomologyPropPK$ShowDataValue("continuous"), whole.eff.ki.2.fltr.extra.table, "bmh", whole.value.names.vec)
continuous.load.acc.whole.future.data.df <- MergeExtraData(load.acc.whole.microhomology.40.40.CalcMicroHomologyPropPK$ShowDataValue("continuous"), whole.acc.ki.2.fltr.extra.table, "bmh", whole.value.names.vec)

continuous.whole.eff.enh.future.data.df <- MergeExtraData(whole.eff.enh.microhomology.40.40.CalcMicroHomologyPropPK$ShowDataValue("continuous"), whole.eff.enh.fltr.extra.table, "bmh", whole.value.names.vec)
continuous.whole.acc.enh.future.data.df <- MergeExtraData(whole.acc.enh.microhomology.40.40.CalcMicroHomologyPropPK$ShowDataValue("continuous"), whole.acc.enh.fltr.extra.table, "bmh", whole.value.names.vec)



# browser()
# make future bool vector
all.left.future.colnames <- (colnames(all.pitch.eff.left.future.data.df) != "sequence") & (colnames(all.pitch.eff.left.future.data.df) != "surrounding.sequence")
all.right.future.colnames <- (colnames(all.pitch.eff.right.future.data.df) != "sequence") & (colnames(all.pitch.eff.right.future.data.df) != "surrounding.sequence")
all.whole.future.colnames <- (colnames(all.pitch.eff.whole.future.data.df) != "sequence") & (colnames(all.pitch.eff.whole.future.data.df) != "surrounding.sequence")
continuous.left.future.colnames <- (colnames(continuous.pitch.eff.left.future.data.df) != "sequence") & (colnames(continuous.pitch.eff.left.future.data.df) != "surrounding.sequence")
continuous.right.future.colnames <- (colnames(continuous.pitch.eff.right.future.data.df) != "sequence") & (colnames(continuous.pitch.eff.right.future.data.df) != "surrounding.sequence")
continuous.whole.future.colnames <- (colnames(continuous.pitch.eff.whole.future.data.df) != "sequence") & (colnames(continuous.pitch.eff.whole.future.data.df) != "surrounding.sequence")
# In this version is "max" and "min" parameter is removed because a diversity of these parameter is too low in the samples.
continuous.small.left.future.colnames <- (colnames(continuous.pitch.eff.left.future.data.df) != "sequence") &
    (colnames(continuous.pitch.eff.left.future.data.df) != "surrounding.sequence") &
    !as.logical(table(factor(c(grep("max", colnames(continuous.pitch.eff.left.future.data.df)), grep("min", colnames(continuous.pitch.eff.left.future.data.df)))
    , levels=1:length(colnames(continuous.pitch.eff.left.future.data.df)))))
continuous.small.right.future.colnames <- (colnames(continuous.pitch.eff.right.future.data.df) != "sequence") &
    (colnames(continuous.pitch.eff.right.future.data.df) != "surrounding.sequence") &
    !as.logical(table(factor(c(grep("max", colnames(continuous.pitch.eff.right.future.data.df)), grep("min", colnames(continuous.pitch.eff.right.future.data.df)))
    , levels=1:length(colnames(continuous.pitch.eff.right.future.data.df)))))
continuous.small.whole.future.colnames <- (colnames(continuous.pitch.eff.whole.future.data.df) != "sequence") &
    (colnames(continuous.pitch.eff.whole.future.data.df) != "surrounding.sequence") &
    !as.logical(table(factor(c(grep("max", colnames(continuous.pitch.eff.whole.future.data.df)), grep("min", colnames(continuous.pitch.eff.whole.future.data.df)))
    , levels=1:length(colnames(continuous.pitch.eff.whole.future.data.df)))))

pitch.eff.left.continuous.df <- data.frame(
    sequence = as.vector(left.eff.ki.1.value.set[,"lmh"]),
    Left.efficiency = as.vector(left.eff.ki.1.value.set[,"Left.efficiency"]),
    stringsAsFactors=FALSE
)
pitch.eff.right.continuous.df <- data.frame(
    sequence = as.vector(right.eff.ki.1.value.set[,"rmh"]),
    Right.efficiency = as.vector(right.eff.ki.1.value.set[,"Right.efficiency"]),
    stringsAsFactors=FALSE
)
pitch.eff.whole.continuous.df <- data.frame(
    sequence = as.vector(whole.eff.ki.1.value.set[,"bmh"]),
    Precise.Knock.in.Efficiency = as.vector(whole.eff.ki.1.value.set[,"Precise.Knock.in.Efficiency"]),
    stringsAsFactors=FALSE
)
load.eff.left.continuous.df <- data.frame(
    sequence = as.vector(left.eff.ki.2.value.set[,"lmh"]),
    Left.efficiency = as.vector(left.eff.ki.2.value.set[,"Left.efficiency"]),
    stringsAsFactors=FALSE
)
load.eff.right.continuous.df <- data.frame(
    sequence = as.vector(right.eff.ki.2.value.set[,"rmh"]),
    Right.efficiency = as.vector(right.eff.ki.2.value.set[,"Right.efficiency"]),
    stringsAsFactors=FALSE
)
load.eff.whole.continuous.df <- data.frame(
    sequence = as.vector(whole.eff.ki.2.value.set[,"bmh"]),
    Precise.Knock.in.Efficiency = as.vector(whole.eff.ki.2.value.set[,"Precise.Knock.in.Efficiency"]),
    stringsAsFactors=FALSE
)
pitch.acc.left.continuous.df <- data.frame(
    sequence = as.vector(left.acc.ki.1.value.set[,"lmh"]),
    Left.accuracy = as.vector(left.acc.ki.1.value.set[,"Left.accuracy"]),
    stringsAsFactors=FALSE
)
pitch.acc.right.continuous.df <- data.frame(
    sequence = as.vector(right.acc.ki.1.value.set[,"rmh"]),
    Right.accuracy = as.vector(right.acc.ki.1.value.set[,"Right.accuracy"]),
    stringsAsFactors=FALSE
)
pitch.acc.whole.continuous.df <- data.frame(
    sequence = as.vector(whole.acc.ki.1.value.set[,"bmh"]),
    Accuracy = as.vector(whole.acc.ki.1.value.set[,"Accuracy"]),
    stringsAsFactors=FALSE
)
load.acc.left.continuous.df <- data.frame(
    sequence = as.vector(left.acc.ki.2.value.set[,"lmh"]),
    Left.accuracy = as.vector(left.acc.ki.2.value.set[,"Left.accuracy"]),
    stringsAsFactors=FALSE
)
load.acc.right.continuous.df <- data.frame(
    sequence = as.vector(right.acc.ki.2.value.set[,"rmh"]),
    Right.accuracy = as.vector(right.acc.ki.2.value.set[,"Right.accuracy"]),
    stringsAsFactors=FALSE
)
load.acc.whole.continuous.df <- data.frame(
    sequence = as.vector(whole.acc.ki.2.value.set[,"bmh"]),
    Accuracy = as.vector(whole.acc.ki.2.value.set[,"Accuracy"]),
    stringsAsFactors=FALSE
)

whole.eff.enh.continuous.df <- data.frame(
    sequence = as.vector(whole.eff.enh.value.set[,"bmh"]),
    Accuracy = as.vector(whole.eff.enh.value.set[,"Efficiency.Enhancement"]),
    stringsAsFactors=FALSE
)
whole.acc.enh.continuous.df <- data.frame(
    sequence = as.vector(whole.acc.enh.value.set[,"bmh"]),
    Accuracy = as.vector(whole.acc.enh.value.set[,"Accuracy.Enhancement"]),
    stringsAsFactors=FALSE
)

# Make scaffold plot
PlotTopRankScattter(continuous.pitch.eff.left.future.data.df,
    pitch.eff.left.continuous.df,
    pitch.eff.left.microhomology.40.40.MeanDecreaseGini.vec$importance,
    length(pitch.eff.left.microhomology.40.40.MeanDecreaseGini.vec$importance),
    "Condition1"
)
PlotTopRankScattter(continuous.pitch.eff.right.future.data.df,
    pitch.eff.right.continuous.df,
    pitch.eff.right.microhomology.40.40.MeanDecreaseGini.vec$importance,
    length(pitch.eff.right.microhomology.40.40.MeanDecreaseGini.vec$importance),
    "Condition1"
)
PlotTopRankScattter(continuous.pitch.eff.whole.future.data.df,
    pitch.eff.whole.continuous.df,
    pitch.eff.whole.microhomology.40.40.MeanDecreaseGini.vec$importance,
    length(pitch.eff.whole.microhomology.40.40.MeanDecreaseGini.vec$importance),
    "Condition1"
)
PlotTopRankScattter(continuous.pitch.acc.left.future.data.df,
    pitch.acc.left.continuous.df,
    pitch.acc.left.microhomology.40.40.MeanDecreaseGini.vec$importance,
    length(pitch.acc.left.microhomology.40.40.MeanDecreaseGini.vec$importance),
    "Condition1"
)
PlotTopRankScattter(continuous.pitch.acc.right.future.data.df,
    pitch.acc.right.continuous.df,
    pitch.acc.right.microhomology.40.40.MeanDecreaseGini.vec$importance,
    length(pitch.acc.right.microhomology.40.40.MeanDecreaseGini.vec$importance),
    "Condition1"
)
PlotTopRankScattter(continuous.pitch.acc.whole.future.data.df,
    pitch.acc.whole.continuous.df,
    pitch.acc.whole.microhomology.40.40.MeanDecreaseGini.vec$importance,
    length(pitch.acc.whole.microhomology.40.40.MeanDecreaseGini.vec$importance),
    "Condition1"
)
PlotTopRankScattter(continuous.load.eff.right.future.data.df,
    load.eff.right.continuous.df,
    load.eff.right.microhomology.40.40.MeanDecreaseGini.vec$importance,
    length(load.eff.right.microhomology.40.40.MeanDecreaseGini.vec$importance),
    "Condition2"
)
PlotTopRankScattter(continuous.load.eff.left.future.data.df,
    load.eff.left.continuous.df,
    load.eff.left.microhomology.40.40.MeanDecreaseGini.vec$importance,
    length(load.eff.left.microhomology.40.40.MeanDecreaseGini.vec$importance),
    "Condition2"
)
PlotTopRankScattter(continuous.load.eff.whole.future.data.df,
    load.eff.whole.continuous.df,
    load.eff.whole.microhomology.40.40.MeanDecreaseGini.vec$importance,
    length(load.eff.whole.microhomology.40.40.MeanDecreaseGini.vec$importance),
    "Condition2"
)
PlotTopRankScattter(continuous.load.acc.right.future.data.df,
    load.acc.right.continuous.df,
    load.acc.right.microhomology.40.40.MeanDecreaseGini.vec$importance,
    length(load.acc.right.microhomology.40.40.MeanDecreaseGini.vec$importance),
    "Condition2"
)
PlotTopRankScattter(continuous.load.acc.left.future.data.df,
    load.acc.left.continuous.df,
    load.acc.left.microhomology.40.40.MeanDecreaseGini.vec$importance,
    length(load.acc.left.microhomology.40.40.MeanDecreaseGini.vec$importance),
    "Condition2"
)
PlotTopRankScattter(continuous.load.acc.whole.future.data.df,
    load.acc.whole.continuous.df,
    load.acc.whole.microhomology.40.40.MeanDecreaseGini.vec$importance,
    length(load.acc.whole.microhomology.40.40.MeanDecreaseGini.vec$importance),
    "Condition2"
)

PlotTopRankScattter(continuous.whole.eff.enh.future.data.df,
    whole.eff.enh.continuous.df,
    whole.eff.enh.microhomology.40.40.MeanDecreaseGini.vec$importance,
    length(whole.eff.enh.microhomology.40.40.MeanDecreaseGini.vec$importance),
    "Efficiency.Enhancement"
)
PlotTopRankScattter(continuous.whole.acc.enh.future.data.df,
    whole.acc.enh.continuous.df,
    whole.acc.enh.microhomology.40.40.MeanDecreaseGini.vec$importance,
    length(whole.acc.enh.microhomology.40.40.MeanDecreaseGini.vec$importance),
    "Accuracy.Enhancement"
)
#############################################
### make colleration box plot of each future set
# common
future.genome.prop.name <- c("Editing.Efficiency", "Chromatin_size", "ExonNumber", "protein_length.aa.", "Position.relative.", "CNV_HMM_293T", "Expression")
future.indelphi.name <- c("Endogenous.InDelphi.MS.Score", "Endogenous.InDelphi.MS.Percentile", "Endogenous.InDelphi.F.Frequency", "Endogenous.InDelphi.F.Percentile")
# special
future.left.accessibility.name <- c("Left_accessibility_20bp", "Left_accessibility_40bp", "Left_accessibility_100bp", "Left_accessibility_500bp", "Left_accessibility_1kbp")
future.right.accessibility.name <- c("Right_accessibility_20bp", "Right_accessibility_40bp", "Right_accessibility_100bp", "Right_accessibility_500bp", "Right_accessibility_1kbp")
future.whole.accessibility.name <- c("Accessibility_20bp", "Accessibility_40bp", "Accessibility_100bp", "Accessibility_500bp", "Accessibility_1kbp")
future.left.MP.name <- c("Endogenous.MP.Score", "Endogenous.MPOF.Score", "LeftRevKI.MP.Score", "LeftRevKI.MPOF.Score", "LeftKI.MP.Score", "LeftKI.MPOF.Score")
future.right.MP.name <- c("Endogenous.MP.Score", "Endogenous.MPOF.Score", "RightRevKI.MP.Score", "RightRevKI.MPOF.Score", "RightKI.MP.Score", "RightKI.MPOF.Score")
future.all.MP.name <- c("Endogenous.MP.Score", "Endogenous.MPOF.Score", "LeftRevKI.MP.Score", "LeftRevKI.MPOF.Score", "RightRevKI.MP.Score", "RightRevKI.MPOF.Score", "LeftKI.MP.Score", "LeftKI.MPOF.Score", "RightKI.MP.Score", "RightKI.MPOF.Score")

base.future.vec.list <- list(
    SgRNA.Align = pitch.left.microhomology.40.40.CalcMicroHomologyPropPK$colnames.align.vec
    , Freq = pitch.left.microhomology.40.40.CalcMicroHomologyPropPK$colnames.freq.vec
    , Packer = pitch.left.microhomology.40.40.CalcMicroHomologyPropPK$colnames.packer.vec
    , Phychem = pitch.left.microhomology.40.40.CalcMicroHomologyPropPK$colnames.phychem.vec
    , Pseknc = pitch.left.microhomology.40.40.CalcMicroHomologyPropPK$colnames.pseknc.vec
    , Thermo = pitch.left.microhomology.40.40.CalcMicroHomologyPropPK$colnames.thermo.vec
    , Genome.Proparty = future.genome.prop.name
    , InDelPhi = future.indelphi.name
)

left.future.vec.list <- c(base.future.vec.list
    , Left.Accessibility = list(future.left.accessibility.name)
    , Left.MP = list(future.left.MP.name))

right.future.vec.list <- c(base.future.vec.list
    , Right.Accessibility = list(future.right.accessibility.name)
    , Right.MP = list(future.right.MP.name))

whole.future.vec.list <- c(base.future.vec.list
    , Accessibility = list(c(future.left.accessibility.name, future.right.accessibility.name, future.whole.accessibility.name))
    , All.MP = list(future.all.MP.name))

left.future.vec.list <- c(left.future.vec.list , All = list(unlist(left.future.vec.list, use.names = FALSE)))
right.future.vec.list <- c(right.future.vec.list , All = list(unlist(right.future.vec.list, use.names = FALSE)))
whole.future.vec.list <- c(whole.future.vec.list , All = list(unlist(whole.future.vec.list, use.names = FALSE)))

print("Left.Efficiency")
for(future.ind in 1:length(left.future.vec.list)){
    MakeCorCompPlot(continuous.pitch.eff.left.future.data.df,
        continuous.load.eff.left.future.data.df,
        pitch.eff.left.continuous.df,
        load.eff.left.continuous.df,
        left.future.vec.list[[future.ind]],
        output.dir, figure.dir,
        names(left.future.vec.list)[future.ind], "Left.Efficiency", c("1.Condition1", "2.Condition2")
        )
}
print("Right.Efficiency")
for(future.ind in 1:length(right.future.vec.list)){
    MakeCorCompPlot(continuous.pitch.eff.right.future.data.df,
        continuous.load.eff.right.future.data.df,
        pitch.eff.right.continuous.df,
        load.eff.right.continuous.df,
        right.future.vec.list[[future.ind]],
        output.dir, figure.dir,
        names(right.future.vec.list)[future.ind], "Right.Efficiency", c("1.Condition1", "2.Condition2")
        )
}
print("Left.Accuracy")
for(future.ind in 1:length(left.future.vec.list)){
    MakeCorCompPlot(continuous.pitch.acc.left.future.data.df,
        continuous.load.acc.left.future.data.df,
        pitch.acc.left.continuous.df,
        load.acc.left.continuous.df,
        left.future.vec.list[[future.ind]],
        output.dir, figure.dir,
        names(left.future.vec.list)[future.ind], "Left.Accuracy", c("1.Condition1", "2.Condition2")
        )
}
print("Right.Accuracy")
for(future.ind in 1:length(right.future.vec.list)){
    MakeCorCompPlot(continuous.pitch.acc.right.future.data.df,
        continuous.load.acc.right.future.data.df,
        pitch.acc.right.continuous.df,
        load.acc.right.continuous.df,
        right.future.vec.list[[future.ind]],
        output.dir, figure.dir,
        names(right.future.vec.list)[future.ind], "Right.Accuracy", c("1.Condition1", "2.Condition2")
        )
}
print("Whole.Efficiency")
for(future.ind in 1:length(whole.future.vec.list)){
    MakeCorCompPlot(continuous.pitch.eff.whole.future.data.df,
        continuous.load.eff.whole.future.data.df,
        pitch.eff.whole.continuous.df,
        load.eff.whole.continuous.df,
        whole.future.vec.list[[future.ind]],
        output.dir, figure.dir,
        names(whole.future.vec.list)[future.ind], "Whole.Efficiency", c("1.Condition1", "2.Condition2")
        )
}
print("Whole.Accuracy")
for(future.ind in 1:length(whole.future.vec.list)){
    MakeCorCompPlot(continuous.pitch.acc.whole.future.data.df,
        continuous.load.acc.whole.future.data.df,
        pitch.acc.whole.continuous.df,
        load.acc.whole.continuous.df,
        whole.future.vec.list[[future.ind]],
        output.dir, figure.dir,
        names(whole.future.vec.list)[future.ind], "Whole.Accuracy", c("1.Condition1", "2.Condition2")
        )
}


message(paste("Ending process record in", output.dir, sep=" "))



