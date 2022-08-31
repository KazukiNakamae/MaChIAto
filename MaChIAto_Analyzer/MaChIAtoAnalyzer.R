debug.flag <- FALSE

message("
                                                                            
,-.-.     ,---.|    |,---.|         ,---.          |                        
| | |,---.|    |---.||---||--- ,---.|---|,---.,---.|    ,   .,---,,---.,---.
| | |,---||    |   |||   ||    |   ||   ||   |,---||    |   | .-' |---'|    
` ' '`---^`---'`   '``   '`---'`---'`   '`   '`---^`---'`---|'---'`---'`    
                                                        `---'               
version 1.1
")

help.message <- "

-------------------------------------------------------------------------------
This is not correct input. The description below would be helpful.
-------------------------------------------------------------------------------


-------------------------------------------------------------------------------
EXAMPLE INPUT
-------------------------------------------------------------------------------
Rscript MaChIAto_Analyzer/MaChIAtoAnalyzer.R <Summary directory> <Output prefix> [<Name of extra data> <Table of extra data> ...]

-------------------------------------------------------------------------------
DESPRIPTION
-------------------------------------------------------------------------------

Summary directory:
The summary directory generated with “collect_MaChIAto_data.py”.

Output prefix:
The directory into which the output directory is saved.

Name of extra data (optional):
The name of feature group that the next argument includes.

Table of extra data (optional):
The pathname of extra data added by the user. The data should be a .csv file.

-------------------------------------------------------------------------------

"

# Get directory path
if(debug.flag){
    script.dir <- "./MaChIAto_Analyzer/"
}else{
    initial.options <- commandArgs(trailingOnly = FALSE)
    file.arg.name <- "--file="
    script.fn <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
    script.dir <- dirname(script.fn)
}
ref.dir <- file.path(script.dir, "ref")
source(file.path(script.dir, "MaChIAtoAnalyzerPackageInstaller.R"))
source(file.path(script.dir, "MaChIAtoAnalyzerFunctions.R"))
source(file.path(script.dir, "ManageAnalysisData-class.R"))
source(file.path(script.dir, "CalcMHScore-class.R"))

# library(randomForest)
library(Biostrings)
library(effsize)
library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)
library(maptools)
library(scales)
# library(ggpubr)

ignore.list.fn <- character(0)
if(debug.flag){
    input.dir <- "/Volumes/databank4/yamamoto_nakamae/yamamoto_summary"
    output.prefix <- "/Volumes/databank4/yamamoto_nakamae"
    #extra.table.vec <- c(
    #    "InDelphi",
    #    "/Volumes/databank2/extra_table/bmh_extra_data_InDelphi.csv",
    #    "Accessibility",
    #    "/Volumes/databank2/extra_table/bmh_extra_data_accessibility.csv",
    #    "Genome_Property",
    #    "/Volumes/databank2/extra_table/bmh_extra_data_genomeprop.csv"
    #)
    extra.table.vec <- character(0)
    target.type <- as.vector(read.csv(file.path(input.dir, "target_type.csv"), header=FALSE , stringsAsFactors=FALSE))[[1]]
}else{
    arg.len <- length(commandArgs(trailingOnly=TRUE))
    if(arg.len < 2){
        message(help.message)
        q()
    }
    # Get commandline arguments
    input.dir <- GetAbsolutePath(commandArgs(trailingOnly=TRUE)[1])
    output.prefix <- GetAbsolutePath(commandArgs(trailingOnly=TRUE)[2])
    if(arg.len > 2){
        extra.table.vec <- commandArgs(trailingOnly=TRUE)[3:arg.len]
    }else{
        extra.table.vec <- character(0)
    }
    target.type <- as.vector(read.csv(file.path(input.dir, "target_type.csv"), header=FALSE , stringsAsFactors=FALSE))[[1]]
}

# Make output directory
time.id <- format(Sys.time(), "%Y%m%d%H%M%S")
output.name <- paste("MaChIAtoAnalyzer_at_", time.id, sep = "")
output.dir <- file.path(output.prefix, output.name)
message(paste("Creating Folder", output.dir, sep=" "))
system(paste("mkdir", output.dir, sep=" "))
message(paste("done"))

# Read data
message(paste0("Read MaChIAto data in ", input.dir))
all.count.table <- read.csv(file.path(input.dir, "sum_ALL_read_countdata.csv"), stringsAsFactors=FALSE)
all.count.table <- cbind(all.count.table,
    sample.target = as.character(sapply(all.count.table$name, function(x){strsplit(x, "-")[[1]]})[1,]),
    sample.label = as.character(sapply(all.count.table$name, function(x){strsplit(x, "-")[[1]]})[2,]),
    stringsAsFactors=FALSE)
all.rate.table <- read.csv(file.path(input.dir, "sum_ALL_read_rates.csv"), stringsAsFactors=FALSE)
label.table <- read.csv(file.path(input.dir, "label_sample_type.csv"), stringsAsFactors=FALSE)
extra.table.list <- list()
temp.extra.feature.name <- ""
if(length(extra.table.vec) > 0){
    for(extra.table.ind in 1:length(extra.table.vec)){
        if(extra.table.ind %% 2 == 1){
            temp.extra.feature.name <- extra.table.vec[extra.table.ind]
        }else if(extra.table.ind %% 2 == 0){
            message(paste0("Read extra data in ", extra.table.vec[extra.table.ind]))
            extra.table.list <- c(extra.table.list, list(read.csv(extra.table.vec[extra.table.ind], stringsAsFactors=FALSE)))
            if(is.null(names(extra.table.list))){
                names(extra.table.list)[1] <- temp.extra.feature.name
            }else{
                temp.name.vec <- names(extra.table.list)[!names(extra.table.list) %in% ""]
                names(extra.table.list) <- c(temp.name.vec, temp.extra.feature.name)
            }
        }
    }
}

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
}else if((n.ko.sample == 1) & (n.ki.sample == 1)){
    message("ANALYSIS TYPE: Single knock-in analysis")
    analysis.mode = 2
}else if((n.ko.sample == 0) & (n.ki.sample == 1)){
    message("ANALYSIS TYPE: Simple knock-in analysis")
    analysis.mode = 3
}else if((n.ko.sample == 2) & (n.ki.sample == 0)){
    message("ANALYSIS TYPE: Double knock-out analysis")
    analysis.mode = 4
}else if((n.ko.sample == 1) & (n.ki.sample == 0)){
    message("ANALYSIS TYPE: Single knock-out analysis")
    analysis.mode = 5
}
message("")
message("---------------------------------------------")
# Check whether machiato_dummy_sample label exists
is.dummy.untreated <- FALSE
if("machiato_dummy_sample" %in% label.table$label){
    is.dummy.untreated <- TRUE
}

################################################################################################################################
### Filtering
################################################################################################################################
message("Filtering data")

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

remove.target.vec <- character(0)

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
message("Check amount of classified reads...")
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
untreated.label <- label.table$label[which(label.table$sample.type == "untreated")]
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

message("Make average table regarding the rate of editing")
# fltr.rate.table %>%
#     group_by_at(c(target.type, "sample.type")) %>%
#     summarise_each(mean, matches("[A-Z]", ignore.case=FALSE)) %>% ungroup() -> mean.rate.table
fltr.rate.table %>%
    group_by_at(c(target.type, "sample.type")) %>%
    summarise(across(matches("[A-Z]", ignore.case=FALSE), mean)) %>% ungroup() -> mean.rate.table

str.col <- c("name", "lmh", "rmh", "bmh", "elmh", "ermh", "ebmh", "protospacer", "scaffold", "fourtybase.region.surrounding.targetseq", "sample.name", "group.name", "sample.type")
add.str.col <- str.col[!(str.col %in% c(target.type, "sample.type"))]
add.str.col.ind <- which(colnames(fltr.rate.table) %in% add.str.col)
for(sample.ind in 1:nrow(mean.rate.table)){
    original.table.ind <- which(fltr.rate.table[, target.type] %in% mean.rate.table[sample.ind, target.type] & fltr.rate.table$sample.type == mean.rate.table$sample.type[sample.ind])[1]
    temp.row <- cbind(
        mean.rate.table[sample.ind, ]
        , fltr.rate.table[original.table.ind, ][1, add.str.col.ind]
    )
    if(sample.ind == 1){
        average.rate.table <- temp.row
    }else{
        average.rate.table <- rbind(average.rate.table, temp.row)
    }
}
write.csv(average.rate.table, file.path(output.dir, "average_read_rates.csv"))
message("The average table is saved into ", file.path(output.dir, "average_read_rates.csv"), ".")

#############################################################################################################################################################################################################

# split dataframe into list
fltr.rate.list <- split(average.rate.table, average.rate.table$group.name)

untreated.value.set <- filter(average.rate.table, sample.type == "untreated")
response.name.vec <- ExtractResopnseName(untreated.value.set)

# Make Value set by sample.type
if(analysis.mode == 1){
    ko.1.value.set <- filter(average.rate.table, sample.type == "knock-out_1")
    ki.1.value.set <- filter(average.rate.table, sample.type == "knock-in_1")
    ki.2.value.set <- filter(average.rate.table, sample.type == "knock-in_2")
    if(nrow(untreated.value.set) == 0){ # if there is no untreated sample, make dummy dataset
        untreated.value.set <- MakeDummyValueSet(ki.1.value.set)
    }
    response.name.vec <- ExtractResopnseName(untreated.value.set)

    analized.1.value.set <- ki.1.value.set
    analized.1.value.name <- "knock-in_1"
    write.csv(analized.1.value.set,
        file.path(output.dir, paste0(analized.1.value.name, "_values.csv"))
        )
    analized.2.value.set <- ki.2.value.set
    analized.2.value.name <- "knock-in_2"
    write.csv(analized.2.value.set,
        file.path(output.dir, paste0(analized.2.value.name, "_values.csv"))
        )

    if(all(analized.1.value.set$group.name == analized.2.value.set$group.name)){ # Point !!!
        can.enhancement.analysis <- TRUE
    }else{
        can.enhancement.analysis <- FALSE
    }
    
    if(can.enhancement.analysis){
        # Make enhancement value set
        analized.enhancement.value.set <- analized.1.value.set
        analized.enhancement.value.set[, response.name.vec] <- analized.2.value.set[, response.name.vec] / (analized.1.value.set[, response.name.vec] + 1e-100) # Point !!!
        enhancement.value.name <- paste0("enhancement_", analized.1.value.name, "-", analized.2.value.name)
        write.csv(analized.enhancement.value.set,
            file.path(output.dir, paste0(enhancement.value.name, "_values.csv"))
            )
    }

}else if(analysis.mode == 2){
    ko.1.value.set <- filter(average.rate.table, sample.type == "knock-out_1")
    ki.1.value.set <- filter(average.rate.table, sample.type == "knock-in_1")
    if(nrow(untreated.value.set) == 0){
        untreated.value.set <- MakeDummyValueSet(ki.1.value.set)
    }
    response.name.vec <- ExtractResopnseName(untreated.value.set)

    analized.1.value.set <- ki.1.value.set
    analized.1.value.name <- "knock-in_1"
    write.csv(analized.1.value.set,
        file.path(output.dir, paste0(analized.1.value.name, "_values.csv"))
        )

    can.enhancement.analysis <- FALSE

}else if(analysis.mode == 3){
    ki.1.value.set <- filter(average.rate.table, sample.type == "knock-in_1")
    if(nrow(untreated.value.set) == 0){
        untreated.value.set <- MakeDummyValueSet(ki.1.value.set)
    }
    response.name.vec <- ExtractResopnseName(untreated.value.set)

    analized.1.value.set <- ki.1.value.set
    analized.1.value.name <- "knock-in_1"
    write.csv(analized.1.value.set,
        file.path(output.dir, paste0(analized.1.value.name, "_values.csv"))
        )

    can.enhancement.analysis <- FALSE

}else if(analysis.mode == 4){
    ko.1.value.set <- filter(average.rate.table, sample.type == "knock-out_1")
    ko.2.value.set <- filter(average.rate.table, sample.type == "knock-out_2")
    if(nrow(untreated.value.set) == 0){
        untreated.value.set <- MakeDummyValueSet(ko.1.value.set)
    }
    response.name.vec <- ExtractResopnseName(untreated.value.set)

    analized.1.value.set <- ko.1.value.set
    analized.1.value.name <- "knock-out_1"
    write.csv(analized.1.value.set,
        file.path(output.dir, paste0(analized.1.value.name, "_values.csv"))
        )
    analized.2.value.set <- ko.2.value.set
    analized.2.value.name <- "knock-out_2"
    write.csv(analized.2.value.set,
        file.path(output.dir, paste0(analized.2.value.name, "_values.csv"))
        )

    if(all(analized.1.value.set$group.name == analized.2.value.set$group.name)){
        can.enhancement.analysis <- TRUE
    }else{
        can.enhancement.analysis <- FALSE
    }
    
    if(can.enhancement.analysis){
        # Make enhancement value set
        analized.enhancement.value.set <- analized.1.value.set
        analized.enhancement.value.set[, response.name.vec] <- analized.2.value.set[, response.name.vec] / (analized.1.value.set[, response.name.vec] + 1e-100) # Point !!!
        enhancement.value.name <- paste0("enhancement_", analized.1.value.name, "-", analized.2.value.name)
        write.csv(analized.enhancement.value.set,
            file.path(output.dir, paste0(enhancement.value.name, "_values.csv"))
            )
    }

}else if(analysis.mode == 5){
    ko.1.value.set <- filter(average.rate.table, sample.type == "knock-out_1")
    if(nrow(untreated.value.set) == 0){
        untreated.value.set <- MakeDummyValueSet(ko.1.value.set)
    }
    response.name.vec <- ExtractResopnseName(untreated.value.set)

    analized.1.value.set <- ko.1.value.set
    analized.1.value.name <- "knock-out_1"
    write.csv(analized.1.value.set,
        file.path(output.dir, paste0(analized.1.value.name, "_values.csv"))
        )

    can.enhancement.analysis <- FALSE

}else{
    stop("Unexpected Error. Please contact the author (Kazuki Nakamae: kazukinakamae@gmail.com)")
}

# make filtered extra table
# Load MHScore
if(analysis.mode != 3){
    message("Load editing activity of knock-out sample")
    ko.editing.eff.table <- ko.1.value.set[, c("Editing.Efficiency", "group.name")]
    colnames(ko.editing.eff.table) <- c("KO.Editing.Efficiency", "group.name")
    extra.table.list <- c(extra.table.list, list(ko.editing.eff.table))
    temp.name.vec <- names(extra.table.list)[!names(extra.table.list) %in% ""]
    names(extra.table.list) <- c(temp.name.vec, "KO.Editing.Efficiency")
}else{
    extra.table.list <- extra.table.list
}
# Calc MHScore
message("Calcurate microhomology score")
# browser()
mh.list <- MakeMHScoreTable(analized.1.value.set[,c("lmh", "rmh", "bmh", "elmh", "ermh", "ebmh", "protospacer", "group.name")], target.type, script.dir, output.dir)
if(length(extra.table.list) != 0){
    extra.table.list <- c(extra.table.list, mh.list)
    temp.name.vec <- names(extra.table.list)[!names(extra.table.list) %in% ""]
    temp.name.vec <- c(temp.name.vec, "Microhomology.Score")
}else{
    extra.table.list <- mh.list
    temp.name.vec <- c("Microhomology.Score")
}
names(extra.table.list) <- temp.name.vec
for(i in 1:length(extra.table.list)){
    write.csv(extra.table.list[[i]],
            file.path(output.dir, paste0(temp.name.vec[i], "_values.csv"))
            )
}

# Calculate base features and perform comparison
for(response.name.ind in 1:length(response.name.vec)){

    response.name <- response.name.vec[response.name.ind]
    # response.name <- "Precise.Knock.in.Efficiency" # debug

    # Make save directory
    analized.1.save.dir <- file.path(output.dir, paste0(analized.1.value.name, "_", response.name))
    system(paste("mkdir", analized.1.save.dir, sep=" "))

    # Copy alignment file
    if(response.name.ind == 1){
        first.analized.1.align.file <- file.path(analized.1.save.dir, paste0(analized.1.value.name, "_", target.type, "_", response.name,  "_align.rds"))
        analized.1.align.file <- first.analized.1.align.file
    }else if(response.name.ind > 1){
        analized.1.align.file <- file.path(analized.1.save.dir, paste0(analized.1.value.name, "_", target.type, "_", response.name,  "_align.rds"))
        CopyFile(first.analized.1.align.file, analized.1.align.file)
    }
    analized.1.ManageAnalysisData <- CreateManageAnalysisData( # analized.1 / not log scale
        value.set = analized.1.value.set,
        response.variable = response.name,
        target.type = target.type,
        should.convert.log = FALSE,
        extra.table.list = extra.table.list,
        save.dir = analized.1.save.dir,
        name = analized.1.value.name
        )
    
    # Make save directory
    log.analized.1.save.dir <- file.path(output.dir, paste0(analized.1.value.name, "_", response.name, "_log"))
    system(paste("mkdir", log.analized.1.save.dir, sep=" "))

    # Copy alignment file
    log.analized.1.align.file <- file.path(log.analized.1.save.dir, paste0(analized.1.value.name, "_", target.type, "_", response.name,  "_align.rds"))
    CopyFile(analized.1.align.file, log.analized.1.align.file)
    log.analized.1.ManageAnalysisData <- CreateManageAnalysisData( # analized.1 / log scale
        value.set = analized.1.value.set,
        response.variable = response.name,
        target.type = target.type,
        should.convert.log = TRUE,
        extra.table.list = extra.table.list,
        save.dir = log.analized.1.save.dir,
        name = analized.1.value.name
        )

    if(analysis.mode == 1 | analysis.mode == 4){ # double sample
        # Make save directory
        analized.2.save.dir <- file.path(output.dir, paste0(analized.2.value.name, "_", response.name))
        system(paste("mkdir", analized.2.save.dir, sep=" "))

        # Copy alignment file
        if(response.name.ind == 1){
            first.analized.2.align.file <- file.path(analized.2.save.dir, paste0(analized.2.value.name, "_", target.type, "_", response.name,  "_align.rds"))
            analized.2.align.file <- first.analized.2.align.file
        }else if(response.name.ind > 1){
            analized.2.align.file <- file.path(analized.2.save.dir, paste0(analized.2.value.name, "_", target.type, "_", response.name,  "_align.rds"))
            CopyFile(first.analized.2.align.file, analized.2.align.file)
        }
        analized.2.ManageAnalysisData <- CreateManageAnalysisData( # analized.2 / not log scale
            value.set = analized.2.value.set,
            response.variable = response.name,
            target.type = target.type,
            should.convert.log = FALSE,
            extra.table.list = extra.table.list,
            save.dir = analized.2.save.dir,
            name = analized.2.value.name
            )

        # Make save directory
        log.analized.2.save.dir <- file.path(output.dir, paste0(analized.2.value.name, "_", response.name, "_log"))
        system(paste("mkdir", log.analized.2.save.dir, sep=" "))

        # Copy alignment file
        log.analized.2.align.file <- file.path(log.analized.2.save.dir, paste0(analized.2.value.name, "_", target.type, "_", response.name,  "_align.rds"))
        CopyFile(analized.2.align.file, log.analized.2.align.file)
        log.analized.2.ManageAnalysisData <- CreateManageAnalysisData( # analized.2 / log scale
            value.set = analized.2.value.set,
            response.variable = response.name,
            target.type = target.type,
            should.convert.log = TRUE,
            extra.table.list = extra.table.list,
            save.dir = log.analized.2.save.dir,
            name = analized.2.value.name
            )

        if(can.enhancement.analysis){
            ### Enhanement
            # Make save directory
            enhancement.save.dir <- file.path(output.dir, paste0(enhancement.value.name, "_", response.name))
            system(paste("mkdir", enhancement.save.dir, sep=" "))

            # Copy alignment file
            if(response.name.ind == 1){
                first.enhancement.align.file <- file.path(enhancement.save.dir, paste0(enhancement.value.name, "_", target.type, "_", response.name,  "_align.rds"))
                enhancement.align.file <- first.enhancement.align.file
            }else if(response.name.ind > 1){
                enhancement.align.file <- file.path(enhancement.save.dir, paste0(enhancement.value.name, "_", target.type, "_", response.name,  "_align.rds"))
                CopyFile(first.enhancement.align.file, enhancement.align.file)
            }

            enhancement.ManageAnalysisData <- CreateManageAnalysisData( # enhancement / not log scale
                value.set = analized.enhancement.value.set,
                response.variable = response.name,
                target.type = target.type,
                should.convert.log = FALSE,
                extra.table.list = extra.table.list,
                save.dir = enhancement.save.dir,
                name = enhancement.value.name
                )
            
            # Make save directory
            log.enhancement.save.dir <- file.path(output.dir, paste0(enhancement.value.name, "_", response.name, "_log"))
            system(paste("mkdir", log.enhancement.save.dir, sep=" "))

            # Copy alignment file
            log.enhancement.align.file <- file.path(log.enhancement.save.dir, paste0(enhancement.value.name, "_", target.type, "_", response.name,  "_align.rds"))
            CopyFile(enhancement.align.file, log.enhancement.align.file)
            log.enhancement.ManageAnalysisData <- CreateManageAnalysisData( # enhancement / log scale
                value.set = analized.enhancement.value.set,
                response.variable = response.name,
                target.type = target.type,
                should.convert.log = TRUE,
                extra.table.list = extra.table.list,
                save.dir = log.enhancement.save.dir,
                name = enhancement.value.name
                )
        }

        ### Comparison
        # Make save directory
        comp.save.dir <- file.path(output.dir, paste0("comparison_", response.name))
        system(paste("mkdir", comp.save.dir, sep=" "))

        # Perform comparison
        if(
            CompDoubleData(
                ManageAnalysisData.1 = analized.1.ManageAnalysisData,
                ManageAnalysisData.2 = analized.2.ManageAnalysisData,
                value.name.1 = analized.1.value.name,
                value.name.2 = analized.2.value.name,
                save.dir = comp.save.dir,
                response.name = response.name,
                should.convert.log = FALSE
            )
        ){
            message("Done")
        }

        # Perform comparison (log scale version)
        if(log.analized.1.ManageAnalysisData$can.convert.log & log.analized.2.ManageAnalysisData$can.convert.log){
            log.comp.save.dir <- file.path(output.dir, paste0("comparison_", response.name, "_log"))
            system(paste("mkdir", log.comp.save.dir, sep=" "))
            if(
                CompDoubleData(
                    ManageAnalysisData.1 = log.analized.1.ManageAnalysisData,
                    ManageAnalysisData.2 = log.analized.2.ManageAnalysisData,
                    value.name.1 = paste0(analized.1.value.name, "_log"),
                    value.name.2 = paste0(analized.2.value.name, "_log"),
                    save.dir = log.comp.save.dir,
                    response.name = response.name,
                    should.convert.log = TRUE
                )
            ){
                message("Done")
            }
        }
    }
    # q()  # debug
}

message(paste0("Analysis is done. Data was saved into ", output.dir))