debug.flag <- FALSE
############################################################################
# NOTICE : The shown total read is not actual number of reads.
# Because MaChIAto Aligner show only read which is avaliable for alignment. we call it as "alignment bias".
# And, Total reads of analysis [1] and [2] is different because the references sequence is not same.
# If reference sequence is same, data will be processed in same alignment bias. This example is comparison between original CRISPResso data and CRISPResso-MaChIAto data.
############################################################################
### TODO
## reverse ki analysis
# Left reverse ki alignment may be less sensitive than Right reverse ki alignment.
# this cause is located in bam alignment setting. but, the exact point remains unclear.
############################################################################
#-- Setting ----------------------------------------------------#

message("
.___  ___.      ___       ______  __    __   __       ___   .___________.  ______   
|   \\/   |     /   \\     /      ||  |  |  | |  |     /   \\  |           | /  __  \\  
|  \\  /  |    /  ^  \\   |  ,----'|  |__|  | |  |    /  ^  \\ `---|  |----`|  |  |  | 
|  |\\/|  |   /  /_\\  \\  |  |     |   __   | |  |   /  /_\\  \\    |  |     |  |  |  | 
|  |  |  |  /  _____  \\ |  `----.|  |  |  | |  |  /  _____  \\   |  |     |  `--'  | 
|__|  |__| /__/     \\__\\ \\______||__|  |__| |__| /__/     \\__\\  |__|      \\______/  
                                                                                    
     ___       __       __    _______ .__   __.  _______ .______      
    /   \\     |  |     |  |  /  _____||  \\ |  | |   ____||   _  \\     
   /  ^  \\    |  |     |  | |  |  __  |   \\|  | |  |__   |  |_)  |    
  /  /_\\  \\   |  |     |  | |  | |_ | |  . `  | |   __|  |      /     
 /  _____  \\  |  `----.|  | |  |__| | |  |\\   | |  |____ |  |\\  \\----.
/__/     \\__\\ |_______||__|  \\______| |__| \\__| |_______|| _| `._____|
                                                                      
Author: Kazuki Nakamae at Sakuma Tetsushi and Takashi Yamamoto lab, Hiroshima Univ, Japan
For support contact kazukinakamae@gmail.com
version 1.0
")

help.message <- "

-------------------------------------------------------------------------------
This is not correct input. The description below would be helpful.
-------------------------------------------------------------------------------


-------------------------------------------------------------------------------
EXAMPLE INPUT
-------------------------------------------------------------------------------
Rscript ./MAChIAto/MaChIAto_Aligner/MaChIAtoAligner.R <Input directory> <Output prefix> <Left extra sequence> <Right extra sequence>

-------------------------------------------------------------------------------
DESPRIPTION
-------------------------------------------------------------------------------

Input directory:
The directory that MaChIAto Classifier generated.

Output prefix:
The directory into which the output directory is saved.

Left extra sequence (optional):
The outside sequence flanking 5'-homology arm of donor. 

Right extra sequence (optional):
The outside sequence flanking 3'-homology arm of donor. 

-------------------------------------------------------------------------------

"

# Get directory path
if(debug.flag == TRUE){
    script.dir <- ""
}else{
    initial.options <- commandArgs(trailingOnly = FALSE)
    file.arg.name <- "--file="
    script.fn <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
    script.dir <- dirname(script.fn)
}
ref.dir <- file.path(script.dir, "ref")



#-- Run Script ----------------------------------------------------#
message("---Load MaChIAto_Aligner programs---")

source(file.path(script.dir, "MaChIAtoAlignerPackageInstaller.R"))
source(file.path(script.dir, "MaChIAtoAlignerFunctions.R"))
source(file.path(script.dir, "RunLocalMutationAnalysis.R"))
source(file.path(script.dir, "RunLocalKIAnalysis.R"))


#-- Load library ----------------------------------------------------#
message("---Load libraries---")

# library(vctrs)
library(grDevices)
library(network)
library(dplyr)
library(reshape2)
library(stringr)
library(Biostrings)
library(Rsamtools)
library(GenomicRanges)
library(CrispRVariants)

#-------------------------------------------------------------------#

message("----------------------------")
message("---Start MaChIAto_Aligner---")
message("----------------------------")



left.extra.seq <- "" # default
right.extra.seq <- "" # default
if(debug.flag == TRUE){
    input.fn <- ""
    output.prefix <- ""
    left.extra.seq <- ""
    right.extra.seq <- ""
}else{
    if(length(commandArgs(trailingOnly=TRUE)) < 2){
      message(help.message)
      q()
    }
    # Get commandline arguments
    input.fn <- GetAbsolutePath(commandArgs(trailingOnly=TRUE)[1])
    output.prefix <- GetAbsolutePath(commandArgs(trailingOnly=TRUE)[2])
    if(length(commandArgs(trailingOnly=TRUE)) > 2){
      left.extra.seq <- commandArgs(trailingOnly=TRUE)[3]
      right.extra.seq <- commandArgs(trailingOnly=TRUE)[4]
    }
}
input.fn <- gsub("/$", "", input.fn)
output.prefix <- gsub("/$", "", output.prefix)

#-- Make output directory ----------------------------------------------------#
message("---Make directory for the analysis---")

#time.id <- "20190404150926"
time.id <- format(Sys.time(), "%Y%m%d%H%M%S")
output.name <- paste0("MaChIAtoAligner_at_", time.id, "_from_", basename(input.fn))
dir.create(file.path(output.prefix), showWarnings=FALSE)
output.dir <- file.path(output.prefix, output.name)

message(paste("Creating Folder", output.dir, sep=" "))

# index
index.dir <- file.path(output.dir, "index")
# data dir
data.dir <- file.path(output.dir, "data")
seq.dir <- file.path(data.dir, "seq")
fasta.dir <- file.path(data.dir, "fasta")
fake_fastq.dir <- file.path(data.dir, "fake_fastq")
bam.dir <- file.path(data.dir, "bam")
# fasta dir
fasta.crispresso.dir <- file.path(fasta.dir, "crispresso")
fasta.machiato.dir <- file.path(fasta.dir, "machiato")
# fake_fastq dir
fake_fastq.crispresso.dir <- file.path(fake_fastq.dir, "crispresso")
fake_fastq.machiato.dir <- file.path(fake_fastq.dir, "machiato")
# bam dir for comparison analysis
bam.crispresso.dir <- file.path(bam.dir, "crispresso")
bam.machiato.dir <- file.path(bam.dir, "machiato")
### result data dir
result.dir <- file.path(output.dir, "result")
## result data dir
comp.with.original.crispresso.data.dir <- file.path(result.dir, "[1]comparison_with_original_crispresso_data")
machiato.local.alignment.dir <- file.path(result.dir, "[2]machiato_local_alignment")
# result data dir for comparison analysis
crispresso.alignment.data.dir <- file.path(comp.with.original.crispresso.data.dir, "[Ⅰ]crispresso_alignment_data")
machiato.alignment.data.dir <- file.path(comp.with.original.crispresso.data.dir, "[Ⅱ]machiato_alignment_data")
# comparison.analysis.data.dir <- file.path(comp.with.original.crispresso.data.dir, "[Ⅲ]comparison_analysis_data")
crispresso.all.dir <- file.path(crispresso.alignment.data.dir, "[Aa]crispresso_all")
crispresso.hdr.dir <- file.path(crispresso.alignment.data.dir, "[Ae]crispresso_hdr")
crispresso.mixed_hdr_nhej.dir <- file.path(crispresso.alignment.data.dir, "[Ad]crispresso_mixed_hdr_nhej")
crispresso.nhej.dir <- file.path(crispresso.alignment.data.dir, "[Ac]crispresso_nhej")
crispresso.unmodified.dir <- file.path(crispresso.alignment.data.dir, "[Ab]crispresso_unmodified")
crispresso.rate.dir <- file.path(crispresso.alignment.data.dir, "[B]crispresso_rate_data")
machiato.all.dir <- file.path(machiato.alignment.data.dir, "[Aa]machiato_all")
machiato.hdr.dir <- file.path(machiato.alignment.data.dir, "[Ae]machiato_hdr")
machiato.mixed_hdr_nhej.dir <- file.path(machiato.alignment.data.dir, "[Ad]machiato_mixed_hdr_nhej")
machiato.nhej.dir <- file.path(machiato.alignment.data.dir, "[Ac]machiato_nhej")
machiato.unmodified.dir <- file.path(machiato.alignment.data.dir, "[Ab]machiato_unmodified")
machiato.unclassified.dir <- file.path(machiato.alignment.data.dir, "[Af]machiato_unclassified")
machiato.rate.dir <- file.path(machiato.alignment.data.dir, "[B]machiato_rate_data")

for(dirname in c(output.dir, index.dir, data.dir, seq.dir, fasta.dir, fake_fastq.dir
  , bam.dir, result.dir, machiato.local.alignment.dir, comp.with.original.crispresso.data.dir
  , crispresso.alignment.data.dir, machiato.alignment.data.dir, fasta.crispresso.dir, fasta.machiato.dir
  , fake_fastq.crispresso.dir, fake_fastq.machiato.dir, bam.crispresso.dir, bam.machiato.dir
  , crispresso.all.dir, crispresso.hdr.dir, crispresso.mixed_hdr_nhej.dir, crispresso.nhej.dir
  , crispresso.unmodified.dir, crispresso.rate.dir, machiato.all.dir, machiato.hdr.dir, machiato.mixed_hdr_nhej.dir
  , machiato.nhej.dir, machiato.unmodified.dir, machiato.unclassified.dir, machiato.rate.dir)){
  dir.create(file.path(dirname), showWarnings=FALSE)
}
# for(data.type in c("ALL", "Unmodified", "NHEJ", "Mixed", "HDR", "Unclassified")){
for(data.type in c("Unmodified", "NHEJ", "Mixed", "HDR")){
  # result file of each data type
  result.datatype.dir <- file.path(machiato.local.alignment.dir, data.type)
  dir.create(file.path(result.datatype.dir), showWarnings=FALSE)
}

#-- make mapped bam file by bwa mem ----------------------------------------------------#
message("---Run bwa mem---")

system(paste("chmod +x", file.path(script.dir, "prepare_ref.sh"), sep=" "))
system(paste(file.path(script.dir, "prepare_ref.sh"), input.fn, index.dir, seq.dir , sep=" "))

mmej.ki.seq <- as.character(readDNAStringSet(file.path(index.dir, "mmej_ki.fa")))
donor.seq <- as.character(readDNAStringSet(file.path(seq.dir, "donor.fa")))
left.homology.seq <- as.character(readDNAStringSet(file.path(seq.dir, "Leftside_homology_arm.fa")))
right.homology.seq <- as.character(readDNAStringSet(file.path(seq.dir, "Rightside_homology_arm.fa")))
nhej.donor.seq <- paste(left.extra.seq, left.homology.seq, donor.seq, right.homology.seq, right.extra.seq, sep="")
nhej.rc.donor.seq <- as.character(reverseComplement(DNAString(nhej.donor.seq)))
nhej.ki.seq <- DNAStringSet(sub(paste0(left.homology.seq, donor.seq), paste0(left.homology.seq, nhej.donor.seq), mmej.ki.seq))
names(nhej.ki.seq) <- c("nhej_ki")
nhej.rc.ki.seq <- DNAStringSet(sub(paste0(left.homology.seq, donor.seq), paste0(left.homology.seq, nhej.rc.donor.seq), mmej.ki.seq))
names(nhej.rc.ki.seq) <- c("nhej_rc_ki")

writeXStringSet(nhej.ki.seq, file.path(index.dir, "nhej_ki.fa"))
writeXStringSet(nhej.rc.ki.seq, file.path(index.dir, "nhej_rc_ki.fa"))
system(paste("chmod +x", file.path(script.dir, "prepare_index.sh"), sep=" "))
system(paste(file.path(script.dir, "prepare_index.sh"), index.dir, sep=" "))
system(paste("chmod +x", file.path(script.dir, "prepare_bam_comparing_CRISPResso.sh"), sep=" "))
system(paste(file.path(script.dir, "prepare_bam_comparing_CRISPResso.sh"), input.fn, index.dir, data.dir, script.dir, sep=" "))

for(data.type in c("ALL", "Unclassified", "Unmodified", "NHEJ", "Mixed", "HDR")){

  # fasta file prefix
  fasta.datatype.dir <- file.path(fasta.dir, data.type)
  # fasta dir
  fasta.muttype.dir <- file.path(fasta.datatype.dir, "mut_type")
  fasta.ki_type.dir <- file.path(fasta.datatype.dir, "ki_type")
  fasta.rc_ki_type.dir <- file.path(fasta.datatype.dir, "rc_ki_type")
  fasta.other_type.dir <- file.path(fasta.datatype.dir, "other_type")
  # fastaq file prefix
  fake_fastq.datatype.dir <- file.path(fake_fastq.dir, data.type)
  # fastq dir
  fake_fastq.muttype.dir <- file.path(fake_fastq.datatype.dir, "mut_type")
  fake_fastq.ki_type.dir <- file.path(fake_fastq.datatype.dir, "ki_type")
  fake_fastq.rc_ki_type.dir <- file.path(fake_fastq.datatype.dir, "rc_ki_type")
  fake_fastq.other_type.dir <- file.path(fake_fastq.datatype.dir, "other_type")
  # bam file prefix
  bam.datatype.dir <- file.path(bam.dir, data.type)
  # bam dir
  bam.muttype.dir <- file.path(bam.datatype.dir, "mut_type")
  bam.ki_type.dir <- file.path(bam.datatype.dir, "ki_type")
  bam.rc_ki_type.dir <- file.path(bam.datatype.dir, "rc_ki_type")
  bam.other_type.dir <- file.path(bam.datatype.dir, "other_type")
  # make bam dir for local alignment analysis
  for(dirname in c(
      fasta.datatype.dir, fasta.muttype.dir, fasta.ki_type.dir, fasta.rc_ki_type.dir, fasta.other_type.dir
      , fake_fastq.datatype.dir, fake_fastq.muttype.dir, fake_fastq.ki_type.dir, fake_fastq.rc_ki_type.dir, fake_fastq.other_type.dir
      , bam.datatype.dir, bam.muttype.dir, bam.ki_type.dir, bam.rc_ki_type.dir, bam.other_type.dir)
    ){
    dir.create(file.path(dirname), showWarnings=FALSE)
  }
}
# make mapped bam file of each machiato output
for(data.type in c("Unclassified", "Unmodified", "NHEJ", "Mixed", "HDR")){
  system(paste("chmod +x", file.path(script.dir, "prepare_bam_analysis_each_outcome_pe.sh"), sep=" "))
  system(paste(file.path(script.dir, "prepare_bam_analysis_each_outcome_pe.sh"), input.fn, index.dir, data.dir, script.dir, data.type, sep=" "))
}

system(paste("chmod +x", file.path(script.dir, "prepare_seq.sh"), sep=" "))
system(paste(file.path(script.dir, "prepare_seq.sh"), input.fn, data.dir, sep=" "))

###################################################################################################################################################
### Load sequence
message("---Load key sequence for analysis---")

wt.seq <- as.character(readDNAStringSet(file.path(index.dir, "wt.fa")))
protospacer.seq <- as.character(readDNAStringSet(file.path(seq.dir, "sgRNA.fa")))
donor.seq <- as.character(readDNAStringSet(file.path(seq.dir, "donor.fa")))
inleft.indicator.seq <- as.character(readDNAStringSet(file.path(seq.dir, "In-Left_indicator_sequence.fa")))
inright.indicator.seq <- as.character(readDNAStringSet(file.path(seq.dir, "In-Right_indicator_sequence.fa")))
left.homology.seq <- as.character(readDNAStringSet(file.path(seq.dir, "Leftside_homology_arm.fa")))
outleft.indicator.seq <- as.character(readDNAStringSet(file.path(seq.dir, "Out-Left_indicator_sequence.fa")))
outright.indicator.seq <- as.character(readDNAStringSet(file.path(seq.dir, "Out-Right_indicator_sequence.fa")))
right.homology.seq <- as.character(readDNAStringSet(file.path(seq.dir, "Rightside_homology_arm.fa")))
untreated.seq <- as.character(readDNAStringSet(file.path(seq.dir, "Untreated_sequence.fa")))

###################################################################################################################################################

### CRISPResso original output
message("---Comparison with original crispresso data---")
###
# ALL
crispresso.all.crispr.set <- MakeWtCrisprSet(
untreated.seq
, "wt"
, wt.seq
, protospacer.seq
, file.path(bam.crispresso.dir, "all.bam")
, "crispresso.all"
, 35)
if(!is.null(crispresso.all.crispr.set)){
  SaveVariantsData(crispresso.all.crispr.set, crispresso.all.dir, mut.type = "mut")
  SaveVariantsPlot(crispresso.all.crispr.set, crispresso.all.dir)
}

###

# HDR
crispresso.hdr.crispr.set <- MakeWtCrisprSet(
untreated.seq
, "wt"
, wt.seq
, protospacer.seq
, file.path(bam.crispresso.dir, "hdr.bam")
, "crispresso.hdr"
, 35)
if(!is.null(crispresso.hdr.crispr.set)){
  SaveVariantsData(crispresso.hdr.crispr.set, crispresso.hdr.dir, mut.type = "mut")
  SaveVariantsPlot(crispresso.hdr.crispr.set, crispresso.hdr.dir)
}

###

# Mixed_HDR_NHEJ
crispresso.mixed_hdr_nhej.crispr.set <- MakeWtCrisprSet(
untreated.seq
, "wt"
, wt.seq
, protospacer.seq
, file.path(bam.crispresso.dir, "mixed_hdr_nhej.bam")
, "crispresso.mixed_hdr_nhej"
, 35)
if(!is.null(crispresso.mixed_hdr_nhej.crispr.set)){
  SaveVariantsData(crispresso.mixed_hdr_nhej.crispr.set, crispresso.mixed_hdr_nhej.dir, mut.type = "mut")
  SaveVariantsPlot(crispresso.mixed_hdr_nhej.crispr.set, crispresso.mixed_hdr_nhej.dir)
}

###

# NHEJ
crispresso.nhej.crispr.set <- MakeWtCrisprSet(
untreated.seq
, "wt"
, wt.seq
, protospacer.seq
, file.path(bam.crispresso.dir, "nhej.bam")
, "crispresso.nhej"
, 35)
if(!is.null(crispresso.nhej.crispr.set)){
  SaveVariantsData(crispresso.nhej.crispr.set, crispresso.nhej.dir, mut.type = "mut")
  SaveVariantsPlot(crispresso.nhej.crispr.set, crispresso.nhej.dir)
}

###

# Unmodified
crispresso.unmodified.crispr.set <- MakeWtCrisprSet(
untreated.seq
, "wt"
, wt.seq
, protospacer.seq
, file.path(bam.crispresso.dir, "unmodified.bam")
, "crispresso.unmodified"
, 35)
if(!is.null(crispresso.unmodified.crispr.set)){
  SaveVariantsData(crispresso.unmodified.crispr.set, crispresso.unmodified.dir, mut.type = "mut")
  SaveVariantsPlot(crispresso.unmodified.crispr.set, crispresso.unmodified.dir)
}

###

### MaChIAto original output
# ALL
machiato.all.crispr.set <- MakeWtCrisprSet(
untreated.seq
, "wt"
, wt.seq
, protospacer.seq
, file.path(bam.machiato.dir, "all.bam")
, "machiato.all"
, 35)
if(!is.null(machiato.all.crispr.set)){
  SaveVariantsData(machiato.all.crispr.set, machiato.all.dir, mut.type = "mut")
  SaveVariantsPlot(machiato.all.crispr.set, machiato.all.dir)
}

###

# HDR
machiato.hdr.crispr.set <- MakeWtCrisprSet(
untreated.seq
, "wt"
, wt.seq
, protospacer.seq
, file.path(bam.machiato.dir, "hdr.bam")
, "machiato.hdr"
, 35)
if(!is.null(machiato.hdr.crispr.set)){
  SaveVariantsData(machiato.hdr.crispr.set, machiato.hdr.dir, mut.type = "mut")
  SaveVariantsPlot(machiato.hdr.crispr.set, machiato.hdr.dir)
}

###

# Mixed_HDR_NHEJ
machiato.mixed_hdr_nhej.crispr.set <- MakeWtCrisprSet(
untreated.seq
, "wt"
, wt.seq
, protospacer.seq
, file.path(bam.machiato.dir, "mixed_hdr_nhej.bam")
, "machiato.mixed_hdr_nhej"
, 35)
if(!is.null(machiato.mixed_hdr_nhej.crispr.set)){
  SaveVariantsData(machiato.mixed_hdr_nhej.crispr.set, machiato.mixed_hdr_nhej.dir, mut.type = "mut")
  SaveVariantsPlot(machiato.mixed_hdr_nhej.crispr.set, machiato.mixed_hdr_nhej.dir)
}

###

# NHEJ
machiato.nhej.crispr.set <- MakeWtCrisprSet(
untreated.seq
, "wt"
, wt.seq
, protospacer.seq
, file.path(bam.machiato.dir, "nhej.bam")
, "machiato.nhej"
, 35)
if(!is.null(machiato.nhej.crispr.set)){
  SaveVariantsData(machiato.nhej.crispr.set, machiato.nhej.dir, mut.type = "mut")
  SaveVariantsPlot(machiato.nhej.crispr.set, machiato.nhej.dir)
}

###

# Unmodified
# browser()
machiato.unmodified.crispr.set <- MakeWtCrisprSet(
untreated.seq
, "wt"
, wt.seq
, protospacer.seq
, file.path(bam.machiato.dir, "unmodified.bam")
, "machiato.unmodified"
, 35)
if(!is.null(machiato.unmodified.crispr.set)){
  SaveVariantsData(machiato.unmodified.crispr.set, machiato.unmodified.dir, mut.type = "mut")
  SaveVariantsPlot(machiato.unmodified.crispr.set, machiato.unmodified.dir)
}

###

# Unclassified
machiato.unclassified.crispr.set <- MakeWtCrisprSet(
untreated.seq
, "wt"
, wt.seq
, protospacer.seq
, file.path(bam.machiato.dir, "unclassified.bam")
, "machiato.unclassified"
, 35)
if(!is.null(machiato.unclassified.crispr.set)){
  SaveVariantsData(machiato.unclassified.crispr.set, machiato.unclassified.dir, mut.type = "mut")
  SaveVariantsPlot(machiato.unclassified.crispr.set, machiato.unclassified.dir)
}

###

if(!is.null(crispresso.all.crispr.set) && !is.null(machiato.all.crispr.set)){
  message("---Make plot each rate of original CRISPResso data---")
  
  crispresso.count.sum <- c(sum(crispresso.unmodified.crispr.set$cigar_freqs[,1])
      , sum(crispresso.nhej.crispr.set$cigar_freqs[,1])
      , sum(crispresso.mixed_hdr_nhej.crispr.set$cigar_freqs[,1])
      , sum(crispresso.hdr.crispr.set$cigar_freqs[,1])
    )
  names(crispresso.count.sum) <- c("Unmodified", "NHEJ", "Mixed HDR-NHEJ", "HDR")
  crispresso.class.count.table <- data.frame(
    Class=names(crispresso.count.sum)
    , Frequency=crispresso.count.sum
    , Percentage=crispresso.count.sum / sum(crispresso.count.sum) * 100
  )
  
  saveRDS(crispresso.class.count.table, file = file.path(crispresso.rate.dir, "[ⅰ]CRISPResso_Unmodified_NHEJ_HDR_table.rds"))
  SavePieChart(crispresso.class.count.table
    , c("#D98F4E", "#A65437", "#F2CA52", "#F2E2CE")
    , "Reads"
    , file.path(crispresso.rate.dir, "[ⅰ]CRISPResso_Unmodified_NHEJ_HDR_pie_chart.png")
    , file.path(crispresso.rate.dir, "[ⅰ]CRISPResso_Unmodified_NHEJ_HDR_table.csv") # TODO : add ID number
  )

  crispresso.pre.bwa.read.sum <- sum(as.numeric(read.csv(file.path(input.fn, "MaChIAto_alignment_pie_chart.csv"), stringsAsFactors=FALSE)[2, ]))
  crispresso.bwa.read.sum <- sum(crispresso.class.count.table$Frequency)
  crispresso.bwa.read.sum.table <- data.frame(
    Position=factor(c("Before Alignment", "After Alignment"), levels = c("Before Alignment", "After Alignment"))
    , Frequency=c(crispresso.pre.bwa.read.sum, crispresso.bwa.read.sum)
  )
  write.table(rbind(as.character(crispresso.bwa.read.sum.table$Position), crispresso.bwa.read.sum.table$Frequency), file = file.path(crispresso.rate.dir, "[i]The_amount_of_alignment_loss_with_CRISPResso-BWAMEM.csv")
    , quote=FALSE, col.names=FALSE, row.names=FALSE, append=TRUE, sep = ",")
  saveRDS(crispresso.bwa.read.sum.table, file = file.path(crispresso.rate.dir, "[i]The_amount_of_alignment_loss_with_CRISPResso-BWAMEM.rds"))
  SaveCompBarPlot(crispresso.bwa.read.sum.table
      , color.code = c("#2768B3", "#FFC152")
      , file.name = file.path(crispresso.rate.dir, "[i]The_amount_of_alignment_loss_with_CRISPResso-BWAMEM.png")
      , xlab = "Change of total reads"
      , ylab = "Survival Reads"
    )
  

  message("---Make plot each rate of original CRISPResso data (including Unclassified class)---")
  
  machiato.including.unclassified.class.count.sum <- c(sum(machiato.unclassified.crispr.set$cigar_freqs[,1])
    , sum(machiato.unmodified.crispr.set$cigar_freqs[,1])
    , sum(machiato.nhej.crispr.set$cigar_freqs[,1])
    , sum(machiato.mixed_hdr_nhej.crispr.set$cigar_freqs[,1])
    , sum(machiato.hdr.crispr.set$cigar_freqs[,1])
  )
  names(machiato.including.unclassified.class.count.sum) <- c("Unclassified", "Unmodified", "NHEJ", "Mixed HDR-NHEJ", "HDR")
  machiato.including.unclassified.class.count.table <- data.frame(
    Class=names(machiato.including.unclassified.class.count.sum)
    , Frequency=machiato.including.unclassified.class.count.sum
    , Percentage=machiato.including.unclassified.class.count.sum / sum(machiato.including.unclassified.class.count.sum) * 100
  )
  
  saveRDS(machiato.including.unclassified.class.count.table, file = file.path(machiato.rate.dir, "[ⅰ]CRISPResso-MaChIAto_Unmodified_NHEJ_HDR_table(including_Unclassified).rds"))
  SavePieChart(machiato.including.unclassified.class.count.table
  , c("#D98F4E",  "#A65437", "#F2CA52", "#F2E2CE", "#ADD4D9")
  , "Reads"
  , file.path(machiato.rate.dir, "[ⅰ]CRISPResso-MaChIAto_Unmodified_NHEJ_HDR_pie_chart(including_Unclassified).png")
  , file.path(machiato.rate.dir, "[ⅰ]CRISPResso-MaChIAto_Unmodified_NHEJ_HDR_table(including_Unclassified).csv")
  )

  message("---Make plot each rate of CRISPResso-MaChIAto data---")
  machiato.class.count.sum <- c(sum(machiato.unmodified.crispr.set$cigar_freqs[,1])
    , sum(machiato.nhej.crispr.set$cigar_freqs[,1])
    , sum(machiato.mixed_hdr_nhej.crispr.set$cigar_freqs[,1])
    , sum(machiato.hdr.crispr.set$cigar_freqs[,1])
  )
  names(machiato.class.count.sum) <- c("Unmodified", "NHEJ", "Mixed HDR-NHEJ", "HDR")
  machiato.class.count.table <- data.frame(
    Class=names(machiato.class.count.sum)
    , Frequency=machiato.class.count.sum
    , Percentage=machiato.class.count.sum / sum(machiato.class.count.sum) * 100
  )
  
  saveRDS(machiato.class.count.table, file = file.path(machiato.rate.dir, "[ⅱ]CRISPResso-MaChIAto_Unmodified_NHEJ_HDR_table.rds"))
  SavePieChart(machiato.class.count.table
  , c("#D98F4E", "#A65437", "#F2CA52", "#F2E2CE")
  , "Reads"
  , file.path(machiato.rate.dir, "[ⅱ]CRISPResso-MaChIAto_Unmodified_NHEJ_HDR_pie_chart.png")
  , file.path(machiato.rate.dir, "[ⅱ]CRISPResso-MaChIAto_Unmodified_NHEJ_HDR_table.csv")
  )

  machiato.pre.bwa.read.sum <- sum(as.numeric(read.csv(file.path(input.fn, "MaChIAto_alignment_pie_chart.csv"), stringsAsFactors=FALSE)[2, ]))
  machiato.bwa.read.sum <- sum(machiato.including.unclassified.class.count.table$Frequency)
  machiato.bwa.read.sum.table <- data.frame(
    Position=factor(c("Before Alignment", "After Alignment"), levels = c("Before Alignment", "After Alignment"))
    , Frequency=c(machiato.pre.bwa.read.sum, machiato.bwa.read.sum)
  )
  write.table(rbind(as.character(machiato.bwa.read.sum.table$Position), machiato.bwa.read.sum.table$Frequency), file = file.path(machiato.rate.dir, "[ii]The_amount_of_alignment_loss_with_MaChIAto-BWAMEM.csv")
    , quote=FALSE, col.names=FALSE, row.names=FALSE, append=TRUE, sep = ",")
  saveRDS(machiato.bwa.read.sum.table, file = file.path(machiato.rate.dir, "[ii]The_amount_of_alignment_loss_with_MaChIAto-BWAMEM.rds"))
  SaveCompBarPlot(machiato.bwa.read.sum.table
      , color.code = c("#2768B3", "#FFC152")
      , file.name = file.path(machiato.rate.dir, "[ii]The_amount_of_alignment_loss_with_MaChIAto-BWAMEM.png")
      , xlab = "Change of total reads"
      , ylab = "Survival Reads"
  )
}
message("---Run MaChIAto local alignment---")

# for(data.type in c("ALL", "Unmodified", "NHEJ", "Mixed", "HDR", "Unclassified")){
for(data.type in c("Unmodified", "NHEJ", "Mixed", "HDR")){
  message(paste("Align ", data.type, " data...", sep=""))
  # bam file prefix
  bam.datatype.dir <- file.path(bam.dir, data.type)
  # result file of each type for local alignment data
  result.datatype.dir <- file.path(machiato.local.alignment.dir, data.type)
  # if(any(data.type %in% c("ALL", "Unclassified", "Unmodified", "NHEJ"))){
  if(any(data.type %in% c("Unmodified", "NHEJ"))){
    # run local mutation analysis
    RunLocalMutationAnalysis(bam.datatype.dir, result.datatype.dir, index.dir, seq.dir)
  }
  # if(any(data.type %in% c("ALL", "Unclassified", "Mixed", "HDR"))){
  if(any(data.type %in% c("Mixed", "HDR"))){
    # run local knock-in analysis
    RunLocalKIAnalysis(bam.datatype.dir, result.datatype.dir, index.dir, seq.dir)
  }
}

message("---The process is completed.---")
quit(save = "no")