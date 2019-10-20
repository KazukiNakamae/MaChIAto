# cd /Users/kazuki/Documents/Lab/MachinelearningForPITChprediction/20190218_development_of_MAChAto;
# conda create -n MaChIAto_Aligner samtools r r-essentials bioconductor-crisprvariants bioconductor-genomicfeatures bioconductor-rsamtools bwa;
# source activate MaChIAto_Aligner;
# source("./MaChIAto_Aligner/MaChIAtoAlignerCore.R")
# quit(save = "no") # for debug
###################################################################################################################################################
############################################################################
# NOTICE : The shown total read is not actual number of reads.
# Because MaChIAto Aligner show only read which is avaliable for alignment. we call it as "alignment bias".
# And, Total reads of analysis [1] and [2] is different because the references sequence is not same.
# If reference sequence is same, data will be processed in same alignment bias. This example is comparison between original CRISPResso data and CRISPResso-MaChIAto data.
############################################################################
### TODO
## reverse ki analysis
# Left reverse ki alignment may be less sensitive than Right reverse ki alignment.
# this cause is located in bam alignment setting. but the exact point remain unclear.
############################################################################
#-- Setting ----------------------------------------------------#

# version 1.3

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
                                                                      
version 1.3.beta
")

help.message <- "

-------------------------------------------------------------------------------
This is not correct input. The description below would be helpful.
-------------------------------------------------------------------------------


-------------------------------------------------------------------------------
EXAMPLE INPUT
-------------------------------------------------------------------------------
Rscript ./MAChIAto/MaChIAto_Analyzer/MaChIAtoAnalyzer.R <summary directory of all results from MaChIAto Classifier> <output directory> <table of extra data> <list of ignore target set>

-------------------------------------------------------------------------------
DESPRIPTION
-------------------------------------------------------------------------------

summary directory of all results from MaChIAto Classifier:
The summary directory of all results from MaChIAto Classifier can be generated using “collect_MaChIAto_data.py” in “./MAChIAto/MaChIAto_Analyzer”.

output directory:
The directory name which is to be put output directory into. The name should include path.

table of extra data:
The table of extra data can contain additional data such as InDelphi score, accessibility, chromatin size, exon number, position of amino acid, protein length, CNV_HMM and expression. “example_data” directory has “extra_data.csv” as example. As the table is described, the format requires value per target. If some values are nothing, the cell should be NaN. The value of example data is calculated with ENCODE database (https://www.encodeproject.org). The detailed script used for calculating the values of example data is described in our GihHub page. Please see the detailed instruction in our GitHub page: https://github.com/Kazuki-Nakamae/MAChIAto.

list of ignore target set:
The list of ignore target set contains target names which are not desired to analyze for some reasons. The data (e.g. DBF4B-A, DBF4B-B, DBF4B-C, DBF4B-D) including target name (e.g. DBF4B) shown in the list is skipped through the process of MaChIAto Analyzer. The format should be comma-separated like “targetA, TargetB, …” “example_data” directory has “ignore_list.csv” as example.

-------------------------------------------------------------------------------

"

message("----------------------------")
message("---Start MaChIAto_Aligner---")
message("----------------------------")

debug.flag <- FALSE

# Get directory path
if(debug.flag == TRUE){
    script.dir <- "./MaChIAto_Aligner"
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

library(CrispRVariants)
library(GenomicFeatures)
library(Biostrings)
library(GenomicRanges)
library(Rsamtools)
library(grDevices)
library(network)
library(dplyr)
library(reshape2)
library(stringr)

left.extra.seq <- "" # default
right.extra.seq <- "" # default
if(debug.flag == TRUE){
    # input.fn <- "./MaChIAto_Aligner_testdata/MaChIAto_from_CRISPResso_at_20190505235404_on_AAVS1-A"
    # input.fn <- "./MaChIAto_Aligner_testdata2/MaChIAto_from_CRISPResso_at_20190508102944_on_DBF4B-C"
    # input.fn <- "/Users/kazuki/Documents/Lab/MachinelearningForPITChprediction/20190218_development_of_MAChAto/MaChIAto_Aligner_testdata2/MaChIAto_from_CRISPResso_at_20190508102945_on_ATP5B-D"
    # input.fn <- "./MaChIAto_Aligner_testdata/MaChIAto_from_CRISPResso_at_20190508102946_on_DBF4B-A"
    # input.fn <- "./MaChIAto_Aligner_testdata2/MaChIAto_from_CRISPResso_at_20190509083822_on_C8orf86-D"
    # input.fn <- "./MaChIAto_Aligner_testdata2/MaChIAto_from_CRISPResso_at_20190508102947_on_ATP5B-C"
    # input.fn <- "/Users/kazuki/Documents/Lab/MachinelearningForPITChprediction/20190621_run_MaChIAto_Aligner/MaChIAto_output_190508/MaChIAto_from_CRISPResso_at_20190508102945_on_CYTH3-D"
    # input.fn <- "/Users/kazuki/Documents/Lab/MachinelearningForPITChprediction/20190621_run_MaChIAto_Aligner/MaChIAto_output_190508/MaChIAto_from_CRISPResso_at_20190508102945_on_FBXO10-D"
    # input.fn <- "/Users/kazuki/Documents/Lab/MachinelearningForPITChprediction/20190621_run_MaChIAto_Aligner/MaChIAto_output_190508/MaChIAto_from_CRISPResso_at_20190508102947_on_ATP5B-C"    
    # input.fn <- "/Users/kazuki/Documents/Lab/MachinelearningForPITChprediction/20190621_run_MaChIAto_Aligner/MaChIAto_output_190508/MaChIAto_from_CRISPResso_at_20190513171739_on_SERPINB5-C"
    # input.fn <- "/Users/kazuki/Documents/Lab/MachinelearningForPITChprediction/20190621_run_MaChIAto_Aligner/MaChIAto_output_190508/MaChIAto_from_CRISPResso_at_20190510050514_on_PQLC1-C"
    # input.fn <- "/Users/kazuki/Documents/Lab/MachinelearningForPITChprediction/20190621_run_MaChIAto_Aligner/MaChIAto_output_190508/MaChIAto_from_CRISPResso_at_20190510050512_on_LGR6-D"
    # input.fn <- "/Users/kazuki/Documents/Lab/MachinelearningForPITChprediction/20190621_run_MaChIAto_Aligner/MaChIAto_v1.1.0_output_190710_2019y50data_th0995/MaChIAto_from_CRISPResso_at_20190716174641_on_DBF4B-C"
    input.fn <- "/Users/kazuki/Documents/Lab/MachinelearningForPITChprediction/20190621_run_MaChIAto_Aligner/MaChIAto2_v1.1.0_output_190719/MaChIAto_from_CRISPResso2_at_20190719181211_on_ATP5B-D"
    # output.prefix <- "./MaChIAtoAligner_output2"
    output.prefix <- "/Volumes/databank1/MaChIAtoAligner_output_190716_from_MaChIAto_output_190508"
    left.extra.seq <- "GTTTGG"
    right.extra.seq <- "CCAAAC"
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

#-- Make output directory ----------------------------------------------------#
message("---Make directory for the analysis---")

#time.id <- "20190404150926"
time.id <- format(Sys.time(), "%Y%m%d%H%M%S")
output.name <- paste("MaChIAtoAligner_at_", time.id, "_from_", basename(input.fn), sep = "")
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

# Start logging
# message(paste("Starting process record in", file.path(output.dir, "MaChIAtoAligner_RUNNING_LOG.txt"), sep=" "))
# if(debug.flag == FALSE){
#    log.fh <- file(, open="wt")
#    sink(file.path(output.dir, "MaChIAtoAligner_RUNNING_LOG.txt"), type="message")
# }

###################################################################################################################################################

#-- make mapped bam file by bwa mem ----------------------------------------------------#
message("---Run bwa mem---")

system(paste("chmod +x", file.path(script.dir, "prepare_ref.sh"), sep=" "))
system(paste(file.path(script.dir, "prepare_ref.sh"), input.fn, index.dir, seq.dir , sep=" "))

mmej.ki.seq <- as.character(readDNAStringSet(file.path(index.dir, "mmej_ki.fa")))
donor.seq <- as.character(readDNAStringSet(file.path(seq.dir, "donor.fa")))
left.homology.seq <- as.character(readDNAStringSet(file.path(seq.dir, "/Leftside_homology_arm.fa")))
right.homology.seq <- as.character(readDNAStringSet(file.path(seq.dir, "Rightside_homology_arm.fa")))
nhej.donor.seq <- paste(left.extra.seq, left.homology.seq, donor.seq, right.homology.seq, right.extra.seq, sep="")
nhej.rc.donor.seq <- as.character(reverseComplement(DNAString(nhej.donor.seq)))
nhej.ki.seq <- DNAStringSet(sub(donor.seq, nhej.donor.seq, mmej.ki.seq))
names(nhej.ki.seq) <- c("nhej_ki")
nhej.rc.ki.seq <- DNAStringSet(sub(donor.seq, nhej.rc.donor.seq, mmej.ki.seq))
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
# make mapped bam file of all machiato output (deprecated in recent and more version)
# system(paste("chmod +x", file.path(script.dir, "prepare_bam_analysis_all_outcome.sh"), sep=" "))
# system(paste(file.path(script.dir, "prepare_bam_analysis_all_outcome.sh"), input.fn, index.dir, data.dir, script.dir, sep=" "))
# make mapped bam file of each machiato output
for(data.type in c("Unclassified", "Unmodified", "NHEJ", "Mixed", "HDR")){
  system(paste("chmod +x", file.path(script.dir, "prepare_bam_analysis_each_outcome.sh"), sep=" "))
  system(paste(file.path(script.dir, "prepare_bam_analysis_each_outcome.sh"), input.fn, index.dir, data.dir, script.dir, data.type, sep=" "))
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
# browser()
# ALL
crispresso.all.crispr.set <- MakeWtCrisprSet(
untreated.seq
, "wt"
, wt.seq
, protospacer.seq
, file.path(bam.crispresso.dir, "all.bam")
, "crispresso.all")
if(!is.null(crispresso.all.crispr.set))SaveVariantsData(crispresso.all.crispr.set, crispresso.all.dir, mut.type = "mut")


# quit(save = "no") # for debug

# HDR
crispresso.hdr.crispr.set <- MakeWtCrisprSet(
untreated.seq
, "wt"
, wt.seq
, protospacer.seq
, file.path(bam.crispresso.dir, "hdr.bam")
, "crispresso.hdr")
if(!is.null(crispresso.hdr.crispr.set))SaveVariantsData(crispresso.hdr.crispr.set, crispresso.hdr.dir, mut.type = "mut")

# Mixed_HDR_NHEJ
crispresso.mixed_hdr_nhej.crispr.set <- MakeWtCrisprSet(
untreated.seq
, "wt"
, wt.seq
, protospacer.seq
, file.path(bam.crispresso.dir, "mixed_hdr_nhej.bam")
, "crispresso.mixed_hdr_nhej")
if(!is.null(crispresso.mixed_hdr_nhej.crispr.set))SaveVariantsData(crispresso.mixed_hdr_nhej.crispr.set, crispresso.mixed_hdr_nhej.dir, mut.type = "mut")

# NHEJ
crispresso.nhej.crispr.set <- MakeWtCrisprSet(
untreated.seq
, "wt"
, wt.seq
, protospacer.seq
, file.path(bam.crispresso.dir, "nhej.bam")
, "crispresso.nhej")
if(!is.null(crispresso.nhej.crispr.set))SaveVariantsData(crispresso.nhej.crispr.set, crispresso.nhej.dir, mut.type = "mut")

# Unmodified
crispresso.unmodified.crispr.set <- MakeWtCrisprSet(
untreated.seq
, "wt"
, wt.seq
, protospacer.seq
, file.path(bam.crispresso.dir, "unmodified.bam")
, "crispresso.unmodified")
if(!is.null(crispresso.unmodified.crispr.set))SaveVariantsData(crispresso.unmodified.crispr.set, crispresso.unmodified.dir, mut.type = "mut")

### MaChIAto original output
# ALL
machiato.all.crispr.set <- MakeWtCrisprSet(
untreated.seq
, "wt"
, wt.seq
, protospacer.seq
, file.path(bam.machiato.dir, "all.bam")
, "machiato.all")
if(!is.null(machiato.all.crispr.set))SaveVariantsData(machiato.all.crispr.set, machiato.all.dir, mut.type = "mut")

# HDR
machiato.hdr.crispr.set <- MakeWtCrisprSet(
untreated.seq
, "wt"
, wt.seq
, protospacer.seq
, file.path(bam.machiato.dir, "hdr.bam")
, "machiato.hdr")
if(!is.null(machiato.hdr.crispr.set))SaveVariantsData(machiato.hdr.crispr.set, machiato.hdr.dir, mut.type = "mut")

# Mixed_HDR_NHEJ
machiato.mixed_hdr_nhej.crispr.set <- MakeWtCrisprSet(
untreated.seq
, "wt"
, wt.seq
, protospacer.seq
, file.path(bam.machiato.dir, "mixed_hdr_nhej.bam")
, "machiato.mixed_hdr_nhej")
if(!is.null(machiato.mixed_hdr_nhej.crispr.set))SaveVariantsData(machiato.mixed_hdr_nhej.crispr.set, machiato.mixed_hdr_nhej.dir, mut.type = "mut")

# NHEJ
machiato.nhej.crispr.set <- MakeWtCrisprSet(
untreated.seq
, "wt"
, wt.seq
, protospacer.seq
, file.path(bam.machiato.dir, "nhej.bam")
, "machiato.nhej")
if(!is.null(machiato.nhej.crispr.set))SaveVariantsData(machiato.nhej.crispr.set, machiato.nhej.dir, mut.type = "mut")

# Unmodified
machiato.unmodified.crispr.set <- MakeWtCrisprSet(
untreated.seq
, "wt"
, wt.seq
, protospacer.seq
, file.path(bam.machiato.dir, "unmodified.bam")
, "machiato.unmodified")
if(!is.null(machiato.unmodified.crispr.set))SaveVariantsData(machiato.unmodified.crispr.set, machiato.unmodified.dir, mut.type = "mut")

# Unclassified
machiato.unclassified.crispr.set <- MakeWtCrisprSet(
untreated.seq
, "wt"
, wt.seq
, protospacer.seq
, file.path(bam.machiato.dir, "unclassified.bam")
, "machiato.unclassified")
if(!is.null(machiato.unclassified.crispr.set))SaveVariantsData(machiato.unclassified.crispr.set, machiato.unclassified.dir, mut.type = "mut")

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
  
  write.table(crispresso.class.count.table, file = file.path(crispresso.rate.dir, "[ⅰ]CRISPResso_Unmodified_NHEJ_HDR_table.txt") # TODO : add ID number
    , quote=FALSE, col.names=FALSE, row.names=FALSE, append=TRUE, sep = ",")
  saveRDS(crispresso.class.count.table, file = file.path(crispresso.rate.dir, "[ⅰ]CRISPResso_Unmodified_NHEJ_HDR_table.rds"))
  SavePieChart(crispresso.class.count.table
    , c("#D98F4E", "#A65437", "#F2CA52", "#F2E2CE")
    , "Reads"
    , file.path(crispresso.rate.dir, "[ⅰ]CRISPResso_Unmodified_NHEJ_HDR_pie_chart.png")
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
  
  write.table(machiato.including.unclassified.class.count.table, file = file.path(machiato.rate.dir, "[ⅰ]CRISPResso-MaChIAto_Unmodified_NHEJ_HDR_table(including_Unclassified).txt")
    , quote=FALSE, col.names=FALSE, row.names=FALSE,append=TRUE, sep = ",")
  saveRDS(machiato.including.unclassified.class.count.table, file = file.path(machiato.rate.dir, "[ⅰ]CRISPResso-MaChIAto_Unmodified_NHEJ_HDR_table(including_Unclassified).rds"))
  SavePieChart(machiato.including.unclassified.class.count.table
  , c("#D98F4E",  "#A65437", "#F2CA52", "#ADD4D9", "#F2E2CE")
  , "Reads"
  , file.path(machiato.rate.dir, "[ⅰ]CRISPResso-MaChIAto_Unmodified_NHEJ_HDR_pie_chart(including_Unclassified).png")
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
  
  write.table(machiato.class.count.table, file = file.path(machiato.rate.dir, "[ⅱ]CRISPResso-MaChIAto_Unmodified_NHEJ_HDR_table.txt")
    , quote=FALSE, col.names=FALSE, row.names=FALSE,append=TRUE, sep = ",")
  saveRDS(machiato.class.count.table, file = file.path(machiato.rate.dir, "[ⅱ]CRISPResso-MaChIAto_Unmodified_NHEJ_HDR_table.rds"))
  SavePieChart(machiato.class.count.table
  , c("#D98F4E", "#A65437", "#F2CA52", "#F2E2CE")
  , "Reads"
  , file.path(machiato.rate.dir, "[ⅱ]CRISPResso-MaChIAto_Unmodified_NHEJ_HDR_pie_chart.png")
  )

  # We have to correct KL divergence calculation function. but this neccesity is not high. so it is suspended.
  # message("---Run comparison analysis between original CRISPResso data and CRISPResso-MaChIAto data---")
  #
  # #CRESPResso_datatype.MaChIAto_datatype.kl
  # 
  # # Prepare matrix contains each KL divergence
  # crispresso.crispr.set.list <- list(
  #   ALL=crispresso.all.crispr.set
  #   , Unmodified=crispresso.unmodified.crispr.set
  #   , NHEJ=crispresso.nhej.crispr.set
  #   , Mixed_HDR_NHEJ=crispresso.mixed_hdr_nhej.crispr.set
  #   , HDR=crispresso.hdr.crispr.set
  # )
  # machiato.crispr.set.list <- list(
  #   ALL=machiato.all.crispr.set
  #   , Unclassified=machiato.unclassified.crispr.set
  #   , Unmodified=machiato.unmodified.crispr.set
  #   , NHEJ=machiato.nhej.crispr.set
  #   , Mixed_HDR_NHEJ=machiato.mixed_hdr_nhej.crispr.set
  #   , HDR=machiato.hdr.crispr.set
  # )
  # crispresso.class.names <- c("ALL", "Unmodified", "NHEJ", "Mixed_HDR_NHEJ", "HDR")
  # machiato.class.names <- rev(c("ALL", "Unmodified", "NHEJ", "Mixed_HDR_NHEJ", "HDR", "Unclassified"))
  # crispresso.machiato.class.kl.matrix <- matrix(
  #   NA
  #   , ncol=length(machiato.class.names)
  #   , nrow=length(crispresso.class.names)
  # )
  # colnames(crispresso.machiato.class.kl.matrix) = machiato.class.names
  # rownames(crispresso.machiato.class.kl.matrix) = crispresso.class.names
  # # Calc KL divergence
  # for(crispresso.datatype in crispresso.class.names){
  #   for(machiato.datatype in machiato.class.names){
  #     crispresso.machiato.class.kl.matrix[crispresso.datatype, machiato.datatype] <- CalcKldivergence(
  #       getElement(machiato.crispr.set.list, machiato.datatype)
  #       , getElement(crispresso.crispr.set.list, crispresso.datatype)
  #     )
  #   }
  # }
  # saveRDS(crispresso.machiato.class.kl.matrix, file = file.path(comparison.analysis.data.dir, "[ⅰ]crispresso.machiato.class.kl.matrix.rds"))
  # # Make heatmap
  # options(warn=-1)
  # crispresso.machiato.class.kl.p <- ggplot(melt(crispresso.machiato.class.kl.matrix)
  #   , aes(x=Var1, y=Var2, fill=value))+
  # geom_raster(aes(fill = value))+
  # labs(title ="Overall Variants Analysis"
  #   , x = "Original CRISPResso"
  #   , y = "CRISPResso-MaChIAto")+
  # theme(axis.ticks = element_blank()
  #   , axis.text.x = element_text(angle = 330, hjust = 0, size=15)
  #   , axis.text.y = element_text(size=15)
  #       , plot.title = element_text(size = 20),  
  #   , axis.title.x = element_text(size = 20),  
  #   , axis.title.y = element_text(size = 20)
  # )+
  # scale_fill_gradient(low = "#F28907", high = "#261201", name = "KL Divergence")
  # ggsave(file = file.path(comparison.analysis.data.dir, "[ⅰ]Overall_Variants_Analysis_heatmap.png")
  #   , plot = crispresso.machiato.class.kl.p, dpi = 300, width = 10, height = 10)
  # options(warn=0)
  # message(paste("Save heatmap in ", file.path(comparison.analysis.data.dir, "[ⅰ]Overall_Variants_Analysis_heatmap.png"), sep=""))

}
############################################################################

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

##################
# save data
# list(untreated.variants.count.data = untreated.variants.count.data
# ,overall_edited_substitution.variants.count.data = overall_edited_substitution.variants.count.data
# ,overall_edited_deletion.variants.count.data = overall_edited_deletion.variants.count.data
# ,overall_edited_insertion.variants.count.data = overall_edited_insertion.variants.count.data)


# End logging
# if(debug.flag == FALSE){
#    sink(type="message")
#    close(log.fh)
# }
# message(paste("Ending process record in", output.dir, sep=" "))

message("---All process is completed.---")
quit(save = "no")