#!/bin/bash

INDIR=$1;
INDEXDIR=$2;
DATADIR=$3;
SCRIPTDIR=$4;

# CRISPResso
awk -F, 'BEGIN{OFS=","}{for (i=1; i<=$10; i++) {if ($10==0 || substr($10,1,1) !~ /^[0-9]/) {print "";exit;} else print ">"$1"\n"$2;}}' $INDIR"/ALL_dataframe.csv" > $DATADIR"/fasta/crispresso/all.fa";
awk -F, 'BEGIN{OFS=","}{for (i=1; i<=$10; i++) {if ($10==0 || substr($10,1,1) !~ /^[0-9]/) {print "";exit;} else print ">"$1"\n"$2;}}' $INDIR"/CRISPResso_HDR_dataframe.csv" > $DATADIR"/fasta/crispresso/hdr.fa";
awk -F, 'BEGIN{OFS=","}{for (i=1; i<=$10; i++) {if ($10==0 || substr($10,1,1) !~ /^[0-9]/) {print "";exit;} else print ">"$1"\n"$2;}}' $INDIR"/CRISPResso_Mixed_HDR_NHEJ_dataframe.csv" > $DATADIR"/fasta/crispresso/mixed_hdr_nhej.fa";
awk -F, 'BEGIN{OFS=","}{for (i=1; i<=$10; i++) {if ($10==0 || substr($10,1,1) !~ /^[0-9]/) {print "";exit;} else print ">"$1"\n"$2;}}' $INDIR"/CRISPResso_NHEJ_dataframe.csv" > $DATADIR"/fasta/crispresso/nhej.fa";
awk -F, 'BEGIN{OFS=","}{for (i=1; i<=$10; i++) {if ($10==0 || substr($10,1,1) !~ /^[0-9]/) {print "";exit;} else print ">"$1"\n"$2;}}' $INDIR"/CRISPResso_UNMODIFIED_dataframe.csv" > $DATADIR"/fasta/crispresso/unmodified.fa";

# MaChIAto
awk -F, 'BEGIN{OFS=","}{ for (i=1; i<=$10; i++) {if ($10==0 || substr($10,1,1) !~ /^[0-9]/) {print "";exit;} else print ">"$1"\n"$2;}}' $INDIR"/ALL_dataframe.csv" > $DATADIR"/fasta/machiato/all.fa";
awk -F, 'BEGIN{OFS=","}{ if ($77 == "HDR") {for (i=1; i<=$10; i++) print ">"$1"\n"$2} }' $INDIR"/ALL_dataframe.csv" > $DATADIR"/fasta/machiato/hdr.fa";
awk -F, 'BEGIN{OFS=","}{ if ($77 == "NHEJ") {for (i=1; i<=$10; i++) print ">"$1"\n"$2} }' $INDIR"/ALL_dataframe.csv" > $DATADIR"/fasta/machiato/nhej.fa";
awk -F, 'BEGIN{OFS=","}{ if ($77 == "Unmodified") {for (i=1; i<=$10; i++) print ">"$1"\n"$2} }' $INDIR"/ALL_dataframe.csv" > $DATADIR"/fasta/machiato/unmodified.fa";
awk -F, 'BEGIN{OFS=","}{ if ($77 == "Unclassified") {for (i=1; i<=$10; i++) print ">"$1"\n"$2} }' $INDIR"/ALL_dataframe.csv" > $DATADIR"/fasta/machiato/unclassified.fa";
awk -F, 'BEGIN{OFS=","}{ gsub ("Mixed HDR-NHEJ", "Mixed", $77); if ($77 == "Mixed") {for (i=1; i<=$10; i++) print ">"$1"\n"$2} }' $INDIR"/ALL_dataframe.csv" > $DATADIR"/fasta/machiato/mixed_hdr_nhej.fa";

# convert to fake fastq
for dir in $(ls -1 $DATADIR"/fasta" | sed -e 's/\..*$//');
  do for file in $(ls -1 $DATADIR"/fasta/"$dir | sed -e 's/\..*$//');
    do python $SCRIPTDIR"/make_fake_fq.py" $DATADIR"/fasta/"$dir"/"$file".fa" $DATADIR"/fake_fastq/"$dir"/"$file".fq";
  done;
done;

# bwa mem mapping
for dir in "crispresso" "machiato";
  do for file in $(ls -1 $DATADIR"/fake_fastq/"$dir | sed -e 's/\..*$//');
    do echo $DATADIR"/fake_fastq/"$dir"/"$file".fq";bwa mem -t 1 $INDEXDIR"/wt" $DATADIR"/fake_fastq/"$dir"/"$file".fq" | samtools view -h | grep -v -e 'XA:Z:' -e 'SA:Z:' | samtools view -b | samtools sort -o $DATADIR"/bam/"$dir"/"$file".bam";# samtools index $DATADIR"/bam/"$dir"/"$file".bam";
  done;
done;

# export from MaChIAto_optimized_param.csv
#awk -F, 'BEGIN{OFS=","}{ if (NR>1) {print ">sgRNA\n"$2} }' $INDIR"/MaChIAto_optimized_param.csv" > $DATADIR"/seq/sgRNA.fa";
# export from MaChIAto_specific_sequences.csv
#awk -F, 'BEGIN{OFS=","}{ if (NR>1) {print ">In-Left_indicator_sequence\n"$1} }' $INDIR"/MaChIAto_specific_sequences.csv" > $DATADIR"/seq/In-Left_indicator_sequence.fa";
#awk -F, 'BEGIN{OFS=","}{ if (NR>1) {print ">In-Right_indicator_sequence\n"$2} }' $INDIR"/MaChIAto_specific_sequences.csv" > $DATADIR"/seq/In-Right_indicator_sequence.fa";
#awk -F, 'BEGIN{OFS=","}{ if (NR>1) {print ">Out-Left_indicator_sequence\n"$7} }' $INDIR"/MaChIAto_specific_sequences.csv" > $DATADIR"/seq/Out-Left_indicator_sequence.fa";
#awk -F, 'BEGIN{OFS=","}{ if (NR>1) {print ">Out-Right_indicator_sequence\n"$8} }' $INDIR"/MaChIAto_specific_sequences.csv" > $DATADIR"/seq/Out-Right_indicator_sequence.fa";
#awk -F, 'BEGIN{OFS=","}{ if (NR>1) {print ">Untreated_sequence\n"$13} }' $INDIR"/MaChIAto_specific_sequences.csv" > $DATADIR"/seq/Untreated_sequence.fa";

echo "done";