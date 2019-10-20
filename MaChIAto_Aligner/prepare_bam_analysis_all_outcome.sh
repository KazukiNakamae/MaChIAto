#!/bin/bash

INDIR=$1;
INDEXDIR=$2;
DATADIR=$3;
SCRIPTDIR=$4;

# Mutation
awk -F, 'BEGIN{OFS=","}{ if ($72 == "unidentified") {for (i=1; i<=$10; i++) print ">"$1"\n"$2} }' $INDIR"/ALL_dataframe.csv" > $DATADIR"/fasta/ALL/other_type/unidentified.fa";
for type in "untreated" "donor_insertion" "unexpected_mutations" "overall_edited_substitution" "overall_edited_insertion" "overall_edited_deletion";
  do awk -F, 'BEGIN{OFS=","}{ if ($72 == type) {for (i=1; i<=$10; i++) print ">"$1"\n"$12} }' type=$type $INDIR"/ALL_dataframe.csv" > $DATADIR"/fasta/ALL/mut_type/"$type".fa";
done;

# whole knock-in
for type in "LEFT_PRECISE_KNOCK_IN" "LEFT_IMPRECISE_KNOCK_IN" "LEFT_IMPRECISE_KNOCK_IN_COMPLEX" "LEFT_OTHER_MUATIONS";
  do awk -F, 'BEGIN{OFS=","}{ if ($73 == type && $60 == "False") {for (i=1; i<=$10; i++) print ">"$1"\n"$12;} }' type=$type $INDIR"/ALL_dataframe.csv" > $DATADIR"/fasta/ALL/ki_type/"$type".fa";
done;
for type in "LEFT_IMPRECISE_KNOCK_IN" "LEFT_IMPRECISE_KNOCK_IN_COMPLEX" "LEFT_OTHER_MUATIONS";
  do awk -F, 'BEGIN{OFS=","}{ if ($73 == type && $60 == "True") {for (i=1; i<=$10; i++) print ">"$1"\n"$12;} }' type=$type $INDIR"/ALL_dataframe.csv" > $DATADIR"/fasta/ALL/rc_ki_type/rc_"$type".fa";
done;
for type in "RIGHT_PRECISE_KNOCK_IN" "RIGHT_IMPRECISE_KNOCK_IN" "RIGHT_IMPRECISE_KNOCK_IN_COMPLEX" "RIGHT_OTHER_MUATIONS";
  do awk -F, 'BEGIN{OFS=","}{ if ($74 == type && $64 == "False") {for (i=1; i<=$10; i++) print ">"$1"\n"$12} }' type=$type $INDIR"/ALL_dataframe.csv" > $DATADIR"/fasta/ALL/ki_type/"$type".fa";
done;
for type in "RIGHT_IMPRECISE_KNOCK_IN" "RIGHT_IMPRECISE_KNOCK_IN_COMPLEX" "RIGHT_OTHER_MUATIONS";
  do awk -F, 'BEGIN{OFS=","}{ if ($74 == type && $64 == "True") {for (i=1; i<=$10; i++) print ">"$1"\n"$12} }' type=$type $INDIR"/ALL_dataframe.csv" > $DATADIR"/fasta/ALL/rc_ki_type/rc_"$type".fa";
done;

# convert to fake fastq
for dir in $(ls -1 $DATADIR"/fasta/ALL" | sed -e 's/\..*$//');
  do for file in $(ls -1 $DATADIR"/fasta/ALL/"$dir | sed -e 's/\..*$//');
    do python $SCRIPTDIR"/make_fake_fq.py" $DATADIR"/fasta/ALL/"$dir"/"$file".fa" > $DATADIR"/fake_fastq/ALL/"$dir"/"$file".fq";
  done;
done;

# bwa mem mapping
for dir in "mut_type" "other_type";
  do for file in $(ls -1 $DATADIR"/fake_fastq/ALL/"$dir | sed -e 's/\..*$//');
    do bwa mem $INDEXDIR"/wt" $DATADIR"/fake_fastq/ALL/"$dir"/"$file".fq" | samtools sort -o $DATADIR"/bam/ALL/"$dir"/"$file".bam";
  done;
done;

for file in "LEFT_PRECISE_KNOCK_IN" "RIGHT_PRECISE_KNOCK_IN";
  do bwa mem $INDEXDIR"/mmej_ki" $DATADIR"/fake_fastq/ALL/ki_type/"$file".fq" | samtools sort -o $DATADIR"/bam/ALL/ki_type/"$file".bam";
done;

for file in "LEFT_IMPRECISE_KNOCK_IN" "LEFT_IMPRECISE_KNOCK_IN_COMPLEX" "LEFT_OTHER_MUATIONS" "RIGHT_IMPRECISE_KNOCK_IN" "RIGHT_IMPRECISE_KNOCK_IN_COMPLEX" "RIGHT_OTHER_MUATIONS";
  do bwa mem $INDEXDIR"/nhej_ki" $DATADIR"/fake_fastq/ALL/ki_type/"$file".fq" | samtools sort -o $DATADIR"/bam/ALL/ki_type/"$file".bam";
done;

for file in "rc_LEFT_IMPRECISE_KNOCK_IN" "rc_LEFT_IMPRECISE_KNOCK_IN_COMPLEX" "rc_LEFT_OTHER_MUATIONS" "rc_RIGHT_IMPRECISE_KNOCK_IN" "rc_RIGHT_IMPRECISE_KNOCK_IN_COMPLEX" "rc_RIGHT_OTHER_MUATIONS";
  do bwa mem $INDEXDIR"/nhej_rc_ki" $DATADIR"/fake_fastq/ALL/rc_ki_type/"$file".fq" | samtools sort -o $DATADIR"/bam/ALL/rc_ki_type/"$file".bam";
done;

