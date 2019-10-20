#!/bin/bash

INDIR=$1;
INDEXDIR=$2;
DATADIR=$3;
SCRIPTDIR=$4;
DATATYPE=$5;
#if [ $# -eq 6 ]
#  then
#    DATATYPE2=$6
#fi

# Mutation
awk -F, 'BEGIN{OFS=","}{ gsub ("Mixed HDR-NHEJ", "Mixed", $76); if ($72 == "unidentified" && $76 == class) {for (i=1; i<=$10; i++) print ">"$1"\n"$2} }' class=$DATATYPE $INDIR"/ALL_dataframe.csv" > $DATADIR"/fasta/"$DATATYPE"/other_type/unidentified.fa";
for type in "untreated" "unexpected_mutations" "overall_edited_substitution" "overall_edited_insertion" "overall_edited_deletion";
  do awk -F, 'BEGIN{OFS=","}{ gsub ("Mixed HDR-NHEJ", "Mixed", $76); if ($72 == type && $76 == class) {for (i=1; i<=$10; i++) print ">"$1"\n"$12} }' type=$type class=$DATATYPE $INDIR"/ALL_dataframe.csv" > $DATADIR"/fasta/"$DATATYPE"/mut_type/"$type".fa";
done;
# We don't analize donor_insertion label because samples which has this label is analized in knock-in analysis below.

# whole knock-in
for type in "LEFT_PRECISE_KNOCK_IN" "LEFT_IMPRECISE_KNOCK_IN" "LEFT_IMPRECISE_KNOCK_IN_COMPLEX" "LEFT_OTHER_MUATIONS";
  do awk -F, 'BEGIN{OFS=","}{ gsub ("Mixed HDR-NHEJ", "Mixed", $76); if ($73 == type && $60 == "False" && $76 == class) {for (i=1; i<=$10; i++) print ">"$1"\n"$12;} }' type=$type class=$DATATYPE $INDIR"/ALL_dataframe.csv" > $DATADIR"/fasta/"$DATATYPE"/ki_type/"$type".fa";
done;
for type in "LEFT_IMPRECISE_KNOCK_IN" "LEFT_IMPRECISE_KNOCK_IN_COMPLEX" "LEFT_OTHER_MUATIONS";
  do awk -F, 'BEGIN{OFS=","}{ gsub ("Mixed HDR-NHEJ", "Mixed", $76); if ($73 == type && $60 == "True" && $76 == class) {for (i=1; i<=$10; i++) print ">"$1"\n"$12;} }' type=$type class=$DATATYPE $INDIR"/ALL_dataframe.csv" > $DATADIR"/fasta/"$DATATYPE"/rc_ki_type/rc_"$type".fa";
done;
for type in "RIGHT_PRECISE_KNOCK_IN" "RIGHT_IMPRECISE_KNOCK_IN" "RIGHT_IMPRECISE_KNOCK_IN_COMPLEX" "RIGHT_OTHER_MUATIONS";
  do awk -F, 'BEGIN{OFS=","}{ gsub ("Mixed HDR-NHEJ", "Mixed", $76); if ($74 == type && $64 == "False" && $76 == class) {for (i=1; i<=$10; i++) print ">"$1"\n"$12} }' type=$type class=$DATATYPE $INDIR"/ALL_dataframe.csv" > $DATADIR"/fasta/"$DATATYPE"/ki_type/"$type".fa";
done;
for type in "RIGHT_IMPRECISE_KNOCK_IN" "RIGHT_IMPRECISE_KNOCK_IN_COMPLEX" "RIGHT_OTHER_MUATIONS";
  do awk -F, 'BEGIN{OFS=","}{ gsub ("Mixed HDR-NHEJ", "Mixed", $76); if ($74 == type && $64 == "True" && $76 == class) {for (i=1; i<=$10; i++) print ">"$1"\n"$12} }' type=$type class=$DATATYPE $INDIR"/ALL_dataframe.csv" > $DATADIR"/fasta/"$DATATYPE"/rc_ki_type/rc_"$type".fa";
done;


# convert to fake fastq
for dir in $(ls -1 $DATADIR"/fasta/"$DATATYPE | sed -e 's/\..*$//');
  do for file in $(ls -1 $DATADIR"/fasta/"$DATATYPE"/"$dir | sed -e 's/\..*$//');
    do python $SCRIPTDIR"/make_fake_fq.py" $DATADIR"/fasta/"$DATATYPE"/"$dir"/"$file".fa" > $DATADIR"/fake_fastq/"$DATATYPE"/"$dir"/"$file".fq";
  done;
done;

# bwa mem mapping
for dir in "mut_type" "other_type";
  do for file in $(ls -1 $DATADIR"/fake_fastq/"$DATATYPE"/"$dir | sed -e 's/\..*$//');
    do bwa mem -t 1 $INDEXDIR"/wt" $DATADIR"/fake_fastq/"$DATATYPE"/"$dir"/"$file".fq" | samtools sort -o $DATADIR"/bam/"$DATATYPE"/"$dir"/"$file".bam";
  done;
done;

for file in "LEFT_PRECISE_KNOCK_IN" "RIGHT_PRECISE_KNOCK_IN";
  do bwa mem -t 1 $INDEXDIR"/mmej_ki" $DATADIR"/fake_fastq/"$DATATYPE"/ki_type/"$file".fq" | samtools sort -o $DATADIR"/bam/"$DATATYPE"/ki_type/"$file".bam";
done;

for file in "LEFT_IMPRECISE_KNOCK_IN" "LEFT_IMPRECISE_KNOCK_IN_COMPLEX" "LEFT_OTHER_MUATIONS" "RIGHT_IMPRECISE_KNOCK_IN" "RIGHT_IMPRECISE_KNOCK_IN_COMPLEX" "RIGHT_OTHER_MUATIONS";
  do bwa mem -t 1 $INDEXDIR"/nhej_ki" $DATADIR"/fake_fastq/"$DATATYPE"/ki_type/"$file".fq" | samtools sort -o $DATADIR"/bam/"$DATATYPE"/ki_type/"$file".bam";
done;

for file in "rc_LEFT_IMPRECISE_KNOCK_IN" "rc_LEFT_IMPRECISE_KNOCK_IN_COMPLEX" "rc_LEFT_OTHER_MUATIONS" "rc_RIGHT_IMPRECISE_KNOCK_IN" "rc_RIGHT_IMPRECISE_KNOCK_IN_COMPLEX" "rc_RIGHT_OTHER_MUATIONS";
  do bwa mem -t 1 $INDEXDIR"/nhej_rc_ki" $DATADIR"/fake_fastq/"$DATATYPE"/rc_ki_type/"$file".fq" | samtools sort -o $DATADIR"/bam/"$DATATYPE"/rc_ki_type/"$file".bam";
done;

