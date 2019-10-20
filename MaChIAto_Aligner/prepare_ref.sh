#!/bin/bash

INDIR=$1;
INDEXDIR=$2;
SEQDIR=$3;

# extract reference sequence
awk -F, 'BEGIN{OFS=","}{ if ($5 <= 1) {print ">wt\n"$1} }' $INDIR"/MaChiIAto_optimized_param.csv" > $INDEXDIR"/wt.fa";
awk -F, 'BEGIN{OFS=","}{ if ($5 <= 1) {print ">mmej_ki\n"$3} }' $INDIR"/MaChiIAto_optimized_param.csv" > $INDEXDIR"/mmej_ki.fa";
awk -F, 'BEGIN{OFS=","}{ if (NR>1) {print ">donor\n"$2} }' $INDIR"/MaChiIAto_optimized_param.csv" > $SEQDIR"/donor.fa";
awk -F, 'BEGIN{OFS=","}{ if (NR>1) {print ">Leftside_homology_arm\n"$5} }' $INDIR"/MaChiIAto_specific_sequences.csv" > $SEQDIR"/Leftside_homology_arm.fa";
awk -F, 'BEGIN{OFS=","}{ if (NR>1) {print ">Rightside_homology_arm\n"$11} }' $INDIR"/MaChiIAto_specific_sequences.csv" > $SEQDIR"/Rightside_homology_arm.fa";
