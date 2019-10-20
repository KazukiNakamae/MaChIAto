#!/bin/bash

INDIR=$1;
DATADIR=$2;

# export from MaChiIAto_optimized_param.csv
awk -F, 'BEGIN{OFS=","}{ if (NR>1) {print ">sgRNA\n"$4} }' $INDIR"/MaChiIAto_optimized_param.csv" > $DATADIR"/seq/sgRNA.fa";
# export from MaChiIAto_specific_sequences.csv
awk -F, 'BEGIN{OFS=","}{ if (NR>1) {print ">In-Left_indicator_sequence\n"$1} }' $INDIR"/MaChiIAto_specific_sequences.csv" > $DATADIR"/seq/In-Left_indicator_sequence.fa";
awk -F, 'BEGIN{OFS=","}{ if (NR>1) {print ">In-Right_indicator_sequence\n"$2} }' $INDIR"/MaChiIAto_specific_sequences.csv" > $DATADIR"/seq/In-Right_indicator_sequence.fa";
awk -F, 'BEGIN{OFS=","}{ if (NR>1) {print ">Out-Left_indicator_sequence\n"$7} }' $INDIR"/MaChiIAto_specific_sequences.csv" > $DATADIR"/seq/Out-Left_indicator_sequence.fa";
awk -F, 'BEGIN{OFS=","}{ if (NR>1) {print ">Out-Right_indicator_sequence\n"$8} }' $INDIR"/MaChiIAto_specific_sequences.csv" > $DATADIR"/seq/Out-Right_indicator_sequence.fa";
awk -F, 'BEGIN{OFS=","}{ if (NR>1) {print ">Untreated_sequence\n"$13} }' $INDIR"/MaChiIAto_specific_sequences.csv" > $DATADIR"/seq/Untreated_sequence.fa";
