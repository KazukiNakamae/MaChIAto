#!/bin/bash

INDEXDIR=$1;

# make index
bwa index -p $INDEXDIR"/wt" $INDEXDIR"/wt.fa";
bwa index -p $INDEXDIR"/mmej_ki" $INDEXDIR"/mmej_ki.fa";
bwa index -p $INDEXDIR"/nhej_ki" $INDEXDIR"/nhej_ki.fa";
bwa index -p $INDEXDIR"/nhej_rc_ki" $INDEXDIR"/nhej_rc_ki.fa";
