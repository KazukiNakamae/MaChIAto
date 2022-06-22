# /bin/bash

function message
{
    printf "\e[34m$1\n\e[m"
}

message "Search each position of sequence using blat.";
cd ./blat/bin;
./blat ../genome/hg19.fa ../fa/$1 ../psl/$1.psl -minMatch=0 -minScore=$2 -noHead;
message "The result was saved as $1.psl";
cd ..;
# 事前にBEDOPSをインストールし、パスを通しておく必要がある。BEDOPS（pre-built packages）：http://bedops.readthedocs.io/en/latest/content/installation.html#installation
message "Convert .psl into .bed";
mkdir ./bed
/usr/local/bin/psl2bed --do-not-sort < ./psl/$1.psl > ./bed/$1.psl.bed;
message "The result was saved as $1.psl.bed";