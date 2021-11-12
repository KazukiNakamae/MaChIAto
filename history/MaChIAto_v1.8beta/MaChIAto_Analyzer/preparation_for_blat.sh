# /bin/bash

function message
{
    printf "\e[34m$1\n\e[m"
}

message "Make directories";
mkdir ./blat;
mkdir ./blat/bin;
mkdir ./blat/genome;
mkdir ./blat/genome/temp;
mkdir ./blat/psl;
mkdir ./blat/fa;

message "Download hg19 (Human Genome version 19)";
for ch in `seq 1 22` X Y;
    do wget -P ./blat/genome/temp http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr$ch.fa.gz;
done;

message "Unzip .gz file";
gunzip -r ./blat/genome/temp;

message "Merge .fa files to make hg19.fa";
cat ./blat/genome/temp/chr*.fa > ./blat/genome/hg19.fa
message "The hg19.fa was saved in ./test/";

message "Remove temp directory";
rm -rf ./blat/genome/temp;

message "Download blat";
if [ "$(uname)" == 'Darwin' ]; then
    message "blat for Mac OSX";
    wget -P ./blat/bin http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/blat/blat;
    chmod +x ./blat/bin/blat;
elif [ "$(expr substr $(uname -s) 1 5)" == 'Linux' ]; then
    message "blat for Linux";
    wget -P ./bin http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/blat/blat;
    chmod +x ./blat/bin/blat;
else
  message "This environment is not supported in this script. Please make sure whether blat is avariable on your system in UCSC web site.";
  message "UCSC web site : https://genome.ucsc.edu/FAQ/FAQblat.html";
  message "If your system is Mac OSX or Linux, Please check value of $(uname). I guess that the variable is not unusual.";
  echo '$(uname)';
  echo $(uname);
  exit 1;
fi

message "Done.";
