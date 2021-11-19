# MaChIAto
##### LATEST VERSION (2019/10/20)

<img src="logo.png" alt="MaChIAto_logo" title="MaChIAto_logo" width="200" height="200">

**MaChIAto** (**M**MEJ-**a**ssisted **Ch**rommosomal **I**ntegration **A**nalysis **to**ols) is a third-party tool of CRISPResso and CRISPResso2.
MaChIAto can perform stricter re-classification of the output of CRISPResso(2). Also, The tools accommodate other functions for profiling a tendency of knock-in reads.

MaChIAto is designed to analyze amplicon sequence including HDR- or MMEJ- assisted knock-in reads. However, the tools are also available for knock-out analysis.

If you want to know more detailed information and instruction, please check our paper (coming soon!). Or, don't hesitate to contact me (kazukinakamae@gmail.com).

---

### Version History

###　[v1.9beta] - 

##### added
- all data that is used for plot is exported to the result in MaChIAto Aligner
- the classification for frameshift variants is available in MaChIAto Aligner

##### fixed
- the most of plot data is exported as .csv format in MaChIAto Aligner

###　[v1.8beta] - 

##### added
- recognition of "Scaffold-incorporated" in Alleles_frequency_table.txt in MaChIAto Classifier.

##### fixed
- table with a certain amplicon name can be recognized when CRSIRPEsso2 format is processed in MaChIAto Classifier.
- some optimization parameters is calculated as float variables in MaChIAto Classifier.
- result of correlation test can be put out in MaChIAto Analyzer.

##### changed
- hybrid-ss-min (OligoArrayAux subset) is installed into conda environment. The path is overwritten in CalcThermo-class.R of MaChIAto Analyzer.

###　[v1.7beta] - 2020-07-27

##### added
- MaChIAto Classifier adapted the format of CRISPResso2 Prime Editing mode.

##### fixed
- Some errors due to the value of zero and Inf in MaChIAto Reviewer.

##### changed
- In knock-out analysis, changes of precise knock-in and impricise knock-in becomes invisible on the box charts in MaChIAto Reviewer.


###　[v1.6beta] - 2020-06-30

##### fixed
- plot window of distribution is not defined in MaChIAto Reviewer.

##### changed
- The class label that has no read becomes invisible on the pie charts in MaChIAto.

###　[v1.5beta] - 2020-06-14

##### fixed
- missing samples included in ignore_list in MaChIAto Reviewer

##### Changed
- installation of oligoarrayaux-3.8 for calculation of MEF in MaChIAto Analyzer 
- order of plyr and dplyr when loding library in MaChIAto Analyzer 
- summarize_each() replaced into summarize(across) in MaChIAto Reviewer

##### Added
- boxplot of each fearure in MaChIAto Analyzer 

###　[v1.5beta] - 2020-05-26

##### Changed
- The script of MaChIAto Analyzer is modified to perform calculation of each rate and sample set that user want to observe.
- The sgRNA and edited homology arm can be used as a calculation target on MaChIAto Analyzer.
- MaChIAto Analyzer can work without the extra table.
- MaChIAto Analyzer and Reviewer can analyze single knock-in sample and single/double knock-out sample.
- The sample with less classified reads will be filtered.

###　[v1.4beta] - 2020-02-25

CrispRVariants instalation with conda did not work. I solved this problem by changing download function.

#### MaChIAto Classifier version 1.7.1
##### Changed
- The input sequecence is converted into upper letter in the parameter checking.

#### MaChIAto Aligner version 1.4
##### Changed
- Pathnames given by user become processed into string with no end-slash
- All R package depends on package intallation program.

###　[1.0alpha-v1.3beta] - 2020-01-30

I developed Prime Editing analysis mode.

#### MaChIAto (Classifer) version 1.7.0
##### Added
- Prime editing analysis mode
- Detection of ambiguous reads
- Piechart of classification
##### Changed
- Default paramters was adjusted for new mode
- Default µHs lengh was changed into 20 bp
- Order of column in xxx_dataframe.csv

#### MaChIAto Aligner version 1.3
##### Changed
- Method to get sequence in shell scripts

#### MaChIAto Analyzer version beta.1.5
##### Changed
- Filter for samples include NA in the value set is removed because it also lose the available values
##### Deleted
- Extra comments

#### MaChIAto Reviewer version beta2.3
##### Deleted
- Codes that only works in local environment of a developer

#### Update commands

```bash
cd /Users/nedo01/nakamae_github;
# git clone https://github.com/KazukiNakamae/MaChIAto.git;
cd MaChIAto;
# 更新ファイルの登録
git add --all;
# コミット
git commit -m "Upload the latest version";
# リモートリポジトリの設定
git remote add origin https://github.com//KazukiNakamae/MaChIAto.git;
# mainへプッシュ
git push origin master;
```
https://github.com/KazukiNakamae/MaChIAto