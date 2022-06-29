# MaChIAto
##### LATEST VERSION (2022/06/28)

<img src="https://github.com/KazukiNakamae/temp/blob/temp-images/logo.png" alt="MaChIAto_logo" title="MaChIAto_logo" width="200" height="200">

MaChIAto (**M**icrohomology-**a**ssociated **Ch**romosomal **I**ntegration/editing **A**nalysis **to**ols); a comprehensive analysis software that can precisely classify, deeply analyze, correctly align, and thoroughly review the targeted amplicon sequencing analysis data obtained by various CRISPR experiments, including template-free gene knock-out, short homology-based gene knock-in, and even a new-class CRISPR methodology, Prime Editing.

MaChIAto would be helpful for people who want to
- get a more accurate classification of knock-out, homology-based knock-in, and Prime Editing.
- quantify the rate of imprecise editing reads, including the reads with donor-reversed integration.
- obtain the local alignment of the target site and knock-in junction.
- analyze the sequence involved with InDels and imprecise editing.
- find the useful features for the efficiency prediction.

Online manual:
https://machiatopage.github.io/

If you have a question and request, don't hesitate to contact me [Kazuki Nakamae, Ph.D.](kazukinakamae@gmail.com).

# Quick Start

We demonstrate a simple example.
Before that, you have to build the environment according to the section of Preparation.

#### 1. Download MaChIAto

Click "Download ZIP".

![download_zip.png](https://github.com/KazukiNakamae/temp/blob/temp-images/download_zip.png)

Alternatively, enter the following command terminal.

```bash
git clone https://github.com/KazukiNakamae/MaChIAto.git;
```

#### 2. Open a terminal and go the directory.

```bash
cd (...)/MaChIAto
```

#### 3. Re-classify the allele frequency table derived from CRISPResso2.

```bash
source activate MaChIAto_env; # *You do not need to enter it again once you did it.
```

The re-classification command of MaChIAto (Classifier) 
```bash
python MaChIAto/MaChIAto.py \
-m CRISPResso2 \
-ccf (Alleles_frequency_table.zip of the CRISPResso2 output) \
-o (output directory) \
-a (wt amplicon sequence) \
-g (protospaser sequence) \
(other parameters);
```

```bash
conda deactivate;
```

The other parameters varies per experiment type {knock-out, homology-based knokc-in, Prime Editing}.

The example is for the knock-out analysis.

##### Case of knock-out
```bash
python MaChIAto/MaChIAto.py \
-ccf (Alleles_frequency_table.zip of the CRISPResso2 output) \
-o (output directory) \
-a (wt amplicon sequence) \
-g (protospaser sequence);
```

*The case of the knock-out and knock-in is shown in the section of MaChIAto (MaChIAto Classifier).

#### 4. Check the classification result.

You can get the classification result of the editing. The result is visualized as the pie chart based on the ALL_dataframe.txt.

We provide the detailed description how to read them in the section of the MaChIAtoClassifier output.

# Quick Start (Optional)

MaChIAto can profile and visualize the MaChIAto output using MaChIAto Aligner, MaChIAto Analyzer, and MaChIAto Reviewer. 

#### 5. Get local alignment and the mutation profiling using MaChIAto Aligner.

```bash
source activate MaChIAto_Aligner_env; # *You do not need to enter it again once you did it.
```

The alignment command of MaChIAto Aligner
```bash
Rscript MaChIAto_Aligner/MaChIAtoAligner.R \
(directory of MaChIAto Classifier) \
(output directory);
```

```bash
conda deactivate;
```

If your donor insert has the extra sequence outside of the homology arm, you should enter it.
Please refer to the section of MaChIAtoAligner.

#### 6. Check the alignment result.

You can get the local alignment of the editing. The result is visualized as the map.

We provide a detailed description of how to read them in the section of the MaChIAtoAligner output.

#### 7. Aggregate the multiple data.

If you have multiple data (n>3) to profile the characteristics, you should aggregate it using the following command.

```bash
source activate MaChIAto_Analyzer_env; # *You do not need to enter it again once you did it.
python MaChIAto_Analyzer/collect_MaChIAto_data.py \
-i (the prefix of MaChIAto Classifier output) \
-o (output directory) \
-ol (knock-out label) \
-t (calculation target) \
```

```bash
conda deactivate;
```

It is used for the knock-out analysis. The calculation target for MaChIAto Analyzer/Reviewer is protospacer.
You can target the homology arm and RT template using another command. Please see the section of MaChIAtoAnalyzer.

#### 8. Investigate the relationship with (epi-)genomic context

You can see the correlation with the >70 (epi-)genomic context using MaChIAto Analyzer.

```bash
Rscript MaChIAto_Analyzer/MaChIAtoAnalyzer.R \
(directory of collect_MaChIAto_data.py) \
(output directory);
```

```bash
conda deactivate;
```

MaChIAtoAnalyzer basically use > 70 genomic context.
If you want to use epigenomic context or other context, please see the section of MaChIAtoAnalyzer.

#### 9. Check the correlation.

You can get the correlation between the editing efficacy and the context. The result is visualized as the scatter plot.

We provide a detailed description of how to read them in the section of the MaChIAtoAnalyzer output.

#### 10. Profile the mutation tendency.


```bash
source activate MaChIAto_Reviewer_env; # *You do not need to enter it again once you did it.
Rscript MaChIAto_Reviewer/MaChIAtoReviewer.R \
(the prefix of MaChIAto Classifier output) \
(the prefix of MaChIAto Aligner output) \
(directory of collect_MaChIAto_data.py) \
(output directory);
```

```bash
conda deactivate;
```

#### 11. Check the mutation profile.

You can get the mutation profile of the editing. The result is visualized as the bar plot and pie chart.

We provide a detailed description of how to read them in the section of the MaChIAtoReviewer output.


# Preparation

### Change default shell

#### 1. Set bash as default shell
```bash
chsh -s /bin/bash;
```

### Install software

#### 1. install miniconda3 from Conda
(https://docs.conda.io/en/latest/miniconda.html)
#### 2. install R (>version R-4.0.1) from the CRAN
(https://cran.ism.ac.jp)
#### 3. install Xcode from the App Store
(https://apps.apple.com/jp/app/xcode/id497799835?mt=12)
```bash
sudo xcodebuild -license;
```
#### 4. install XQuartz from the XQuartz project
(xquartz.macosforge.org)

### Build the environment using conda

#### 1. Environment for MaChIAto_(MaChIAto Classifier)
```bash
conda create --name MaChIAto_env;
source activate MaChIAto_env;
conda install -c anaconda python=3.8;
pip install --upgrade pip;
pip install regex tqdm argparse biopython numpy matplotlib GPy gpyopt datetime pandas;
conda deactivate;
```

#### 2. Environment for MaChIAtoAligner
```bash
conda create --name MaChIAto_Aligner_env;
source activate MaChIAto_Aligner_env;
conda install -c bioconda python=3.8 bwa=0.7.17 samtools=1.9;
conda deactivate;
```

#### 3. Environment for MaChIAtoAnalyzer
```bash
conda create --name MaChIAto_Analyzer_env;
source activate MaChIAto_Analyzer_env;
conda install -c anaconda python=3.8 wget;
conda install -c bioconda emboss;
conda install -c bioconda oligoarrayaux;
pip install --upgrade pip;
pip install numpy regex tqdm pandas;
conda deactivate;
```

#### 4. Environment for MaChIAtoReviewer
```bash
conda create -n MaChIAto_Reviewer_env --clone MaChIAto_Aligner_env;
source activate MaChIAto_Reviewer_env;
conda deactivate;
```

# Usage

## Usage of MaChIAto Classifier

### Parameter list

**-m** or **--mode** {str}: This parameter allows for the specification of the type of analysis: “CRISPResso” and “CRISPResso2” is allowed in the latest version

**-a** or **--amplicon_seq** {str}: This parameter allows the user to enter the amplicon sequence used for the CRISPResso. The length should be >105bp due to the setting for the other parameter.

**-g** or **--guide_seq** {str}: This parameter allows for the specification of the sgRNA sequence used for the CRISPResso. The length of the sequence should be 20nt without PAM. The MaChIAto convention is to depict the expected cleavage position using the value of the parameter three nt 3' from the end of the guide.", required=True)

**-cf** or **--crispreeso_file** {str}: This parameter allows for the specification of the “Alleles_frequency_table.txt” from CRISPResso. When this parameter is used, “CRISPResso” should be entered as -m parameter.', (default: "./Alleles_frequency_table.txt") (optional)

**-ccf** or **--crispreeso2_file** {str}: This parameter allows for the specification of the “Alleles_frequency_table.zip” from CRISPResso2. When this parameter is used, “CRISPResso2” should be entered as -m parameter.', (default: "./Alleles_frequency_table.zip") (optional)

**-d** or **--donor_seq** {str}:This parameter allows for the specification of the expected HDR amplicon used for the CRISPResso. The length of sequence should be >12bp in knock-out/knock-in analysis. In knock-out analysis, this parameter should not be entered, and then the value will be given fake parameter (“TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT”), and some internal settings are changed for knock-out analysis. However, the fake parameter will be poly-C|G|A if amplicon sequence contains poly-T sequence.', (default: "") (optional)

**-e** or **--expected_ki_amplicon_seq** {str}: This parameter allows for the specification of the expected knock-in amplicon sequence used for the CRISPResso after HDR. The length of the sequence should be >12bp. In knock-out analysis, this parameter should not be entered, and then the value will be given fake parameter including fake donor sequence, and some internal settings are changed for knock-out analysis.', (default: "") (optional)

**-o** or **--output_folder** {str}: This parameter allows for the specification of the output directory to use for the analysis (default: current directory) (default: "./") (optional)

**-lh** or **--length_left_homologyarm** {int}: This parameter allows for the specification of the length of 5’ homology arm (default: 20). The length of the sequence should be >17bp, and the flanking sequence with the homology arm needs the length of >24bp in the expected amplicon sequence. In the prime editing analysis, the 5’ homology arm is considered as the RT template (default: 20) (optional)

**-rh** or **--length_right_homologyarm** {int}: This parameter allows for the specification of the length of 3’ homology arm (default: 20). The length of the sequence should be >17bp, and the flanking sequence with the homology arm needs the length of >24bp in the expected amplicon sequence. In the prime editing analysis, the 3’ homology arm is considered as the prime binding site (default: 20) (optional)

**-cn** or **--location_comp_nick** {int}: This parameter allows for the specification of the complementary strand nick location [3prime direction is +] (default: 90). This parameter is used in the prime editing and should be over the homology arm's length to which the nickase is adjacent' (default: 90) (optional)

**-n** or **--name**  : This parameter allows for the specification of the name, which will be included output directory (default: “untitled-X”). If MaChIAto Analyzer and MaChIAto Reviewer will be used in the following analysis, the value should be “{target_name}-{sample label }” (e.g., DBF4B-C) and underbar "_" should not be used', (default: "untitled-X") (optional)

**--primeediting_analysis** : Re-classify the data as prime editing analysis. This option forces the setting to change for prime editing analysis.' (optional)

#### Advanced parameter (for developers)

**--force_knockout_analysis** : Usually, MaChIAto re-classify the data as knock-in analysis. This option forces the setting to change for knock-out analysis. Under this mode, the length of indicator on the knock-in donor is the maximum value, and the threshold value for alignment of the knock-in sequence is 1.0. If this parameter is not entered, MaChIAto can automatically set up this mode by finding some characteristics of the knock-out sample. For example, MaChIAto checks that there is no donor sequence or expected knock-in sequence as input, and there are less three kinds HDR variants among input data. (optional)

**--skip_optimization** : Usually, MaChIAto runs Bayesian optimization for finding the optimized setting. This option allows MaChIAto to skip the process of optimization. The option is made for debugging. So, the option should not be used in the usual analysis. However, if the optimization process disturbs an accurate analysis, this option might be useful. (optional)

**--copy_optimization** : Usually, MaChIAto runs Bayesian optimization for finding the optimized setting. This option allows MaChIAto to use the provided optimization data instead of the optimization. If you analyze NGS data with the setting of previously analyzed data, this option is useful. Especially in substitution editing, we can recommend applying the option to the negative/positive control sample. However, you should understand that a classification error might frequently occur. (optional)

**--provided_optimization_file** : This parameter allows for the specification of the “MaChIAto_optimized_param.csv” from the analyzed MaChIAto folder. When this parameter is used, --copy_optimization parameter is required.', (default: "./MaChIAto_optimized_param.csv") (optional)

### Template command

##### Case of knock-out
```bash
python MaChIAto/MaChIAto.py \
-m CRISPResso2 \
-ccf (Alleles_frequency_table.zip of the CRISPResso2 output) \
-o (output directory) \
-a (wt amplicon sequence) \
-g (protospaser sequence) \
-n (sample name)-(label name);
```

##### Case of homology-based knock-in
```bash
python MaChIAto/MaChIAto.py \
-m CRISPResso2 \
-ccf (Alleles_frequency_table.zip of the CRISPResso2 output) \
-o (output directory) \
-a (wt amplicon sequence) \
-g (protospaser sequence) \
-e (expected editing amplicon seqeunce) \
-d (donor insert sequence) \
-lh (length of 5\' homology arm) \
-rh (length of 3\' homology arm) \
-n (sample name)-(label name);
```

##### Case of Prime Editing (Substitution/Deletion editing)
```bash
python MaChIAto/MaChIAto.py \
-m CRISPResso2 \
-ccf (Alleles_frequency_table.zip of the CRISPResso2 output) \
-o (output directory) \
-a (wt amplicon sequence) \
-g (protospaser sequence) \
-e (expected editing amplicon seqeunce) \
-lh (length of prime binding site) \
-rh (length of RT template) \
-cn (distance between two nick sites) \
--primeediting_analysis \
-n (sample name)-(label name);
```

##### Case of Prime Editing (Insertion editing)
```bash
python MaChIAto/MaChIAto.py \
-m CRISPResso2 \
-ccf (Alleles_frequency_table.zip of the CRISPResso2 output) \
-o (output directory) \
-a (wt amplicon sequence) \
-g (protospaser sequence) \
-e (expected editing amplicon seqeunce) \
-d (donor insert sequence) \
-lh (length of prime binding site) \
-rh (length of RT template) \
-cn (distance between two nick sites) \
--primeediting_analysis \
-n (sample name)-(label name);
```

If you want to use the result of CRISPResso version1, the -ccf should be replased into -cf.
Example for CRISPResso version1
```bash
python MaChIAto/MaChIAto.py \
-m CRISPResso \
-cf (Alleles_frequency_table.txt of the CRISPResso output) \
-o (output directory) \
-a (wt amplicon sequence) \
-g (protospaser sequence) \
-n (sample name)-(label name);
```

*The **sample name** and **label name** can be **arbitrary**. If you analyze the multiple experiment (e.g. knock-out and knock-in), the name should be different from others.


# MaChIAtoClassifier output

The classification result is visualized in the "MaChIAto_alignment_pie_chart.png"
"MaChIAto_alignment_pie_chart.png" has the following labels:
- Unmodified: **The unedited amplicon sequence**
- Perfect Editing: **The edited sequence** that a user intended to generate (e.g., precise knock-in sequence, expected substitution/deletion/insertion...)
- Imperfect Editing: The sequence that is simular to "Perfect Editing" sequence. However, **the junction sequence is imperfect**. 
*The class is specilized for the homology-based knock-in. In knock-out and Prime Editing, the "Imperfect Editing" reads can be condidered as just InDels.
InDels: The sequence that has a trace of editing (e.g., deletion/insertion/frequent substitutions) around the target site.
- Deletion: The InDel sequence whose length is **less** than that of the unedited amplicon sequence
- Insertion: The InDel sequence whose length is **more** than that of the unedited amplicon sequence
- Frequent substitutions: The InDel sequence whose length is **equal** to that of the unedited amplicon sequence
UC: Unclear sequence that has the insufficient information to be recognized as the focusing sequence.
- UC (Ambiguous Editing): There is the expected sequence (e.g., knock-in cassette and Reverse transcript). However, the information of the end of amplicon is unclear. For example, The edited reads, which do not have the site of complementary strand nicking, are classified into this class.
- UC (Others): There is not enough sequence informaton for the classification. The unedited amplicon read with sequence errors and read amplified from the other locus can be classified into the this class.



# MaChIAto Aligner

### Parameter list

Input directory:
The directory that MaChIAto Classifier generated.

Output prefix:
The directory into which the output directory is saved.

Left extra sequence (optional):
The outside sequence flanking 5'-homology arm of the donor. 

Right extra sequence (optional):
The outside sequence flanking 3'-homology arm of the donor. 

### Template command

```bash
Rscript MaChIAto_Aligner/MaChIAtoAligner.R \
(Input directory) \
(Output prefix) \
(Left extra sequence) \
(Right extra sequence);
```

*when you use external storage for save output, the process can be too slow. We recommend that you use the tool in the local storage.


# MaChIAto Analyzer

## collect_MaChIAto_data.py

### Parameter list

**-i** or **--indir** {str}: input directory which includes the outputs of MaChIAto Classifier

**-o** or **--outdir** {str}: output directory

**-ul** or **--untreated_label** {str}: untreated label (This must be one.)' (default: "machiato_dummy_sample") (optional)

**-ol** or **--knock_out_label** {str}: knock-out label' (default: []) (optional)

**-il** or **--knock_in_label** {str}: knock-in label' (default: []) (optional)

**-sc** or **--scaffold_seq** {str}: scaffold sequence of sgRNA (default: "gttttagagctagaaatagcaagttaaaataaggctagtccgttatcaacttgaaaaagtggcaccgagtcggtgc") (optional)

**-t** or **--target_type** {str}: target sequence type {"lmh" | "rmh" | "bmh" | "elmh" | "ermh" | "ebmh" | "protospacer"} (default: "bmh") (optional)

**--ignore_list** {str}: The list of ignore target set contains target names which are not desired to analyze for some reasons. The data (e.g. DBF4B-A, DBF4B-B, DBF4B-C, DBF4B-D) including target name (e.g. DBF4B) shown in the list is skipped through the process of MaChIAto Analyzer. The format should be comma-separated like “TargetA, TargetB, …” “example_data” directory has “ignore_list.csv” as example. (default: "") (optional)

### Template command

##### Case of Double knock-in analysis (*ADVANCE: The analysis includes the comparison between two knock-in methods.)

```bash
python MaChIAto_Analyzer/collect_MaChIAto_data.py \
-i (the prefix of MaChIAto Classifier output) \
-o (output directory) \
-sc (scaffold sequence of sgRNA) \
-ul (untreated label) \
-ol (knock-out label) \
-il (knock-in label 1) (knock-in label 2) \
-t (calculation target) \
--ignore_list (list of samples ignored);
```

##### Case of Single knock-in analysis (*STANDARD)

```bash
python MaChIAto_Analyzer/collect_MaChIAto_data.py \
-i (the prefix of MaChIAto Classifier output) \
-o (output directory) \
-sc (scaffold sequence of sgRNA) \
-ul (untreated label) \
-ol (knock-out label) \
-il (knock-in label) \
-t (calculation target) \
--ignore_list (list of samples ignored);
```

##### Case of Simple knock-in analysis (*SIMPLE: The analysis can be applied when there is no knock-out sample used as control.)

```bash
python MaChIAto_Analyzer/collect_MaChIAto_data.py \
-i (the prefix of MaChIAto Classifier output) \
-o (output directory) \
-sc (scaffold sequence of sgRNA) \
-il (knock-in label) \
-t (calculation target) \
--ignore_list (list of samples ignored);
```

##### Case of Double knock-out analysis (*ADVANCE: The analysis includes the comparison between two knock-out methods.)

```bash
python MaChIAto_Analyzer/collect_MaChIAto_data.py \
-i (the prefix of MaChIAto Classifier output) \
-o (output directory) \
-sc (scaffold sequence of sgRNA) \
-ul (untreated label) \
-ol (knock-out label 1) (knock-out label 2) \
-t (calculation target) \
--ignore_list (list of samples ignored);
```

##### Case of Single knock-out analysis (*STANDARD)

```bash
python MaChIAto_Analyzer/collect_MaChIAto_data.py \
-i (the prefix of MaChIAto Classifier output) \
-o (output directory) \
-sc (scaffold sequence of sgRNA) \
-ul (untreated label) \
-ol (knock-out label) \
-t (calculation target) \
--ignore_list (list of samples ignored);
```

*If you do not enter "-ul (untreated label)", the process can work. However, some filtering process will be skipped.


## MaChIAto Analyzer

### Parameter list

Summary directory:
The summary directory generated with “collect_MaChIAto_data.py”.

Output prefix:
The directory into which the output directory is saved.

Name of extra data (optional):
The name of the feature group that the next argument includes.

Table of extra data (optional):
The pathname of extra data added by the user. The data should be a .csv file.

### Template command

```bash
Rscript MaChIAto_Analyzer/MaChIAtoAnalyzer.R \
(Summary directory) \
(Output prefix) \
(Name of extra data) \
(Table of extra data);
```


# MaChIAto Reviewer

### Parameter list

Summary directory:
The summary directory generated with “collect_MaChIAto_data.py”.

Output prefix:
The directory into which the output directory is saved.

### Template command

Rscript MaChIAto_Reviewer/MaChIAtoReviewer.R \
(Summary directory) \
(Output prefix);

