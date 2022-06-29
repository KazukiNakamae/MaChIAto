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

We discribe the template command the the example in the [Online manual](https://machiatopage.github.io/).

**Let's Check it!!!**
https://machiatopage.github.io/


If you have a question and request, don't hesitate to contact me [Kazuki Nakamae, Ph.D.](kazukinakamae@gmail.com).

# Download MaChIAto

Click [HERE](https://github.com/KazukiNakamae/MaChIAto/archive/refs/heads/master.zip).

Alternatively, enter the following command terminal.

```bash
git clone https://github.com/KazukiNakamae/MaChIAto.git;
```

You can check the [Quick Start](https://machiatopage.github.io/2100/01/02/Quick-Start/#more) to start using MaChIAto.

# Environment

MaChIAto needs miniconda3, R (>version R-4.0.1), Xcode, and XQuartz.
All program can run under the conda environment. The detailed setup procedure is [HERE](https://machiatopage.github.io/2100/01/01/Preparation/#more).

**We are preparing the Docker container for MaChIAto. Please look forward to it**😁

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



# MaChIAto Analyzer

You have to execute the collect_MaChIAto_data.py before using MaChIAto Analyzer.

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



# MaChIAto Reviewer

### Parameter list

Summary directory:
The summary directory generated with “collect_MaChIAto_data.py”.

Output prefix:
The directory into which the output directory is saved.



# LICENCE

MaChIAto is released under the MIT Licence. You can use our software for any purpose (including commercial use). Please enjoy MaChIAto!!!
