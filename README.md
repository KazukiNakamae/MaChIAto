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

If you have a question and request, don't hesitate to contact me.

Kazuki Nakamae, Ph.D.
- kazukinakamae[at mark]gmail.com

## References

[1] Nakamae et. al, Detailed profiling with MaChIAto reveals various genomic and epigenomic features affecting the efficacy of knock-out, short homology-based knock-in and Prime Editing. bioRxiv. 2022 https://doi.org/10.1101/2022.06.27.496697

## Installation of MaChIAto

The MaChIAto consists of python and R with various packages other developers established. So the installation process is complex. However, you can quickly install MaChIAto using the docker container.

## Install MaChIAto using the docker (Strongly Recommended)

1. Download and install Docker Desktop: https://docs.docker.com/engine/install/#desktop

2. Enter Docker settings menu to adjust the memory allocation (~10G of memory is recommmended)

3. Type the following commands in the terminal.

```
# Download Docker image for MaChIAto Classifier
docker run --rm -v `pwd`:/DATA -w /DATA -i kazukinakamae/machiato:latest python ../MaChIAto.py -h;
# Download Docker image for MaChIAto Aligner
docker run --rm -v `pwd`:/DATA -w /DATA -i kazukinakamae/machiato_aligner:latest Rscript ../MaChIAtoAligner.R
# Download Docker image for collect_MaChIAto_data.py
docker run --rm -v `pwd`:/DATA -w /DATA -i kazukinakamae/collect_machiato_data:latest python ../collect_MaChIAto_data.py;
# Download Docker image for MaChIAto Analyzer
docker run --rm -v `pwd`:/DATA -w /DATA -i kazukinakamae/machiato_analyzer:latest Rscript ../MaChIAtoAnalyzer.R;
# Download Docker image for MaChIAto Reviewer
docker run --rm -v `pwd`:/DATA -w /DATA -i kazukinakamae/machiato_reviewer:latest Rscript ../MaChIAtoReviewer.R;
```

Then, if you wish, you can follow the manual below.

- Quick Start: https://machiatopage.github.io/2100/01/02/Quick-Start/
- Usage of MaChIAto Classifier: https://machiatopage.github.io/2022/06/22/Usage-of-MaChIAto-Classifier/
- Usage of MaChIAto Aligner: https://machiatopage.github.io/2022/06/22/Usage-of-MaChIAto-Aligner/
- Usage of collect_MaChIAto_data.py: https://machiatopage.github.io/2022/07/27/Usage-of-collect-MaChIAto-data-py/
- Usage of MaChIAto Analyzer: https://machiatopage.github.io/2022/06/22/Usage-of-MaChIAto-Analyzer/
- Usage of MaChIAto Reviewer: https://machiatopage.github.io/2022/06/22/Usage-of-MaChIAto-Reviewer/

## Time

The installation time of each MaChIAto function is approximately <10 minutes.

- MaChIAto Classifier
- MaChIAto Aligner
- collect_MaChIAto_data.py
- MaChIAto Reviewer

generally takes < 1 hour to complete the run.

- MaChIAto Analyzer

generally takes <1 day to complete run because the calculation uses numerous variants.

## Install MaChIAto within the conda environment (Deprecated)

The installation procedure is deprecated. We are **not supporting** the procedure because the environment is not stable and is difficult to reproduce.

### MacOSX

1. Click [HERE](https://github.com/KazukiNakamae/MaChIAto/archive/refs/heads/master.zip) to download MaChIAto.

Alternatively, enter the following command terminal.

```bash
git clone https://github.com/KazukiNakamae/MaChIAto.git;
```

2. Install miniconda3, R (>version R-4.0.1), Xcode, and XQuartz.

All program can run under the conda environment.The detailed setup procedure is [HERE](https://machiatopage.github.io/2100/01/01/Preparation/#more).

### Linux

1. Click [HERE](https://github.com/KazukiNakamae/MaChIAto/archive/refs/heads/master.zip) to download MaChIAto.

Alternatively, enter the following command terminal.

```bash
git clone https://github.com/KazukiNakamae/MaChIAto.git;
```

2. Install micromamba (version 0.13.1). All packages can be installed using mamba and apt-get in the same way as generating the docker images.
Please refer to the dockerfile of each directory (MaChIAto, MaChIAto_Aligner, collect_MaChIAto_data, MaChIAto_Analyzer, MaChIAto_Reviewer).


# LICENCE

MaChIAto Classifier is released under the MIT Licence. You can use our software for any purpose including commercial use.
MaChIAto Aligner, MaChIAto Analyzer, and MaChIAto Reviewer are freely available to use for academic and commercial use.
Please enjoy MaChIAto!!!
