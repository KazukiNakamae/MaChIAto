############################################################
# Dockerfile to build MaChIAto Aligner
############################################################

FROM mambaorg/micromamba:0.13.1

# File Author / Maintainer
MAINTAINER Kazuki Nakamae
RUN conda install -c anaconda xz=5.2.4 --yes \
  && conda install -c conda-forge ncurses=6.1 perl bzip2 curl=7.64.1 libdeflate=1.3 --yes \
  && conda install -c bioconda bwa=0.7.17 samtools=1.9 --yes \
  && conda install -c anaconda python --yes \
  && conda clean --all --yes \
  && apt-cache search "^r-.*" | sort \
  && apt-get update \
  && apt-get install r-base r-base-dev xorg openbox -y --no-install-recommends \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/* \
  && rm -rf /usr/share/man/* \
  && rm -rf /usr/share/doc/*

# install MaChIAto
COPY . .
WORKDIR .

RUN Rscript MaChIAtoAligner.R

