############################################################
# Dockerfile to build MaChIAto Analyzer
############################################################

FROM mambaorg/micromamba:0.13.1

# File Author / Maintainer
MAINTAINER Kazuki Nakamae
RUN apt-cache search "^r-.*" | sort \
  && apt-get update \
  && apt-get install r-base r-base-dev r-cran-rprojroot xorg openbox cmake r-cran-seqinr r-cran-ape r-cran-tcltk2 r-bioc-biostrings r-cran-ggplot2 r-cran-reshape2 r-cran-maptools r-cran-scales r-cran-nloptr r-cran-rcppeigen -y --no-install-recommends \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/* \
  && rm -rf /usr/share/man/* \
  && rm -rf /usr/share/doc/* \
  && conda install -c anaconda python wget --yes \
  && conda install -c conda-forge libstdcxx-ng=7.3.0 --yes \
  && conda install -c bioconda emboss --yes \
  && conda clean --all --yes \
  && wget http://www.unafold.org/download/oligoarrayaux-3.8.tar.gz \
  && tar -zxvf oligoarrayaux-3.8.tar.gz \
  && cd oligoarrayaux-3.8 \
  && ./configure \
  && make \
  && make install \
  && cd ..

# install MaChIAto
COPY . .
WORKDIR .

RUN chmod +x trf409.linux64
RUN Rscript MaChIAtoAnalyzer.R

