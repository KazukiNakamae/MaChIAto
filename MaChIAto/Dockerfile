############################################################
# Dockerfile to build MaChIAto
############################################################

FROM mambaorg/micromamba:0.13.1

# File Author / Maintainer
MAINTAINER Kazuki Nakamae
RUN conda install -c anaconda python=3.8 --yes && pip install --upgrade pip \
  && pip install regex tqdm argparse biopython numpy matplotlib GPy gpyopt datetime pandas \
  && conda clean --all --yes

# install MaChIAto
COPY . .
WORKDIR .

RUN python MaChIAto.py -h

