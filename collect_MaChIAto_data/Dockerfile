############################################################
# Dockerfile to build collect_MaChIAto_data.py
############################################################

FROM mambaorg/micromamba:0.13.1

# File Author / Maintainer
MAINTAINER Kazuki Nakamae
RUN conda install -c anaconda python=3.8 --yes && pip install --upgrade pip \
  && pip install regex tqdm argparse pandas \
  && conda clean --all --yes

# install MaChIAto
COPY . .
WORKDIR .

RUN python collect_MaChIAto_data.py -h

