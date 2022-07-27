############################################################
# Dockerfile to build MaChIAto Reviewer
############################################################

FROM kazukinakamae/machiato_aligner:lastest

# File Author / Maintainer
MAINTAINER Kazuki Nakamae

# install MaChIAto
COPY . .
WORKDIR .

RUN Rscript MaChIAtoReviewer.R

