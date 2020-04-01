
FROM continuumio/miniconda3:4.8.2

LABEL org.label-schema.license="GPL-2.0" \
      org.label-schema.vcs-url="https://github.com/rsettlag" \
      maintainer="Robert Settlage <rsettlag@vt.edu>"

RUN apt-get update -y \
  && apt install -y git wget
RUN conda create -y -n viralrecall python=3.7.4
RUN conda install -y -n viralrecall biopython \
  && conda install -y -n viralrecall hmmer prodigal -c bioconda \
  && conda install -y -n viralrecall pandas matplotlib
RUN git clone https://github.com/faylward/viralrecall

