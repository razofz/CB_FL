FROM continuumio/miniconda3:4.10.3p1

RUN apt update -y
RUN apt install -y build-essential

COPY *.yaml /home/

RUN conda install mamba -c conda-forge -y
RUN mamba env create -f /home/DESeq2.yaml
RUN mamba env create -f /home/notebooks.yaml
RUN mamba env create -f /home/seurat.yaml

RUN apt install -y direnv

EXPOSE 8888 8888
LABEL author=https://github.com/razofz
LABEL maintainer=razofz
LABEL preprint=https://doi.org/10.1101/2021.12.29.474425

RUN mkdir /fl
ENV PROJECT_PATH=/fl
WORKDIR /fl
