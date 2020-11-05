FROM continuumio/miniconda3:4.8.2

LABEL authors="Gaofei Zhao"

# Updates and install necessary packages
RUN apt-get update
RUN apt-get install -y libgsl-dev unzip

# Install R packages
RUN conda install -y -c r r-base r-essentials r-devtools

# Install project packages
RUN ln -s /bin/tar /bin/gtar
RUN R -e 'options(unzip = "internal")'
RUN R -e 'devtools::install_github("mskcc/tempoSig")'

# Copy and set work directory
COPY . ./tempoSig
WORKDIR ./tempoSig

RUN mkdir input
RUN mkdir output
