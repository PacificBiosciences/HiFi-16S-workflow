FROM continuumio/miniconda3

WORKDIR /app

COPY pb-16s-vis-conda.yml .
RUN rm -rf /var/lib/apt/lists/*
RUN apt-get update
RUN apt-get install -y procps
RUN conda install -n base -c conda-forge mamba
RUN mamba env create -n pb-16S-vis -f pb-16s-vis-conda.yml
RUN echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc
RUN echo "conda activate pb-16S-vis" >> ~/.bashrc
ENV PATH /opt/conda/envs/pb-16S-vis/bin:$PATH
