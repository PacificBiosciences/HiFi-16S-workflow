FROM continuumio/miniconda3

WORKDIR /app

COPY qiime2-2022.2-py38-linux-conda.yml .
RUN rm -rf /var/lib/apt/lists/*
RUN apt-get update
RUN apt-get install -y procps bc
RUN conda install -n base -c conda-forge mamba
RUN mamba env create -n qime2-2022.2 -f qiime2-2022.2-py38-linux-conda.yml
RUN echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc
RUN echo "conda activate qime2-2022.2" >> ~/.bashrc
SHELL ["conda", "run", "-n", "qime2-2022.2", "/bin/bash", "-c"]
RUN pip install git+https://github.com/kaanb93/q2-krona.git
ENV PATH /opt/conda/envs/qime2-2022.2/bin:$PATH
