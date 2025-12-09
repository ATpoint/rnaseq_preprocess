FROM condaforge/mambaforge:latest

COPY ["environment.yml", "./"]

RUN apt update && \
    DEBIAN_FRONTEND=noninteractive apt install -y --no-install-recommends tzdata procps nano && apt autoclean

RUN mamba env update --name base --file environment.yml
