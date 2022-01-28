FROM condaforge/mambaforge:4.11.0-0

COPY ["environment.yml", "./"]

RUN apt update && \
    DEBIAN_FRONTEND=noninteractive apt install -y --no-install-recommends tzdata procps nano && \
    apt-get clean && \
    wget -q https://github.com/COMBINE-lab/salmon/releases/download/v1.6.0/salmon-1.6.0_linux_x86_64.tar.gz && \
    tar zxf salmon-1.6.0_linux_x86_64.tar.gz 

ENV PATH=/salmon-1.6.0_linux_x86_64/bin:$PATH

RUN mamba env update --name base --file environment.yml
            
