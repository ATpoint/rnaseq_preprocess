FROM r-base:4.1.2

RUN apt-get update && \
    apt-get install -y ca-certificates wget nano procps && \
    wget -q https://github.com/COMBINE-lab/salmon/releases/download/v1.6.0/salmon-1.6.0_linux_x86_64.tar.gz && \
    tar zxf salmon-1.6.0_linux_x86_64.tar.gz

ENV PATH=/salmon-1.6.0_linux_x86_64/bin:$PATH
    