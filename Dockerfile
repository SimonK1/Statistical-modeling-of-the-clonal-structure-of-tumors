FROM ubuntu:23.04

COPY app /app

RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    r-base \
    r-cran-rcurl \
    r-base-dev \
    git \
    wget \
    libcurl4-openssl-dev \
    libxml2-dev \
    libssl-dev 

RUN pip3 install -r reqs.txt --break-system-packages

RUN pip3 install git+https://github.com/Roth-Lab/pyclone-vi.git --break-system-packages

RUN wget https://cran.r-project.org/src/contrib/Archive/NORMT3/NORMT3_1.0.4.tar.gz

RUN R -e "install.packages('NORMT3_1.0.4.tar.gz', repos = NULL, type = 'source')"

RUN Rscript installation.R

WORKDIR /app

EXPOSE 8888

CMD ["/bin/bash", "-c", "jupyter lab --ip='0.0.0.0' --port=8888 --no-browser --allow-root"]