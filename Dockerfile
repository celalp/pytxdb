FROM ubuntu:20.04


RUN apt-get update -y
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y tzdata
RUN apt-get install -y software-properties-common curl libpq-dev \
    python3.11-dev python3.10-dev python3.9-dev python3.8-dev git

RUN curl -sS https://downloads.mariadb.com/MariaDB/mariadb_repo_setup | bash \
    apt-get update -y

RUN add-apt-repository ppa:deadsnakes/ppa
RUN apt install -y wget python3.8 python3.9 python3.10 python3.11 \
    samtools gcc zlib1g-dev liblzma-dev libcurl4-openssl-dev \
    libmariadb3 libmariadb-dev \
    wget https://bootstrap.pypa.io/get-pip.py \
    python3.11 get-pip.py \
    python3.10 get-pip.py \
    python3.9 get-pip.py \
    python3.8 get-pip.py

RUN git clone https://github.com/celalp/pytxdb.git \
    cd pytxdb \
    python3.11 -m pip install . \
    python3.10 -m pip install . \
    python3.9 -m pip install . \
    python3.8 -m pip install .





