From debian:stretch
MAINTAINER HSPH QBRC <ddeconti@hsph.harvard.edu>

RUN apt-get update && \
  apt-get install -y build-essential \
    wget \
    zlib1g-dev \
    libcurl4-openssl-dev \
    libncurses-dev \
    libbz2-dev \
    liblzma-dev \
    libxml2-dev \
    zip \
    python3-dev \
    python3-pip \
    openjdk-8-jdk \
    pandoc \
    git

# Install all software under /opt/software:
RUN mkdir -p /opt/software 

# Install fastQC:
RUN cd /opt/software && \
  wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip && \
  unzip fastqc_v0.11.8.zip && \
  chmod +x /opt/software/FastQC/fastqc
ENV PATH="/opt/software/FastQC:${PATH}"

# Install libStatGen and fastQValidator:
RUN cd /opt/software && \
  git clone https://github.com/statgen/libStatGen.git && \
  cd libStatGen && \
  make && \
  cd .. && \
  git clone https://github.com/statgen/fastQValidator.git && \
  cd fastQValidator && \
  make && \
make install

# Install samtools:
RUN cd /opt/software && \
  wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
  tar -xjf samtools-1.9.tar.bz2 && \
  cd samtools-1.9 && \
  ./configure --prefix=/opt/software/samtools && \
  make && \
  make install && \
  rm /opt/software/samtools-1.9.tar.bz2
ENV PATH="/opt/software/samtools/bin:${PATH}"

# Get the Picard JAR:
RUN mkdir -p /opt/software/picard && \
  cd /opt/software/picard && \
  wget https://github.com/broadinstitute/picard/releases/download/2.18.26/picard.jar
ENV PICARD_JAR=/opt/software/picard/picard.jar

# Install BWA
RUN mkdir -p /opt/software/ && \
  cd /opt/software/ && \
  wget https://github.com/lh3/bwa/archive/v0.7.17.tar.gz && \
  tar xzf v0.7.17.tar.gz && \
  cd bwa-0.7.17 && \
  make
ENV PATH="/opt/software/bwa-0.7.17:${PATH}"

# Get the GATK JAR:
RUN mkdir -p /opt/software/gatk && \
  cd /opt/software/gatk && \
  wget -O gatk-4.1.1.0.zip  https://github.com/broadinstitute/gatk/releases/download/4.1.1.0/gatk-4.1.1.0.zip && \
  unzip gatk-4.1.1.0.zip
ENV GATK_JAR=/opt/software/gatk/gatk-4.1.1.0/gatk-package-4.1.1.0-local.jar

# Install some Python3 libraries:
ADD requirements.txt /opt/software/
RUN pip3 install -r /opt/software/requirements.txt

# Install multiqc:
RUN pip3 install multiqc

# Install Jinja2 for templating.  Likely installed with other python packages, but just to be sure:
RUN pip3 install Jinja2

# for multiqc to read commandline:
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

# The scripts for performing pre-workflow checks:
RUN mkdir /opt/software/precheck
ADD check_fastq.py /opt/software/precheck/

# The script for generating the markdown report:
ADD generate_report.py /usr/local/bin/
RUN chmod +x /usr/local/bin/generate_report.py

# Add the report template and styling:
RUN mkdir /opt/report
ADD report.md /opt/report/
ADD report.css /opt/report/

ENTRYPOINT ["/bin/bash"]