FROM ubuntu:latest

MAINTAINER Eugene de Beste
RUN groupadd -r plink && useradd -r -g plink plink

# Install the packages needed to download and extract plink
RUN apt-get update && apt-get install -y \
    unzip\
    wget

# Download and extract binary to /usr/bin

RUN wget https://www.cog-genomics.org/static/bin/plink160816/plink_linux_x86_64.zip && \
    unzip plink_linux_x86_64.zip -d /usr/bin/

RUN rm -rf plink_linux_x86_64.zip

#RUN mkdir /home/plink

#WORKDIR /home/plink

#USER plink
