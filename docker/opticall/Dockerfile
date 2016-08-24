FROM ubuntu:latest

MAINTAINER Eugene de Beste

RUN apt-get update -y && apt-get install wget build-essential -y

WORKDIR comp

RUN wget http://bitbucket.org/tss101/opticall/get/tip.tar.gz
RUN tar -zxf tip.tar.gz && \
    rm -rf tip.tar.gz && \
    mv * opticall && \
    cd opticall && cd opticall && \
    make && \
    mv opticall /usr/bin && \
    rm -rf /comp

WORKDIR /
