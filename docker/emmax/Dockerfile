FROM ubuntu:latest

MAINTAINER Eugene de Beste

ENV pkg "emmax-beta-07Mar2010"
ENV ext "tar.gz"
ENV url "http://genetics.cs.ucla.edu/emmax/${pkg}.${ext}"

#RUN groupadd -r emmax && useradd -r -g emmax emmax

# Install the packages needed to download and extract plink
RUN apt-get update && apt-get install -y \
    wget

# Download and extract binary to /usr/bin

RUN wget $url && \
    tar zxf $pkg.$ext && \
    mv $pkg/* /usr/bin

RUN rm -rf $pkg && \
    rm -rf $pkg.$ext
