# use the ubuntu base image
FROM ubuntu:14.04

MAINTAINER Tobias Rausch rausch@embl.de

# install packages
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    g++ \
    gfortran \
    git \
    hdf5-tools \
    libboost-date-time-dev \
    libboost-program-options-dev \
    libboost-system-dev \
    libboost-filesystem-dev \
    libboost-iostreams-dev \
    libbz2-dev \
    libhdf5-dev \
    libncurses-dev \
    liblzma-dev \
    zlib1g-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# set environment
ENV BOOST_ROOT /usr

# install dicey
RUN cd /opt \
    && git clone --recursive https://github.com/gear-genomics/dicey.git \
    && cd /opt/dicey/ \
    && make STATIC=1 all \
    && make install

# Multi-stage build
FROM alpine:latest
RUN mkdir -p /opt/dicey/bin
WORKDIR /opt/dicey/bin
COPY --from=0 /opt/dicey/bin/dicey .

# Workdir
WORKDIR /root/

# Add Dicey to PATH
ENV PATH="/opt/dicey/bin:${PATH}"

# by default /bin/sh is executed
CMD ["/bin/sh"]
