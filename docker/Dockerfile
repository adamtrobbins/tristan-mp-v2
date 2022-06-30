FROM ubuntu:20.04
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    git zsh neovim wget curl \
    python3 python3-pip build-essential libopenmpi-dev \
    libhdf5-openmpi-dev hdf5-tools

RUN ln -s $(which python3) $(dirname $(which python3))/python
RUN ln -s $(which neovim) $(dirname $(which neovim))/vim
RUN ln -s $(which neovim) $(dirname $(which neovim))/vi

RUN chsh -s /bin/zsh

RUN apt-get clean && \
    apt-get autoclean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/cache/* && \
    rm -rf /var/lib/log/*

ENTRYPOINT ["tail", "-f", "/dev/null"]
