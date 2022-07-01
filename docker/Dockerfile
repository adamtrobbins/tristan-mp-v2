FROM ubuntu:20.04
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    git zsh vim wget curl \
    python3 python3-pip build-essential libopenmpi-dev \
    libhdf5-openmpi-dev hdf5-tools

RUN ln -s /bin/python3 /bin/python

RUN chsh -s /bin/zsh

RUN git clone https://github.com/tarjoilija/zgen.git "${HOME}/.zgen"
RUN echo 'source "${HOME}/.zgen/zgen.zsh"\n\
zgen load denysdovhan/spaceship-zsh-theme spaceship\n\
alias ls="ls --color=auto"' >> /root/.zshrc

RUN pip install autopep8 fprettify json5 packaging pycodestyle pyparsing ConfigArgParse

RUN apt-get clean && \
    apt-get autoclean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/cache/* && \
    rm -rf /var/lib/log/*

ENTRYPOINT ["tail", "-f", "/dev/null"]
