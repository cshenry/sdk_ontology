FROM kbase/kbase:sdkbase.latest
MAINTAINER KBase Developer
# -----------------------------------------

# Insert apt-get instructions here to install
# any required dependencies for your module.

# RUN apt-get update

# -----------------------------------------

RUN cpanm -i Config::IniFiles
RUN apt-get install libffi-dev libssl-dev
RUN pip install --upgrade requests[security]
RUN apt-add-repository ppa:webupd8team/java
RUN apt-get update
RUN apt-get -q install -y oracle-java8-installer

# Install InterProScan
RUN \
    echo 'Downloading big interproscan file...' && \
    wget -nv ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.18-57.0/interproscan-5.18-57.0-64-bit.tar.gz && \
    tar xf interproscan-5.18-57.0-64-bit.tar.gz && \
    mv interproscan-5.18-57.0 /kb/deployment/interproscan && \
    echo 'export INTERPROSCAN_INSTALL=/kb/deployment/interproscan' >> /kb/deployment/user-env.sh && \
    echo 'export PATH=$PATH:$INTERPROSCAN_INSTALL' >> /kb/deployment/user-env.sh

# Copy local wrapper files, and build

COPY ./ /kb/module
RUN mkdir -p /kb/module/work

WORKDIR /kb/module

RUN make

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
