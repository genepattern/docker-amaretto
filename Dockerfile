# copyright 2017-2018 Regents of the University of California and the Broad Institute. All rights reserved.

FROM r-base:3.5.0

RUN mkdir /build

RUN echo "==================== 1 ============================"
RUN R --version
RUN apt-mark hold r-base

RUN apt-mark hold r-base && \
    apt-get update && apt-get autoclean && \
    apt-get install -t unstable libssl-dev  --yes && \
    apt-get install libxml2-dev --yes && \
    apt-get install libcurl4-gnutls-dev --yes && \
    apt-get install mesa-common-dev --yes && \ 
    apt-get update && apt-get install -y --no-install-recommends apt-utils && \
    apt-get install aptitude -y && \
    apt-get install libxml2-dev -y && \
    aptitude install libglib2.0-dev -y && \
    apt-get install libcairo2-dev -y && \
    apt-get install  xvfb xauth xfonts-base libxt-dev -y && \
    apt-get install pandoc -y  -t unstable git && \
    rm -rf /var/lib/apt/lists/*

RUN echo "==================== 2 ============================"
RUN R --version

COPY sources.list /etc/apt/sources.list
COPY Rprofile.gp.site ~/.Rprofile
COPY Rprofile.gp.site /usr/lib/R/etc/Rprofile.site
ENV R_LIBS=/usr/local/lib/R/site-library

RUN echo "==================== 3 ============================"
RUN R --version

RUN apt-mark hold r-base && \
   apt-get update  && \
   apt-get install -t unstable  libv8-3.14-dev --yes 

COPY src/install_stuff.R /usr/local/bin/amaretto/install_stuff.R
COPY src/callr.R /usr/local/bin/amaretto/callr.R
RUN echo "=================== 4 ============================="

RUN R --version


RUN   mkdir /source1 && \
   mv /source1 /source && \
   cd /source && \
   echo "i" && \
   wget https://install-github.me/r-lib/callr && \
   mv callr callr.R && \
   git clone https://github.com/gevaertlab/AMARETTO.git && \
   cd AMARETTO && \
   git checkout master
#   git checkout develop 
RUN echo "================ 5 ================================"

RUN R --version

RUN Rscript /usr/local/bin/amaretto/install_stuff.R

RUN echo "================= 6 ==============================="

RUN R --version

# now copy the rest of the files in case of a change
# and then the second install

COPY src/ /usr/local/bin/amaretto/
COPY R_environ ~/.Renviron
RUN echo "================ 7 ================================"
RUN R --version

# this is where amaretto is actually installed.
RUN Rscript /usr/local/bin/amaretto/install_amaretto.R


CMD ["Rscript", "/usr/local/bin/cogaps/run_amaretto_module.R" ]

