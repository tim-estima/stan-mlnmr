FROM nicholasjclark/brms

# system libraries of general use
RUN apt-get update && apt-get install -y \
    sudo \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev

# system library dependency for the euler app
# the R package Rmpfr requires the system library libmpfr-dev to be available
# note that if you don't need to use the Rmpfr package, you can delete this line
RUN apt-get update && apt-get install -y \
    libmpfr-dev \
    libgsl0-dev 

# install basic shiny functionality to R
RUN R -e "install.packages(c('tidyverse', 'shinystan', 'broom', 'randtoolbox', 'copula', \
'xtable', 'boot', 'paralell', 'logitnorm', 'broom.mixed', 'Rccp'))"

RUN R -e "install.packages(c('tidyverse', 'Rccp'))"

RUN mkdir /analysis
COPY  analysis /analysis

WORKDIR /analysis
RUN ls -1 
CMD Rscript psoriasis_run_analysis.R
