## In R:
 install.packages(c("MASS", 
                   "ggplot2",
                   "plyr",
                   "reshape2",
                   "knitr",
                   "Cairo",
                   "gridExtra",
                   "markdown",
                   "matrixStats",
                   "multcomp",
                   "lme4",
                   "parallel",
                   "fastICA",
                   "quadprog",
                   "betareg"))
 source("http://bioconductor.org/biocLite.R")
 biocLite(c("illuminaio",
           "limma",
           "sva",
           "DNAcopy"))

## From the command line:
 wget https://github.com/perishky/meffil/archive/master.zip
 unzip master.zip
 mv meffil-master meffil
 R CMD INSTALL meffil

