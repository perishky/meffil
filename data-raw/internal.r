library(devtools)

library(meffil)
options(mc.cores=10)

probe.info <- meffil:::collate.probe.info()
## I tried to make probe.info a global variable that
## is computed only the first time after meffil.probe.info()
## is called after the package is loaded.
## I ran into two problems:
## 
## 1. Global variables in packages are frozen after the package is loaded.
## The solution here is to create a global environment when the package loads.
## Just add the following function to the code in the R directory:
## 
## .onLoad <- function(libname, pkgname) {
##     assign("pkg.globals", new.env(), envir=parent.env(environment()))
## }
##
## Global variables can be created and assigned using "assign('var-name', value, pkg.globals)"
## and retrieved using "get('var-name', pkg.globals)".
##
## 2. If variable is assigned or created within a function called by mclapply(),
## then the change is reversed after mclapply exits.
## Parallelization is one of the points of the package so another solution was required.
## I then decided to precompute the object and load it when referenced
## (LazyData is set to 'yes' in the package DESCRIPTION file).

gse35069.references <- meffil:::create.gse35069.references(probes=probe.info, verbose=T)

use_data(gse35069.references, probe.info, internal=TRUE, overwrite=TRUE)

