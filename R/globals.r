.onLoad <- function(libname, pkgname) {
    assign("pkg.globals", new.env(), envir=parent.env(environment()))

    probe.info.filename <- system.file("probe-info.rda", package="meffil")
    if (file.exists(probe.info.filename))
        load(probe.info.filename, pkg.globals)

    assign("reference.globals", new.env(), envir=parent.env(environment()))
    
    references.filename <- system.file("gse35069-references.rda", package="meffil")
    if (file.exists(references.filename))
        load(references.filename, reference.globals)
}

## To create a global variable:
## assign("name", value, pkg.globals)
## To retrieve a global variable:
## get("name", pkg.globals)
## To see if it exists:
## exists("name", pkg.globals)

