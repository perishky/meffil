.onLoad <- function(libname, pkgname) {
    assign("probe.globals", new.env(), envir=parent.env(environment()))
    assign("featureset.globals", new.env(), envir=parent.env(environment()))
    assign("reference.globals", new.env(), envir=parent.env(environment()))
    assign("cnv.globals", new.env(), envir=parent.env(environment()))
    assign("features.globals", new.env(), envir=parent.env(environment()))
    load.globals()
}

load.globals <- function() {
    load.env <- function(filename, env) {
        if (file.exists(filename))
            load(filename, env)
    }

    load.env(system.file("probes.rda", package="meffil"), probe.globals)
    load.env(system.file("featuresets.rda", package="meffil"), featureset.globals)
    load.env(system.file("references.rda", package="meffil"), reference.globals)
    load.env(system.file("cnv.rda", package="meffil"), cnv.globals)
    assign("features", get.all.features(), features.globals)
}

# called by ../data-raw/globals.r to save generated global variables to Rdata files
# for loading whenever the package is loaded.
save.globals <- function(dir) {
    require(remotes)
    require(pkgload)
    save.env <- function(filename, env) {
        save(list=ls(env),
             file=filename,
             envir=env)
        file.copy(filename, inst("meffil"), overwrite=T)
    }

    save.env(file.path(dir, "probes.rda"), probe.globals)
    save.env(file.path(dir, "featuresets.rda"), featureset.globals)
    save.env(file.path(dir, "references.rda"), reference.globals)
    save.env(file.path(dir, "cnv.rda"), cnv.globals)
}
