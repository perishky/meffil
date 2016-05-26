#' List of available copy number references
#'
#' @export
meffil.list.cnv.references <- function() {
    ls(cnv.globals)
}
 
get.cnv.reference <- function(name) {
    stopifnot(is.character(name) && name %in% meffil.list.cnv.references())
    get(name, cnv.globals)
}

#' Create a copy number reference object
#'
#' Create a copy number reference object for estimating copy number variation
#' with the Infinium HumanMethylation450 BeadChip.
#'
#' @param name Character string providing the name of the reference.
#' @param M Matrix of methylated probe intensities (rows=CpG sites, columns=samples).
#' @param U Matrix of unmethylatd probe intensities (rows=CpG sites, columns=samples).
#' @param chip Name returned by \code{\link{meffil.list.chips()}} (Default: NA).
#' @param featureset Name returned by \code{\link{meffil.list.featuresets()}}
#' (Default: \code{chip}).
#' @param object A previously created copy number reference object created by this function.
#' If not \code{NULL}, then this reference is added with the given \code{name}
#' and all other function arguments are ignored (Default: NULL).
#' @param verbose If \code{TRUE}, then status messages are printed during execution
#' (Default: \code{FALSE}).
#' @return A list specifying a copy number reference object that can be used by
#' \code{\link{meffil.calculate.cnv}()} to estimate copy number variation in another dataset.
#' @export
meffil.add.cnv.reference <- function(name, M, U, chip=NA, featureset=chip, object=NULL, verbose=T) {
    if (is.null(object)) {
        object <- create.cnv.reference(M, U, chip, featureset, verbose)
        object$name <- name
    }
    else
        stopifnot(is.cnv.reference(object))
    assign(object$name, object, cnv.globals)
    invisible(object)
}

is.cnv.reference <- function(object) {
    "class" %in% names(object) && object[[class]] == "cnv.reference"
}

#' Create copy number references from CopyNumber450kData
#'
#' Two copy number references are created using data from the
#' Bioconductor CopyNumber450kData R package.
#' Reference "copynumber450k" is created using the "450k" feature set,
#' and reference "copynumber450k-common" is created using the "common" feature set
#' so it can be used with datasets with mixed 450K and EPIC chips.
#' 
#' @export
meffil.add.copynumber450k.references <- function(verbose=T) {
    require(minfi)
    require(CopyNumber450kData)

    data(RGcontrolSetEx)
    rgset <- preprocessIllumina(RGcontrolSetEx, bg.correct=TRUE, normalize=NULL)
    M <- getMeth(rgset)
    U <- getUnmeth(rgset)

    list("copynumber450k"=meffil.add.cnv.reference("copynumber450k", M, U,
                             chip="450k", featureset="450k", verbose=verbose),
         "copynumber450k-common"=meffil.add.cnv.reference("copynumber450k-common", M, U,
             chip="450k", featureset="common", verbose=verbose))
}



create.cnv.reference <- function(M, U, chip=NA, featureset=chip, verbose=T) {
    chip <- guess.chip(M, chip)
    
    if (is.na(featureset))
        featureset <- chip

    stopifnot(is.compatible.chip(featureset, chip))

    stopifnot(nrow(M) == nrow(U))
    stopifnot(ncol(M) == ncol(M))
    stopifnot(all(rownames(M) == rownames(U)))

    cn <- M + U

    msg("Predicting sex", verbose=verbose)
    sex <- apply(cn, 2, meffil:::cnv.predict.sex, featureset=featureset)

    if (length(table(sex)) < 2)
        warning("Only one sex is represented among the CNV reference samples")

    msg("Performing quantile normalisation", verbose=verbose)
    sites.sex <- c(meffil.get.x.sites(featureset),
                   meffil.get.y.sites(featureset))
    sites.aut <- meffil.get.autosomal.sites(featureset)

    int.sex.m <- cn[sites.sex,sex=="M"]
    int.sex.f <- cn[sites.sex,sex=="F"]
    int.aut <- cn[sites.aut,]

    nom.sex.m <- dimnames(int.sex.m)
    nom.sex.f <- dimnames(int.sex.f)
    nom.aut <- dimnames(int.aut)

    if (ncol(int.sex.m) > 0)
        int.sex.m <- preprocessCore::normalize.quantiles(int.sex.m)
    if (ncol(int.sex.f) > 0)
        int.sex.f <- preprocessCore::normalize.quantiles(int.sex.f)
    int.aut <- preprocessCore::normalize.quantiles(int.aut)

    dimnames(int.sex.m) <- nom.sex.m
    dimnames(int.sex.f) <- nom.sex.f
    dimnames(int.aut) <- nom.aut
    
    msg("Calculating control medians", verbose=verbose)
    control.medians.sex.m <- rowMedians(int.sex.m, na.rm=T)
    control.medians.sex.f <- rowMedians(int.sex.f, na.rm=T)
    control.medians.aut <- rowMedians(int.aut, na.rm=T)

    return(list(
        class="cnv.reference",
        version=packageVersion("meffil"),
        featureset=featureset,
        chip=chip,
        intensity.sex=list(M=int.sex.m, F=int.sex.f), 
        intensity.aut=int.aut,
        sex=sex,
        control.medians.sex=list(M=control.medians.sex.m, F=control.medians.sex.f),
        control.medians.aut=control.medians.aut
	))
}
