#' Matrix of SNP 'beta' values
#'
#' @param qc.objects List of objects obtained from \code{\link{meffil.qc}()}
#' or \code{\link{meffil.create.qc.object}()}.
#' 
#' @export 
meffil.snp.betas <- function(qc.objects) {
    stopifnot(sapply(qc.objects, is.qc.object))
    sapply(qc.objects, function(object) object$snp.betas)
}

#' Obtain the list of identifiers for the SNPs on the microarray.
#'
#' @param featureset Name from \code{\link{meffil.list.featuresets}()} (Default: "450k").
#' @export
meffil.snp.names <- function(featureset="450k") {
    features <- meffil.featureset(featureset)
    features$name[which(features$target == "snp")]
}

extract.snp.betas <- function(rg, probes, verbose=F) {
    msg(verbose=verbose)
    stopifnot(is.rg(rg))
    
    probes.M.R <- probes[which(probes$target == "M-snp" & probes$dye == "R"),]
    probes.M.G <- probes[which(probes$target == "M-snp" & probes$dye == "G"),]
    probes.U.R <- probes[which(probes$target == "U-snp" & probes$dye == "R"),]
    probes.U.G <- probes[which(probes$target == "U-snp" & probes$dye == "G"),]

    stopifnot(all(c(probes.M.R$address, probes.U.R$address) %in% rownames(rg$R)))
    stopifnot(all(c(probes.M.G$address, probes.U.G$address) %in% rownames(rg$G)))
    
    M <- c(rg$R[probes.M.R$address,"Mean"],
           rg$G[probes.M.G$address,"Mean"])
    U <- c(rg$R[probes.U.R$address,"Mean"],
           rg$G[probes.U.G$address,"Mean"])

    names(M) <- c(probes.M.R$name, probes.M.G$name)
    names(U) <- c(probes.U.R$name, probes.U.G$name)

    get.beta(M,U[names(M)])
}

