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

extract.snp.probe.betas <- function(rg, verbose=F) {
    msg(verbose=verbose)

    probes <- meffil.probe.info()

    probes.M.R <- probes[which(probes$target == "MG" & probes$dye == "R"),]
    probes.M.G <- probes[which(probes$target == "MG" & probes$dye == "G"),]
    probes.U.R <- probes[which(probes$target == "UG" & probes$dye == "R"),]
    probes.U.G <- probes[which(probes$target == "UG" & probes$dye == "G"),]

    M <- c(rg$R[probes.M.R$address,"Mean"],
           rg$G[probes.M.G$address,"Mean"])
    U <- c(rg$R[probes.U.R$address,"Mean"],
           rg$G[probes.U.G$address,"Mean"])

    names(M) <- c(probes.M.R$name, probes.M.G$name)
    names(U) <- c(probes.U.R$name, probes.U.G$name)

    get.beta(M,U[names(M)])
}
