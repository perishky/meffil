# Dear Matt
# Please implement this to extract detectionp values AND bead number 
# Lots of love
# Gib and Josine
# P.S. Here is the minfi code that does it.

detectionP <- function(rgSet, type = "m+u") {
    locusNames <- getManifestInfo(rgSet, "locusNames")
    detP <- matrix(NA_real_, ncol = ncol(rgSet), nrow = length(locusNames),
                   dimnames = list(locusNames, sampleNames(rgSet)))

    controlIdx <- getControlAddress(rgSet, controlType = "NEGATIVE")   
    r <- getRed(rgSet)
    rBg <- r[controlIdx,]
    rMu <- matrixStats::colMedians(rBg)
    rSd <- matrixStats::colMads(rBg)

    g <- getGreen(rgSet)
    gBg <- g[controlIdx,]
    gMu <- matrixStats::colMedians(gBg)
    gSd <- matrixStats::colMads(gBg)

    TypeII <- getProbeInfo(rgSet, type = "II")
    TypeI.Red <- getProbeInfo(rgSet, type = "I-Red")
    TypeI.Green <- getProbeInfo(rgSet, type = "I-Green")
    for (i in 1:ncol(rgSet)) {   
        ## Type I Red
        intensity <- r[TypeI.Red$AddressA, i] + r[TypeI.Red$AddressB, i]
        detP[TypeI.Red$Name, i] <- 1-pnorm(intensity, mean=rMu[i]*2, sd=rSd[i]*2)
        ## Type I Green
        intensity <- g[TypeI.Green$AddressA, i] + g[TypeI.Green$AddressB, i]
        detP[TypeI.Green$Name, i] <- 1-pnorm(intensity, mean=gMu[i]*2, sd=gSd[i]*2)
        ## Type II
        intensity <- r[TypeII$AddressA, i] + g[TypeII$AddressA, i]
        detP[TypeII$Name, i] <- 1-pnorm(intensity, mean=rMu[i]+gMu[i], sd=rSd[i]+gSd[i])
    }
    detP
}



#' Normalization object
#'
#' Create a normalization object for a given Infinium HumanMethylation450 BeadChip.
#'
#' @param basename IDAT file basename (see \code{\link{meffil.basenames}}).
#' @param number.quantiles Number of quantiles to compute for probe subset (Default: 500).
#' @param dye.intensity Reference intensity for scaling each color channel (Default: 5000).
#' @param probes Probe annotation used to construct the control matrix
#' (Default: \code{\link{meffil.probe.info}()}).
#' @param verbose If \code{TRUE}, then status messages are printed during execution (Default: \code{FALSE}).
#' @param pvalue.threshold Default value = 0.01. All probes ABOVE this detection threhsold are returned
#' @param bead.threshold Default value = 3. All probes with less than this number of beads detected.
#' @return List containing control probe information, probe summaries
#' and quantiles.
#'
#' @export
meffil.compute.normalization.object <- function(basename,
                                                number.quantiles=500,
                                                dye.intensity=5000,
                                                probes=meffil.probe.info(),
                                                verbose=F,
                                                pvalue.threshold=0.01,
                                                bead.threshold=3) {
    stopifnot(number.quantiles >= 100)
    stopifnot(dye.intensity >= 100)
    stopifnot(nrow(probes) > 100000)

    rg <- meffil.read.rg(basename, probes, verbose=verbose)

    controls <- extract.controls(rg, probes, verbose=verbose)

    rg.correct <- meffil.background.correct(rg, probes, verbose=verbose)

    intensity.R <- calculate.intensity.R(rg.correct, probes)
    intensity.G <- calculate.intensity.G(rg.correct, probes)

    rg.correct <- meffil.dye.bias.correct(rg.correct, dye.intensity, probes, verbose=verbose)

    mu <- meffil.rg.to.mu(rg.correct, probes, verbose=verbose)

    probes.x <- probes$name[which(probes$chr == "chrX")]
    x.signal <- median(log(mu$M[probes.x] + mu$U[probes.x], 2), na.rm=T)

    probes.y <- probes$name[which(probes$chr == "chrY")]
    y.signal <- median(log(mu$M[probes.y] + mu$U[probes.y], 2), na.rm=T)

    probs <- seq(0,1,length.out=number.quantiles)

    quantiles <- lapply(get.quantile.probe.subsets(probes), function(sets) {
        list(M=unname(quantile(mu$M[sets$M], probs=probs,na.rm=T)),
             U=unname(quantile(mu$U[sets$U], probs=probs,na.rm=T)))
    })


    bad.probes.detectionp=array()
    bad.probes.beadnum=array()

    #extract SNP probes
    snp.probes=array()


    list(origin="meffil.compute.normalization.object",
         basename=basename,
         controls=controls,
         quantiles=quantiles,
         dye.intensity=dye.intensity,
         intensity.R=intensity.R,
         intensity.G=intensity.G,
         x.signal=x.signal,
         y.signal=y.signal,
         median.m.signal=median(mu$M),
         median.u.signal=median(mu$U),
         bad.probes.detectionp=bad.probes.detectionp,
         bad.probes.beadnum=bad.probes.beadnum,
         snp.probes=snp.probes
         )
}

is.normalization.object <- function(object) {
    (all(c("quantiles","dye.intensity","origin","basename","x.signal","y.signal","controls",
           "intensity.R","intensity.G")
         %in% names(object))
     && object$origin == "meffil.compute.normalization.object")
}

get.quantile.probe.subsets <- function(probes=meffil.probe.info()) {
    rm.na <- function(x) {
        x[which(is.na(x))] <- F
        x
    }

    is.iG <- rm.na(probes$type == "i" & probes$dye == "G")
    is.iR <- rm.na(probes$type == "i" & probes$dye == "R")
    is.ii <- rm.na(probes$type == "ii")
    is.genomic <- !is.na(probes$chr)
    is.sex <- rm.na(is.genomic & probes$chr %in% c("chrX","chrY"))
    is.x <- rm.na(is.genomic & probes$chr == "chrX")
    is.y <- rm.na(is.genomic & probes$chr == "chrY")
    is.autosomal <- rm.na(is.genomic & !is.sex)
    is.not.y <- rm.na(is.genomic & probes$chr != "chrY")

    get.probe.subsets <- function(in.subset) {
        list(M=probes$name[which(probes$target == "M" & in.subset)],
             U=probes$name[which(probes$target == "U" & in.subset)])
    }

    list(genomic.iG = get.probe.subsets(is.iG & is.genomic),
         genomic.iR = get.probe.subsets(is.iR & is.genomic),
         genomic.ii = get.probe.subsets(is.ii & is.genomic),
         autosomal.iG = get.probe.subsets(is.iG & is.autosomal),
         autosomal.iR = get.probe.subsets(is.iR & is.autosomal),
         autosomal.ii = get.probe.subsets(is.ii & is.autosomal),
         not.y.iG = get.probe.subsets(is.iG & is.not.y),
         not.y.iR = get.probe.subsets(is.iR & is.not.y),
         not.y.ii = get.probe.subsets(is.ii & is.not.y),
         sex = get.probe.subsets(is.sex),
         chry = get.probe.subsets(is.y),
         chrx = get.probe.subsets(is.x))
}

sex.specific.quantile.probe.subsets <- function() {
    c("genomic.iG",
      "genomic.iR",
      "genomic.ii",
      "not.y.iG",
      "not.y.iR",
      "not.y.ii",
      "sex",
      "chry",
      "chrx")
}

applicable.quantile.probe.subsets <- function(sex, both.sexes) {
    if (both.sexes && sex == "M") return(c("autosomal.iG","autosomal.iR","autosomal.ii","sex"))
    if (both.sexes && sex == "F") return(c("autosomal.iG","autosomal.iR","autosomal.ii","chrx","chry"))
    if (!both.sexes && sex == "M") return(c("genomic.iG", "genomic.iR", "genomic.ii"))
    if (!both.sexes && sex == "F") return(c("not.y.iG", "not.y.iR", "not.y.ii","chry"))
    stop("invalid input", "sex =", sex, "both.sexes =", both.sexes)
}



extract.controls <- function(rg, probes=meffil.probe.info(), verbose=F) {
    stopifnot(is.rg(rg))

    x.mean <- function(x, na.rm=T) {
        stopifnot(length(x) > 1)
        mean(x,na.rm=na.rm)
    }
    x.which <- function(x) {
        i <- which(x)
        stopifnot(length(i) > 0)
        i
    }

    msg(verbose=verbose)
    probes.G <- probes[x.which(probes$dye == "G"),]
    probes.R <- probes[x.which(probes$dye == "R"),]
    probes.G <- probes.G[match(names(rg$G), probes.G$address),]
    probes.R <- probes.R[match(names(rg$R), probes.R$address),]

    bisulfite2 <- x.mean(rg$R[x.which(probes.R$target == "BISULFITE CONVERSION II")])

    bisulfite1.G <- rg$G[x.which(probes.G$target == "BISULFITE CONVERSION I"
                               & probes.G$ext
                               %in% sprintf("BS Conversion I%sC%s", c(" ", "-", "-"), 1:3))]
    bisulfite1.R <- rg$R[x.which(probes.R$target == "BISULFITE CONVERSION I"
                               & probes.R$ext %in% sprintf("BS Conversion I-C%s", 4:6))]
    bisulfite1 <- x.mean(bisulfite1.G + bisulfite1.R)

    stain.G <- rg$G[x.which(probes.G$target == "STAINING" & probes.G$ext == "Biotin (High)")]

    stain.R <- rg$R[x.which(probes.R$target == "STAINING" & probes.R$ext == "DNP (High)")]

    extension.R <- rg$R[x.which(probes.R$target == "EXTENSION"
                              & probes.R$ext %in% sprintf("Extension (%s)", c("A", "T")))]
    extension.G <- rg$G[x.which(probes.G$target == "EXTENSION"
                              & probes.G$ext %in% sprintf("Extension (%s)", c("C", "G")))]

    hybe <- rg$G[x.which(probes.G$target == "HYBRIDIZATION")]

    targetrem <- rg$G[x.which(probes.G$target %in% "TARGET REMOVAL")]

    nonpoly.R <- rg$R[x.which(probes.R$target == "NON-POLYMORPHIC"
                            & probes.R$ext %in% sprintf("NP (%s)", c("A", "T")))]

    nonpoly.G <- rg$G[x.which(probes.G$target == "NON-POLYMORPHIC"
                            & probes.G$ext %in% sprintf("NP (%s)", c("C", "G")))]

    spec2.G <- rg$G[x.which(probes.G$target == "SPECIFICITY II")]
    spec2.R <- rg$R[x.which(probes.R$target == "SPECIFICITY II")]
    spec2.ratio <- x.mean(spec2.G,na.rm=T)/x.mean(spec2.R,na.rm=T)

    ext <- sprintf("GT Mismatch %s (PM)", 1:3)
    spec1.G <- rg$G[x.which(probes.G$target == "SPECIFICITY I" & probes.G$ext %in% ext)]
    spec1.Rp <- rg$R[x.which(probes.R$target == "SPECIFICITY I" & probes.R$ext %in% ext)]
    spec1.ratio1 <- x.mean(spec1.Rp,na.rm=T)/x.mean(spec1.G,na.rm=T)

    ext <- sprintf("GT Mismatch %s (PM)", 4:6)
    spec1.Gp <- rg$G[x.which(probes.G$target == "SPECIFICITY I" & probes.G$ext %in% ext)]
    spec1.R <- rg$R[x.which(probes.R$target == "SPECIFICITY I" & probes.R$ext %in% ext)]
    spec1.ratio2 <- x.mean(spec1.Gp,na.rm=T)/x.mean(spec1.R,na.rm=T)

    spec1.ratio <- (spec1.ratio1 + spec1.ratio2)/2

    normA <- x.mean(rg$R[x.which(probes.R$target == "NORM_A")], na.rm = TRUE)
    normT <- x.mean(rg$R[x.which(probes.R$target == "NORM_T")], na.rm = TRUE)
    normC <- x.mean(rg$G[x.which(probes.G$target == "NORM_C")], na.rm = TRUE)
    normG <- x.mean(rg$G[x.which(probes.G$target == "NORM_G")], na.rm = TRUE)

    dye.bias <- (normC + normG)/(normA + normT)

    probs <- c(0.01, 0.5, 0.99)
    oob.G <- quantile(rg$G[with(probes.G, x.which(target == "OOB"))], na.rm=T, probs=probs)
    oob.R <- quantile(rg$R[with(probes.R, x.which(target == "OOB"))], na.rm=T, probs=probs)
    oob.ratio <- oob.G[["50%"]]/oob.R[["50%"]]

    c(bisulfite1=bisulfite1,
      bisulfite2=bisulfite2,
      extension.G=extension.G,
      extension.R=extension.R,
      hybe=hybe,
      stain.G=stain.G,
      stain.R=stain.R,
      nonpoly.G=nonpoly.G,
      nonpoly.R=nonpoly.R,
      targetrem=targetrem,
      spec1.G=spec1.G,
      spec1.R=spec1.R,
      spec2.G=spec2.G,
      spec2.R=spec2.R,
      spec1.ratio1=spec1.ratio1,
      spec1.ratio=spec1.ratio,
      spec2.ratio=spec2.ratio,
      spec1.ratio2=spec1.ratio2,
      normA=normA,
      normC=normC,
      normT=normT,
      normG=normG,
      dye.bias=dye.bias,
      oob.G=oob.G,
      oob.ratio=oob.ratio)
}

