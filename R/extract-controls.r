extract.controls <- function(rg, probes, verbose=F) {
    stopifnot(is.rg(rg))

    x.mean <- function(x, na.rm=T) {
        if (length(x) <= 1)
            stop("It seems that the IDAT files do not match the supplied chip annotation.")
        mean(x,na.rm=na.rm)
    }
    x.which <- function(x) {
        i <- which(x)
        if (length(i) == 0)
            stop("It seems that the IDAT files do not match the supplied chip annotation")
        i
    }

    msg(verbose=verbose)

    probes.G <- probes[x.which(probes$dye == "G"),]
    probes.R <- probes[x.which(probes$dye == "R"),]
    rg$R <- rg$R[match(probes.R$address, rownames(rg$R)),]
    rg$G <- rg$G[match(probes.G$address, rownames(rg$G)),]

    bisulfite2 <- x.mean(rg$R[x.which(probes.R$target == "BISULFITE CONVERSION II"), "Mean"])

    bisulfite1.G <- rg$G[x.which(probes.G$target == "BISULFITE CONVERSION I"
                               & probes.G$name
                               %in% sprintf("BS Conversion I%sC%s", c(" ", "-", "-"), 1:3)),"Mean"]
    bisulfite1.R <- rg$R[x.which(probes.R$target == "BISULFITE CONVERSION I"
                               & probes.R$name %in% sprintf("BS Conversion I-C%s", 4:6)),"Mean"]
    bisulfite1 <- x.mean(bisulfite1.G + bisulfite1.R)

    stain.G <- rg$G[x.which(probes.G$target == "STAINING" & probes.G$name == "Biotin (High)"),"Mean"]

    stain.R <- rg$R[x.which(probes.R$target == "STAINING" & probes.R$name == "DNP (High)"),"Mean"]

    extension.R <- rg$R[x.which(probes.R$target == "EXTENSION"
                              & probes.R$name %in% sprintf("Extension (%s)", c("A", "T"))),"Mean"]
    extension.G <- rg$G[x.which(probes.G$target == "EXTENSION"
                              & probes.G$name %in% sprintf("Extension (%s)", c("C", "G"))),"Mean"]

    hybe <- rg$G[x.which(probes.G$target == "HYBRIDIZATION"),"Mean"]

    targetrem <- rg$G[x.which(probes.G$target %in% "TARGET REMOVAL"),"Mean"]

    nonpoly.R <- rg$R[x.which(probes.R$target == "NON-POLYMORPHIC"
                            & probes.R$name %in% sprintf("NP (%s)", c("A", "T"))),"Mean"]

    nonpoly.G <- rg$G[x.which(probes.G$target == "NON-POLYMORPHIC"
                            & probes.G$name %in% sprintf("NP (%s)", c("C", "G"))),"Mean"]

    spec2.G <- rg$G[x.which(probes.G$target == "SPECIFICITY II"),"Mean"]
    spec2.R <- rg$R[x.which(probes.R$target == "SPECIFICITY II"),"Mean"]
    spec2.ratio <- x.mean(spec2.G,na.rm=T)/x.mean(spec2.R,na.rm=T)

    name <- sprintf("GT Mismatch %s (PM)", 1:3)
    spec1.G <- rg$G[x.which(probes.G$target == "SPECIFICITY I" & probes.G$name %in% name),"Mean"]
    spec1.Rp <- rg$R[x.which(probes.R$target == "SPECIFICITY I" & probes.R$name %in% name),"Mean"]
    spec1.ratio1 <- x.mean(spec1.Rp,na.rm=T)/x.mean(spec1.G,na.rm=T)

    name <- sprintf("GT Mismatch %s (PM)", 4:6)
    spec1.Gp <- rg$G[x.which(probes.G$target == "SPECIFICITY I" & probes.G$name %in% name),"Mean"]
    spec1.R <- rg$R[x.which(probes.R$target == "SPECIFICITY I" & probes.R$name %in% name),"Mean"]
    spec1.ratio2 <- x.mean(spec1.Gp,na.rm=T)/x.mean(spec1.R,na.rm=T)

    spec1.ratio <- (spec1.ratio1 + spec1.ratio2)/2

    normA <- x.mean(rg$R[x.which(probes.R$target == "NORM_A"),"Mean"], na.rm = TRUE)
    normT <- x.mean(rg$R[x.which(probes.R$target == "NORM_T"),"Mean"], na.rm = TRUE)
    normC <- x.mean(rg$G[x.which(probes.G$target == "NORM_C"),"Mean"], na.rm = TRUE)
    normG <- x.mean(rg$G[x.which(probes.G$target == "NORM_G"),"Mean"], na.rm = TRUE)

    dye.bias <- (normC + normG)/(normA + normT)

    probs <- c(0.01, 0.5, 0.99)
    oob.G <- quantile(rg$G[with(probes.G, x.which(target == "OOB")),"Mean"], na.rm=T, probs=probs)
    oob.R <- quantile(rg$R[with(probes.R, x.which(target == "OOB")),"Mean"], na.rm=T, probs=probs)
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


