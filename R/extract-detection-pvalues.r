extract.detection.pvalues <- function(rg, probes, verbose=F) {
    msg(verbose=verbose)

    dyes <- c(R="R", G="G")
    negative.intensities <- lapply(dyes, function(color) {
        addresses <- with(probes, address[which(dye == color & target == "NEGATIVE")])
        rg[[color]][which(rownames(rg[[color]]) %in% addresses),"Mean"]
    })
    negative.med <- lapply(negative.intensities, median, na.rm=T)
    negative.sd  <- lapply(negative.intensities, mad, na.rm=T)

    features <- unique(na.omit(probes$name))

    detection.stats <- lapply(c(M="M", U="U"), function(target) {
        probes <- probes[which(probes$target == target),]
        probes <- probes[match(features, probes$name),]
        r.idx <- which(probes$dye == "R")
        g.idx <- which(probes$dye == "G")
        
        intensity <- rep(NA, nrow(probes))

        addresses <- probes$address[r.idx]
        if(!all(addresses %in% rownames(rg$R)))
            stop("It seems that the IDAT files do not match the supplied chip annotation")
        intensity[r.idx] <- rg$R[match(probes$address[r.idx], rownames(rg$R)),"Mean"]

        addresses <- probes$address[g.idx]
        if(!all(addresses %in% rownames(rg$G)))
            stop("It seems that the IDAT files do not match the supplied chip annotation")
        intensity[g.idx] <- rg$G[match(probes$address[g.idx], rownames(rg$G)),"Mean"]

        med <- rep(NA, nrow(probes))
        med[r.idx] <- negative.med$R
        med[g.idx] <- negative.med$G

        sd <- rep(NA, nrow(probes))
        sd[r.idx] <- negative.sd$R
        sd[g.idx] <- negative.sd$G

        data.frame(intensity=intensity, med=med, sd=sd)
    })

    p.value <- with(detection.stats, 1-pnorm(M$intensity + U$intensity,
                                             mean=M$med + U$med,
                                             sd=M$sd + U$sd))
    names(p.value) <- features
    p.value
}
