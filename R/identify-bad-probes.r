identify.bad.probes.beadnum <- function(rg, threshold=3, verbose=F) {
    msg(verbose=verbose)

    probes <- meffil.probe.info()
    
    cpg.sites <- unique(na.omit(probes$name))
 
    targets <- c(M="M", U="U")
    n.beads <- sapply(targets, function(target) {
        probes <- probes[which(probes$target == target),]
        probes <- probes[match(cpg.sites, probes$name),]
        r.idx <- which(probes$dye == "R")
        g.idx <- which(probes$dye == "G")        
        n.beads <- rep(NA, nrow(probes))
        n.beads[r.idx] <- rg$R[match(probes$address[r.idx], rownames(rg$R)), "NBeads"]
        n.beads[g.idx] <- rg$G[match(probes$address[g.idx], rownames(rg$G)), "NBeads"]
        n.beads
    })
    n.beads <- pmin(n.beads[,"U"], n.beads[,"M"])
    names(n.beads) <- cpg.sites
    n.beads[which(n.beads < threshold)]
}

identify.bad.probes.detectionp <- function(rg, threshold=0.01, verbose=F) {
    msg(verbose=verbose)

    probes <- meffil.probe.info()
    
    dyes <- c(R="R", G="G")
    negative.intensities <- lapply(dyes, function(color) {
        addresses <- with(probes, address[which(dye == color & target == "NEGATIVE")])
        rg[[color]][which(rownames(rg[[color]]) %in% addresses),"Mean"]
    })
    negative.med <- lapply(negative.intensities, median, na.rm=T)
    negative.sd  <- lapply(negative.intensities, mad, na.rm=T)

    cpg.sites <- unique(na.omit(probes$name))

    targets <- c(M="M", U="U")
    detection.stats <- lapply(targets, function(target) {
        probes <- probes[which(probes$target == target),]
        probes <- probes[match(cpg.sites, probes$name),]
        r.idx <- which(probes$dye == "R")
        g.idx <- which(probes$dye == "G")

        intensity <- rep(NA, nrow(probes))
        intensity[r.idx] <- rg$R[match(probes$address[r.idx], rownames(rg$R)),"Mean"]
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
    names(p.value) <- cpg.sites
    p.value[which(p.value > threshold)]
}


