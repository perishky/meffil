extract.beadnum <- function(rg, probes, verbose=F) {
    msg(verbose=verbose)

    features <- unique(na.omit(probes$name))
 
    n.beads <- sapply(c(M="M", U="U"), function(target) {
        probes <- probes[which(probes$target == target),]
        probes <- probes[match(features, probes$name),]
        r.idx <- which(probes$dye == "R")
        g.idx <- which(probes$dye == "G")        
        n.beads <- rep(NA, nrow(probes))

        addresses <- probes$address[r.idx]
        if(!all(addresses %in% rownames(rg$R)))
            stop("It seems that the IDAT files do not match the supplied chip annotation")
        n.beads[r.idx] <- rg$R[match(addresses, rownames(rg$R)), "NBeads"]

        addresses <- probes$address[g.idx]
        if(!all(addresses %in% rownames(rg$G)))
            stop("It seems that the IDAT files do not match the supplied chip annotation")
        n.beads[g.idx] <- rg$G[addresses, "NBeads"]
        
        n.beads
    })
    n.beads <- pmin(n.beads[,"U"], n.beads[,"M"])
    names(n.beads) <- features
    n.beads
}

