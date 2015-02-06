require(illuminaio) ## for readIDAT
require(IlluminaHumanMethylation450kmanifest)
require(MASS) ## for huber
require(limma) ## lm.fit


meffil.extract.control.matrix <- function(basenames, probes=probe.info()) {
    basenames <- sub("_Grn\\.idat$", "", basenames)
    basenames <- sub("_Red\\.idat$", "", basenames)
    basenames <- unique(basenames)
    
    msg("Samples \n", paste(basenames, collapse="\n "))
    
    sapply(basenames, function(basename) {
        rg <- read.rg(basename)
        collect.controls(rg, probes)
    })
}


meffil.compute.normalization.object <- function(basenames, control.matrix,
                                                number.pcs=2, number.quantiles=500,
                                                probes=probe.info()) {
    stopifnot(length(basenames) == ncol(control.matrix))
    stopifnot(number.pcs >= 2)
    stopifnot(number.quantiles >= 100)

    colnames(control.matrix) <- basenames

    msg("calculating reference dye intensity") ## must happen before the control matrix gets scaled
    dye.reference <- which.min(abs(control.matrix["dye.bias",]-1))
    dye.intensity <- mean(control.matrix[c("normA","normT","normC","normG"), dye.reference])
    
    msg("cleaning up the control matrix")
    control.matrix <- impute.matrix(control.matrix)
    control.matrix <- scale(t(control.matrix))
    control.matrix[control.matrix > 3] <- 3
    control.matrix[control.matrix < -3] <- -3
    control.matrix <- t(scale(control.matrix))

    quantile.sets <- expand.grid(target=c("M","U"), type=c("iR","iG","ii"),
                                        stringsAsFactors=F)
    
    msg("compute probe quantiles")
    quantiles <- sapply(basenames, function(basename) {
        rg <- read.rg(basename)
        rg.correct <- background.correct(rg, probes)        
        rg.correct <- dye.bias.correct(rg.correct, dye.intensity, probes)
        mu <- rg.to.mu(rg.correct, probes)

        sapply(1:nrow(quantile.sets), function(i) {
            target <- quantile.sets$target[i]
            type <- quantile.sets$type[i]
            msg("computing quantiles of", target, type, "probes")
            probe.names <- probes$name[which(probes$target == target
                                             & probes$type3 == type
                                             & probes$chr.type == "autosomal")]
            probs <- seq(0,1,length.out=number.quantiles)
            quantile(mu[[target]][probe.names], probs=probs, na.rm=T)
        })
    }, simplify=F)

    msg("merging corresponding sample quantiles into matrices")
    quantile.sets$quantiles <- lapply(1:nrow(quantile.sets), function(i) 
                                   sapply(quantiles, function(sample) sample[,i]))
    
    msg("normalizing quantiles")
    quantile.sets$norm <- lapply(quantile.sets$quantiles, function(original) { 
        normalize.quantiles(original, control.matrix, number.pcs)
    })
    
    list(origin="meffil.compute.normalization.object",
         basenames=basenames,
         control.matrix=control.matrix,
         quantile.sets=quantile.sets,
         dye=list(reference=names(dye.reference), intensity=dye.intensity))
}

is.normalization.object <- function(object) {
    (all(c("control.matrix", "quantile.sets","dye","origin","basenames") %in% names(object))
     && all(colnames(object$control.matrix) == object$basenames)
     && object$origin == "meffil.compute.normalization.object")
}

meffil.normalize <- function(object,probes=probe.info()) {
    stopifnot(is.normalization.object(object))


    probe.names <- unique(na.omit(probes$name))

    U <- M <- matrix(NA_integer_,
                     ncol=length(object$basenames),
                     nrow=length(probe.names),
                     dimnames=list(probe.names,object$basenames))
    
    x.ignore <- sapply(object$basenames, function(basename) {
        sample.idx <- which(object$basenames == basename)
        rg <- read.rg(basename)
        rg.correct <- background.correct(rg, probes)
        rg.correct <- dye.bias.correct(rg.correct, object$dye$intensity, probes)
        mu <- rg.to.mu(rg.correct, probes)
        
        mu$M <- mu$M[probe.names]
        mu$U <- mu$U[probe.names]

        for (i in 1:nrow(object$quantile.sets)) {
            target <- object$quantile.sets$target[i]
            type <- object$quantile.sets$type[i]
            msg("normalizing", target, type, "probes for", basename)

            probe.idx <- which(names(mu[[target]])
                                     %in% probes$name[which(probes$target == target
                                                            & probes$type3 == type)])

            is.autosomal <- (names(mu[[target]])[probe.idx]
                             %in% probes$name[which(probes$chr.type=="autosomal")])
            
            mu[[target]][probe.idx] <- normalize.sample(mu[[target]][probe.idx],
                                                        object$quantile.sets$quantiles[[i]][,sample.idx],
                                                        object$quantile.sets$norm[[i]][,sample.idx],
                                                        is.autosomal)
        }
        U[names(mu$U),basename] <<- mu$U
        M[names(mu$M),basename] <<- mu$M
        NULL
    })

    list(M=M,U=U)

}


msg <- function(..., verbose=T) {
    x <- paste(list(...))
    name <- sys.call(sys.parent(1))[[1]]
    cat(paste("[", name, "]", sep=""), date(), x, "\n")
}

read.rg <- function(basename) {
    rg <- list(G=read.idat(paste(basename, "_Grn.idat", sep = "")),
               R=read.idat(paste(basename, "_Red.idat", sep="")))
}

probe.info <- function() {
    probe.locations <- function(array="IlluminaHumanMethylation450k",annotation="ilmn12.hg19") {
        annotation <- paste(array, "anno.", annotation, sep="")

        msg("loading probe genomic location annotation", annotation)
            
        require(annotation,character.only=T)
        data(list=annotation)
        as.data.frame(get(annotation)@data$Locations)
    }
    probe.characteristics <- function(type) {
        msg("extracting", type)
        getProbeInfo(IlluminaHumanMethylation450kmanifest, type=type)
    }
    
    type1.R <- probe.characteristics("I-Red")
    type1.G <- probe.characteristics("I-Green")
    type2 <- probe.characteristics("II")
    controls <- probe.characteristics("Control")

    msg("reorganizing type information")
    ret <- rbind(data.frame(type="i",target="M", dye="R", address=type1.R$AddressB, name=type1.R$Name,ext=NA),
                 data.frame(type="i",target="M", dye="G", address=type1.G$AddressB, name=type1.G$Name,ext=NA),
                 data.frame(type="ii",target="M", dye="G", address=type2$AddressA, name=type2$Name,ext=NA),
                 
                 data.frame(type="i",target="U", dye="R", address=type1.R$AddressA, name=type1.R$Name,ext=NA),
                 data.frame(type="i",target="U", dye="G", address=type1.G$AddressA, name=type1.G$Name,ext=NA),
                 data.frame(type="ii",target="U", dye="R", address=type2$AddressA, name=type2$Name,ext=NA),
                 
                 data.frame(type="i",target="OOB", dye="G", address=type1.R$AddressA, name=NA,ext=NA),
                 data.frame(type="i",target="OOB", dye="G", address=type1.R$AddressB, name=NA,ext=NA),
                 data.frame(type="i",target="OOB", dye="R", address=type1.G$AddressA, name=NA,ext=NA),
                 data.frame(type="i",target="OOB", dye="R", address=type1.G$AddressB, name=NA,ext=NA),
                 
                 data.frame(type="control",target=controls$Type,dye="R",address=controls$Address, name=NA,ext=controls$ExtendedType),
                 data.frame(type="control",target=controls$Type,dye="G",address=controls$Address, name=NA,ext=controls$ExtendedType))

    for (col in setdiff(colnames(ret), "pos")) ret[,col] <- as.character(ret[,col])

    locations <- probe.locations()
    ret <- cbind(ret, locations[match(ret$name, rownames(locations)),])
    
    ret$chr.type <- ret$chr
    ret$chr.type[with(ret,which(chr.type !="chrX" & chr.type !="chrY"))]<-"autosomal"
    
    ret$type3 <- ret$type
    ret$type3[which(ret$type == "i" & ret$dye == "R")] <- "iR"
    ret$type3[which(ret$type == "i" & ret$dye == "G")] <- "iG"

    for (col in setdiff(colnames(ret), "pos")) ret[,col] <- as.character(ret[,col])
    ret
}


read.idat <- function(filename) {
    msg("Reading", filename)
    
    if (!file.exists(filename))
        stop("Filename does not exist:", filename)
    readIDAT(filename)$Quants[,"Mean"]
}

is.rg <- function(rg) {
    (all(c("R","G") %in% names(rg))
     && is.vector(rg$R) && is.vector(rg$G)
     #&& length(rg$R) == length(rg$G)
     && length(names(rg$G)) == length(rg$G)
     && length(names(rg$R)) == length(rg$R))
     ##&& all(names(rg$R) == names(rg$G)))
}

rg.to.mu <- function(rg, probes=probe.info()) {
    stopifnot(is.rg(rg))

    msg("converting red/green to methylated/unmethylated signal")
    probes.M.R <- probes[which(probes$target == "M" & probes$dye == "R"),]
    probes.M.G <- probes[which(probes$target == "M" & probes$dye == "G"),]
    probes.U.R <- probes[which(probes$target == "U" & probes$dye == "R"),]
    probes.U.G <- probes[which(probes$target == "U" & probes$dye == "G"),]
    
    M <- c(rg$R[probes.M.R$address], rg$G[probes.M.G$address])
    U <- c(rg$R[probes.U.R$address], rg$G[probes.U.G$address])
    
    names(M) <- c(probes.M.R$name, probes.M.G$name)
    names(U) <- c(probes.U.R$name, probes.U.G$name)
    
    U <- U[names(M)]
    list(M=M,U=U)
}

background.correct <- function(rg, probes=probe.info(), offset=15) {
    stopifnot(is.rg(rg))
    
    lapply(c(R="R",G="G"), function(dye) {
        msg("background correction for dye =", dye)
        addresses <- probes$address[which(probes$target != "OOB" & probes$dye == dye)]
        xf <- rg[[dye]][addresses]
        xf[which(xf <= 0)] <- 1
        
        addresses <- probes$address[which(probes$target == "OOB" & probes$dye == dye)]
        oob <- rg[[dye]][addresses]
        
        ests <- MASS::huber(oob) 
        mu <- ests$mu
        sigma <- log(ests$s)
        alpha <- log(max(MASS::huber(xf)$mu - mu, 10))
        xf.bkg <- limma::normexp.signal(as.numeric(c(mu,sigma,alpha)), xf) + offset
        names(xf.bkg) <- names(xf)
        xf.bkg
    })
}

dye.bias.correct <- function(rg, reference, probes=probe.info()) {
    stopifnot(is.rg(rg))
    msg()
    addresses <- probes$address[which(probes$target %in% c("NORM_A","NORM_T")
                                      & probes$dye == "R")]
    factor.R <- mean(rg$R[addresses],na.rm=T)
    
    addresses <- probes$address[which(probes$target %in% c("NORM_C","NORM_G")
                                      & probes$dye == "G")]
    factor.G <- mean(rg$G[addresses],na.rm=T)

    rg$R <- rg$R * reference/factor.R
    rg$G <- rg$G * reference/factor.G
    rg
}

collect.controls <- function(rg, probes=probe.info()) {
    stopifnot(is.rg(rg))

    msg()
    
    probes.G <- probes[which(probes$dye == "G"),]
    probes.R <- probes[which(probes$dye == "R"),]
    probes.G <- probes.G[match(names(rg$G), probes.G$address),]
    probes.R <- probes.R[match(names(rg$R), probes.R$address),]
    
    bisulfite2 <- mean(rg$R[which(probes.R$target == "BISULFITE CONVERSION II")], na.rm=T)
    
    bisulfite1.G <- rg$G[which(probes.G$target == "BISULFITE CONVERSION I"
                               & probes.G$ext
                               %in% sprintf("BS Conversion I%sC%s", c(" ", "-", "-"), 1:3))]
    bisulfite1.R <- rg$R[which(probes.R$target == "BISULFITE CONVERSION I"
                               & probes.R$ext %in% sprintf("BS Conversion I-C%s", 4:6))]
    bisulfite1 <- mean(bisulfite1.G + bisulfite1.R, na.rm=T)
    
    stain.G <- rg$G[which(probes.G$target == "STAINING" & probes.G$ext == "Biotin (High)")]
    
    stain.R <- rg$R[which(probes.R$target == "STAINING" & probes.R$ext == "DNP (High)")]
    
    extension.R <- rg$R[which(probes.R$target == "EXTENSION"
                              & probes.R$ExternalType %in% sprintf("Extension (%s)", c("A", "T")))]
    extension.G <- rg$G[which(probes.G$target == "EXTENSION"
                              & probes.G$ExternalType %in% sprintf("Extension (%s)", c("C", "G")))]
    
    hybe <- rg$G[which(probes.G$target == "HYBRIDIZATION")]
    
    targetrem <- rg$G[which(probes.G$target %in% "TARGET REMOVAL")]
    
    nonpoly.R <- rg$R[which(probes.R$target == "NON-POLYMORPHIC"
                            & probes.R$ext %in% sprintf("NP (%s)", c("A", "T")))]
    
    nonpoly.G <- rg$G[which(probes.G$target == "NON-POLYMORPHIC"
                            & probes.G$ext %in% sprintf("NP (%s)", c("C", "G")))]
    
    spec2.G <- rg$G[which(probes.G$target == "SPECIFICITY II")]
    spec2.R <- rg$R[which(probes.R$target == "SPECIFICITY II")]
    spec2.ratio <- mean(spec2.G,na.rm=T)/mean(spec2.R,na.rm=T)
    
    ext <- sprintf("GT Mismatch %s (PM)", 1:3)
    spec1.G <- rg$G[which(probes.G$target == "SPECIFICITY I" & probes.G$ext %in% ext)]
    spec1.R <- rg$R[which(probes.R$target == "SPECIFICITY I" & probes.R$ext %in% ext)]
    spec1.ratio1 <- mean(spec1.R,na.rm=T)/mean(spec2.G,na.rm=T)
    
    ext <- sprintf("GT Mismatch %s (PM)", 4:6)
    spec1.G <- rg$G[which(probes.G$target == "SPECIFICITY I" & probes.G$ext %in% ext)]
    spec1.R <- rg$R[which(probes.R$target == "SPECIFICITY I" & probes.R$ext %in% ext)]
    spec1.ratio2 <- mean(spec1.R,na.rm=T)/mean(spec2.G,na.rm=T)
    
    spec1.ratio <- (spec1.ratio1 + spec1.ratio2)/2
    
    normA <- mean(rg$R[which(probes.R$target == "NORM_A")], na.rm = TRUE)
    normT <- mean(rg$R[which(probes.R$target == "NORM_T")], na.rm = TRUE)
    normC <- mean(rg$G[which(probes.G$target == "NORM_C")], na.rm = TRUE)
    normG <- mean(rg$G[which(probes.G$target == "NORM_G")], na.rm = TRUE)

    dye.bias <- (normA + normT)/(normC + normG)
    
    probs <- c(0.01, 0.5, 0.99)
    oob.G <- quantile(rg$G[with(probes.G, which(target == "OOB" & dye == "G"))], na.rm=T, probs=probs)
    oob.R <- quantile(rg$R[with(probes.R, which(target == "OOB" & dye == "R"))], na.rm=T, probs=probs)
    oob.ratio <- oob.G[["50%"]]/oob.R[["50%"]]
    
    model.matrix <- c(bisulfite1=bisulfite1,
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




impute.matrix <- function(x, FUN=function(x) mean(x, na.rm=T)) {
    idx <- which(is.na(x), arr.ind=T)
    if (length(idx) > 0) {
        na.rows <- unique(idx[,"row"])
        v <- apply(x[na.rows,],1,FUN)
        v[which(is.na(v))] <- FUN(v) ## if any row imputation is NA ...
        x[idx] <- v[match(idx[,"row"],na.rows)]
    }
    x
}


normalize.quantiles <- function(quantiles, control.matrix, number.pcs) {
    stopifnot(is.matrix(quantiles))
    stopifnot(is.matrix(control.matrix))
    stopifnot(ncol(quantiles) == ncol(control.matrix))
    stopifnot(number.pcs >= 2)
    
    quantiles[1,] <- 0
    ## quantiles[nrow(quantiles),] <- quantiles[nrow(quantiles)-1,] + 1000
    
    mean.quantiles <- rowMeans(quantiles)
    control.components <- prcomp(t(control.matrix))$x[,1:number.pcs,drop=F]
    design <- model.matrix(~control.components-1)
    fits <- lm.fit(x=design, y=t(quantiles - mean.quantiles))
    mean.quantiles - t(residuals(fits))
}


compute.quantiles.target <- function(quantiles) {
    n <- length(quantiles)
    unlist(lapply(1:(n-1), function(j) {
        start <- quantiles[j]
        end <- quantiles[j+1]
        seq(start,end,(end-start)/n)[-n]
    }))
 }   

normalize.sample <- function(orig.signal, orig.quantiles, norm.quantiles, primary) {
    stopifnot(length(orig.signal) == length(primary))
    stopifnot(length(orig.quantiles) == length(norm.quantiles))
    msg()
    
    norm.target <- compute.quantiles.target(norm.quantiles)

    norm.signal <- orig.signal
    norm.signal[primary] <- preprocessCore::normalize.quantiles.use.target(matrix(orig.signal[primary]), norm.target)
    if (sum(!primary) > 0) {
        orig.target <- compute.quantiles.target(orig.quantiles)
        orig.intervals <- findInterval(orig.signal[!primary], orig.target)
        norm.signal[!primary] <- norm.target[orig.intervals]
    }
    norm.signal
}

        


