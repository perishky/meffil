require(illuminaio) ## for readIDAT
require(IlluminaHumanMethylation450kmanifest)
require(MASS) ## for huber
require(limma) ## lm.fit

meffil.extract.controls <- function(basename, probes=meffil.probe.info()) {
    msg("sample file", basename)
    rg <- read.rg(basename)
    extract.controls(rg, probes)
}

meffil.compute.normalization.object <- function(basename, control.matrix,
                                                number.quantiles=500,
                                                probes=meffil.probe.info()) {
    sample.idx <- match(basename, colnames(control.matrix))
    stopifnot(!is.na(sample.idx))
    dye.bias.factors <- calculate.dye.bias.factors(control.matrix, sample.idx)
        
    rg <- read.rg(basename)
    rg.correct <- background.correct(rg, probes)
    rg.correct <- dye.bias.correct(rg.correct, dye.bias.factors$R, dye.bias.factors$G)
    mu <- rg.to.mu(rg.correct, probes)

    probes.x <- probes[which(probes$chr == "chrX"),]
    probes.y <- probes[which(probes$chr == "chrY"),]
    x.signal <- median(log(mu$M[probes.x$name] + mu$U[probes.x$name], 2), na.rm=T)
    y.signal <- median(log(mu$M[probes.y$name] + mu$U[probes.y$name], 2), na.rm=T)

    probs <- seq(0,1,length.out=number.quantiles)
    quantile.sets <- define.quantile.probe.sets(probes)
    quantile.sets$names <- get.quantile.probe.sets(quantile.sets)
    quantile.sets$quantiles <- lapply(1:nrow(quantile.sets), function(i) {
        probe.names <- quantile.sets$names[[i]]
        target <- quantile.sets$target[i]
        quantile(mu[[target]][probe.names], probs=probs, na.rm=T)   
    })
    quantile.sets$names <- NULL

    list(origin="meffil.compute.normalization.object",
         basename=basename,
         quantile.sets=quantile.sets,
         dye.bias.factors=dye.bias.factors,
         x.signal=x.signal,
         y.signal=y.signal)         
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

extract.controls <- function(rg, probes=meffil.probe.info()) {
    stopifnot(is.rg(rg))

    msg()
    probes.G <- probes[which(probes$dye == "G"),]
    probes.R <- probes[which(probes$dye == "R"),]
    probes.G <- probes.G[match(names(rg$G), probes.G$address),]
    probes.R <- probes.R[match(names(rg$R), probes.R$address),]
    
    bisulfite2 <- mean(rg$R[which(probes.R$target == "BISULFITE CONVERSION II")], na.rm=T)
    
    bisulfite1.G <- rg$R[which(probes.R$target == "BISULFITE CONVERSION I"
                               & probes.R$ext
                               %in% sprintf("BS Conversion I%sC%s", c(" ", "-", "-"), 1:3))]
    ## minfi does this. shouldn't it be green, not red??????
    bisulfite1.R <- rg$R[which(probes.R$target == "BISULFITE CONVERSION I"
                               & probes.R$ext %in% sprintf("BS Conversion I-C%s", 4:6))]
    bisulfite1 <- mean(bisulfite1.G + bisulfite1.R, na.rm=T)
    
    stain.G <- rg$G[which(probes.G$target == "STAINING" & probes.G$ext == "Biotin (High)")]
    
    stain.R <- rg$R[which(probes.R$target == "STAINING" & probes.R$ext == "DNP (High)")]
    
    extension.R <- rg$R[which(probes.R$target == "EXTENSION"
                              & probes.R$ext %in% sprintf("Extension (%s)", c("A", "T")))]
    extension.G <- rg$G[which(probes.G$target == "EXTENSION"
                              & probes.G$ext %in% sprintf("Extension (%s)", c("C", "G")))]
    
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
    spec1.Rp <- rg$R[which(probes.R$target == "SPECIFICITY I" & probes.R$ext %in% ext)]
    spec1.ratio1 <- mean(spec1.Rp,na.rm=T)/mean(spec1.G,na.rm=T)
    
    ext <- sprintf("GT Mismatch %s (PM)", 4:6)
    spec1.Gp <- rg$G[which(probes.G$target == "SPECIFICITY I" & probes.G$ext %in% ext)]
    spec1.R <- rg$R[which(probes.R$target == "SPECIFICITY I" & probes.R$ext %in% ext)]
    spec1.ratio2 <- mean(spec1.Gp,na.rm=T)/mean(spec1.R,na.rm=T)
    ## color swap vs spec1.ratio1??? just following minfi but that seems weird
    
    spec1.ratio <- (spec1.ratio1 + spec1.ratio2)/2
    
    normA <- mean(rg$R[which(probes.R$target == "NORM_A")], na.rm = TRUE)
    normT <- mean(rg$R[which(probes.R$target == "NORM_T")], na.rm = TRUE)
    normC <- mean(rg$G[which(probes.G$target == "NORM_C")], na.rm = TRUE)
    normG <- mean(rg$G[which(probes.G$target == "NORM_G")], na.rm = TRUE)

    dye.bias <- (normC + normG)/(normA + normT)

    rg.bg <- background.correct(rg, probes)
    addresses <- probes.R$address[which(probes.R$target %in% c("NORM_A", "NORM_T"))]
    intensity.bc.R <- mean(rg.bg$R[addresses], na.rm = TRUE)
    addresses <- probes.G$address[which(probes.G$target %in% c("NORM_G", "NORM_C"))]
    intensity.bc.G <- mean(rg.bg$G[addresses], na.rm = TRUE)
    
    probs <- c(0.01, 0.5, 0.99)
    oob.G <- quantile(rg$G[with(probes.G, which(target == "OOB" & dye == "G"))], na.rm=T, probs=probs)
    oob.R <- quantile(rg$R[with(probes.R, which(target == "OOB" & dye == "R"))], na.rm=T, probs=probs)
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
      oob.ratio=oob.ratio,
      intensity.bc.G=intensity.bc.G,
      intensity.bc.R=intensity.bc.R)
}

calculate.dye.bias.factors <- function(control.matrix,sample.idx) {
    ratios <- control.matrix["intensity.bc.G",]/control.matrix["intensity.bc.R",]
    reference <- which.min(abs(ratios-1))
    intensity <- (control.matrix["intensity.bc.G",] + control.matrix["intensity.bc.R",])[reference]/2
    list(R=intensity/control.matrix["intensity.bc.R",sample.idx],
         G=intensity/control.matrix["intensity.bc.G",sample.idx])
}


meffil.probe.info <- function() {
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

    ret$type3 <- ret$type
    ret$type3[which(ret$type == "i" & ret$dye == "R")] <- "iR"
    ret$type3[which(ret$type == "i" & ret$dye == "G")] <- "iG"

    ret$chr.type <- ifelse(is.na(ret$chr), NA, "autosomal")
    ret$chr.type[which(ret$chr %in% c("chrX","chrY"))] <- "sex"

    for (col in setdiff(colnames(ret), "pos")) ret[,col] <- as.character(ret[,col])
    ret
}

define.quantile.probe.sets <- function(probes=meffil.probe.info()) {
    cbind(expand.grid(target=c("M","U"),
                      type3=c("iG","iR","ii"),
                      chr=NA,
                      chr.type=c(NA,"autosomal"),
                      stringsAsFactors=F),
          expand.grid(target=c("M","U"),
                      type3=NA,
                      chr="chrX",
                      chr.type="sex",
                      stringsAsFactors=F),
          expand.grid(target=c("M","U"),
                      type3=NA,
                      chr=NA,
                      chr.type="sex",
                      stringsAsFactors=F))
}

eq.wild <- function(x,y) {
    is.na(y) | x == y
}

get.quantile.probe.sets <- function(quantile.sets) {
    lapply(1:nrow(quantile.sets), function(i) {
        probes$name[which(eq.wild(probes$target, quantile.sets$target[i])
                          & eq.wild(probes$type3, quantile.sets$type3[i])
                          & eq.wild(probes$chr, quantile.sets$chr[i])
                          & eq.wild(probes$chr.type, quantile.sets$chr.type[i]))]
    })
}

select.normalization.subsets <- function(quantile.sets, sex="M", mixture=T) {
    (mixture & eq.wild("autosomal", quantile.sets$chr.type)
     | mixture & sex == "M" & eq.wild("sex", quantile.sets$chr.type)
     | mixture & sex == "F" & eq.wild("chrX", quantile.sets$chr)     
     | !mixture & sex == "M" & is.na(quantile.sets$chr.type)
     #| !mixture & sex == "M" & eq.wild("autosomal",quantile.sets$chr.type)
     #| !mixture & sex == "M" & eq.wild("sex",quantile.sets$chr.type)
     | !mixture & sex == "F" & eq.wild("autosomal", quantile.sets$chr.type)
     | !mixture & sex == "F" & eq.wild("chrX", quantile.sets$chr))
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
     && length(names(rg$G)) == length(rg$G)
     && length(names(rg$R)) == length(rg$R))
}

rg.to.mu <- function(rg, probes=meffil.probe.info()) {
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

background.correct <- function(rg, probes=meffil.probe.info(), offset=15) {
    stopifnot(is.rg(rg))
    
    lapply(c(R="R",G="G"), function(dye) {
        msg("background correction for dye =", dye)
        addresses <- probes$address[which(probes$target %in% c("M","U") & probes$dye == dye)]
        xf <- rg[[dye]][addresses]
        xf[which(xf <= 0)] <- 1

        addresses <- probes$address[which(probes$type == "control" & probes$dye == dye)]
        xc <- rg[[dye]][addresses]
        xc[which(xc <= 0)] <- 1
        
        addresses <- probes$address[which(probes$target == "OOB" & probes$dye == dye)]
        oob <- rg[[dye]][addresses]
        
        ests <- MASS::huber(oob) 
        mu <- ests$mu
        sigma <- log(ests$s)
        alpha <- log(max(MASS::huber(xf)$mu - mu, 10))
        bg <- limma::normexp.signal(as.numeric(c(mu,sigma,alpha)), c(xf,xc)) + offset
        names(bg) <- c(names(xf), names(xc))
        bg
    })
}

dye.bias.correct <- function(rg, factor.R, factor.G) {
    rg$R <- rg$R * factor.R
    rg$G <- rg$G * factor.G
    rg
}

is.normalization.object <- function(object) {
    (all(c("quantile.sets","dye.bias.factors","origin","basename","x.signal","y.signal")
         %in% names(object))
     && object$origin == "meffil.compute.normalization.object")
}

meffil.normalize.objects <- function(objects, control.matrix, 
                                    number.pcs=2, sex.cutoff=-2, sex=NULL,
                                    probes=meffil.probe.info()) {
    stopifnot(length(objects) == ncol(control.matrix))
    stopifnot(is.null(sex) || length(sex) == length(objects) && all(sex %in% c("F","M")))
    stopifnot(number.pcs >= 2)

    msg("cleaning up the control matrix")
    control.matrix <- control.matrix[-grep("^intensity.bc", rownames(control.matrix)),]
    control.matrix <- impute.matrix(control.matrix)
    control.matrix <- scale(t(control.matrix))
    control.matrix[control.matrix > 3] <- 3
    control.matrix[control.matrix < -3] <- -3
    control.matrix <- t(scale(control.matrix))

    if (is.null(sex)) {
        msg("predicting sex")
        x.signal <- sapply(objects, function(obj) obj$x.signal)
        y.signal <- sapply(objects, function(obj) obj$y.signal)
        xy.diff <- y.signal-x.signal
        sex <- ifelse(xy.diff < sex.cutoff, "F","M")
    }
    
    msg("normalizing quantiles")
    quantile.sets <- define.quantile.probe.sets(probes)
    quantile.sets$sex.diff <- (length(unique(sex)) >= 2
                               & with(quantile.sets, !is.na(chr.type) & chr.type != "autosomal"))
    normalized.quantiles <- lapply(1:nrow(quantile.sets), function(i) {
        original <- sapply(objects, function(obj) obj$quantile.sets$quantiles[[i]])        
        if (quantile.sets$sex.diff[i]) {
            norm <- original
            for (sex.value in unique(na.omit(sex))) {
                sample.idx <- which(sex == sex.value)
                norm[,sample.idx] <- normalize.quantiles(original[,sample.idx],
                                                         control.matrix[,sample.idx], number.pcs)
            }
            norm
        }
        else            
            normalize.quantiles(original, control.matrix, number.pcs)            
    })
    
    for (i in 1:length(objects)) {
        objects[[i]]$sex.cutoff <- sex.cutoff
        objects[[i]]$xy.diff <- xy.diff[i]
        objects[[i]]$sex <- sex[i]
        objects[[i]]$quantile.sets$sex.diff <- quantile.sets$sex.diff
        objects[[i]]$quantile.sets$norm <- lapply(normalized.quantiles,
                                                  function(sample.quantiles) sample.quantiles[,i])
    }
    objects
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
    quantiles[nrow(quantiles),] <- quantiles[nrow(quantiles)-1,] + 1000
    
    mean.quantiles <- rowMeans(quantiles)
    control.components <- prcomp(t(control.matrix))$x[,1:number.pcs,drop=F]
    design <- model.matrix(~control.components-1)
    fits <- lm.fit(x=design, y=t(quantiles - mean.quantiles))
    mean.quantiles - t(residuals(fits))
}


meffil.normalize.sample <- function(object, probes=meffil.probe.info()) {
    stopifnot(is.normalization.object(object))

    probe.names <- unique(na.omit(probes$name))

    U <- M <- rep(NA_integer_, length(probe.names))
    names(U) <- names(M) <- probe.names

    rg <- read.rg(object$basename)
    rg.correct <- background.correct(rg, probes)
    rg.correct <- dye.bias.correct(rg.correct, object$dye.bias.factors$R, object$dye.bias.factors$G)
    mu <- rg.to.mu(rg.correct, probes)

    mu$M <- mu$M[probe.names]
    mu$U <- mu$U[probe.names]

    object$quantile.sets$names <- get.quantile.probe.sets(object$quantile.sets)
    mixture <- sum(object$quantile.sets$sex.diff) == 0
    object$quantile.sets$apply <-select.normalization.subsets(object$quantile.sets,object$sex,mixture)

    for (i in which(object$quantile.sets$apply)) {
        target <- object$quantile.sets$target[i]
        probe.idx <- which(names(mu[[target]]) %in% object$quantile$names[[i]])

        orig.signal <- mu[[target]][probe.idx]
        norm.target <- compute.quantiles.target(object$quantile.sets$norm[[i]])
        norm.signal <- preprocessCore::normalize.quantiles.use.target(matrix(orig.signal),
                                                                      norm.target)
        mu[[target]][probe.idx] <- norm.signal
    }
    mu
}

meffil.normalize.samples <- function(objects, probes=meffil.probe.info()) {
    M <- U <- NA
    for (i in 1:length(objects)) {
        msg(i)
        mu <- meffil.normalize.sample(objects[[i]], probes)
        if (i == 1) {
            U <- M <- matrix(NA_integer_,
                             nrow=length(mu$M), ncol=length(objects),
                             dimnames=list(names(mu$M), names(objects)))
        }
        M[,i] <- mu$M
        U[,i] <- mu$U
    }
    list(M=M,U=U)
}

meffil.get.beta <- function(mu) mu$M/(mu$M+mu$U+100)


compute.quantiles.target <- function(quantiles) {
    n <- length(quantiles)
    unlist(lapply(1:(n-1), function(j) {
        start <- quantiles[j]
        end <- quantiles[j+1]
        seq(start,end,(end-start)/n)[-n]
    }))
}   


#############.....................################
check.sex <- function(xy.diff) {
    fit <- kmeans(xy.diff, centers=range(xy.diff))
    sex.kmeans <- ifelse(fit$cluster == which.min(fit$centers), "F", "M")        
}

