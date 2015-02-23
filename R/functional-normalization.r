#library(illuminaio) ## for readIDAT()
#library(IlluminaHumanMethylation450kmanifest) ## for getProbeInfo()
#library(MASS) ## for huber()
#library(IlluminaHumanMethylation450kanno.ilmn12.hg19) ## for probe locations



meffil.basenames <- function(path,recursive=FALSE) {    
    grn.files <- list.files(path, pattern = "_Grn.idat$", recursive = recursive, 
                            ignore.case = TRUE, full.names = TRUE)
    red.files <- list.files(path, pattern = "_Red.idat$", recursive = recursive, 
                            ignore.case = TRUE, full.names = TRUE)
    intersect(sub("_Grn.idat$", "", grn.files), 
              sub("_Red.idat$", "", red.files))
}

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

meffil.probe.info <- function(array="IlluminaHumanMethylation450k",annotation="ilmn12.hg19") {    
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

    locations <- probe.locations(array, annotation)
    ret <- cbind(ret, locations[match(ret$name, rownames(locations)),])

    for (col in setdiff(colnames(ret), "pos")) ret[,col] <- as.character(ret[,col])
    ret
}

meffil.read.rg <- function(basename, probes=meffil.probe.info()) {
    rg <- list(G=read.idat(paste(basename, "_Grn.idat", sep = "")),
               R=read.idat(paste(basename, "_Red.idat", sep="")))
    rg$R <- rg$R[which(names(rg$R) %in% probes$address[which(probes$dye == "R")])]
    rg$G <- rg$G[which(names(rg$G) %in% probes$address[which(probes$dye == "G")])]
    rg
}

meffil.rg.to.mu <- function(rg, probes=meffil.probe.info()) {
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

meffil.background.correct <- function(rg, probes=meffil.probe.info(), offset=15) {
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

calculate.intensity.R <- function(rg, probes=meffil.probe.info()) {
    addresses <- probes$address[which(probes$target %in% c("NORM_A", "NORM_T")
                                      & probes$dye == "R")]
    mean(rg$R[match(addresses, names(rg$R))], na.rm=T)
}

calculate.intensity.G <- function(rg, probes=meffil.probe.info()) {
    addresses <- probes$address[which(probes$target %in% c("NORM_G", "NORM_C")
                                      & probes$dye == "G")]
    mean(rg$G[match(addresses, names(rg$G))], na.rm=T)
}

meffil.dye.bias.correct <- function(rg, intensity=5000, probes=meffil.probe.info()) {
    msg()
    stopifnot(is.rg(rg))
    
    rg$R <- rg$R * intensity/calculate.intensity.R(rg, probes)
    rg$G <- rg$G * intensity/calculate.intensity.G(rg, probes)
    rg
}


meffil.compute.normalization.object <- function(basename, 
                                                number.quantiles=500,
                                                dye.intensity=5000,
                                                probes=meffil.probe.info()) {        
    rg <- meffil.read.rg(basename, probes)

    controls <- extract.controls(rg, probes)
 
    rg.correct <- meffil.background.correct(rg, probes)

    intensity.R <- calculate.intensity.R(rg.correct, probes)
    intensity.G <- calculate.intensity.G(rg.correct, probes)

    rg.correct <- meffil.dye.bias.correct(rg.correct, dye.intensity, probes)
    
    mu <- meffil.rg.to.mu(rg.correct, probes)

    probes.x <- probes$name[which(probes$chr == "chrX")]
    x.signal <- median(log(mu$M[probes.x] + mu$U[probes.x], 2), na.rm=T)

    probes.y <- probes$name[which(probes$chr == "chrY")]
    y.signal <- median(log(mu$M[probes.y] + mu$U[probes.y], 2), na.rm=T)
 
    probs <- seq(0,1,length.out=number.quantiles)
    
    quantiles <- lapply(get.quantile.probe.subsets(probes), function(sets) {
        list(M=quantile(mu$M[sets$M], probs=probs,na.rm=T),
             U=quantile(mu$U[sets$U], probs=probs,na.rm=T))
    })

    list(origin="meffil.compute.normalization.object",
         basename=basename,
         controls=controls,
         quantiles=quantiles,
         dye.intensity=dye.intensity,
         intensity.R=intensity.R,
         intensity.G=intensity.G,
         x.signal=x.signal,
         y.signal=y.signal)         
}

meffil.preprocess.control.matrix <- function(control.matrix) {
    control.matrix <- impute.matrix(control.matrix)
    control.matrix <- scale(t(control.matrix))
    control.matrix[control.matrix > 3] <- 3
    control.matrix[control.matrix < -3] <- -3
    t(scale(control.matrix))
}

meffil.normalize.objects <- function(objects, 
                                    number.pcs=2, sex.cutoff=-2, sex=NULL) {
    stopifnot(is.null(sex) || length(sex) == length(objects) && all(sex %in% c("F","M")))
    stopifnot(number.pcs >= 2)

    msg("preprocessing the control matrix")
    control.matrix <- sapply(objects, function(object) object$controls)
    control.matrix <- meffil.preprocess.control.matrix(control.matrix)

    msg("selecting dye correction reference")
    intensity.R <- sapply(objects, function(object) object$intensity.R)
    intensity.G <- sapply(objects, function(object) object$intensity.G)
    reference.idx <- which.min(abs(intensity.R/intensity.G-1))
    dye.intensity <- (intensity.R + intensity.G)[reference.idx]/2
                          
    msg("predicting sex")
    x.signal <- sapply(objects, function(obj) obj$x.signal)
    y.signal <- sapply(objects, function(obj) obj$y.signal)
    xy.diff <- y.signal-x.signal
    predicted.sex <- ifelse(xy.diff < sex.cutoff, "F","M")
    if (is.null(sex))
        sex <- predicted.sex
    
    sex.summary <- table(sex)
    has.both.sexes <- length(sex.summary) >= 2 & min(sex.summary) > 1

    male.idx <- which(sex == "M")
    female.idx <- which(sex == "F")

    msg("normalizing quantiles")
    subset.names <- names(objects[[1]]$quantiles)
    normalized.quantiles <- sapply(subset.names, function(name) {
        sapply(c("M","U"), function(target) {
            msg(name, target)
            original <- sapply(objects, function(object) {
                object$quantiles[[name]][[target]] * dye.intensity/object$dye.intensity
            })

            if (name %in% sex.specific.quantile.probe.subsets() && has.both.sexes) {
                norm.male <- normalize.quantiles(original[,male.idx,drop=F],
                                                 control.matrix[,male.idx,drop=F],
                                                 number.pcs)
                norm.female <- normalize.quantiles(original[,female.idx,drop=F],
                                                   control.matrix[,female.idx,drop=F],
                                                   number.pcs)
                norm <- original
                norm[,male.idx] <- norm.male
                norm[,female.idx] <- norm.female
            }
            else
                norm <- normalize.quantiles(original, control.matrix, number.pcs)
            norm
        }, simplify=F)
    }, simplify=F)
        
    lapply(1:length(objects), function(i) {
        object <- objects[[i]]
        object$sex.cutoff <- sex.cutoff
        object$xy.diff <- xy.diff[i]
        object$predicted.sex <- predicted.sex[i]
        object$sex <- sex[i]
        object$reference.intensity <- dye.intensity

        subset.names <- applicable.quantile.probe.subsets(object$sex, has.both.sexes)
        object$norm <- sapply(subset.names, function(subset.name) {
            list(M=normalized.quantiles[[subset.name]]$M[,i],
                 U=normalized.quantiles[[subset.name]]$U[,i])
        },simplify=F)
        
        object
    })
}

meffil.normalize.sample <- function(object, probes=meffil.probe.info()) {
    stopifnot(is.normalization.object(object))

    probe.names <- unique(na.omit(probes$name))

    rg <- meffil.read.rg(object$basename, probes)
    rg.correct <- meffil.background.correct(rg, probes)
    rg.correct <- meffil.dye.bias.correct(rg.correct, object$reference.intensity, probes)
    mu <- meffil.rg.to.mu(rg.correct, probes)

    mu$M <- mu$M[probe.names]
    mu$U <- mu$U[probe.names]

    msg("Normalizing methylated and unmethylated signals.")
    probe.subsets <- get.quantile.probe.subsets(probes)
    for (name in names(object$norm)) {
        for (target in names(object$norm[[name]])) {
            probe.idx <- which(names(mu[[target]]) %in% probe.subsets[[name]][[target]])
            orig.signal <- mu[[target]][probe.idx]
            norm.target <- compute.quantiles.target(object$norm[[name]][[target]])            
            norm.signal <- preprocessCore::normalize.quantiles.use.target(matrix(orig.signal),
                                                                          norm.target)
            mu[[target]][probe.idx] <- norm.signal
        }
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

meffil.get.beta <- function(mu, pseudo=100) {
    mu$M/(mu$M+mu$U+pseudo)
}


msg <- function(..., verbose=T) {
    x <- paste(list(...))
    name <- sys.call(sys.parent(1))[[1]]
    cat(paste("[", name, "]", sep=""), date(), x, "\n")
}


extract.controls <- function(rg, probes=meffil.probe.info()) {
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
    
    spec1.ratio <- (spec1.ratio1 + spec1.ratio2)/2
    
    normA <- mean(rg$R[which(probes.R$target == "NORM_A")], na.rm = TRUE)
    normT <- mean(rg$R[which(probes.R$target == "NORM_T")], na.rm = TRUE)
    normC <- mean(rg$G[which(probes.G$target == "NORM_C")], na.rm = TRUE)
    normG <- mean(rg$G[which(probes.G$target == "NORM_G")], na.rm = TRUE)

    dye.bias <- (normC + normG)/(normA + normT)
    
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
      oob.ratio=oob.ratio)
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

is.normalization.object <- function(object) {
    (all(c("quantiles","dye.intensity","origin","basename","x.signal","y.signal","controls",
           "intensity.R","intensity.G")
         %in% names(object))
     && object$origin == "meffil.compute.normalization.object")
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
    safe.increment <- 1000 
    quantiles[nrow(quantiles),] <- quantiles[nrow(quantiles)-1,] + safe.increment
    
    mean.quantiles <- rowMeans(quantiles)
    control.components <- prcomp(t(control.matrix))$x[,1:number.pcs,drop=F]
    design <- model.matrix(~control.components-1)
    fits <- lm.fit(x=design, y=t(quantiles - mean.quantiles))
    mean.quantiles + t(residuals(fits))
}

compute.quantiles.target <- function(quantiles) {
    n <- length(quantiles)
    unlist(lapply(1:(n-1), function(j) {
        start <- quantiles[j]
        end <- quantiles[j+1]
        seq(start,end,(end-start)/n)[-n]
    }))
}   


meffil.normalize.dataset <- function(basenames, data.dir, recursive=F, number.pcs=2,
                                     sex=NULL, probes=meffil.probe.info()) {
    stopifnot(missing(basenames) && missing(data.dir))
    
    if (missing(basenames))
        basenames <- meffil.basenames(data.dir, recursive)
    
    norm.objects <- lapply(basenames, function(basename) {
        meffil.compute.normalization.object(basename, probes=probes)
    })

    norm.objects <- meffil.normalize.objects(norm.objects, number.pcs=number.pcs, sex=sex)

    sapply(norm.objects, function(object) {
        meffil.get.beta(meffil.normalize.sample(object, probes)) 
    })
}
