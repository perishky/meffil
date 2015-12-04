## library(plyr)
## library(meffil)
## library(splitstackshape)
## library(matrixStats)

#' Calculate CNVs from IDAT files
#'
#' Based on the algorithm developed in \code{R/CopyNumber450k} bioconductor package
#' 
#' @param samplesheet Output from \code{meffil.create.samplesheet}
#' @param cnv.reference Name returned by \code{\link{meffil.list.cnv.references}()}.
#' @param chip Name returned by \code{\link{meffil.list.chips()}} (Default: \code{NULL}).
#' @param verbose Default = FALSE
#' @param ... Extra parameters to be passed to \code{DNAcopy} for segmentation. See details.
#' 
#' @details
#' The following default values are being used:
#' - trim = 0.1
#' - min.width = 5
#' - nperm = 10000
#' - alpha = 0.001
#' - undo.splits = "sdundo"
#' - undo.SD = 2
#'
#' @export
#' @return Dataframe of segmented results
meffil.calculate.cnv <- function(samplesheet, cnv.reference, chip=NULL, verbose=FALSE, ...) {
    cnv.reference <- meffil:::get.cnv.reference(cnv.reference)
    l1 <- meffil:::get.index.list(nrow(samplesheet), options("mc.cores")[[1]])
    
    l <- lapply(l1, function(x) {		
        l2 <- mclapply(x, function(i) {
            calculate.cnv(samplesheet$Basename[i],
                          samplesheet$Sample_Name[i],
                          cnv.reference=cnv.reference,
                          chip=chip,
                          verbose=verbose, ...)
        })
        names(l2) <- samplesheet$Sample_Name[x]
        return(l2) 
    })
    nom <- unlist(sapply(l, names))
    l <- unlist(l, recursive=FALSE)
    names(l) <- nom
    return(l) # list of data frames (output of calculate.cnv), one per sample
}

calculate.cnv <- function(bname, samplename=basename(bname), cnv.reference,
                          chip=NULL, 
                          trim=0.1, min.width = 5, nperm = 10000, alpha = 0.001,
                          undo.splits = "sdundo", undo.SD = 2, verbose = TRUE, smoothing=TRUE) {
    msg("Reading idat file for", bname, verbose=verbose)
    rg <- meffil:::read.rg(bname, verbose=verbose)
    
    chip <- meffil:::guess.chip(rg, chip)
    
    featureset <- cnv.reference$featureset
    
    probes <- meffil.probe.info(chip, featureset)
    case <- meffil:::mu.to.cn(meffil:::rg.to.mu(rg, probes))
    
    msg("Predicting sex", verbose=verbose)
    sex <- meffil:::cnv.predict.sex(case, featureset)

    if (ncol(cnv.reference$intensity.sex[[sex]]) == 0)
        stop(paste("Sample is predicted to be", 
                   ifelse(sex == "M", "male", "female"),
                   "but none found in CNV reference."))
    
    msg("Normalisting against controls", verbose=verbose)
    int.sex <- case[rownames(cnv.reference$intensity.sex[[sex]])]
    int.sex[order(int.sex)] <- sort(cnv.reference$intensity.sex[[sex]][,1])

    int.aut <- case[rownames(cnv.reference$intensity.aut)]
    int.aut[order(int.aut)] <- sort(cnv.reference$intensity.aut[,1])
    
    msg("Estimating CNVs", verbose=verbose)
    int.sex.prop <- log2(int.sex / cnv.reference$control.medians.sex[[sex]])
    int.aut.prop <- log2(int.aut / cnv.reference$control.medians.aut)

    features <- meffil:::cnv.features(featureset)
    
    int.sex.prop <- int.sex.prop[names(int.sex.prop) %in% features$name]
    int.aut.prop <- int.aut.prop[names(int.aut.prop) %in% features$name]

    p.sex <- features[match(names(int.sex.prop), features$name),c("name","chromosome","position")]
    p.aut <- features[match(names(int.aut.prop), features$name),c("name","chromosome","position")]

    cna.sex <- CNA(int.sex.prop, chrom=p.sex$chromosome, maploc=p.sex$position,
                   data.type="logratio", sampleid=samplename)
    cna.aut <- CNA(int.aut.prop, chrom=p.aut$chromosome, maploc=p.aut$position,
                   data.type="logratio", sampleid=samplename)
    
    if(smoothing) {
        cna.sex <- smooth.CNA(cna.sex, trim=trim)
        cna.aut <- smooth.CNA(cna.aut, trim=trim)
    }
    
    segment.cna.sex <- segment(cna.sex, min.width=min.width, verbose = verbose, nperm = nperm,
                               alpha = alpha, undo.splits = undo.splits, undo.SD = undo.SD,
                               trim = trim)
    segment.cna.aut <- segment(cna.aut, min.width=min.width, verbose = verbose, nperm = nperm,
                               alpha = alpha, undo.splits = undo.splits, undo.SD = undo.SD,
                               trim = trim)
    out <- rbind(segment.cna.aut$output, segment.cna.sex$output)
    out$chrom <- ordered(out$chrom, levels = paste("chr", c(1:22, "X","Y"), sep=""))
    out$ID <- samplename

    msg("Calculating p-values", verbose=verbose)
    out <- meffil:::cnv.pvalue(out, int.aut, int.sex, sex, cnv.reference)
    return(out)
}

get.index.list <- function(n, mc.cores) {
    mc.cores <- ifelse(is.null(mc.cores) || mc.cores < 1, 1, min(mc.cores, n))
    div <- floor(n / mc.cores)
    rem <- n %% mc.cores
    l1 <- lapply(1:div, function(x) (x-1) * mc.cores + 1:mc.cores)
    if(rem != 0) l1[[div+1]] <- l1[[div]][mc.cores] + 1:rem
    return(l1)
}

mu.to.cn <- function(mu) {
    mu$M + mu$U
}

cnv.predict.sex <- function(cn, featureset, sex.cutoff=-2) {
    log.cn <- log2(cn)
    
    sites.x <- meffil.get.x.sites(featureset)
    stopifnot(all(sites.x %in% names(log.cn)))
    x.signal <- median(log.cn[sites.x], na.rm=T)
    
    sites.y <- meffil.get.y.sites(featureset)    
    stopifnot(all(sites.y %in% names(log.cn)))
    y.signal <- median(log.cn[sites.y], na.rm=T)
    
    xy.diff <- y.signal-x.signal
    ifelse(xy.diff < sex.cutoff, "F", "M")
}

cnv.pvalue <- function(out, int.aut, int.sex, sex, cnv.reference) {
    ## Get the sum of the intensities for the probes in each of the controls
    ## Get the mean sum of intensities
    ## Get the SD of the sum of intensities
    ## Get the sum of intensities of probes in the target
    features <- cnv.features(cnv.reference$featureset)
    for(i in 1:nrow(out)) {
        feature.sub <- subset(features, (chromosome == out$chrom[i]
                                         & position <= out$loc.end[i]
                                         & position >= out$loc.start[i]))$name
        if(out$chrom[i] %in% c("chrX", "chrY")) {
            control.int.sum <- colSums(cnv.reference$intensity.sex[[sex]][feature.sub, , drop=FALSE])
            sample.sum <- sum(int.sex[feature.sub])
        } else {
            control.int.sum <- colSums(cnv.reference$intensity.aut[feature.sub, , drop=FALSE])
            sample.sum <- sum(int.aut[feature.sub])
        }
        control.mean <- mean(control.int.sum)
        control.sd <- sd(control.int.sum)
        
        z.score <- (sample.sum - control.mean)/control.sd
        pval <- 2 * pnorm(-abs(z.score))
        
        ## Compute p-value (2 sided t-test)
        out$nprobe[i] <- length(feature.sub)
        out$sample.sum[i] <- sample.sum
        out$control.mean[i] <- control.mean
        out$control.sd[i] <- control.sd
        out$z.score[i] <- z.score
        out$pvalue[i] <- pval
    }
    out$adjusted.pvalue <- p.adjust(out$pvalue, method="bonferroni")
    return(out)
}

cnv.features <- function(featureset) {
    features <- meffil.get.features(featureset)
    sites <- features[which(!is.na(features$chromosome)
                            & features$target == "methylation"
                            & !features$snp.exclude),]
    sites$chromosome <- ordered(sites$chromosome, levels=paste("chr", c(1:22, "X", "Y"), sep=""))
    sites <- sites[with(sites, order(chromosome, position, decreasing=F)),]
    sites
}
