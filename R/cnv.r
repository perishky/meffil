## library(CopyNumber450k)
## library(CopyNumber450kData)
## library(plyr)
## library(meffil)
## library(splitstackshape)


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
  

#' Create controls for calculating CNVs
#' @param featureset Name returned by \code{\link{meffil.list.featuresets()}} (Default: \code{"450k"}).
#' Must be compatible with the "450k" chip.
#' @param verbose (Default: FALSE)
#' @export
#' @return List with pre-calculated CNV data
meffil.cnv.controls <- function(featureset=NULL, verbose=FALSE)
{
        chip <- "450k" ## because the control data below comes from 450K microarrays

        if (is.null(featureset))
            featureset <- chip
        
        if (!is.compatible.chip(featureset, chip))
            stop(paste("feature set '", featureset, "' is not compatible with the control chip '", chip, "'", sep=""))

	require(CopyNumber450kData)
            
	data(RGcontrolSetEx)
	x <- preprocessRaw(RGcontrolSetEx)
	x <- getMeth(x) + getUnmeth(x)

	msg("Predicting sex", verbose=verbose)
	s <- apply(x, 2, cnv.predict.sex, featureset=featureset)

	# Quantile normalize controls

	msg("Performing quantile normalisation", verbose=verbose)
	sites.sex <- c(meffil.get.x.sites(featureset),
                        meffil.get.y.sites(featureset))
	sites.aut <- meffil.get.autosomal.sites(featureset)

	int_sex_m <- x[sites.sex,s=="M"]
	int_sex_f <- x[sites.sex,s=="F"]
	int_aut <- x[sites.aut,]

	nom_sex_m <- list(rownames(int_sex_m), colnames(int_sex_m))
	nom_sex_f <- list(rownames(int_sex_f), colnames(int_sex_f))
	nom_aut <- list(rownames(int_aut), colnames(int_aut))

	int_sex_m <- preprocessCore::normalize.quantiles(int_sex_m)
	int_sex_f <- preprocessCore::normalize.quantiles(int_sex_f)
	int_aut <- preprocessCore::normalize.quantiles(int_aut)

	rownames(int_sex_m) <- nom_sex_m[[1]]
	colnames(int_sex_m) <- nom_sex_m[[2]]
	rownames(int_sex_f) <- nom_sex_f[[1]]
	colnames(int_sex_f) <- nom_sex_f[[2]]
	rownames(int_aut) <- nom_aut[[1]]
	colnames(int_aut) <- nom_aut[[2]]

	# Get medians for controls

	msg("Calculating control medians", verbose=verbose)
	control_medians_sex_m <- apply(int_sex_m, 1, median, na.rm=T)
	control_medians_sex_f <- apply(int_sex_f, 1, median, na.rm=T)
	control_medians_aut <- apply(int_aut, 1, median, na.rm=T)

	return(list(
                featureset=featureset,
                chip=chip,
		intensity_sex=list(M=int_sex_m, F=int_sex_f), 
		intensity_aut=int_aut, 
		sex=s, 
		control_medians_sex=list(M=control_medians_sex_m, F=control_medians_sex_f),
		control_medians_aut=control_medians_aut
	))
}


quantile.normalize.from.reference <- function(x, ref)
{
	quan <- sort(ref[,1])
	ord <- order(x)
	x[ord] <- quan
	return(x)
}


calculate.cnv <- function(bname, samplename=basename(bname), controls, chip=NULL, featureset=NULL, trim=0.1, min.width = 5, nperm = 10000, alpha = 0.001,
                          undo.splits = "sdundo", undo.SD = 2, verbose = TRUE, smoothing=TRUE)
{
	require(DNAcopy)
	msg("Reading idat file for", bname, verbose=verbose)
        rg <- read.rg(bname, verbose=verbose)

        if (is.null(chip))
            chip <- guess.chip(rg)
        if (is.null(featureset))
            featureset <- controls$featureset
        
        probes <- meffil.probe.info(chip, featureset)
        case <- mu.to.cn(rg.to.mu(rg, probes))

	msg("Predicting sex", verbose=verbose)
	sex <- cnv.predict.sex(case, featureset)

	msg("Normalisting against controls", verbose=verbose)
	sites.sex <- c(meffil.get.x.sites(featureset),
                        meffil.get.y.sites(featureset))
	sites.aut <- meffil.get.autosomal.sites(featureset)

	int_sex <- quantile.normalize.from.reference(case[sites.sex], controls$intensity_sex[[sex]])
	int_aut <- quantile.normalize.from.reference(case[sites.aut], controls$intensity_aut)

	msg("Estimating CNVs", verbose=verbose)
	int_sex <- log2(int_sex / controls$control_medians_sex[[sex]])
	int_aut <- log2(int_aut / controls$control_medians_aut)

        features <- meffil.get.features(featureset)
        sites <- features[which(!is.na(features$chromosome)
                                & features$target == "methylation"
                                & !features$snp.exclude),]

	int_sex <- int_sex[names(int_sex) %in% sites$name]
	int_aut <- int_aut[names(int_aut) %in% sites$name]
	sites$chromosome <- ordered(sites$chromosome,
                                    levels = c(paste("chr", 1:22, sep = ""), "chrX", "chrY"))
	p_sex <- sites[match(names(int_sex), sites$name), c("name", "chromosome","position")]
	p_aut <- sites[match(names(int_aut), sites$name), c("name", "chromosome","position")]


	cna_sex <- CNA(int_sex, chrom=p_sex$chromosome, maploc=p_sex$position, data.type="logratio", sampleid=samplename)
	cna_aut <- CNA(int_aut, chrom=p_aut$chromosome, maploc=p_aut$position, data.type="logratio", sampleid=samplename)

	if(smoothing)
	{
		cna_sex <- smooth.CNA(cna_sex, trim=trim)
		cna_aut <- smooth.CNA(cna_aut, trim=trim)
	}

	segment_cna_sex <- segment(cna_sex, min.width=min.width, verbose = verbose, nperm = nperm, alpha = alpha, undo.splits = undo.splits, undo.SD = undo.SD, trim = trim)
	segment_cna_aut <- segment(cna_aut, min.width=min.width, verbose = verbose, nperm = nperm, alpha = alpha, undo.splits = undo.splits, undo.SD = undo.SD, trim = trim)
	out <- rbind(segment_cna_aut$output, segment_cna_sex$output)
	out$chrom <- ordered(out$chrom, levels = c(paste("chr", 1:22, sep = ""), "chrX", "chrY"))
	out$ID <- samplename
	return(out)
}

#' Calculate CNVs from IDAT files
#'
#' Based on the algorithm developed in \code{R/CopyNumber450k} bioconductor package
#' 
#' @param samplesheet Output from \code{meffil.create.samplesheet}
#' @param chip Name returned by \code{\link{meffil.list.chips()}} (Default: \code{NULL}).
#' @param featureset Name returned by \code{\link{meffil.list.featuresets()}} (Default: \code{chip}).
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
meffil.calculate.cnv <- function(samplesheet, chip=NULL, featureset=chip, verbose=FALSE, ...)
{
	controls <- meffil.cnv.controls(featureset, verbose=verbose)

	l1 <- get.index.list(nrow(samplesheet), options("mc.cores")[[1]])

	l <- lapply(l1, function(x)
	{		
		l2 <- mclapply(x, function(i)
		{
			calculate.cnv(samplesheet$Basename[i], samplesheet$Sample_Name[i], controls, chip, featureset, verbose=verbose, ...)
		})
		names(l2) <- samplesheet$Sample_Name[x]
		return(l2) 
	})
	nom <- unlist(sapply(l, names))
	l <- unlist(l, recursive=FALSE)
	names(l) <- nom
	return(l) # list of data frames (output of calculate.cnv), one per sample
}

get.index.list <- function(n, mc.cores)
{
	mc.cores <- ifelse(is.null(mc.cores) || mc.cores < 1, 1, min(mc.cores, n))
	div <- floor(n / mc.cores)
	rem <- n %% mc.cores
	l1 <- lapply(1:div, function(x) (x-1) * mc.cores + 1:mc.cores)
	if(rem != 0) l1[[div+1]] <- l1[[div]][mc.cores] + 1:rem
	return(l1)
}


#' Create matrix of CNV values
#'
#' @param cnv Output from \code{\link{meffil.calculate.cnv}()}.
#' @param featureset Name from \code{\link{meffil.list.featuresets}()}.
#' @return Matrix of ncpg x nsample
#' @export
meffil.cnv.matrix <- function(cnv, featureset="450k") {
    features <- meffil.get.features(featureset)
    sites <- features[which(!is.na(features$chromosome)
                            & features$target == "methylation"
                            & !features$snp.exclude),]

    sites$chromosome <- as.character(sites$chromosome)
    sites <- sites[with(sites, order(chromosome, position, decreasing=F)),]
    sites$id <- with(sites, paste(chromosome, position))

    cnv.matrix <- sapply(cnv, function(segments) {
        segments$chrom <- as.character(segments$chrom)
        segments$loc.start <- as.integer(segments$loc.start)
        segments <- segments[with(segments, order(chrom, loc.start, decreasing=F)),]
        segments$id <- with(segments, paste(chrom, loc.start))

        for (missing.idx in which(!(segments$id %in% sites$id))) {
            chrom.idx <- which(sites$chromosome == segments$chrom[missing.idx])
            idx <- findInterval(segments$loc.start[missing.idx], sites$position[chrom.idx])
            if (idx == 0) idx <- 1
            segments$id[missing.idx] <- sites$id[chrom.idx[idx]]                
        }
        
        first <- match(segments$id, sites$id)
        num.sites <- c(first[-1], nrow(sites)+1) - first
        rep(as.numeric(segments$seg.mean), num.sites)
    })
    rownames(cnv.matrix) <- sites$name
    cnv.matrix
}
           
