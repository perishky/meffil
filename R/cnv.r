# library(CopyNumber450k)
# library(CopyNumber450kData)
# library(plyr)
# library(meffil)
# library(splitstackshape)



read.mu <- function(basename, verbose=F) {
  rg <- read.rg(basename, verbose=verbose)
  rg.to.mu(rg)
}

mu.to.cn <- function(mu) {
  mu$M + mu$U
}

read.cn <- function(basename, verbose=F) {
  mu <- read.mu(basename, verbose=verbose)
  mu.to.cn(mu)
}

cnv.predict.sex <- function(cn, sex.cutoff=-2) {
    probes <- meffil.probe.info()
  probes.x <- meffil.get.x.sites()
  x.signal <- median(cn[probes.x], na.rm=T)
  probes.y <- meffil.get.y.sites()
  y.signal <- median(cn[probes.y], na.rm=T)
   xy.diff <- y.signal-x.signal
   ifelse(xy.diff < sex.cutoff, "F", "M")
}
  

#' Create controls for calculating CNVs
#' @param verbose Default = FALSE>
#' @export
#' @return List with pre-calculated CNV data
meffil.cnv.controls <- function(verbose=FALSE)
{
	require(CopyNumber450kData)
	data(RGcontrolSetEx)
	x <- preprocessRaw(RGcontrolSetEx)
	x <- getMeth(x) + getUnmeth(x)

	msg("Predicting sex", verbose=verbose)
	s <- apply(log2(x), 2, cnv.predict.sex)

	# Quantile normalize controls

	msg("Performing quantile normalisation", verbose=verbose)
	probes.sex <- c(meffil.get.x.sites(), meffil.get.y.sites())
	probes.aut <- meffil.get.autosomal.sites()

	int_sex_m <- x[probes.sex,s=="M"]
	int_sex_f <- x[probes.sex,s=="F"]
	int_aut <- x[probes.aut,]

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


calculate.cnv <- function(bname, samplename=basename(bname), controls, trim=0.1, min.width = 5, nperm = 10000, alpha = 0.001, undo.splits = "sdundo", undo.SD = 2, verbose = TRUE, smoothing=TRUE)
{
	require(DNAcopy)
	msg("Reading idat file for", bname, verbose=verbose)
	case <- read.cn(bname)
	msg("Predicting sex", verbose=verbose)
	sex <- cnv.predict.sex(log2(case))

	msg("Normalisting against controls", verbose=verbose)
	probes.sex <- c(meffil.get.x.sites(), meffil.get.y.sites())
	probes.aut <- meffil.get.autosomal.sites()

	int_sex <- quantile.normalize.from.reference(case[probes.sex], controls$intensity_sex[[sex]])
	int_aut <- quantile.normalize.from.reference(case[probes.aut], controls$intensity_aut)

	msg("Estimating CNVs", verbose=verbose)
	int_sex <- log2(int_sex / controls$control_medians_sex[[sex]])
	int_aut <- log2(int_aut / controls$control_medians_aut)

	probes <- meffil.probe.info()
	int_sex <- int_sex[names(int_sex) %in% probes$name[!probes$snp_exclude]]
	int_aut <- int_aut[names(int_aut) %in% probes$name[!probes$snp_exclude]]
	probes$chr <- ordered(probes$chr, levels = c(paste("chr", 1:22, sep = ""), "chrX", "chrY"))
	p_sex <- probes[match(names(int_sex), probes$name), c("name", "chr","pos")]
	p_aut <- probes[match(names(int_aut), probes$name), c("name", "chr","pos")]


	cna_sex <- CNA(int_sex, chrom=p_sex$chr, maploc=p_sex$pos, data.type="logratio", sampleid=samplename)
	cna_aut <- CNA(int_aut, chrom=p_aut$chr, maploc=p_aut$pos, data.type="logratio", sampleid=samplename)

	if(smoothing)
	{
		cna_sex <- smooth.CNA(cna_sex, trim=trim)
		cna_aut <- smooth.CNA(cna_aut, trim=trim)
	}

	segment_cna_sex <- segment(cna_sex, min.width=min.width, verbose = verbose, nperm = nperm, alpha = alpha, undo.splits = undo.splits, undo.SD = undo.SD, trim = trim)
	segment_cna_aut <- segment(cna_aut, min.width=min.width, verbose = verbose, nperm = nperm, alpha = alpha, undo.splits = undo.splits, undo.SD = undo.SD, trim = trim)
	out <- rbind(segment_cna_aut$output, segment_cna_sex$output)
	out$chrom <- ordered(out$chrom, levels = c(paste("chr", 1:22, sep = ""), "chrX", "chrY"))
	return(out)
}

#' Calculate CNVs from IDAT files
#'
#' Based on the algorithm developed in \code{R/CopyNumber450k} bioconductor package
#' 
#' @param samplesheet Output from \code{meffil.create.samplesheet}
#' @param verbose Default = FALSE
#' @param ... Extra parameters to be passed to \code{CNAcopy} for segmentation. See details.
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
meffil.calculate.cnv <- function(samplesheet, verbose=FALSE, ...)
{
	controls <- meffil.cnv.controls(TRUE)
	l <- mclapply(1:nrow(samplesheet), function(i)
	{
		calculate.cnv(samplesheet$Basename[i], samplesheet$Sample_Name[i], controls, verbose=verbose, ...)
	})
	names(l) <- samplesheet$Sample_Name
	return(l)
}

#' Create matrix of CNV values
#'
#' @param cnv Output from \code{meffil.calculate.cnv}
#' @export
#' @return Matrix of ncpg x nsample
meffil.cnv.matrix <- function(cnv)
{
	probenames <- meffil.get.sites()
	probeinfo <- subset(meffil.probe.info(), name %in% probenames & target == "M" & ! snp_exclude, select=c(name, chr, pos))
	probeinfo$chr <- ordered(probeinfo$chr, levels=paste("chr", c(1:22, "X", "Y"), sep=""))
	probeinfo <- probeinfo[order(probeinfo$chr, probeinfo$pos), ]
	probeinfo$chr <- as.character(probeinfo$chr)

	res <- mcsapply.safe(seq_along(cnv), function(X)
	{
		x <- cnv[[X]]
		message(X, " of ", length(cnv))
		res <- lapply(1:nrow(x), function(i)
		{
			index <- with(probeinfo, chr == x$chrom[i] & pos >= x$loc.start[i] & pos <= x$loc.end[i])
			p <- rep(x$seg.mean[i], sum(index))
			names(p) <- probeinfo$name[index]
			return(p)
		})
		res <- unlist(res)
		return(res)
	})
	colnames(res) <- names(cnv)
	return(res)
}
