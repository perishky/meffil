# library(CopyNumber450k)
# library(CopyNumber450kData)
# library(plyr)
# library(meffil)
# library(splitstackshape)

# library(matrixStats)



#' Calculate CNVs from IDAT files
#'
#' Based on the algorithm developed in \code{R/CopyNumber450k} bioconductor package
#' 
#' @param samplesheet Output from \code{meffil.create.samplesheet}
#' @param controls Output from \code{meffil.cnv.controls}
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
meffil.calculate.cnv <- function(samplesheet, controls, verbose=FALSE, ...)
{

	l1 <- get.index.list(nrow(samplesheet), options("mc.cores")[[1]])

	l <- lapply(l1, function(x)
	{		
		l2 <- mclapply(x, function(i)
		{
			calculate.cnv(samplesheet$Basename[i], samplesheet$Sample_Name[i], controls, verbose=verbose, ...)
		})
		names(l2) <- samplesheet$Sample_Name[x]
		return(l2)
	})
	nom <- unlist(sapply(l, names))
	l <- unlist(l, recursive=FALSE)
	names(l) <- nom
	return(l)
}

#' Create controls for calculating CNVs
#'
#' @param rgset An rgset of samples to be used as controls
#' @param verbose Default = FALSE
#' @export
#' @return List with intensities, sex, pre-calculated median values for use in CNV calling
meffil.cnv.controls <- function(rgset, verbose=FALSE)
{
	msg("Normalising", verbose=verbose)
	x <- preprocessIllumina(rgset, bg.correct=TRUE, normalize=NULL)
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
	# control_medians_sex_m <- apply(int_sex_m, 1, median, na.rm=T)
	# control_medians_sex_f <- apply(int_sex_f, 1, median, na.rm=T)
	# control_medians_aut <- apply(int_aut, 1, median, na.rm=T)

	control_medians_sex_m <- rowMedians(int_sex_m, na.rm=T)
	control_medians_sex_f <- rowMedians(int_sex_f, na.rm=T)
	control_medians_aut <- rowMedians(int_aut, na.rm=T)


	return(list(
		intensity_sex=list(M=int_sex_m, F=int_sex_f), 
		intensity_aut=int_aut, 
		sex=s, 
		control_medians_sex=list(M=control_medians_sex_m, F=control_medians_sex_f),
		control_medians_aut=control_medians_aut
	))
}


#' Create matrix of CNV values
#'
#' @param cnv Output from \code{meffil.calculate.cnv}
#' @param abs.thresh Seg.mean threshold for CNV to be called
#' @param p.thresh Bonferroni adjusted p-value threshold for CNV to be called
#' @export
#' @return Matrix of ncpg x nsample
meffil.cnv.matrix <- function(cnv, abs.thresh=0.25, p.thresh=0.05)
{
	probeinfo <- cnv.probeinfo()

	l1 <- get.index.list(length(cnv), options("mc.cores")[[1]])
	l <- lapply(l1, function(ii)
	{
		res <- mclapply(ii, function(X)
		{
			x <- cnv[[X]]
			x$seg.mean[x$pval_adj > p.thresh | abs(x$seg.mean) < abs.thresh] <- 0
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
		names(res) <- names(cnv)[ii]
		return(res)
	})
	nom <- unlist(lapply(l, names))
	l <- do.call(cbind, unlist(l, recursive=FALSE))
	colnames(l) <- nom
	return(l)
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
	int_sex_prop <- log2(int_sex / controls$control_medians_sex[[sex]])
	int_aut_prop <- log2(int_aut / controls$control_medians_aut)

	probes <- meffil.probe.info()
	int_sex_prop <- int_sex_prop[names(int_sex_prop) %in% probes$name[!probes$snp.exclude]]
	int_aut_prop <- int_aut_prop[names(int_aut_prop) %in% probes$name[!probes$snp.exclude]]
	probes$chr <- ordered(probes$chr, levels = c(paste("chr", 1:22, sep = ""), "chrX", "chrY"))
	p_sex <- probes[match(names(int_sex_prop), probes$name), c("name", "chr","pos")]
	p_aut <- probes[match(names(int_aut_prop), probes$name), c("name", "chr","pos")]


	cna_sex <- CNA(int_sex_prop, chrom=p_sex$chr, maploc=p_sex$pos, data.type="logratio", sampleid=samplename)
	cna_aut <- CNA(int_aut_prop, chrom=p_aut$chr, maploc=p_aut$pos, data.type="logratio", sampleid=samplename)

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

	msg("Calculating p-values", verbose=verbose)
	out <- cnv.pvalue(out, int_aut, int_sex, sex, controls)
	return(out)
}


cnv.pvalue <- function(out, int_aut, int_sex, sex, controls)
{
	# Get the sum of the intensities for the probes in each of the controls
	# Get the mean sum of intensities
	# Get the SD of the sum of intensities
	# Get the sum of intensities of probes in the target
	probeinfo <- cnv.probeinfo()
	for(i in 1:nrow(out))
	{
		probes <- subset(probeinfo, chr == out$chrom[i] & pos <= out$loc.end[i] & pos >= out$loc.start[i])$name
		if(out$chrom[i] %in% c("chrX", "chrY"))
		{
			control_int_sum <- colSums(controls$intensity_sex[[sex]][probes, , drop=FALSE])
			sample_sum <- sum(int_sex[probes])
		} else {
			control_int_sum <- colSums(controls$intensity_aut[probes, , drop=FALSE])
			sample_sum <- sum(int_aut[probes])
		}
		control_mean <- mean(control_int_sum)
		control_sd <- sd(control_int_sum)

		z_score <- (sample_sum - control_mean)/control_sd
		pval <- 2 * pnorm(-abs(z_score))

		# Compute p-value (2 sided t-test)
		out$nprobe[i] <- length(probes)
		out$sample_sum[i] <- sample_sum
		out$control_mean[i] <- control_mean
		out$control_sd[i] <- control_sd
		out$z_score[i] <- z_score
		out$pval[i] <- pval
	}
	out$pval_adj <- p.adjust(out$pval, method="bonferroni")
	return(out)
}


cnv.probeinfo <- function()
{
	probenames <- meffil.get.sites()
	probeinfo <- subset(meffil.probe.info(), name %in% probenames & target == "M" & ! snp.exclude, select=c(name, chr, pos))
	probeinfo$chr <- ordered(probeinfo$chr, levels=paste("chr", c(1:22, "X", "Y"), sep=""))
	probeinfo <- probeinfo[order(probeinfo$chr, probeinfo$pos), ]
	probeinfo$chr <- as.character(probeinfo$chr)
	rownames(probeinfo) <- NULL
	return(probeinfo)
}

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
  

get.index.list <- function(n, mc.cores)
{
	mc.cores <- ifelse(is.null(mc.cores) || mc.cores < 1, 1, mc.cores)
	div <- floor(n / mc.cores)
	rem <- n %% mc.cores
	l1 <- lapply(1:div, function(x) (x-1) * mc.cores + 1:mc.cores)
	if(rem != 0) l1[[div+1]] <- l1[[div]][mc.cores] + 1:rem
	return(l1)
}



quantile.normalize.from.reference <- function(x, ref)
{
	quan <- sort(ref[,1])
	ord <- order(x)
	x[ord] <- quan
	return(x)
}
