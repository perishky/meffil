#' Create sample sheet if an Illumina one isn't available
#'
#' If necessary generates two columns necessary for some functions: \code{Sample_Name} and \code{Sex}
#' @param basenames Output from \code{\link{meffil.basenames}}
#' @param  Sample_Name Array of unique sample IDs
#' @param  Sex Array of values denoting sex for each sample, must be "M", "F" or NA
#' @param  delim Optional delim character to separate \code{Sample_Name} into multiple columns. Default: "_"
#' @export
#' @return Sample sheet data frame
meffil.create.samplesheet <- function(path, Sample_Name = NULL, Sex = NULL, delim = "_")
{
	basenames <- meffil.basenames(path)
	if(is.null(Sex)) Sex <- rep(NA, length(basenames))
	if(any(!Sex %in% c("M", "F", NA))) stop("Sex column must only contain 'M', 'F' or NA values")
	if(is.null(Sample_Name)) Sample_Name <- make.samplename.from.basename(basenames)
	stopifnot(length(Sample_Name) == length(basenames))
	stopifnot(length(Sex) == length(basenames))

	dat <- data.frame(do.call(rbind, strsplit(basename(basenames), split=delim)), stringsAsFactors=FALSE)
	idcol <- which(apply(dat, 2, function(x) all(!duplicated(x))))
	if(length(idcol) == 1)
	{
		Sample_Name <- dat[,idcol]
		dat <- dat[,-idcol]
	}
	sentrixpos <- grep("^R[0-9][0-9]C[0-9][0-9]$", as.character(unlist(dat[1,])))
	if(length(sentrixpos)==1)
	{
		temp <- do.call(rbind, strsplit(as.character(dat[,sentrixpos]), split="C"))
		dat$sentrix_row <- gsub("R", "", temp[,1])
		dat$sentrix_col <- temp[,2]
		dat <- dat[,-sentrixpos]
	}
	samplesheet <- data.frame(Sample_Name = Sample_Name, Sex = Sex, dat, Basename = basenames, stringsAsFactors=FALSE)
	return(samplesheet)
}


make.samplename.from.basename <- function(basenames)
{
	Sample_Name <- basename(basenames)
	if(any(duplicated(Sample_Name)))
	{
		warning("Some duplicated Sample_Name entries")
		Sample_Name <- make.unique(Sample_Name)
	}
	return(Sample_Name)
}


