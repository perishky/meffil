#' @title Function to read Illumina "Sample Sheet" taken from R/minfi
#'
#' @description
#' Reading an Illumina methylation sample sheet, containing pheno-data 
#' information for the samples in an experiment.
#'
#' @param base The base directory from which the search is started.
#' @param  pattern = "csv$" What pattern is used to identify a sample sheet file, see \code{list.files}
#' @param  ignore.case = TRUE Should the file search be case sensitive?
#' @param  recursive = TRUE Should the file search be recursive, see \code{list.files}?
#' @param  verbose = TRUE Should the function be verbose?
#'
#' @details
#' This function search the directory \code{base} (possibly including
#' subdirectories depending on the argument \code{recursive} for
#' \dQuote{sample sheet} files (see below).  These files are identified
#' solely on the base of their filename given by the arguments
#' \code{pattern} and \code{ignore.case} (note the use of a dollarsign to
#' mean end of file name).#' 
#'
#' In case multiple sheet files are found, they are all read and the
#' return object will contain the concatenation of the files.
#' 
#' A sample sheet file is essentially a CSV (comma-separated) file
#' containing one line per sample, with a number of columns describing
#' pheno-data or other important information about the sample.  The file
#' may contain a header, in which case it is assumed that all lines up to
#' and including a line starting with \code{\[Data\]} should be dropped.
#' This is modelled after a sample sheet file Illumina provides.  It is
#' also very similar to the \code{targets} file made used by the popular
#' \code{limma} package (see the extensive package vignette).#' 
#'
#' An attempt at guessing the file path to the IDAT files represented in
#' the sheet is made.  This should be doublechecked and might need to
#' manually changed.
#'
#' @export
#'
#' @return
#' A \code{data.frame} containing the columns of all the sample sheets. 
#' As described in details, a column named \code{Sentrix_Position} is renamed 
#' to \code{Array} and \code{Sentrix_ID} is renamed to \code{Slide}.  In addition 
#' the \code{data.frame} will contain a column named \code{Basename}.
#'
read.450k.sheet <- function(base, pattern = "csv$", ignore.case = TRUE, recursive = TRUE, verbose = TRUE) {
	readSheet <- function(file) {
		dataheader <- grep("^\\[DATA\\]", readLines(file), ignore.case = TRUE)
		if(length(dataheader) == 0)
			dataheader <- 0
		df <- read.csv(file, stringsAsFactor = FALSE, skip = dataheader)
		if(length(nam <- grep("Sentrix_Position", names(df), ignore.case = TRUE, value = TRUE)) == 1) {
			df$Array <- df[, nam]
			df[, nam] <- NULL
		}
		if(length(nam <- grep("Array[\\._]ID", names(df), ignore.case = TRUE, value = TRUE)) == 1) {
			df$Array <- df[, nam]
			df[, nam] <- NULL
		}
		if(! "Array" %in% names(df))
			warning(sprintf("Could not infer array name for file: %s", file))
		if(length(nam <- grep("Sentrix_ID", names(df), ignore.case = TRUE, value = TRUE)) == 1) {
			df$Slide <- df[, nam]
			df[, nam] <- NULL
		}
		if(length(nam <- grep("Slide[\\._]ID", names(df), ignore.case = TRUE, value = TRUE)) == 1) {
			df$Slide <- df[, nam]
			df[, nam] <- NULL
		}
		if(! "Slide" %in% names(df))
			warning(sprintf("Could not infer slide name for file: %s", file))
		else
			df[, "Slide"] <- as.character(df[, "Slide"])
		if(length(nam <- grep("Plate[\\._]ID", names(df), ignore.case = TRUE, value = TRUE)) == 1) {
			df$Plate <- df[, nam]
			df[, nam] <- NULL
		}
		if(!is.null(df$Array)) {
			patterns <- sprintf("%s_%s_Grn.idat", df$Slide, df$Array)
			allfiles <- list.files(dirname(file), recursive = recursive, full.names = TRUE)
			basenames <- sapply(patterns, function(xx) grep(xx, allfiles, value = TRUE))
			names(basenames) <- NULL
			basenames <- sub("_Grn\\.idat", "", basenames, ignore.case = TRUE)
			df$Basename <- basenames
		}
		df
	}
	if(!all(file.exists(base)))
		stop("'base' does not exists")
	info <- file.info(base)
	if(!all(info$isdir) && !all(!info$isdir))
		stop("'base needs to be either directories or files")
	if(all(info$isdir)) {
		csvfiles <- list.files(base, recursive = recursive, pattern = pattern,
							   ignore.case = ignore.case, full.names = TRUE)
		if(verbose) {
			cat("[read.450k.sheet] Found the following CSV files:\n")
			print(csvfiles)
		}
	} else
		csvfiles <- list.files(base, full.names = TRUE)
	dfs <- lapply(csvfiles, readSheet)
	namesUnion <- Reduce(union, lapply(dfs, names))
	df <- do.call(rbind, lapply(dfs, function(df) {
		newnames <- setdiff(namesUnion, names(df))
		newdf <- matrix(NA, ncol = length(newnames), nrow = nrow(df), dimnames = list(NULL, newnames))
		cbind(df, as.data.frame(newdf))
	}))
	df
}




