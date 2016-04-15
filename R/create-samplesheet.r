#' Create sample sheet if an Illumina one isn't available
#'
#' If necessary generates two columns necessary for some functions: \code{Sample_Name} and \code{Sex}
#' @param basenames Output from \code{\link{meffil.basenames}}
#' @param  delim Optional delim character to separate \code{Sample_Name} into multiple columns. Default: "_"
#' @return Sample sheet data frame
#'
#' @export
meffil.create.samplesheet <- function(path, basenames=meffil.basenames(path), delim = "_") {
    if (length(basenames) == 0) {
        warning("No idat files found.")
        return(NULL)
    }
        
    dat <- data.frame(do.call(rbind, strsplit(basename(basenames), split=delim)), stringsAsFactors=FALSE)

    if (ncol(dat) < 2)
        warning("The basenames in ", path, " do not appear to correspond to idat files")
    
    sentrixpos <- grep("^R[0-9][0-9]C[0-9][0-9]$", as.character(unlist(dat[1,])))
    if(length(sentrixpos)==1) {
        temp <- do.call(rbind, strsplit(as.character(dat[,sentrixpos]), split="C"))
        dat$sentrix_row <- gsub("R", "", temp[,1])
        dat$sentrix_col <- temp[,2]
        dat <- dat[,-sentrixpos]
    }
    
    slidecol <- grep("^[0-9]{9}[0-9]*$", as.character(unlist(dat[1,])))
    if (length(slidecol) == 1) {
        colnames(dat)[slidecol] <- "Slide"
    }
    
    idcol <- which(apply(dat, 2, function(x) all(!duplicated(x))))
    if(length(idcol) >= 1) {
        Sample_Name <- dat[,idcol[1]]
        dat <- dat[,-idcol[1],drop=F]
    }
    else {
        Sample_Name <- make.samplename.from.basename(basenames)
    }
    
    samplesheet <- data.frame(Sample_Name = Sample_Name, Sex = NA, dat, Basename = basenames, stringsAsFactors=FALSE)
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


