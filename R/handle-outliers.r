#' Handle outliers in a methylation matrix
#'
#' @param beta Methylation matrix (rows=CpG sites, columns=samples, values=methylation levels).
#' @param winsorize.pct Apply all regression models to methylation levels
#' winsorized to the given level. Set to NA to avoid winsorizing (Default: 0.05).
#' @param outlier.iqr.factor For each CpG site, prior to fitting regression models,
#' set methylation levels less than
#' \code{Q1 - outlier.iqr.factor * IQR} or more than
#' \code{Q3 + outlier.iqr.factor * IQR} to NA. Here IQR is the inter-quartile
#' range of the methylation levels at the CpG site, i.e. Q3-Q1.
#' Set to NA to skip this step (Default: NA).
#' @return \code{beta} after winsorizing and outliers set to NA.
#'
#' @export
meffil.handle.outliers <- function(beta, winsorize.pct=0.05, outlier.iqr.factor=NA) {
    if (is.numeric(winsorize.pct))  
        beta <- meffil:::winsorize(beta, pct=winsorize.pct)
    
    if (is.numeric(outlier.iqr.factor)) {
        if (is.matrix(beta)) {        
            q <- rowQuantiles(beta, probs = c(0.25, 0.75), na.rm = T)
            iqr <- q[,2] - q[,1]
            too.hi <- which(beta > q[,2] + outlier.iqr.factor * iqr, arr.ind=T)
            too.lo <- which(beta < q[,1] - outlier.iqr.factor * iqr, arr.ind=T)
            if (nrow(too.hi) > 0) beta[too.hi] <- NA
            if (nrow(too.lo) > 0) beta[too.lo] <- NA
        }
        else {
            q <- quantile(beta, probs = c(0.25, 0.75), na.rm = T)
            iqr <- q[2] - q[1]
            too.hi <- which(beta > q[2] + outlier.iqr.factor * iqr)
            too.lo <- which(beta < q[1] - outlier.iqr.factor * iqr)
            if (length(too.hi) > 0) beta[too.hi] <- NA
            if (length(too.lo) > 0) beta[too.lo] <- NA
        }
    }

    beta
}
