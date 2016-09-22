#' Number of control matrix principal components
#'
#' Fits probe intensities to principal components of the microarray control matrix
#' and calculates the resulting mean squared residuals for different
#' numbers of principal components.
#' 
#' @param qc.objects A list of outputs from \code{\link{meffil.create.qc.object}()}.
#' @param number.pcs Number of principal components to include in the design matrix (Default: all).
#' @param fixed.effects Names of columns in samplesheet that should be included as fixed effects
#' along with control matrix principal components (Default: NULL).
#' @param random.effects Names of columns in samplesheet that should be included as random effects
#' (Default: NULL).
#' @return A list containing a data frame with the mean squared residuals for different numbers of principal components
#' and a plot of these residuals.
#'
#' @export
meffil.plot.pc.fit <- function(qc.objects, fixed.effects=NULL, random.effects=NULL, n.cross=10, name="autosomal.ii") {
    stopifnot(is.valid.site.subset(name))    
    stopifnot(all(sapply(qc.objects, is.qc.object)))
    
    if (2*n.cross > length(qc.objects)) 
        n.cross <- floor(length(qc.objects)/2)
    
    n.quantiles <- length(qc.objects[[1]]$quantiles[[1]]$M)
    max.pcs <- min(ncol(meffil.control.matrix(qc.objects)),
                   length(qc.objects) - ceiling(length(qc.objects)/n.cross))
    
    stats <- mclapply(1:max.pcs, function(number.pcs) {
        residuals <- list(M=matrix(NA,nrow=n.quantiles,ncol=length(qc.objects)),
                          U=matrix(NA,nrow=n.quantiles,ncol=length(qc.objects)))
        group <- sample(rep(1:n.cross, length.out=length(qc.objects)), length(qc.objects), replace=F)
        
        for (test.group in unique(group)) {
            msg("pcs", number.pcs, "group", test.group)
            test.idx <- which(group == test.group)
            design.matrix <- predict.design.matrix(qc.objects, number.pcs, test.idx,
                                                   fixed.effects=fixed.effects,
                                                   random.effects=random.effects)

            intensity.R <- sapply(qc.objects, function(object) object$intensity.R)
            intensity.G <- sapply(qc.objects, function(object) object$intensity.G)
            valid.idx <- which(intensity.R + intensity.G > 200)
            if(length(valid.idx) == 0) {
                valid.idx <- 1:length(intensity.R)
                warning("Most or all of the microarrays have very low intensity.")
            }
            reference.idx <- valid.idx[which.min(abs(intensity.R/intensity.G-1)[valid.idx])]
            dye.intensity <- (intensity.R + intensity.G)[reference.idx]/2
            for (target in c("M","U")) {
                original <- sapply(qc.objects[test.idx], function(object) {
                    object$quantiles[[name]][[target]] * dye.intensity/object$dye.intensity
                })
                residuals[[target]][,test.idx] <- (normalize.quantiles(original, design.matrix)
                                                   - rowMeans(original))
            }
        }
        
        c(n=number.pcs,M=mean((residuals$M[,-1])^2), U=mean((residuals$U[,-1])^2))
    })

    stats <- as.data.frame(do.call(rbind, stats))

    list(data=stats,
         plot=(ggplot(stats, aes(x=n)) +
               geom_line(aes(y=M, colour="M")) +
               geom_line(aes(y=U, colour="U")) +
               ggtitle("Fit residuals for different numbers of PCs") +
               labs(x="number of PCs", y="Mean squared residuals") +
               scale_x_continuous(breaks=seq(0,max(stats$n),by=5 )) +
               theme(legend.title=element_blank())))
}
