#' Cell count estimate quality plot
#'
#' @param count.objects A list of objects each obtained from
#' \code{\link{meffil.estimate.cell.counts}()}.
#' @return Two \link{ggplot2} boxplot objects:
#' - \code{betas} Contains one box per sample or reference cell type
#' representing the distribution of methylation levels for the CpG sites used to
#' estimate cell counts.
#' - \code{counts} Contains one box per reference cell type
#' representing the distribution of cell count estimates across the samples.
#' 
#' @export 
meffil.cell.count.qc.plots <- function(count.objects) {
    if (is.null(names(count.objects)))
        names(count.objects) <- paste0("s", 1:length(count.objects))
    
    reference <- count.objects[[1]]$reference
    reference.object <- get.cell.type.reference.object(reference)
    
    beta <- cbind(reference.object$beta, sapply(count.objects, function(object) object$beta))

    dat <- sapply(colnames(beta), function(sample)
                  data.frame(sample=sample, beta=beta[,sample]), simplify=F)
    dat <- dat[order(sapply(dat, function(x) mean(x$beta)))]
    dat <- do.call(rbind, dat)
    dat$col <- dat$sample %in% colnames(reference.object$beta)
    
    beta.plot <- (ggplot(dat, aes(x=sample, y=beta, fill=col)) +
                  geom_boxplot() +
                  guides(fill=FALSE) +
                  stat_summary(fun.y=mean, geom="point", shape=5, size=3) +
                  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
                  ggtitle("Distribution of beta values per sample"))

    counts <- t(sapply(count.objects, function(object) object$counts))
    dat <- sapply(colnames(counts), function(cell.type)
                  data.frame(cell.type=cell.type, count=counts[,cell.type]), simplify=F)
    dat <- dat[order(sapply(dat, function(x) mean(x$count)))]
    dat <- do.call(rbind, dat)
    
    count.plot <- (ggplot(dat, aes(x=cell.type, y=count)) +
                   geom_boxplot() +
                   guides(fill=FALSE) +
                   stat_summary(fun.y=mean, geom="point", shape=5, size=3) +
                   theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
                   ggtitle("Distribution of cell count estimates per cell type"))
    
    list(betas=beta.plot, counts=count.plot, reference=reference)
}

