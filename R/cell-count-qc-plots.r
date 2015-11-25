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
    stopifnot(is.list(count.objects) && length(count.objects) > 0)
    stopifnot(sapply(count.objects, is.cell.count.object))
    stopifnot(length(unique(sapply(count.objects, function(object) object$reference))) == 1)
    
    if (is.null(names(count.objects)))
        names(count.objects) <- paste0("s", 1:length(count.objects))
    
    reference <- count.objects[[1]]$reference
    reference.object <- get.cell.type.reference(reference)
    
    beta <- cbind(reference.object$beta, sapply(count.objects, function(object) object$beta))

    if (ncol(beta) <= 200) {
        dat <- sapply(colnames(beta), function(sample)
                      data.frame(sample=sample, beta=beta[,sample]), simplify=F)
        dat <- dat[order(sapply(dat, function(x) mean(x$beta)))]
        dat <- do.call(rbind, dat)
        dat$is.reference <- dat$sample %in% colnames(reference.object$beta)        
        beta.plot <- (ggplot(dat, aes(x=sample, y=beta, fill=is.reference)) +
                      geom_boxplot() +
                      guides(fill=FALSE) +
                      stat_summary(fun.y=mean, geom="point", shape=5, size=3) +
                      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,
                                                       face="bold", size=rel(2))) +
                      ggtitle("Distribution of beta values per sample"))
    }
    else {
        dat <- data.frame(sample=colnames(beta), beta=colMeans(beta), na.rm=T)
        dat <- dat[order(dat$beta, decreasing=F),]
        dat$x <- 1:nrow(dat)
        dat$is.reference <- dat$sample %in% colnames(reference.object$beta)
        beta.plot <- (ggplot(dat, aes(x=x, y=beta)) +
                      geom_point(size=3, colour="skyblue3") +
                      geom_text(data=dat[dat$is.reference,],
                                aes(x=x, y=beta, label=sample),
                                vjust=0.5, hjust=0.5,
                                angle=0, fontface="bold") +
                      labs(y="Mean methylation level",
                           x="Samples (sorted by mean methylation level)") +
                      ggtitle("Mean beta values per sample") +
                      theme(axis.text.x = element_blank()))
    }

    counts <- t(sapply(count.objects, function(object) object$counts))
    dat <- sapply(colnames(counts), function(cell.type)
                  data.frame(cell.type=cell.type, count=counts[,cell.type]), simplify=F)
    dat <- dat[order(sapply(dat, function(x) mean(x$count)))]
    dat <- do.call(rbind, dat)
    
    count.plot <- (ggplot(dat, aes(x=cell.type, y=count)) +
                   geom_boxplot() +
                   guides(fill=FALSE) +
                   stat_summary(fun.y=mean, geom="point", shape=5, size=3) +
                   theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1,
                                                    face="bold", size=rel(2))) +
                   ggtitle("Distribution of cell count estimates per cell type"))
    
    list(betas=beta.plot, counts=count.plot, reference=reference)
}

