library(ggplot2)
library(plyr)
library(reshape2)

# ## pre normalisation data

# sample missmatches concordance with genotype data
# - select snps with 99% concordance and then samples with 80% concordance

# gender
# - check predicted sex against manifest files

#' Plot predicted gender
#'
#' @param samplesheet From \code{read.450k.sheet}
#' @param  norm.objects From \code{meffil.normalize.objects}
#' @export
#' @return Data frame of results plus plot
#' @examples \dontrun{
#'
#'}
meffil.plot.gender <- function(samplesheet, norm.objects)
{
    dat <- data.frame(
        Sample_Name = samplesheet$Sample_Name,
        xy.diff = sapply(norm.objects, function(x) x$xy.diff),
        predicted.sex = sapply(norm.objects, function(x) x$predicted.sex),
        declared.sex = samplesheet$sex
    )
    dat$agree <- dat$declared.sex == dat$predicted.sex
    p1 <- ggplot(dat, aes(y=1, x=xy.diff)) +
        geom_jitter(aes(shape=predicted.sex, colour=agree), size=3) +
        scale_colour_manual(values=c("black", "red")) +
        labs(shape="Predicted sex", x="XY diff", y="", colour="Correct\nprediction") +
        theme_bw() +
        theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.line.y=element_blank())
    return(list(graph=p1, tab=dat))
}


#' Plot average methylated vs unmethylated levels for each individuals
#'
#' plot raw control probes and fit linear regression, remove samples that have sd(y - yhat) > mean*3
#'
#' @param samplesheet From \code{read.450k.sheet}
#' @param  norm.objects From \code{meffil.normalize.objects}
#' @param  outlier.sd Cut off for declaring outliers. Default = 3
#' @param  colour.code Array of length n samples to colour code points. Defaults to NULL
#' @export
#' @return Data frame of results plus plot
#' @examples \dontrun{
#'
#'}
meffil.plot.meth.unmeth <- function(samplesheet, norm.objects, outlier.sd=3, colour.code = NULL)
{
    dat <- data.frame(
        Sample_Name = samplesheet$Sample_Name,
        methylated = sapply(norm.objects, function(x) x$median.m.signal),
        unmethylated = sapply(norm.objects, function(x) x$median.u.signal)
    )
    if(!is.null(colour.code))
    {
        dat$colour.code <- samplesheet[[colour.code]]
        g <- "legend"
    } else {
        dat$colour.code <- 1
        g <- FALSE
    }
    dat$resids <- residuals(lm(methylated ~ unmethylated, dat))
    dat$outliers <- dat$resids > mean(dat$resids)+outlier.sd*sd(dat$resids) | dat$resids < mean(dat$resids)-outlier.sd*sd(dat$resids)
    p1 <- ggplot(dat, aes(y=methylated, x=unmethylated)) +
        geom_point(aes(colour=colour.code)) +
        guides(colour = g) +
        geom_point(data=subset(dat, outliers), shape=1, size=3.5) +
        labs(y = "Median methylated signal", x = "Median unmethylated signal", colour = colour.code) +
        stat_smooth(method="lm", se=FALSE)
    return(list(graph=p1, tab=dat))
}

# meffil.plot.meth.unmeth(samplesheet, norm.objects, outlier.sd=2, colour.code = "sex")


# calculate means for each sample from control probes of each
# - dyebias
# - oob ratio
# - bisulphide conversion
# - spec2red2
# - hybe3
# - remaining 37

# plot all of them - sample against mean ,  colour by plate

#' Plot the means of control probes for each sample and for each control probe type
#'
#' @param samplesheet From \code{read.450k.sheet}
#' @param  norm.objects From \code{meffil.normalize.objects}
#' @param  control.categories Which control probe categories to plot. Defaults to all available
#' @param  colour.code Array of length n samples to colour code points. Defaults to NULL
#' @param  outlier.sd Cut off for declaring outliers. Default = 5
#' @export
#' @return Data frame of results plus plot
#' @examples \dontrun{
#'
#'}
meffil.plot.controlmeans <- function(samplesheet, norm.objects, control.categories=names(norm.objects[[1]]$controls), colour.code = NULL, outlier.sd=5)
{
    dat <- data.frame(Sample_Name = samplesheet$Sample_Name)
    if(!is.null(colour.code))
    {
        dat$colour.code <- samplesheet[[colour.code]]
        g <- "legend"
    } else {
        dat$colour.code <- 1
        g <- FALSE
    }

    d <- data.frame(t(sapply(norm.objects, function(x) x$controls)))
    names(d) <- names(norm.objects[[1]]$controls)
    dat <- data.frame(dat, subset(d, select=control.categories))
    names(dat) <- c("Sample_Name", "colour.code", control.categories)
    dat <- dat[order(dat$colour.code), ]
    dat$colour.code <- as.character(dat$colour.code)
    dat$id <- 1:nrow(dat)
    dat <- reshape2::melt(dat, id.vars=c("Sample_Name", "colour.code", "id"))
    dat <- ddply(dat, .(variable), function(x)
    {
        x <- mutate(x)
        x$outliers <- x$value > mean(x$value) + outlier.sd * sd(x$value) | x$value < mean(x$value) - outlier.sd * sd(x$value)
        return(x)
    })
    p1 <- ggplot(dat, aes(x=id, y=value)) +
        geom_point(aes(colour=colour.code)) +
        geom_point(data=subset(dat, outliers), shape=1, size=3.5) +
        guides(colour=g) +
        facet_wrap(~ variable, scales="free_y") +
        labs(y="Mean signal", x="ID", colour=colour.code)
    return(list(graph=p1, tab=dat))
}

# meffil.plot.controlmeans(samplesheet, norm.objects, colour.code="sex", outlier.sd=5)




#' Plot detection p values from idat files
#'
#' @param samplesheet From \code{read.450k.sheet}
#' @param  norm.objects From \code{meffil.normalize.objects}
#' @export
#' @return Data frame of results plus plot
#' @examples \dontrun{
#'
#'}
meffil.plot.detectionp <- function(samplesheet, norm.objects, threshold = 0.05, colour.code=NULL)
{
    nprobe <- length(unique(meffil.probe.info()$name))
    dat <- data.frame(
        Sample_Name = samplesheet$Sample_Name,
        prop.badprobes = sapply(norm.objects, function(x) length(x$bad.probes.detectionp) / nprobe)
    )
    if(!is.null(colour.code))
    {
        dat$colour.code <- samplesheet[[colour.code]]
        g <- "legend"
    } else {
        dat$colour.code <- 1
        g <- FALSE
    }
    dat <- dat[order(dat$colour.code), ]
    dat$id <- 1:nrow(dat)
    p1 <- ggplot(dat, aes(y=prop.badprobes, x=id)) +
        geom_point(aes(colour=colour.code)) +
        guides(colour = g) +
        geom_hline(yintercept=threshold) +
        labs(y = paste("Proportion CpG sites with p >", norm.objects$bad.probes.detectionp.threshold), x = "Sample ID", colour = colour.code)
    return(list(graph=p1, tab=dat))
}

# meffil.plot.detectionp(samplesheet, norm.objects)

#' Manhattan plot of detection pval per probe - percentage with pvalue < 0.01
#'
#' @param samplesheet From \code{read.450k.sheet}
#' @param  norm.objects From \code{meffil.normalize.objects}
#' @export
#' @return Data frame of results plus plot
#' @examples \dontrun{
#'
#'}
meffil.plot.detectionp.manhattan <- function(samplesheet, norm.objects)
{

}

#' Plot number of beads per sample
#'
#' @param samplesheet From \code{read.450k.sheet}
#' @param  norm.objects From \code{meffil.normalize.objects}
#' @export
#' @return Data frame of results plus plot
#' @examples \dontrun{
#'
#'}
meffil.plot.beadnum <- function(samplesheet, norm.objects)
{
    
}

#' Manhattan plot of number of beads by probe - percentage of probes with beads < 3 for each sample
#'
#' @param samplesheet From \code{read.450k.sheet}
#' @param  norm.objects From \code{meffil.normalize.objects}
#' @export
#' @return Data frame of results plus plot
#' @examples \dontrun{
#'
#'}
meffil.plot.beadnum.manhattan <- function(samplesheet, norm.objects)
{
    
}




