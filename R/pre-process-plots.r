library(ggplot2)
library(plyr)
library(reshape2)


#' Plot predicted sex
#'
#' @param samplesheet From \code{\link{read.450k.sheet}()}
#' @param  norm.objects From \code{\link{meffil.normalize.objects}()}
#' @export
#' @return Data frame of results plus plot
#' @examples \dontrun{
#'
#'}
meffil.plot.sex <- function(samplesheet, norm.objects, outlier.sd=3)
{
    dat <- data.frame(
        Sample_Name = samplesheet$Sample_Name,
        xy.diff = sapply(norm.objects, function(x) x$xy.diff),
        predicted.sex = sapply(norm.objects, function(x) x$predicted.sex),
        declared.sex = samplesheet$sex
    )
    dat <- ddply(dat, .(predicted.sex), function(x)
    {
        x <- mutate(x)
        x$outliers <- with(x, xy.diff > mean(xy.diff, na.rm=T) + outlier.sd * sd(xy.diff, na.rm=T) | xy.diff < mean(xy.diff, na.rm=T) - outlier.sd * sd(xy.diff, na.rm=T))
        return(x)
    })
    dat$sex.mismatch <- dat$declared.sex != dat$predicted.sex
    p1 <- ggplot(dat, aes(y=1, x=xy.diff)) +
        geom_jitter(aes(shape=predicted.sex, colour=sex.mismatch), size=3) +
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
#' @param samplesheet From \code{\link{read.450k.sheet}()}
#' @param  norm.objects From \code{\link{meffil.normalize.objects}()}
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


#' Plot the means of control probes for each sample and for each control probe type
#'
#' @param samplesheet From \code{\link{read.450k.sheet}()}
#' @param  norm.objects From \code{\link{meffil.normalize.objects}()}
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


#' Plot detection p values from idat files
#'
#' @param samplesheet From \code{\link{read.450k.sheet}()}
#' @param  norm.objects From \code{\link{meffil.normalize.objects}()}
#' @param  threshold Cut off value for proportion of CpGs with poor detection p values. Default 0.05
#' @param  colour.code Array of length n samples to colour code points. Defaults to NULL
#' @export
#' @return Data frame of results plus plot
#' @examples \dontrun{
#'
#'}
meffil.plot.detectionp.samples <- function(samplesheet, norm.objects, threshold = 0.05, colour.code=NULL)
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
    dat$outliers <- dat$prop.badprobes > threshold
    p1 <- ggplot(dat, aes(y=prop.badprobes, x=id)) +
        geom_point(aes(colour=colour.code)) +
        guides(colour = g) +
        geom_hline(yintercept=threshold) +
        labs(y = paste("Proportion CpG sites with p >", norm.objects[[1]]$bad.probes.detectionp.threshold), x = "Sample ID", colour = colour.code)
    return(list(graph=p1, tab=dat))
}

#' Manhattan plot of detection pval per probe - percentage with pvalue < 0.01
#'
#' @param samplesheet From \code{\link{read.450k.sheet}()}
#' @param  norm.objects From \code{\link{meffil.normalize.objects}()}
#' @param  threshold Cut off value for proportion of samples with poor detection p values. Default 0.05.
#' @export
#' @return Data frame of results plus plot
#' @examples \dontrun{
#'
#'}
meffil.plot.detectionp.cpgs <- function(samplesheet, norm.objects, threshold=0.05)
{
    n.badprobes = as.data.frame(table(unlist(sapply(norm.objects, function(x) names(x$bad.probes.detectionp)))))
    names(n.badprobes) <- c("name", "n")
    probe.info <- subset(meffil.probe.info(), !duplicated(name) & chr %in% paste("chr", c(1:22, "X", "Y"), sep=""))
    probe.info <- merge(probe.info, n.badprobes, by="name")
    probe.info$n[is.na(probe.info$n)] <- 0
    probe.info$n <- probe.info$n / length(norm.objects)
    probe.info$chr <- factor(gsub("chr", "", probe.info$chr), levels=c(1:22, "X", "Y"))
    probe.info$chr.colour <- 0
    probe.info$chr.colour[probe.info$chr %in% c(seq(1,22,2), "X")] <- 1
    probe.info <- subset(probe.info, select=c(name, chr, pos, n, chr.colour))
    probe.info$outliers <- probe.info$n > threshold
    p1 <- ggplot(probe.info, aes(x=pos, y=n)) +
        geom_point(aes(colour=chr.colour)) +
        facet_grid(. ~ chr, space="free_x", scales="free_x") +
        guides(colour=FALSE) +
        labs(x="Position", y=paste("Proportion samples with p >", norm.objects[[1]]$bad.probes.detectionp.threshold)) +
        geom_hline(yintercept=threshold) +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    return(list(graph=p1, tab=subset(probe.info, select=-c(chr.colour))))
}

#' Plot number of beads per sample
#'
#' @param samplesheet From \code{\link{read.450k.sheet}()}
#' @param  norm.objects From \code{\link{meffil.normalize.objects}()}
#' @param  threshold Cut off value for proportion of CpGs with low bead numbers. Default 0.05
#' @param  colour.code Array of length n samples to colour code points. Defaults to NULL
#' @export
#' @return Data frame of results plus plot
#' @examples \dontrun{
#'
#'}
meffil.plot.beadnum.samples <- function(samplesheet, norm.objects, threshold = 0.05, colour.code=NULL)
{
    nprobe <- length(unique(meffil.probe.info()$name))
    dat <- data.frame(
        Sample_Name = samplesheet$Sample_Name,
        prop.badprobes = sapply(norm.objects, function(x) length(x$bad.probes.beadnum) / nprobe)
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
    dat$outliers <- dat$prop.badprobes > threshold
    p1 <- ggplot(dat, aes(y=prop.badprobes, x=id)) +
        geom_point(aes(colour=colour.code)) +
        guides(colour = g) +
        geom_hline(yintercept=threshold) +
        labs(y = paste("Proportion CpG sites with bead number < ", norm.objects[[1]]$bad.probes.beadnum.threshold), x = "Sample ID", colour = colour.code)
    return(list(graph=p1, tab=dat))
}

#' Manhattan plot of number of beads by probe - percentage of probes with beads < 3 for each sample
#'
#' @param samplesheet From \code{\link{read.450k.sheet}()}
#' @param  norm.objects From \code{\link{meffil.normalize.objects}()}
#' @param  threshold Cut off value for proportion of samples with poor detection p values. Default 0.05.
#' @export
#' @return Data frame of results plus plot
#' @examples \dontrun{
#'
#'}
meffil.plot.beadnum.cpgs <- function(samplesheet, norm.objects, threshold = 0.05)
{
    n.badprobes = as.data.frame(table(unlist(sapply(norm.objects, function(x) names(x$bad.probes.beadnum)))))
    names(n.badprobes) <- c("name", "n")
    probe.info <- subset(meffil.probe.info(), !duplicated(name) & chr %in% paste("chr", c(1:22, "X", "Y"), sep=""))
    probe.info <- merge(probe.info, n.badprobes, by="name")
    probe.info$n[is.na(probe.info$n)] <- 0
    probe.info$n <- probe.info$n / length(norm.objects)
    probe.info$chr <- factor(gsub("chr", "", probe.info$chr), levels=c(1:22, "X", "Y"))
    probe.info$chr.colour <- 0
    probe.info$chr.colour[probe.info$chr %in% c(seq(1,22,2), "X")] <- 1
    probe.info <- subset(probe.info, select=c(name, chr, pos, n, chr.colour))
    probe.info$outliers <- probe.info$n > threshold
    p1 <- ggplot(probe.info, aes(x=pos, y=n)) +
        geom_point(aes(colour=chr.colour)) +
        facet_grid(. ~ chr, space="free_x", scales="free_x") +
        guides(colour=FALSE) +
        labs(x="Position", y=paste("Proportion samples with bead number <", norm.objects[[1]]$bad.probes.beadnum.threshold)) +
        geom_hline(yintercept=threshold) +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    return(list(graph=p1, tab=subset(probe.info, select=-c(chr.colour))))
}



# Bad quality probes
# Bad quality samples
# Gender mixups
# Sample mixups

#' Generate lists of bad probes and bad samples
#'
#' @param samplesheet From \code{\link{read.450k.sheet}()}
#' @param  norm.objects From \code{\link{meffil.normalize.objects}()}
#' @param  colour.code See above
#' @param  control.categories See above
#' @param  sex.outlier.sd See above
#' @param  meth.unmeth.outlier.sd See above
#' @param  control.means.outlier.sd See above
#' @param  detectionp.samples.threshold See above
#' @param  beadnum.samples.threshold See above
#' @param  detectionp.cpgs.threshold See above
#' @param  beadnum.cpgs.threshold See above
#' @export
#' @return List
#' @examples \dontrun{
#'
#'}
meffil.pre.processing <- function(samplesheet, norm.objects, colour.code, control.categories=names(norm.objects[[1]]$controls), sex.outlier.sd = 3, meth.unmeth.outlier.sd = 3, control.means.outlier.sd = 5, detectionp.samples.threshold = 0.05, beadnum.samples.threshold = 0.05, detectionp.cpgs.threshold = 0.05, beadnum.cpgs.threshold = 0.05)
{

    y.probes <- with(meffil.probe.info(), name[which(chr == "chrY")])
    for (i in 1:length(norm.objects)) 
        if (norm.objects[[i]]$predicted.sex == "F") {
            norm.objects[[i]]$bad.probes.detectionp <- setdiff(norm.objects[[i]]$bad.probes.detectionp, y.probes)
            norm.objects[[i]]$bad.probes.beadnum <- setdiff(norm.objects[[i]]$bad.probes.beadnum, y.probes)
        }       
    
    p1 <- meffil.plot.sex(samplesheet, norm.objects, outlier.sd=sex.outlier.sd)
    p2 <- meffil.plot.meth.unmeth(samplesheet, norm.objects, colour.code=colour.code, outlier.sd=meth.unmeth.outlier.sd)
    p3 <- meffil.plot.controlmeans(samplesheet, norm.objects, control.categories=control.categories, colour.code=colour.code, outlier.sd=control.means.outlier.sd)
    p4 <- meffil.plot.detectionp.samples(samplesheet, norm.objects, colour.code=colour.code, threshold=detectionp.samples.threshold)
    p5 <- meffil.plot.detectionp.cpgs(samplesheet, norm.objects, threshold=detectionp.cpgs.threshold)
    p6 <- meffil.plot.beadnum.samples(samplesheet, norm.objects, colour.code=colour.code, threshold=beadnum.samples.threshold)
    p7 <- meffil.plot.beadnum.cpgs(samplesheet, norm.objects, threshold=beadnum.cpgs.threshold)


    # Sex mismatches
    sex <- subset(p1$tab, sex.mismatch, select=c(Sample_Name, predicted.sex, declared.sex))


    # Bad quality samples
    sexo <- subset(p1$tab, outliers, select=c(Sample_Name))
    if(nrow(sexo) > 0) sexo$issue <- "X-Y ratio outlier"
    methunmeth <- subset(p2$tab, outliers, select=c(Sample_Name))
    if(nrow(methunmeth) > 0) methunmeth$issue <- "Methylated vs Unmethylated"
    controlmeans <- subset(p3$tab, outliers, select=c(Sample_Name, variable))
    names(controlmeans) <- c("Sample_Name", "issue")
    if(nrow(controlmeans) > 0) controlmeans$issue <- paste("Control probe (", controlmeans$issue, ")", sep="")
    detectionp <- subset(p4$tab, outliers, select=c(Sample_Name))
    if(nrow(detectionp) > 0) detectionp$issue <- "Detection p-value"
    beadnum <- subset(p6$tab, outliers, select=c(Sample_Name))
    if(nrow(beadnum) > 0) beadnum$issue <- "Low bead numbers"

    removeids <- rbind(methunmeth, controlmeans, detectionp, beadnum, sexo)
    removeids <- removeids[order(removeids$Sample_Name),]


    # Bad quality probes
    detectionp.cpg <- subset(p5$tab, outliers, select=c(name))
    if(nrow(detectionp.cpg) > 0) detectionp.cpg$issue <- "Detection p-value"
    beadnum.cpg <- subset(p7$tab, outliers, select=c(name))
    if(nrow(beadnum.cpg) > 0) beadnum.cpg$issue <- "Low bead number"

    removecpgs <- rbind(detectionp.cpg, beadnum.cpg)
    removecpgs <- ddply(removecpgs, .(name), summarise, issue = paste(issue, collapse=", "))

    return(list(sex.mismatch = sex, bad.ids = removeids, bad.cpgs = removecpgs, p1, p2, p3, p4, p5, p6, p7))
}




# Still to do:
# - Sample mismatches



