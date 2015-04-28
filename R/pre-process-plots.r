library(ggplot2)
library(reshape2)
library(knitr)
library(rmarkdown)
library(pander)

#' Plot predicted sex
#'
#' @param  qc.objects From \code{meffil.normalize.objects}
#' @export
#' @return Data frame of results plus plot
#' @examples \dontrun{
#'
#'}
meffil.plot.sex <- function(qc.objects, outlier.sd=3)
{
    dat <- data.frame(
        Sample_Name = sapply(qc.objects, function(x) x$Sample_Name),
        xy.diff = sapply(qc.objects, function(x) x$xy.diff),
        predicted.sex = sapply(qc.objects, function(x) x$predicted.sex),
        declared.sex = sapply(qc.objects, function(x) x$sex)
    )
    dat <- ddply(dat, .(predicted.sex), function(x)
    {
        x <- mutate(x)
        x$outliers <- with(x, xy.diff > mean(xy.diff, na.rm=T) + outlier.sd * sd(xy.diff, na.rm=T) | xy.diff < mean(xy.diff, na.rm=T) - outlier.sd * sd(xy.diff, na.rm=T))
        return(x)
    })
    dat$sex.mismatch <- dat$declared.sex != dat$predicted.sex
    dat$sex.mismatch[is.na(dat$sex.mismatch)] <- "Sex not specified"
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
#' @param  qc.objects From \code{meffil.normalize.objects}
#' @param  outlier.sd Cut off for declaring outliers. Default = 3
#' @param  colour.code Array of length n samples to colour code points. Defaults to NULL
#' @export
#' @return Data frame of results plus plot
#' @examples \dontrun{
#'
#'}
meffil.plot.meth.unmeth <- function(qc.objects, outlier.sd=3, colour.code = NULL)
{
    dat <- data.frame(
        Sample_Name = sapply(qc.objects, function(x) x$Sample_Name),
        methylated = sapply(qc.objects, function(x) x$median.m.signal),
        unmethylated = sapply(qc.objects, function(x) x$median.u.signal)
    )
    if(length(colour.code) == nrow(dat)) {
        dat$colour.code <- colour.code
        g <- "legend"
    } else if(length(colour.code) == 1) {
        if(colour.code %in% names(qc.objects[[1]]$samplesheet)) {
            dat$colour.code <- sapply(qc.objects, function(x) x$samplesheet[[colour.code]])
            g <- "legend"
        } else {
            stop("colour.code unknown")
        }
    } else if(is.null(colour.code)) {
        dat$colour.code <- 1
        g <- FALSE
    } else {
        stop("colour.code unknown")
    }
    dat$resids <- residuals(lm(methylated ~ unmethylated, dat))
    dat$outliers <- dat$resids > mean(dat$resids)+outlier.sd*sd(dat$resids) | dat$resids < mean(dat$resids)-outlier.sd*sd(dat$resids)

    p1 <- ggplot(dat, aes(y=methylated, x=unmethylated)) +
        geom_point(aes(colour=colour.code)) +
        guides(colour = g) +
        geom_point(data=subset(dat, outliers), shape=1, size=3.5) +
        labs(y = "Median methylated signal", x = "Median unmethylated signal", colour = colour.code) +
        stat_smooth(method="lm", se=FALSE, colour="red") +
        geom_smooth()
    return(list(graph=p1, tab=dat))
}


#' Plot the means of control probes for each sample and for each control probe type
#'
#' @param  qc.objects From \code{meffil.normalize.objects}
#' @param  control.categories Which control probe categories to plot. Defaults to all available
#' @param  colour.code Array of length n samples to colour code points. Defaults to NULL
#' @param  outlier.sd Cut off for declaring outliers. Default = 5
#' @export
#' @return Data frame of results plus plot
#' @examples \dontrun{
#'
#'}
meffil.plot.controlmeans <- function(qc.objects, control.categories=names(qc.objects[[1]]$controls), colour.code = NULL, outlier.sd=5)
{
    dat <- data.frame(Sample_Name = sapply(qc.objects, function(x) x$Sample_Name))
    if(length(colour.code) == nrow(dat)) {
        dat$colour.code <- colour.code
        g <- "legend"
    } else if(length(colour.code) == 1) {
        if(colour.code %in% names(qc.objects[[1]]$samplesheet)) {
            dat$colour.code <- sapply(qc.objects, function(x) x$samplesheet[[colour.code]])
            g <- "legend"
        } else {
            stop("colour.code unknown")
        }
    } else if(is.null(colour.code)) {
        dat$colour.code <- 1
        g <- FALSE
    } else {
        stop("colour.code unknown")
    }

    d <- data.frame(t(sapply(qc.objects, function(x) x$controls)))
    names(d) <- names(qc.objects[[1]]$controls)
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
#' @param  qc.objects From \code{meffil.normalize.objects}
#' @param  threshold Cut off value for proportion of CpGs with poor detection p values. Default 0.05
#' @param  colour.code Array of length n samples to colour code points. Defaults to NULL
#' @export
#' @return Data frame of results plus plot
#' @examples \dontrun{
#'
#'}
meffil.plot.detectionp.samples <- function(qc.objects, threshold = 0.05, colour.code=NULL)
{
    nprobe <- length(unique(meffil.probe.info()$name))
    dat <- data.frame(
        Sample_Name = sapply(qc.objects, function(x) x$Sample_Name),
        prop.badprobes = sapply(qc.objects, function(x) length(x$bad.probes.detectionp) / nprobe)
    )
    if(length(colour.code) == nrow(dat)) {
        dat$colour.code <- colour.code
        g <- "legend"
    } else if(length(colour.code) == 1) {
        if(colour.code %in% names(qc.objects[[1]]$samplesheet)) {
            dat$colour.code <- sapply(qc.objects, function(x) x$samplesheet[[colour.code]])
            g <- "legend"
        } else {
            stop("colour.code unknown")
        }
    } else if(is.null(colour.code)) {
        dat$colour.code <- 1
        g <- FALSE
    } else {
        stop("colour.code unknown")
    }
    dat <- dat[order(dat$colour.code), ]
    dat$id <- 1:nrow(dat)
    dat$outliers <- dat$prop.badprobes > threshold
    p1 <- ggplot(dat, aes(y=prop.badprobes, x=id)) +
        geom_point(aes(colour=colour.code)) +
        guides(colour = g) +
        geom_hline(yintercept=threshold) +
        labs(y = paste("Proportion CpG sites with p >", qc.objects[[1]]$bad.probes.detectionp.threshold), x = "Sample ID", colour = colour.code)
    return(list(graph=p1, tab=dat))
}

#' Manhattan plot of detection pval per probe - percentage with pvalue < 0.01
#'
#' @param  qc.objects From \code{meffil.normalize.objects}
#' @param  threshold Cut off value for proportion of samples with poor detection p values. Default 0.05.
#' @export
#' @return Data frame of results plus plot
#' @examples \dontrun{
#'
#'}
meffil.plot.detectionp.cpgs <- function(qc.objects, threshold=0.05)
{
    n.badprobes = as.data.frame(table(unlist(sapply(qc.objects, function(x) names(x$bad.probes.detectionp)))))
    names(n.badprobes) <- c("name", "n")
    probe.info <- subset(meffil.probe.info(), !duplicated(name) & chr %in% paste("chr", c(1:22, "X", "Y"), sep=""))
    probe.info <- merge(probe.info, n.badprobes, by="name")
    probe.info$n[is.na(probe.info$n)] <- 0
    probe.info$n <- probe.info$n / length(qc.objects)
    probe.info$chr <- factor(gsub("chr", "", probe.info$chr), levels=c(1:22, "X", "Y"))
    probe.info$chr.colour <- 0
    probe.info$chr.colour[probe.info$chr %in% c(seq(1,22,2), "X")] <- 1
    probe.info <- subset(probe.info, select=c(name, chr, pos, n, chr.colour))
    probe.info$outliers <- probe.info$n > threshold
    p1 <- ggplot(probe.info, aes(x=pos, y=n)) +
        geom_point(aes(colour=chr.colour)) +
        facet_grid(. ~ chr, space="free_x", scales="free_x") +
        guides(colour=FALSE) +
        labs(x="Position", y=paste("Proportion samples with p >", qc.objects[[1]]$bad.probes.detectionp.threshold)) +
        geom_hline(yintercept=threshold) +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    return(list(graph=p1, tab=subset(probe.info, select=-c(chr.colour))))
}

#' Plot number of beads per sample
#'
#' @param  qc.objects From \code{meffil.normalize.objects}
#' @param  threshold Cut off value for proportion of CpGs with low bead numbers. Default 0.05
#' @param  colour.code Array of length n samples to colour code points. Defaults to NULL
#' @export
#' @return Data frame of results plus plot
#' @examples \dontrun{
#'
#'}
meffil.plot.beadnum.samples <- function(qc.objects, threshold = 0.05, colour.code=NULL)
{
    nprobe <- length(unique(meffil.probe.info()$name))
    dat <- data.frame(
        Sample_Name = sapply(qc.objects, function(x) x$Sample_Name),
        prop.badprobes = sapply(qc.objects, function(x) length(x$bad.probes.beadnum) / nprobe)
    )
    if(length(colour.code) == nrow(dat)) {
        dat$colour.code <- colour.code
        g <- "legend"
    } else if(length(colour.code) == 1) {
        if(colour.code %in% names(qc.objects[[1]]$samplesheet)) {
            dat$colour.code <- sapply(qc.objects, function(x) x$samplesheet[[colour.code]])
            g <- "legend"
        } else {
            stop("colour.code unknown")
        }
    } else if(is.null(colour.code)) {
        dat$colour.code <- 1
        g <- FALSE
    } else {
        stop("colour.code unknown")
    }
    dat <- dat[order(dat$colour.code), ]
    dat$id <- 1:nrow(dat)
    dat$outliers <- dat$prop.badprobes > threshold
    p1 <- ggplot(dat, aes(y=prop.badprobes, x=id)) +
        geom_point(aes(colour=colour.code)) +
        guides(colour = g) +
        geom_hline(yintercept=threshold) +
        labs(y = paste("Proportion CpG sites with bead number < ", qc.objects[[1]]$bad.probes.beadnum.threshold), x = "Sample ID", colour = colour.code)
    return(list(graph=p1, tab=dat))
}

#' Manhattan plot of number of beads by probe - percentage of probes with beads < 3 for each sample
#'
#' @param  qc.objects From \code{meffil.normalize.objects}
#' @param  threshold Cut off value for proportion of samples with poor detection p values. Default 0.05.
#' @export
#' @return Data frame of results plus plot
#' @examples \dontrun{
#'
#'}
meffil.plot.beadnum.cpgs <- function(qc.objects, threshold = 0.05)
{
    n.badprobes = as.data.frame(table(unlist(sapply(qc.objects, function(x) names(x$bad.probes.beadnum)))))
    names(n.badprobes) <- c("name", "n")
    probe.info <- subset(meffil.probe.info(), !duplicated(name) & chr %in% paste("chr", c(1:22, "X", "Y"), sep=""))
    probe.info <- merge(probe.info, n.badprobes, by="name")
    probe.info$n[is.na(probe.info$n)] <- 0
    probe.info$n <- probe.info$n / length(qc.objects)
    probe.info$chr <- factor(gsub("chr", "", probe.info$chr), levels=c(1:22, "X", "Y"))
    probe.info$chr.colour <- 0
    probe.info$chr.colour[probe.info$chr %in% c(seq(1,22,2), "X")] <- 1
    probe.info <- subset(probe.info, select=c(name, chr, pos, n, chr.colour))
    probe.info$outliers <- probe.info$n > threshold
    p1 <- ggplot(probe.info, aes(x=pos, y=n)) +
        geom_point(aes(colour=chr.colour)) +
        facet_grid(. ~ chr, space="free_x", scales="free_x") +
        guides(colour=FALSE) +
        labs(x="Position", y=paste("Proportion samples with bead number <", qc.objects[[1]]$bad.probes.beadnum.threshold)) +
        geom_hline(yintercept=threshold) +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    return(list(graph=p1, tab=subset(probe.info, select=-c(chr.colour))))
}


control.probe.categories <- function()
{
    return(c("bisulfite1", "bisulfite2", "extension.G.31698466", "extension.G.74666473", "extension.R.47640365", "extension.R.63642461", "hybe.21771417", "hybe.26772442", "hybe.28684356", "stain.G", "stain.R", "nonpoly.G.23663352", "nonpoly.G.70645401", "nonpoly.R.18773482", "nonpoly.R.24701411", "targetrem.13643320", "targetrem.42790394", "spec1.G.10673427", "spec1.G.23777311", "spec1.G.51804467", "spec1.R.46779338", "spec1.R.53740460", "spec1.R.59783305", "spec2.G.17661470", "spec2.G.29662396", "spec2.G.34730329", "spec2.R.17661470", "spec2.R.29662396", "spec2.R.34730329", "spec1.ratio1", "spec1.ratio", "spec2.ratio", "spec1.ratio2", "normA", "normC", "normT", "normG", "dye.bias", "oob.G.1%", "oob.G.50%", "oob.G.99%", "oob.ratio"))
}


#' Specify parameters for QC
#'
#' 
#' @param  colour.code Default value = NULL <what param does>
#' @param  control.categories Default value = control.probe.categories() <what param does>
#' @param  sex.outlier.sd Default value = 3 <what param does>
#' @param  meth.unmeth.outlier.sd Default value = 3 <what param does>
#' @param  control.means.outlier.sd Default value = 5 <what param does>
#' @param  detectionp.samples.threshold Default value = 0.05 <what param does>
#' @param  beadnum.samples.threshold Default value = 0.05 <what param does>
#' @param  detectionp.cpgs.threshold Default value = 0.05 <what param does>
#' @param  beadnum.cpgs.threshold Default value = 0.05 <what param does>
#' @export
#' @return List of parameter values
#' @examples \dontrun{
#'
#'}
meffil.qc.parameters <- function(colour.code = NULL, control.categories = control.probe.categories(), sex.outlier.sd = 3, meth.unmeth.outlier.sd = 3, control.means.outlier.sd = 5, detectionp.samples.threshold = 0.05, beadnum.samples.threshold = 0.05, detectionp.cpgs.threshold = 0.05, beadnum.cpgs.threshold = 0.05) {
    parameters <- list(
        colour.code = colour.code, 
        control.categories = control.categories,
        sex.outlier.sd = sex.outlier.sd,
        meth.unmeth.outlier.sd = meth.unmeth.outlier.sd,
        control.means.outlier.sd = control.means.outlier.sd,
        detectionp.samples.threshold = detectionp.samples.threshold,
        beadnum.samples.threshold = beadnum.samples.threshold,
        detectionp.cpgs.threshold = detectionp.cpgs.threshold,
        beadnum.cpgs.threshold = beadnum.cpgs.threshold
    )
    return(parameters)
}


#' Perform QC analysis on idat files
#'
#' Performs a number of QC analyses including checking for sex differences, methylated vs unmethylated levels, 
#' deviation from control probe means, detection p-values and bead numbers per sample and probe.
#'
#' Also returns list of sample IDs and CPGs that are low quality.
#'
#' @param  qc.objects From \code{meffil.normalize.objects}
#' @param  parameters Default = meffil.qc.parameters(). List of parameter values. See \code{\link{meffil.qc.parameters}}
#' @export
#' @return List
#' @examples \dontrun{
#'
#'}
meffil.qc.summary <- function(qc.objects, parameters = meffil.qc.parameters(), verbose=TRUE) {
    msg("Sex summary", verbose)
    sex.summary <- meffil.plot.sex(
        qc.objects,
        outlier.sd=parameters$sex.outlier.sd
    )
    msg("Meth vs unmeth summary", verbose=verbose)
    meth.unmeth.summary <- meffil.plot.meth.unmeth(
        qc.objects, 
        colour.code=parameters$colour.code, 
        outlier.sd=parameters$meth.unmeth.outlier.sd
    )
    msg("Control means summary", verbose=verbose)
    controlmeans.summary <- meffil.plot.controlmeans(
        qc.objects,
        control.categories=parameters$control.categories,
        colour.code=parameters$colour.code,
        outlier.sd=parameters$control.means.outlier.sd
    )
    msg("Sample detection summary", verbose=verbose)
    sample.detectionp.summary <- meffil.plot.detectionp.samples(
        qc.objects,
        colour.code=parameters$colour.code,
        threshold=parameters$detectionp.samples.threshold
    )
    msg("CpG detection summary", verbose=verbose)
    cpg.detectionp.summary <- meffil.plot.detectionp.cpgs(
        qc.objects,
        threshold=parameters$detectionp.cpgs.threshold
    )
    msg("Sample bead numbers summary", verbose=verbose)
    sample.beadnum.summary <- meffil.plot.beadnum.samples(
        qc.objects,
        colour.code=parameters$colour.code,
        threshold=parameters$beadnum.samples.threshold
    )
    msg("CpG bead numbers summary", verbose=verbose)
    cpg.beadnum.summary <- meffil.plot.beadnum.cpgs(
        qc.objects,
        threshold=parameters$beadnum.cpgs.threshold
    )


    # Sex mismatches
    sex <- subset(sex.summary$tab, sex.mismatch == "FALSE" | outliers, select=c(Sample_Name, predicted.sex, declared.sex, xy.diff))


    # Bad quality samples
    sexo <- subset(sex.summary$tab, outliers, select=c(Sample_Name))
    if(nrow(sexo) > 0) sexo$issue <- "X-Y ratio outlier"
    methunmeth <- subset(meth.unmeth.summary$tab, outliers, select=c(Sample_Name))
    if(nrow(methunmeth) > 0) methunmeth$issue <- "Methylated vs Unmethylated"
    controlmeans <- subset(controlmeans.summary$tab, outliers, select=c(Sample_Name, variable))
    names(controlmeans) <- c("Sample_Name", "issue")
    if(nrow(controlmeans) > 0) controlmeans$issue <- paste("Control probe (", controlmeans$issue, ")", sep="")
    detectionp <- subset(sample.detectionp.summary$tab, outliers, select=c(Sample_Name))
    if(nrow(detectionp) > 0) detectionp$issue <- "Detection p-value"
    beadnum <- subset(sample.beadnum.summary$tab, outliers, select=c(Sample_Name))
    if(nrow(beadnum) > 0) beadnum$issue <- "Low bead numbers"

    removeids <- rbind(methunmeth, controlmeans, detectionp, beadnum, sexo)
    removeids <- removeids[order(removeids$Sample_Name),]


    # Bad quality probes
    detectionp.cpg <- subset(cpg.detectionp.summary$tab, outliers, select=c(name))
    if(nrow(detectionp.cpg) > 0) detectionp.cpg$issue <- "Detection p-value"
    beadnum.cpg <- subset(cpg.beadnum.summary$tab, outliers, select=c(name))
    if(nrow(beadnum.cpg) > 0) beadnum.cpg$issue <- "Low bead number"

    removecpgs <- rbind(detectionp.cpg, beadnum.cpg)
    removecpgs <- ddply(removecpgs, .(name), summarise, issue = paste(issue, collapse=", "))


    return(list(
        sex.check = sex,
        bad.ids = removeids,
        bad.cpgs = removecpgs,
        parameters = parameters,
        sex.summary = sex.summary,
        meth.unmeth.summary = meth.unmeth.summary,
        controlmeans.summary = controlmeans.summary,
        sample.detectionp.summary = sample.detectionp.summary,
        cpg.detectionp.summary = cpg.detectionp.summary,
        sample.beadnum.summary = sample.beadnum.summary,
        cpg.beadnum.summary = cpg.beadnum.summary
    ))
}


#' Generate QC report
#'
#' Generate HTML file that summarises the QC. 
#'
#' @param  qc.summary Output from \code{meffil.qc.summary}.
#' @param  output.file Default = "meffil.qc.report.html"
#' If specified then a html report will be generated summarising the QC.
#' @param  author Default = "Analyst". Author name to be specified on report.
#' @param  studyname Default = "IlluminaHuman450 data". Study name to be specified on report.
#' @param  ... Arguments to be passed to \code{\link{rmarkdown::render}}
#' @export
#' @return List of tables and graphs describing QC.
#' @examples \dontrun{
#'
#'}
meffil.qc.report <- function(
    qc.summary,
    output.file = "meffil.qc.report.html",
    author = "Analyst",
    studyname = "IlluminaHuman450 data",
    ...
) {
    cat("Writing report as html file to", output.file, "\n")
    save(qc.summary,
        author,
        studyname,
        file = file.path(tempdir(), "meffil.qc.report.rdata")
    )
    output.dir <- ifelse(dirname(output.file) == ".", getwd(), dirname(output.file))
    rmarkdown::render(system.file("reports", "meffil.qc.report.rmd", package="meffil"), output_file=basename(output.file), output_dir=output.dir, ...)
}



#' Remove samples from QC objects
#'
#'
#' @param  qc.objects Output from \code{\link{meffil.qc}}
#' @param  idlist.remove Array of Sample_Name IDs to be removed
#' @export
#' @return qc.objects with samples removed
#' @examples \dontrun{
#'
#'}
meffil.remove.ids <- function(qc.objects, idlist.remove)
{
    stopifnot(all(idlist.remove %in% names(qc.objects)))
    qc.objects <- qc.objects[!names(qc.objects) %in% idlist.remove]
    return(qc.objects)
}

# Still to do:
# - Sample mismatches
