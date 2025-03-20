#' Generate QC report
#'
#' Generate HTML file that summarises the QC. 
#'
#' @param  qc.summary Output from \code{meffil.qc.summary}.
#' @param  output.file Default = "meffil-qc-report.html".
#' If the file extension is not .htm, .html, .HTM or .HTML then
#' output will be in markdown format.
#' @param  author Default = "Analyst". Author name to be specified on report.
#' @param  study Default = "Illumina methylation data". Study name to be specified on report.
#' @param  ... Arguments to be passed to \code{\link{knitr::knit}}
#' @export
#' @return NULL
#' @examples \dontrun{
#'
#'}
meffil.qc.report <- function(
    qc.summary,
    output.file = "qc-report.html",
    author = "Analyst",
    study = "Illumina methylation data",
    ...
) {
    msg("Writing report as html file to", output.file)
    report.path <- system.file("reports", package="meffil")
    require(knitr)
    require(Cairo)
    require(gridExtra)
    opts <- opts_chunk$get()
    on.exit(opts_chunk$set(opts))
    opts_chunk$set(warning=FALSE, echo=FALSE, message=FALSE, results="asis", fig.width=10, fig.height=10, dev="CairoPNG")
    knit.report(file.path(report.path, "qc-report.rmd"), output.file, ...)
}


#' Perform QC analysis on idat files
#'
#' Performs a number of QC analyses including checking for sex differences, methylated vs unmethylated levels, 
#' deviation from control probe means, detection p-values and bead numbers per sample and probe.
#'
#' Also returns list of sample IDs and CPGs that are low quality.
#'
#' @param  qc.objects From \code{meffil.qc}
#' @param genotypes Optional output from \code{\link{meffil.extract.genotypes}()}.
#' Sample genotypes are matched to sample qc.objects using
#' \code{colnames(genotypes)} and \code{names(qc.objects)}.
#' @param  parameters Default = meffil.qc.parameters(). List of parameter values. See \code{\link{meffil.qc.parameters}}
#' @export
#' @return List
#' @examples \dontrun{
#'
#'}
meffil.qc.summary <- function(qc.objects, genotypes = NULL,
                              parameters = meffil.qc.parameters(), verbose=TRUE) {
    stopifnot(sapply(qc.objects, is.qc.object))
    
    parameters$detection.threshold <- qc.objects[[1]]$bad.probes.detectionp.threshold
    parameters$bead.threshold <- qc.objects[[1]]$bad.probes.beadnum.threshold
    parameters$sex.cutoff <- qc.objects[[1]]$sex.cutoff
    
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
    msg("Cell count summary", verbose=verbose)
    cell.counts.summary <- meffil.plot.cell.counts(
        qc.objects
    )
    msg("Genotype concordance", verbose=verbose)
    genotype.summary <- meffil.plot.genotypes(
        qc.objects,
        genotypes=genotypes,
        snp.threshold=parameters$snp.concordance.threshold,
        sample.threshold=parameters$sample.genotype.concordance.threshold)
    
    # Sex mismatches
    sex.check <- subset(sex.summary$tab, as.character(sex.mismatch) == "TRUE" | outliers,
                  select=c(sample.name, predicted.sex, declared.sex, xy.diff, status))
    sex.check <- sex.check[with(sex.check, order(status,predicted.sex,declared.sex,xy.diff)),]
  
    # Bad quality samples
    ## 1a. X-Y ratio outliers
    sexo <- subset(sex.summary$tab, outliers, select=c(sample.name))
    if(nrow(sexo) > 0) 
        sexo$issue <- "X-Y ratio outlier"
    ## 1b. sex mismatches
    sexm <- subset(sex.summary$tab, as.character(sex.mismatch) == "TRUE", select=c(sample.name))
    if (nrow(sexm) > 0)
        sexm$issue <- "Sex mismatch"
    ## 2. meth/unmeth outliers
    methunmeth <- subset(meth.unmeth.summary$tab, outliers, select=c(sample.name))
    if(nrow(methunmeth) > 0) methunmeth$issue <- "Methylated vs Unmethylated"
    ## 3. control means outliers
    controlmeans <- subset(controlmeans.summary$tab, outliers, select=c(sample.name, variable))
    names(controlmeans) <- c("sample.name", "issue")
    if(nrow(controlmeans) > 0) controlmeans$issue <- paste("Control probe (", controlmeans$issue, ")", sep="")
    ## 4. too many undetected probes
    detectionp <- subset(sample.detectionp.summary$tab, outliers, select=c(sample.name))
    if(nrow(detectionp) > 0) detectionp$issue <- "Detection p-value"
    ## 5. too many low bead number probes
    beadnum <- subset(sample.beadnum.summary$tab, outliers, select=c(sample.name))
    if(nrow(beadnum) > 0) beadnum$issue <- "Low bead numbers"
    ## 6. genotype mismatches
    geno <- NULL
    if (!is.null(genotype.summary$graphs$snp.concordance)) {
        geno <- subset(genotype.summary$tabs$samples, !is.concordant, select=c(sample.name))
        if (nrow(geno) > 0)
            geno$issue <- "Genotype mismatch"
    }
    
    bad.samples <- rbind(methunmeth, controlmeans, detectionp, beadnum, sexo, sexm, geno)
    bad.samples <- bad.samples[order(bad.samples$sample.name),,drop=F]

    # Bad quality probes
    detectionp.cpg <- subset(cpg.detectionp.summary$tab, outliers, select=c(name))
    if(nrow(detectionp.cpg) > 0) detectionp.cpg$issue <- "Detection p-value"
    beadnum.cpg <- subset(cpg.beadnum.summary$tab, outliers, select=c(name))
    if(nrow(beadnum.cpg) > 0) beadnum.cpg$issue <- "Low bead number"

    bad.cpgs <- rbind(detectionp.cpg, beadnum.cpg)
    bad.cpgs <- ddply(bad.cpgs, .(name), summarise, issue = paste(issue, collapse=", "))


    return(list(
        sex.check = sex.check,
        bad.samples = bad.samples,
        bad.cpgs = bad.cpgs,
        parameters = parameters,
        sex.summary = sex.summary,
        meth.unmeth.summary = meth.unmeth.summary,
        controlmeans.summary = controlmeans.summary,
        sample.detectionp.summary = sample.detectionp.summary,
        cpg.detectionp.summary = cpg.detectionp.summary,
        sample.beadnum.summary = sample.beadnum.summary,
        cpg.beadnum.summary = cpg.beadnum.summary,
        cell.counts.summary = cell.counts.summary,
        genotype.summary = genotype.summary
    ))
}


#' Plot predicted sex
#'
#' @param  qc.objects From \code{meffil.qc}
#' @export
#' @return Data frame of results plus plot
#' @examples \dontrun{
#'
#'}
meffil.plot.sex <- function(qc.objects, outlier.sd=3)
{
    stopifnot(sapply(qc.objects, is.qc.object))
     
    dat <- data.frame(
        sample.name = sapply(qc.objects, function(x) x$sample.name),
        xy.diff = sapply(qc.objects, function(x) x$xy.diff),
        predicted.sex = sapply(qc.objects, function(x) x$predicted.sex),
        declared.sex = sapply(qc.objects, function(x) x$sex),
        stringsAsFactor=FALSE
    )
    
    dat <- ddply(dat, .(predicted.sex), function(x)
    {
        x <- mutate(x)
        interval.size <- outlier.sd * sd(x$xy.diff, na.rm=T)
        x$upper <- mean(x$xy.diff, na.rm=T) + interval.size
        x$lower <- mean(x$xy.diff, na.rm=T) - interval.size
        x$outliers <- with(x, xy.diff > upper | xy.diff < lower)
        return(x)
    })
    
    dat$sex.mismatch <- as.character(dat$declared.sex) != as.character(dat$predicted.sex)
    dat$sex.mismatch[is.na(dat$sex.mismatch)] <- "Unspecified sex"
    p1 <- (ggplot(dat, aes(y=1, x=xy.diff)) +
           geom_jitter(aes(shape=predicted.sex, colour=sex.mismatch), size=3) +
           scale_colour_manual(values=c("black", "red", "blue")) + ## alphabetic order: FALSE=black, TRUE=red, Unspecified sex=blue
           labs(shape="Predicted sex", x="XY diff", y="", colour="Incorrect\nprediction") +
           theme_bw() +
           theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.line.y=element_blank()) +
           geom_vline(xintercept=unique(c(dat$upper, dat$lower)), linetype="dashed", colour="purple"))

    dat$status <- "good"
    dat$status[which(dat$outliers)] <- "outlier"
    dat$status[which(as.character(dat$sex.mismatch) == "TRUE")] <- "mismatched"
    return(list(graph=p1, tab=dat))
}


#' Plot average methylated vs unmethylated levels for each individuals
#'
#' plot raw control probes and fit linear regression, remove samples that have sd(y - yhat) > mean*3
#'
#' @param  qc.objects From \code{meffil.qc}
#' @param  outlier.sd Cut off for declaring outliers. Default = 3
#' @param  colour.code Array of length n samples to colour code points. Defaults to NULL
#' @export
#' @return Data frame of results plus plot
#' @examples \dontrun{
#'
#'}
meffil.plot.meth.unmeth <- function(qc.objects, outlier.sd=3, colour.code = NULL)
{
    stopifnot(sapply(qc.objects, is.qc.object))

    dat <- data.frame(
        sample.name = sapply(qc.objects, function(x) x$sample.name),
        methylated = sapply(qc.objects, function(x) x$median.m.signal),
        unmethylated = sapply(qc.objects, function(x) x$median.u.signal),
        stringsAsFactors=F
    )

    fit <- lm(methylated ~ unmethylated, dat)
    dat$resids <- residuals(fit)
    dat$methylated.lm <- predict(fit)
    interval.size <- outlier.sd*sd(dat$resids)
  
    dat$upper.lm <- dat$methylated.lm + interval.size
    dat$lower.lm <- dat$methylated.lm - interval.size

    dat$outliers <- (dat$resids > mean(dat$resids) + interval.size
                     | dat$resids < mean(dat$resids) - interval.size)

    if (!is.null(colour.code)) {
        if(length(colour.code) == nrow(dat)) {
            dat$colour.code <- colour.code
            g <- "legend"
        } else {
            if(colour.code %in% names(qc.objects[[1]]$samplesheet)) {
                dat$colour.code <- sapply(qc.objects, function(x) x$samplesheet[[colour.code]])
                g <- "legend"
            } else {
                stop("colour.code unknown")
            }
        }

        p1 <- (ggplot(dat, aes(y=methylated, x=unmethylated)) +
               geom_point(aes(colour=colour.code)) +
               guides(colour = g))
        if (any(dat$outliers))
            p1 <- p1 + geom_point(data=subset(dat, outliers), aes(colour=colour.code), size=3.5)
    } else {
        dat$slide <- as.character(sapply(qc.objects, function(object) object$samplesheet$Slide))
        outlier.counts <- table(factor(dat$slide)[dat$outliers])
        dat$slide.outlier.count <- outlier.counts[dat$slide]
        dat$outlier.slide <- "no outlier"
        outlier.slide.idx <- which(dat$slide.outlier.count > 0)
        dat$outlier.slide[outlier.slide.idx] <- dat$slide[outlier.slide.idx]
        dat$outlier.slide <- paste(dat$outlier.slide,
                                   dat$slide.outlier.count,
                                   c("array","arrays")[sign(dat$slide.outlier.count > 1) + 1])
        
        if (length(unique(dat$outlier.slide)) > 1) 
            g <- "legend"
        else
            g <- FALSE

        p1 <- (ggplot(dat, aes(y=methylated, x=unmethylated)) +
               geom_point(colour="black") +
               guides(colour=g))
        if (any(dat$outliers))
            p1 <- p1 + geom_point(data=subset(dat, outliers), aes(colour=outlier.slide), size=3.5)
        dat$slide <- dat$slide.outlier.count <- dat$outlier.slide <- NULL
    }

    p1 <- (p1 +
           labs(y = "Median methylated signal", x = "Median unmethylated signal") + 
           geom_smooth() +
           geom_line(aes(y=methylated.lm), col="red") +
           geom_line(aes(y=upper.lm), col="red", linetype="dashed") +
           geom_line(aes(y=lower.lm), col="red", linetype="dashed"))

    return(list(graph=p1, tab=dat))
}


#' Plot the means of control probes for each sample and for each control probe type
#'
#' @param  qc.objects From \code{meffil.qc}
#' @param  control.categories Which control probe categories to plot. Defaults to all available.
#' @param  colour.code Array of length n samples to colour code points. Defaults to NULL
#' @param  outlier.sd Cut off for declaring outliers. Default = 5
#' @export
#' @return Data frame of results plus plot
#' @examples \dontrun{
#'
#'}
meffil.plot.controlmeans <- function(qc.objects, control.categories=NULL, colour.code = NULL, outlier.sd=5)
{
    stopifnot(sapply(qc.objects, is.qc.object))

    dat <- data.frame(sample.name = sapply(qc.objects, function(x) x$sample.name),
                      stringsAsFactors=F)

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

    d <- data.frame(t(sapply(qc.objects, function(x) x$controls)), stringsAsFactors=F)
    names(d) <- names(qc.objects[[1]]$controls)

    if (is.null(control.categories))
        control.categories <- names(qc.objects[[1]]$controls)
    
    dat <- data.frame(dat, subset(d, select=control.categories),  stringsAsFactors=F)
    names(dat) <- c("sample.name", "colour.code", control.categories)
    dat <- dat[order(dat$colour.code), ]
    dat$colour.code <- as.character(dat$colour.code)
    dat$id <- 1:nrow(dat)
    dat <- reshape2::melt(dat, id.vars=c("sample.name", "colour.code", "id"))
    dat <- ddply(dat, .(variable), function(x)
    {
        x <- mutate(x)
        x$outliers <- x$value > mean(x$value) + outlier.sd * sd(x$value) | x$value < mean(x$value) - outlier.sd * sd(x$value)
        return(x)
    })
    p1 <- (ggplot(dat, aes(x=id, y=value)) +
           geom_point(aes(colour=colour.code)) +
           guides(colour=g) +
           facet_wrap(~ variable, scales="free_y", ncol=4) +
           labs(y="Mean signal", x="ID", colour=colour.code))
    if (any(dat$outliers))
        p1 <- p1 + geom_point(data=subset(dat, outliers), shape=1, size=3.5)
    return(list(graph=p1, tab=dat))
}


#' Plot detection p values from idat files
#'
#' @param  qc.objects From \code{meffil.qc}
#' @param  threshold Cut off value for proportion of CpGs with poor detection p values. Default 0.05
#' @param  colour.code Array of length n samples to colour code points. Defaults to NULL
#' @export
#' @return Data frame of results plus plot
#' @examples \dontrun{
#'
#'}
meffil.plot.detectionp.samples <- function(qc.objects, threshold = 0.05, colour.code=NULL)
{
    stopifnot(sapply(qc.objects, is.qc.object))
    featureset <- qc.objects[[1]]$featureset

    if (is.null(featureset)) { ## backwards compatibility
        featureset <- "450k"
    }   

    probes <- meffil.get.features(featureset)$name
    y.probes <- meffil.get.y.sites(featureset)
    not.y.probes <- setdiff(probes, y.probes)
        
    dat <- data.frame(
        sample.name = sapply(qc.objects, function(x) x$sample.name),
        prop.badprobes = sapply(qc.objects, function(x) {
            bad.probes <- intersect(probes, names(x$bad.probes.detectionp))
            if (as.character(x$predicted.sex) == "F") {
                bad.probes <- setdiff(bad.probes, y.probes)
                probes <- not.y.probes
            }
            length(bad.probes)/length(probes)
        }),
        stringsAsFactors=F
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
    p1 <- (ggplot(dat, aes(y=prop.badprobes, x=id)) +
           geom_point(aes(colour=colour.code)) +
           guides(colour = g) +
           geom_hline(yintercept=threshold) +
           labs(y = paste("Proportion CpG sites with p >",
                    qc.objects[[1]]$bad.probes.detectionp.threshold),
                x = "Sample ID",
                colour = colour.code))
    return(list(graph=p1, tab=dat))
}

#' Manhattan plot of detection pval per probe - percentage with pvalue < 0.01
#'
#' @param  qc.objects From \code{meffil.qc}
#' @param  threshold Cut off value for proportion of samples with poor detection p values. Default 0.05.
#' @export
#' @return Data frame of results plus plot
#' @examples \dontrun{
#'
#'}
meffil.plot.detectionp.cpgs <- function(qc.objects, threshold=0.05)
{
    stopifnot(sapply(qc.objects, is.qc.object))

    featureset <- qc.objects[[1]]$featureset

    if (is.null(featureset)) { ## backwards compatibility
        featureset <- "450k"
    }   
    
    y.probes <- meffil.get.y.sites(featureset)
    bad.probes <- unlist(lapply(qc.objects, function(x) {
        bad.probes <- names(x$bad.probes.detectionp)
        if (as.character(x$predicted.sex) == "F")
            bad.probes <- setdiff(bad.probes, y.probes)
        bad.probes
    }))
    if (length(bad.probes) == 0)
        return(list(graph=NULL, tab=NULL))
    
    n.badprobes <- as.data.frame(table(bad.probes), stringsAsFactors=F)
    names(n.badprobes) <- c("name", "n")
    
    probes <- meffil.get.features(featureset)
    probes <- merge(probes, n.badprobes, by="name")
    probes$n[is.na(probes$n)] <- 0
    
    n.males <- sum(sapply(qc.objects, function(x) as.character(x$predicted.sex) == "M"))
    probes$n.samples <- length(qc.objects)
    probes$n.samples[which(probes$name %in% y.probes)] <- n.males
    
    probes$n <- probes$n / probes$n.samples
    probes$chromosome <- factor(gsub("chr", "", probes$chromosome), levels=c(1:22, "X", "Y"))
    probes$chr.colour <- 0
    probes$chr.colour[probes$chromosome %in% c(seq(1,22,2), "X")] <- 1
    probes <- subset(probes, select=c(name, chromosome, position, n, chr.colour))
    probes$outliers <- probes$n > threshold
    p1 <- (ggplot(probes, aes(x=position, y=n)) +
           geom_point(aes(colour=chr.colour)) +
           facet_grid(. ~ chromosome, space="free_x", scales="free_x") +
           guides(colour=FALSE) +
           labs(x="Position", y=paste("Proportion samples with p >", qc.objects[[1]]$bad.probes.detectionp.threshold)) +
           geom_hline(yintercept=threshold) +
           theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()))
    return(list(graph=p1, tab=subset(probes, select=-c(chr.colour))))
}

#' Plot number of beads per sample
#'
#' @param  qc.objects From \code{meffil.qc}
#' @param  threshold Cut off value for proportion of CpGs with low bead numbers. Default 0.05
#' @param  colour.code Array of length n samples to colour code points. Defaults to NULL
#' @export
#' @return Data frame of results plus plot
#' @examples \dontrun{
#'
#'}
meffil.plot.beadnum.samples <- function(qc.objects, threshold = 0.05, colour.code=NULL)
{
    stopifnot(sapply(qc.objects, is.qc.object))

    featureset <- qc.objects[[1]]$featureset

    if (is.null(featureset)) { ## backwards compatibility
        featureset <- "450k"
    }   
    
    probes <- meffil.get.features(featureset)$name
    y.probes <- meffil.get.y.sites(featureset)
    not.y.probes <- setdiff(probes, y.probes)

    dat <- data.frame(
        sample.name = sapply(qc.objects, function(x) x$sample.name),
        prop.badprobes = sapply(qc.objects, function(x) {
            bad.probes.beadnum <- intersect(names(x$bad.probes.beadnum), probes)
            if (as.character(x$predicted.sex) == "F") {
                bad.probes.beadnum <- setdiff(bad.probes.beadnum, y.probes)
                probes <- not.y.probes
            }
            length(bad.probes.beadnum)/length(probes)
        }),
	stringsAsFactors=F
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
    p1 <- (ggplot(dat, aes(y=prop.badprobes, x=id)) +
           geom_point(aes(colour=colour.code)) +
           guides(colour = g) +
           geom_hline(yintercept=threshold) +
           labs(y = paste("Proportion CpG sites with bead number < ", qc.objects[[1]]$bad.probes.beadnum.threshold),
                x = "Sample ID",
                colour = colour.code))
    return(list(graph=p1, tab=dat))
}

#' Manhattan plot of number of beads by probe - percentage of probes with beads < 3 for each sample
#'
#' @param  qc.objects From \code{meffil.qc}
#' @param  threshold Cut off value for proportion of samples with poor detection p values. Default 0.05.
#' @export
#' @return Data frame of results plus plot
#' @examples \dontrun{
#'
#'}
meffil.plot.beadnum.cpgs <- function(qc.objects, threshold = 0.05)
{
    stopifnot(sapply(qc.objects, is.qc.object))

    featureset <- qc.objects[[1]]$featureset

    if (is.null(featureset)) { ## backwards compatibility
        featureset <- "450k"
    }   

    y.probes <- meffil.get.y.sites(featureset)
    bad.probes.beadnum <- unlist(lapply(qc.objects, function(x) {
        bad.probes.beadnum <- names(x$bad.probes.beadnum)
        if (as.character(x$predicted.sex) == "F")
            bad.probes.beadnum <- setdiff(bad.probes.beadnum, y.probes)
        bad.probes.beadnum
    }))
    if (length(bad.probes.beadnum) == 0)
        return(list(graph=NULL, tab=NULL))

    n.badprobes <- as.data.frame(table(bad.probes.beadnum), stringsAsFactors=F)
    names(n.badprobes) <- c("name", "n")
    probes <- meffil.get.features(featureset)
    probes <- merge(probes, n.badprobes, by="name")
    probes$n[is.na(probes$n)] <- 0
    probes$n <- probes$n / length(qc.objects)
    probes$chromosome <- factor(gsub("chr", "", probes$chromosome), levels=c(1:22, "X", "Y"))
    probes$chr.colour <- 0
    probes$chr.colour[probes$chromosome %in% c(seq(1,22,2), "X")] <- 1
    probes <- subset(probes, select=c(name, chromosome, position, n, chr.colour))
    probes$outliers <- probes$n > threshold
    p1 <- (ggplot(probes, aes(x=position, y=n)) +
           geom_point(aes(colour=chr.colour)) +
           facet_grid(. ~ chromosome, space="free_x", scales="free_x") +
           guides(colour=FALSE) +
           labs(x="Position", y=paste("Proportion samples with bead number <", qc.objects[[1]]$bad.probes.beadnum.threshold)) +
           geom_hline(yintercept=threshold) +
           theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()))
    return(list(graph=p1, tab=subset(probes, select=-c(chr.colour))))
}


#' Cell count estimate quality plot
#'
#' @param qc.objects Output from \code{\link{meffil.qc}()}.
#' @param reference Object describing methylation profiles of purified cell populations
#' obtained from \code{\link{meffil.add.cell.type.reference}()}.
#' @return Two \link{ggplot2} boxplot objects:
#' - \code{betas} Contains one box per sample or reference cell type
#' representing the distribution of methylation levels for the CpG sites used to
#' estimate cell counts.
#' - \code{counts} Contains one box per reference cell type
#' representing the distribution of cell count estimates across the samples.
#' 
#' @export
meffil.plot.cell.counts <- function(qc.objects) {
    stopifnot(sapply(qc.objects, is.qc.object))
    
    count.objects <- lapply(qc.objects, function(object) object$cell.counts)
    not.null.idx <- which(sapply(count.objects, function(object) !is.null(object)))
    if (length(not.null.idx) > 0) 
        meffil.cell.count.qc.plots(count.objects[not.null.idx])
    else
        NULL
}


#' Plot SNP beta and sample genotype concordances
#'
#' @param qc.objects Output from \code{\link{meffil.qc}()}.
#' @param genotypes Optional output from \code{\link{meffil.extract.genotypes}()}.
#' Sample genotypes are matched to sample qc.objects using
#' \code{colnames(genotypes)} and \code{names(qc.objects)}.
#' @param sample.threshold Concordance threshold below which the Illumina 450K
#' and genetic profiles for a sample are deemed a mismatch (Default: 0.9).
#' @param snp.threshold Concordance threshold below which the Illumina 450K
#' and genetic profiles for a SNP are deemed a mismatch (Default: 0.99).
#' @return A list consisting of:
#' - \code{graphs} A list of \link{ggplot2} objects.  The first \code{snp.betas}
#' plots the beta distributions of each SNP probe in the microarray.
#' The second and third plots are added only if the \code{genotypes} parameter
#' is not \code{NULL}. The second plot shows the distribution of SNP concordances,
#' and the third plot shows the distribution of sample concordances.
#' - \code{tabs} Contains two data frames if the
#' \code{genotypes} parameter is not \code{NULL}.  The first \code{samples}
#' lists the concordances of each sample, the second \code{snps} lists the
#' concordances of each SNP.
#' 
#' @export
meffil.plot.genotypes <- function(qc.objects, genotypes=NULL,
                                  sample.threshold=0.9, snp.threshold=0.99) {
    stopifnot(sapply(qc.objects, is.qc.object))
    
    graphs <- list()
    tabs <- list()
    
    snp.betas <- meffil.snp.betas(qc.objects)

    show.snps <- rownames(snp.betas)[1:min(65,nrow(snp.betas))]
    data <- lapply(show.snps, function(snp.name) 
                   data.frame(snp=snp.name, beta=snp.betas[snp.name,],stringsAsFactors=F))
    data <- do.call(rbind, data)
    graphs$snp.beta <- (ggplot(data=data, aes(beta)) +
                        geom_histogram() +
                        facet_wrap(~snp))

    if (!is.null(genotypes)) {
        genotypes <- genotypes[,which(colSums(is.na(genotypes)) < nrow(genotypes)), drop=F]
        genotypes <- genotypes[which(rowSums(is.na(genotypes)) < ncol(genotypes)),, drop=F]
        
        common.samples <- intersect(colnames(snp.betas), colnames(genotypes))
        common.snps <- intersect(rownames(snp.betas), rownames(genotypes))

        if (length(common.snps) == 0)
            stop("SNPs in genotypes agument are not measured on these microarrays")
        if (length(common.samples) == 0)
            stop("genotype samples do not match samples with QC objects")
        
        concordance <- meffil.snp.concordance(snp.betas[common.snps, common.samples, drop=F],
                                              genotypes[common.snps, common.samples, drop=F],
                                              snp.threshold=snp.threshold,
                                              sample.threshold=sample.threshold)
        
        graphs$snp.concordance <- (ggplot(data=data.frame(concordance=concordance$snp),
                                          aes(concordance)) +
                                   geom_histogram())
        
        graphs$sample.concordance <- (ggplot(data=data.frame(concordance=concordance$sample),
                                             aes(concordance)) +
                                      geom_histogram())

        tabs$samples <- with(concordance, data.frame(sample.name=names(sample),concordance=sample))
        tabs$samples$is.concordant <- tabs$samples$concordance > sample.threshold
        
        tabs$snps <- with(concordance, data.frame(snp.name=names(snp),concordance=snp))
        tabs$snps$is.concordant <- tabs$snps$concordance > snp.threshold        
    }
    
    list(graphs=graphs,
         tabs=tabs)                     
}



#' Specify parameters for QC
#'
#' 
#' @param  colour.code Default value = NULL 
#' @param  control.categories Default value = control.probe.categories()
#' @param  sex.outlier.sd Sets the standard deviation multiple at which sex outliers are identified. Default value = 3.
#' @param  meth.unmeth.outlier.sd Sets the standard deviation multiple at which methylated/unmethylated signal outliers are identified. Default value = 3.
#' @param  control.means.outlier.sd Sets the standard deviation multiple at which control probe signals are identified as outliers. Default value = 5 
#' @param  beadnum.samples.threshold Maximum threshold on the fraction of probes with too few detected beads (minimum number of detected beads is defined by setting the `beads.threshold` parameter of `meffil.qc()` or `meffil.normalize.dataset()`). Samples with probe fractions above this will be excluded from the final dataset. Default value = 0.2
#' @param  detectionp.samples.threshold Maximum threshold on the fraction of undetected probes (probe detection is defined by setting the maximum probe detection p-value threshold parameter `detection.threshold` of `meffil.qc()` or `meffil.normalize.dataset()`).  Samples with probe fractions above this will be excluded from the final dataset. Default value = 0.2
#' @param  beadnum.cpgs.threshold Same as `beadnum.samples.threshold` but used to identify poor quality probes in terms of the fraction of samples in which the probe has too few detected beads. Default value = 0.2
#' @param  detectionp.cpgs.threshold Same as `detectionp.cpgs.threshold` but used to identify poor quality probes in terms of the fraction of samples in which the probe is undetected. Default value = 0.2
#' @param  snp.concordance.threshold Minimum required concordance between supplied genotypes and genotypes estimated from a SNP probe. Default value = 0.99
#' @param  sample.genotype.concordance.threshold Minimum required concordance between supplied genotypes and genotypes estimated from SNP probes for a given individual. Default value = 0.9 <what param does>
#' @export
#' @return List of parameter values
#' @examples \dontrun{
#'
#'}
meffil.qc.parameters <- function(colour.code = NULL, control.categories = NULL, sex.outlier.sd = 3,
                                 meth.unmeth.outlier.sd = 3, control.means.outlier.sd = 5,
                                 detectionp.samples.threshold = 0.2,
                                 beadnum.samples.threshold = 0.2, detectionp.cpgs.threshold = 0.2,
                                 beadnum.cpgs.threshold = 0.2,
                                 snp.concordance.threshold = 0.9,
                                 sample.genotype.concordance.threshold = 0.9
                                 ) {
    parameters <- list(
        colour.code = colour.code, 
        control.categories = control.categories,
        sex.outlier.sd = sex.outlier.sd,
        meth.unmeth.outlier.sd = meth.unmeth.outlier.sd,
        control.means.outlier.sd = control.means.outlier.sd,
        detectionp.samples.threshold = detectionp.samples.threshold,
        beadnum.samples.threshold = beadnum.samples.threshold,
        detectionp.cpgs.threshold = detectionp.cpgs.threshold,
        beadnum.cpgs.threshold = beadnum.cpgs.threshold,        
        snp.concordance.threshold = snp.concordance.threshold,
        sample.genotype.concordance.threshold = sample.genotype.concordance.threshold

    )
    return(parameters)
}




#' Remove samples from QC objects
#'
#'
#' @param  qc.objects Output from \code{\link{meffil.qc}}
#' @param  sample.ids Array of sample.name IDs to be removed
#' @export
#' @return qc.objects with samples removed
#' @examples \dontrun{
#'
#'}
meffil.remove.samples <- function(qc.objects, sample.ids)
{
    stopifnot(sapply(qc.objects, is.qc.object))
    stopifnot(all(sample.ids %in% names(qc.objects)))
    qc.objects <- qc.objects[!names(qc.objects) %in% sample.ids]
    return(qc.objects)
}




