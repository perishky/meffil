#' Remove points from a scatter plot where density is really high
#' @param x x-coordinates vector
#' @param y y-coordinates vector
#' @param resolution number of partitions for the x and y-dimensions.
#' @param max.per.cell maximum number of points per x-y partition.
#' @return index into the points that omits points from x-y partitions
#' so that each has at most \code{max.per.cell} points.
scatter.thinning <- function(x,y,resolution=100,max.per.cell=100) {
    x.cell <- floor((resolution-1)*(x - min(x,na.rm=T))/diff(range(x,na.rm=T))) + 1
    y.cell <- floor((resolution-1)*(y - min(y,na.rm=T))/diff(range(y,na.rm=T))) + 1
    z.cell <- x.cell * resolution + y.cell
    frequency.table <- table(z.cell)
    frequency <- rep(0,max(z.cell, na.rm=T))
    frequency[as.integer(names(frequency.table))] <- frequency.table
    f.cell <- frequency[z.cell]
    
    big.cells <- length(which(frequency > max.per.cell))
    sort(c(which(f.cell <= max.per.cell),
           sample(which(f.cell > max.per.cell),
                  size=big.cells * max.per.cell, replace=F)),
         decreasing=F)
}

#' QQ plot for EWAS
#'
#' @param ewas.object Return object from \code{\link{meffil.ewas()}}.
#' @param sig.threshold P-value threshold for significance (Default: 1e-7).
#' @param sig.color Color for points corresponding to significant tests (Default: "red").
#' @param title Title for the plot (Default: "QQ plot").
#' @param xlab Label for the x-axis (Default: -log_10(expected p-values)).
#' @param ylab Label for the y-axis (Default: -log_10(observed p-values)).
#' @return \code{\link{ggplot}} showing the QQ plot. 
#' @export
meffil.ewas.qq.plot <- function(ewas.object,
                           sig.threshold=1e-7,
                           sig.color="red",
                           title="QQ plot",
                           xlab=bquote(-log[10]("expected p-values")),
                           ylab=bquote(-log[10]("observed p-values"))) {
    stopifnot(is.ewas.object(ewas.object))
    
    sapply(names(ewas.object$analyses), function(name) {
        p.values <- sort(ewas.object$analyses[[name]]$table$p.value, decreasing=T)
        stats <- data.frame(is.sig=p.values < sig.threshold,
                            expected=-log((length(p.values):1 - 0.5)/length(p.values),10),
                            observed=-log(p.values, 10))
        lambda <- median(qchisq(p.values,1,lower.tail=FALSE), na.rm=T)/qchisq(0.5,1)

        label.x <- min(stats$expected) + diff(range(stats$expected))*0.25
        label.y <- min(stats$expected) + diff(range(stats$observed))*0.75

        selection.idx <- scatter.thinning(stats$observed, stats$expected,
                                          resolution=100, max.per.cell=100)

        (ggplot(stats[selection.idx,], aes(x=expected, y=observed)) + 
         geom_abline(intercept = 0, slope = 1, colour="black") +              
         geom_point(aes(colour=factor(sign(is.sig)))) +
         scale_colour_manual(values=c("black", "red"),
                             name="Significant",
                             breaks=c("0","1"),
                             labels=c(paste("p-value >", sig.threshold),
                                 paste("p-value <", sig.threshold))) +
         annotate(geom="text", x=label.x, y=label.y,
                  label=paste("lambda ==", format(lambda, digits=3)), parse=T) +
         xlab(xlab) + ylab(ylab) +
         ggtitle(paste(title, ": ", name, sep="")))
     }, simplify=F)     
}

#' Manhattan plot for EWAS
#'
#' @param ewas.object Return object from \code{\link{meffil.ewas()}}.
#' @param sig.threshold P-value threshold for significance (Default: 1e-7).
#' @param title Title for the plot (Default: "Manhattan plot").
#' @return \code{\link{ggplot}} showing the Manhattan plot. 
#' @export
meffil.ewas.manhattan.plot <- function(ewas.object, sig.threshold=1e-7,
                                       title="Manhattan plot") {
    stopifnot(is.ewas.object(ewas.object))
    
    chromosomes <- paste("chr", c(1:22, "X","Y"), sep="")
    sapply(names(ewas.object$analyses), function(name) {
        stats <- ewas.object$analyses[[name]]$table
        stats$chromosome <- factor(as.character(stats$chromosome), levels=chromosomes)
        stats$chr.colour <- 0
        stats$chr.colour[stats$chromosomes %in% chromosomes[seq(1,length(chromosomes),2)]] <- 1
        stats$stat <- -log(stats$p.value,10) * sign(stats$coefficient)

        stats <- stats[order(stats$stat, decreasing=T),]

        chromosome.lengths <- sapply(chromosomes, function(chromosome)
                                     max(stats$position[which(stats$chromosome == chromosome)]))
        chromosome.lengths <- as.numeric(chromosome.lengths)
        chromosome.starts <- c(1,cumsum(chromosome.lengths)+1)
        names(chromosome.starts) <- c(chromosomes, "NA")
        stats$global <- stats$position + chromosome.starts[stats$chromosome] - 1

        selection.idx <- scatter.thinning(stats$global, stats$stat,
                                          resolution=100, max.per.cell=100)
        
        (ggplot(stats[selection.idx,], aes(x=position, y=stat)) +
         geom_point(aes(colour=chr.colour)) +
         facet_grid(. ~ chromosome, space="free_x", scales="free_x") +
         theme(strip.text.x = element_text(angle = 90)) +
         guides(colour=FALSE) +
         labs(x="Position",
              y=bquote(-log[10]("p-value") * sign(beta))) +             
         geom_hline(yintercept=log(sig.threshold,10), colour="red") +
         geom_hline(yintercept=-log(sig.threshold,10), colour="red") +
         theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
         ggtitle(paste(title, ": ", name, sep=""))) 
    }, simplify=F)        
}


#' Scatter plots for a CpG site in an EWAS
#'
#' @param ewas.object Return object from \code{\link{meffil.ewas()}}.
#' @param cpg CpG site to plot.
#' @param title Title of the plot (Default: \code{cpg}).
#' @param beta Matrix of methylation levels used to create the \code{ewas.object}.
#' @param \code{\link{ggplot}} object showing the scatterplots of DNA methylation vs the variable of interest
#' in the EWAS.  Each plot corresponds to a covariate set.
#' Methylation levels are in fact residuals from fitting a model with DNA methylation and the covariates.
#' 
#' @export
meffil.ewas.cpg.plot <- function(ewas.object, cpg, beta, title=cpg) {
    stopifnot(is.ewas.object(ewas.object))
    stopifnot(is.matrix(beta) && cpg %in% rownames(beta))
    
    variable <- ewas.object$variable
    
    lapply(names(ewas.object$analyses), function(name) {
        ewas <- ewas.object$analyses[[name]]

        if (!all(rownames(ewas$design) %in% colnames(beta)))
            stop("EWAS samples do not match those in the beta argument (methylation matrix)")
        methylation <- beta[cpg,rownames(ewas$design)]

        covariates <- subset(data.frame(ewas$design), select=c(-variable,-intercept))
        if (ncol(covariates) == 0)
            covariates <- NULL
        cpg.plot(methylation, variable, covariates, title=paste(name, ": ", title, sep=""))
    })
}

cpg.plot <- function(methylation, variable, covariates=NULL, title="") {
    ## linear model fit
    if (is.null(covariates)) {
        fit <- lm(methylation ~ variable)
        base <- lm(methylation ~ 1)
    }
    else {
        fit <- lm(methylation ~ variable + ., data=covariates)
        base <- lm(methylation ~ ., data=covariates)
    }
    p.value.lm <- anova(fit,base)[2,"Pr(>F)"]

    stats.desc <- paste("variable\np[lm]= ", format(p.value.lm, digits=3), sep="")
                        
    has.betareg <- all(c("lmtest", "betareg") %in% rownames(installed.packages()))
    if (has.betareg) {
        require("betareg")
        require("lmtest")
        ## beta regression model fit
        if (is.null(covariates)) {
            fit <- betareg(methylation ~ variable)
            base <- betareg(methylation ~ 1)
        }
        else {
            fit <- betareg(methylation ~ variable + ., data=covariates)
            base <- betareg(methylation ~ ., data=covariates)
        }
        p.value.beta <- lrtest(fit, base)[2,"Pr(>Chisq)"]

        stats.desc <- paste(stats.desc, "; p[beta] = ", format(p.value.beta, digits=3), sep="")
    }
    if (!is.null(covariates))
        methylation <- residuals(base)

    ## plot
    data <- data.frame(methylation=methylation, variable=variable)
    if (is.factor(variable) || length(unique(variable)) <= 20) {
        data$variable <- as.factor(data$variable)
        p <- (ggplot(data, aes(x=variable, y=methylation)) +
              geom_boxplot())
    } else {
        p <- (ggplot(data, aes(x=variable, y=methylation)) +
              geom_point() + geom_smooth(method=lm))
    }

    (p + ggtitle(title) +
     xlab(stats.desc) + ylab("DNA methylation"))
}
