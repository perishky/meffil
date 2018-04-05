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
#' @param lambda.method Method for calculating genomic inflation lambda.
#' Valid values are "median", "regression", or "robust" (Default: "median").
#' @return List of \code{\link{ggplot}} for each analysis in \code{ewas.object}.
#' @export
meffil.ewas.qq.plot <- function(ewas.object,
                                sig.threshold=1e-7,
                                sig.color="red",
                                title="QQ plot",
                                xlab=bquote(-log[10]("expected p-values")),
                                ylab=bquote(-log[10]("observed p-values")),
                                lambda.method="median") {
    stopifnot(is.ewas.object(ewas.object))
    
    sapply(names(ewas.object$analyses), function(name) {
        p.values <- sort(ewas.object$analyses[[name]]$table$p.value, decreasing=T)
        p.values[which(p.values < .Machine$double.xmin)] <- .Machine$double.xmin
        stats <- data.frame(is.sig=p.values < sig.threshold,
                            expected=-log(sort(ppoints(p.values),decreasing=T),10),
                            observed=-log(p.values, 10))
        lambda <- qq.lambda(p.values[which(p.values > sig.threshold)],
                            method=lambda.method)

        label.x <- min(stats$expected) + diff(range(stats$expected))*0.1
        label.y <- min(stats$expected) + diff(range(stats$observed))*0.9

        lambda.label <- paste("lambda == ", format(lambda$estimate,digits=3),
                              "%+-%", format(lambda$se, digits=3),
                              "~(", lambda.method, ")", sep="")

        selection.idx <- scatter.thinning(stats$observed, stats$expected,
                                          resolution=100, max.per.cell=100)

        lim <- range(c(0, stats$expected, stats$observed))
        sig.threshold <- format(sig.threshold, digits=3)
        
        (ggplot(stats[selection.idx,], aes(x=expected, y=observed)) + 
         geom_abline(intercept = 0, slope = 1, colour="black") +              
         geom_point(aes(colour=factor(sign(is.sig)))) +
         scale_colour_manual(values=c("black", "red"),
                             name="Significant",
                             breaks=c("0","1"),
                             labels=c(paste("p-value >", sig.threshold),
                                 paste("p-value <", sig.threshold))) +
         annotate(geom="text", x=label.x, y=label.y, hjust=0,
                  label=lambda.label,
                  parse=T) +
         xlim(lim) + ylim(lim) + 
         xlab(xlab) + ylab(ylab) +
         coord_fixed() +
         ggtitle(paste(title, ": ", name, sep="")))
     }, simplify=F)     
}

qq.lambda <- function(p.values, method="median", B=100) {
    stopifnot(method %in% c("median","regression","robust"))
    p.values <- na.omit(p.values)
    observed <- qchisq(p.values, df=1, lower.tail = FALSE)
    observed <- sort(observed)
    expected <- qchisq(ppoints(length(observed)), df=1, lower.tail=FALSE)
    expected <- sort(expected)

    lambda <- se <- NA
    if (method == "median")  {
        lambda <- median(observed)/qchisq(0.5, df=1)
        boot.medians <- sapply(1:B, function(i) median(sample(observed, replace=T)))
        se <- sd(boot.medians/qchisq(0.5,df=1))
    } else if (method %in% c("regression","robust")) {
        if (method == "regression")
            coef.table <- summary(lm(observed ~ 0 + expected))$coeff
        else
            coef.table <- summary(rlm(observed ~ 0 + expected))$coef
        lambda <- coef.table["expected",1]
        se <- coef.table["expected", "Std. Error"]
    }
    list(method=method, estimate=lambda, se=se)
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
        p.values <- stats$p.value
        p.values[which(p.values < .Machine$double.xmin)] <- .Machine$double.xmin
        stats$stat <- -log(p.values,10) * sign(stats$coefficient)

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
         geom_point(aes(colour=chr.colour,size=stat,fill=chr.colour)) + ##!! edit - scaling points by -log10(p) !!##
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
    ## remove missing values
    idx <- which(!is.na(methylation) & !is.na(variable))
    if (length(idx) < 3) {
        warning("Not enough data points to plot CpG methylation.")
        return(NULL)
    }
    methylation <- methylation[idx]
    variable <- variable[idx]
    
    ## linear model fit
    if (is.null(covariates)) {
        fit <- lm(methylation ~ variable)
        base <- lm(methylation ~ 1)
    }
    else {
        covariates <- covariates[idx,,drop=F]
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
        min.value <- 1e-20
        meth <- methylation
        meth[which(meth < min.value)] <- min.value
        meth[which(meth > 1-min.value)] <- 1-min.value

        if (is.null(covariates)) {
            fit <- betareg(meth ~ variable)
            base <- betareg(meth ~ 1)
        }
        else {
            fit <- betareg(meth ~ variable + ., data=covariates)
            base <- betareg(meth ~ ., data=covariates)
        }
        p.value.beta <- lrtest(fit, base)[2,"Pr(>Chisq)"]

        stats.desc <- paste(stats.desc, "; p[beta] = ", format(p.value.beta, digits=3), sep="")
    }
    y.axis.label <- "DNA methylation"
    if (!is.null(covariates)) {
        methylation <- residuals(base)
        y.axis.label <- "adjusted DNA methylation"
    }

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
     xlab(stats.desc) + ylab(y.axis.label))
}
