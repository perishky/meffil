#' Describe EWAS samples using the variable of interest and covariates.
#'
#' @param ewas.object Output of \code{\link{meffil.ewas}()}.
#' @return A data frame with one row for each continuous or ordinal variable
#' and one row for each level of each categorical variable.
#' In the first case, each row provides the name, mean value and standard deviation of each variable.
#' In the second case (categorical), each row provides the name of the variable level and
#' the number of cases and percentage of cases at that level.
#' 
#' @export
meffil.ewas.sample.characteristics <- function(ewas.object) {
    stopifnot(is.ewas.object(ewas.object))

    msg("summarizing variables")
    summarize.variable <- function(name, variable) {
        msg(name)
        if (is.character(variable)) variable <- as.factor(variable)
        if (is.factor(variable)) {
            n <- sapply(levels(variable), function(level) sum(variable == level,na.rm=T))
            data.frame(variable=name,
                       value=names(n),
                       mean=format(n),
                       var=format(n/length(variable)*100,digits=3,nsmall=1),
                       row.names=NULL)                   
        }
        else {
            data.frame(variable=name,
                       value="",
                       mean=format(mean(variable,na.rm=T)),
                       var=format(sd(variable,na.rm=T)),
                       row.names=NULL)
        }
    }

    var.summary <- summarize.variable("variable of interest", ewas.object$variable)
    cov.summary <- NULL
    if (!is.null(ewas.object$covariates)) {
        cov.summary <- lapply(1:ncol(ewas.object$covariates), function(i) {
            summarize.variable(colnames(ewas.object$covariates)[i],
                               ewas.object$covariates[,i])
        })
        cov.summary <- do.call(rbind, cov.summary)
    }
    rbind(var.summary, cov.summary)
}

#' Describe associations between EWAS covariates and the variable of interest.
#'
#' @param ewas.object Output of \code{\link{meffil.ewas}()}.
#' @return A data frame with one or more rows for each covariate.
#'
#' If both the variable of interest and covariate are continuous or ordinal,
#' then the covariate uses one row showing the name, mean
#' and standard deviation of the covariate following the significance
#' of the association between the covariate and the variable of interest.
#'
#' If the covariate is categorical, then there is additionally one row for
#' each level showing the mean and standard deviation of the variable
#' of interest for samples at that covariate level.
#'
#' If the variable of interest is categorical but the covariate is not,
#' then there is one row for each variable level showing the
#' mean and standard deviation of the covariate at the given level.
#'
#' If both the variable of interest and covariate are categorical,
#' then mean is replaced with the number of samples at each pair of
#' variable/categorical levels and standard deviation with the percentage.
#' P-values indicate the significance of association using Fisher's exact test.
#' 
#' @export
meffil.ewas.covariate.associations <- function(ewas.object) {
    stopifnot(is.ewas.object(ewas.object))
    if (is.null(ewas.object$covariates)) return(NULL)

    msg("covariate associations")
    ret <- lapply(1:ncol(ewas.object$covariates), function(i) {
        msg(colnames(ewas.object$covariates)[i])
        vars <- data.frame(variable=ewas.object$variable,
                           covariate=ewas.object$covariates[,i])
        for (j in colnames(vars)) {
            if (!is.numeric(vars[[j]])
                && length(unique(na.omit(vars[[j]]))) <= 2)
                vars[[j]] <- as.factor(vars[[j]])
        }
        colnames(vars)[2] <- colnames(ewas.object$covariates)[i]
        meffil.summarize.relationship(vars)
    })
    names(ret) <- colnames(ewas.object$covariates)
    ret
}

#' Describe the relationship between two variables.
#'
#' @param vars A data frame with at least two columns.
#' The first two columns will be compared.
#' @return A list whose elements depends on the types of the two variables.
#' In each case, the list contains the following elements:
#' \describe{
#' \item{var1}{Name of the first variable, i.e. colnames(vars)[1].}
#' \item{var2}{Name of the second variable.}
#' \item{r}{Spearman's correlation between the variables. This may be meaningless if one variable is an unordered factor.}
#' \item{r.p}{P-value corresponding to the correlation between the variables.}
#' \item{output}{The contents of the list formatted to be printed as markdown text.}
#' \item{plot}{A plot (ggplot2) visualizing the relationship.}
#' }
#'
#' If both variables are factors, then the list will include a frequency table (\code{freq})
#' and a corresponding matrix of p-values (\code{p.values}) obtained using Fisher's test to test for
#' enrichment in each cell of the frequency table.
#'
#' If one variable is a factor and the other numeric, then list will include
#' the F-statistic (\code{F.stat}) and p-value (\code{p.value}) from one-way analysis of variance.
#' There will also be a data frame (\code{cases}) with each row providing statistics from a t-test
#' comparing the numeric variable within and without each level of the factor variable.
#'
#' If both variables are numeric, then the list will include
#' the F-statistic (\code{F.stat}) and p-value (\code{p.value})
#' from the linear model fitting the variables.
#'
#' @export
meffil.summarize.relationship <- function(vars) {
    stopifnot(is.data.frame(vars))
    values <- lapply(vars, function(var) if (!is.factor(var)) NA else levels(var))
    samples <- do.call(expand.grid, c(values, list(stringsAsFactors=F)))
    samples$idx <- lapply(1:nrow(samples), function(i)
                          which((is.na(samples[[1]][i]) | vars[[1]] == samples[[1]][i]) &
                                (is.na(samples[[2]][i]) | vars[[2]] == samples[[2]][i])))
    if (is.factor(vars[[1]]) && is.factor(vars[[2]])) {        
        freq <- do.call(table, vars)
        ret <- list(freq=freq)
        fit <- cor.test(as.integer(vars[[1]]), as.integer(vars[[2]]), method="s")
        ret$r <- fit$estimate[["rho"]]
        ret$r.p <- fit$p.value
        ret$p.values <- sapply(values[[1]], function(level1)
                               sapply(values[[2]], function(level2) {
                                   is.value <- mapply(function(var,value) var == value,
                                                      vars,
                                                      list(level1, level2))
                                   freq <- do.call(table, as.data.frame(is.value))
                                   p <- NA
                                   try(p <- fisher.test(freq, alternative="greater")$p.value, silent=T)
                                   p
                               }))
        ret$p.values <- t(ret$p.values)
        ret$sig <- ret$p.values
        ret$sig[] <- ""
        ret$sig[which(ret$p.values < 0.05)] <- "+"
        ret$sig[which(ret$p.values < 0.01)] <- "*"
        ret$sig[which(ret$p.values < 0.001)] <- "**"
        ret$sig[which(ret$p.values < 0.0001)] <- "***"
        ret$sig[which(ret$p.values < 0.00001)] <- "****"
    } else {
        var1 <- sign(is.factor(vars[[1]])) + 1 ## variable that is not a factor
        var2 <- 3-var1 ## the other variable
        fit <- lm(vars[[var1]] ~ vars[[var2]])        
        f.stat <- summary(fit)$fstatistic
        p.value <- 1 - pf(f.stat["value"],f.stat["numdf"],f.stat["dendf"])[["value"]]

        ret <- list(f.stat=f.stat[["value"]],
                    p.value=p.value,
                    n=sum(!is.na(vars[[var1]]) & !is.na(vars[[var2]])))
        
        if (!is.factor(vars[[var1]]) && !is.factor(vars[[var2]])) {
            fit <- cor.test(vars[[1]], vars[[2]], method="s")
            ret$r <- fit$estimate[["rho"]]
            ret$r.p <- fit$p.value
        } else { ## vars[[var2]] is a factor
            fit <- cor.test(vars[[var1]], as.integer(vars[[var2]]), method="s")
            ret$r <- fit$estimate[["rho"]]
            ret$r.p <- fit$p.value
                         
            ret$cases <- do.call(rbind, lapply(1:nrow(samples), function(i) {
                level <- samples[[var2]][i]
                is.level <- vars[[var2]] == level

                if (sum(is.level, na.rm=T) > 1 && sum(!is.level, na.rm=T) > 1)
                    fit <- t.test(vars[[var1]][which(is.level)], 
                                  vars[[var1]][which(!is.level)])
                else
                    fit <- list(p.value=NA,statistic=list(t=NA))

                data.frame(level=level,
                           mean=mean(vars[[var1]][which(is.level)], na.rm=T),
                           var=var(vars[[var1]][which(is.level)], na.rm=T),
                           n=sum(is.level, na.rm=T),
                           t.stat=fit$statistic[["t"]],
                           p.value=fit$p.value)
            }))
            names(ret$cases)[1] <- names(vars)[var2]
            rownames(ret$cases) <- NULL
        }
    }

    ret$var1 <- names(vars)[1]
    ret$var2 <- names(vars)[2]

    ret$output <- format.relationship(ret)
    ret$plot <- visualize.relationship(vars, ret)
    
    ret
}


format.relationship <- function(ret) {
    if ("p.values" %in% names(ret)) {        
        list("statistics"=kable(data.frame(var1=ret$var1, var2=ret$var2,
                 R=ret$r, "p-value"=ret$r.p, check.names=F)),
             "frequencies"=kable(ret$freq),
             "enrichment p-values"=kable(ret$p.values))
    }
    else if ("cases" %in% names(ret)) {
        list("statistics"=kable(data.frame(var1=ret$var1, var2=ret$var2,
                              F=ret$f.stat, "p-value"=ret$p.value,
                              R=ret$r, "p-value"=ret$r.p, check.names=F)),
             "cases"=kable(ret$cases))
    }
    else {
        list("statistics"=kable(data.frame(var1=ret$var1, var2=ret$var2,
                 F=ret$f.stat, "p-value"=ret$p.value,
                 R=ret$r, "p-value"=ret$r.p, check.names=F)))
    }
}

visualize.relationship <- function(vars, ret) {
    if ("p.values" %in% names(ret)) {
        data <- melt(-log(ret$p.values,10))
        colnames(data) <- c("var1","var2","p.value")
        data$freq <- melt(ret$freq)$value
        (ggplot(data, aes(x=as.factor(var2), y=as.factor(var1)))
         + geom_tile(aes(fill=p.value))
         + geom_text(aes(label=freq))
         + labs(x=names(vars)[2], y=names(vars)[1], fill="-log10(p-value)")
         + theme(axis.text.x = element_text(angle = 90, hjust = 1),
                 aspect.ratio=1)
         + ggtitle(""))
    }
    else {
        if (!is.factor(vars[[1]]) && !is.factor(vars[[2]])) {
            n.levels <- sapply(vars, function(var) length(unique(var)))
            if (n.levels[1] <= n.levels[2] && n.levels[1] < 20)
                vars[[1]] <- factor(vars[[1]], levels=sort(unique(vars[[1]])), ordered=T)
            if (n.levels[2] < n.levels[1] && n.levels[2] < 20)
                vars[[2]] <- factor(vars[[2]], levels=sort(unique(vars[[2]])), ordered=T)
        }
        if (is.factor(vars[[2]])) {
            vars[1:2] <- vars[2:1]
            names(vars)[1:2] <- names(vars)[2:1]
        } ## if one is a factor, then it is the first one

        data <- data.frame(x=vars[[1]], y=vars[[2]])
        if (!is.factor(data$x) && !is.factor(data$y)) {        
            (ggplot(data, aes(x=x, y=y))
             + geom_point()
             + geom_smooth(method=lm)
             + xlab(names(vars)[1])
             + ylab(names(vars)[2])
             + ggtitle(""))
        } else {            
            (ggplot(data, aes(x=x,y=y))
             + geom_boxplot()
             + xlab(names(vars)[1])
             + ylab(names(vars)[2])
             + theme(axis.text.x = element_text(angle = 90, hjust = 1))
             + ggtitle(""))
        } 
    }
}
        
