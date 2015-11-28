#' Describe EWAS samples using variable of interest and covariates.
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
    
    summarize.variable <- function(name, variable) {
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

    rbind(summarize.variable("variable of interest", ewas.object$variable),
          do.call(rbind, lapply(1:ncol(ewas.object$covariates), function(i) {
              summarize.variable(colnames(ewas.object$covariates)[i],
                                 ewas.object$covariates[,i])
          })))
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

    summarize.covariate <- function(name,covariate, variable) {
        samples <- expand.grid(variable=if (!is.factor(variable)) "" else levels(variable),
                               covariate=if (!is.factor(covariate)) "" else levels(covariate),
                               stringsAsFactors=F)
        samples$idx <- lapply(1:nrow(samples), function(i)
                              which((is.na(samples$variable[i])
                                     | variable == samples$variable[i]) &
                                    (is.na(samples$covariate[i])
                                     | covariate == samples$covariate[i])))
    
        if (is.factor(variable) && is.factor(covariate)) {
            contingency.table <- table(variable,covariate)
            all.ret <- data.frame(covariate=name,
                                  value="",
                                  variable="",
                                  mean="",
                                  var="",
                                  p.value=format(fisher.test(contingency.table, simulate.p.value=T)$p.value,digits=3,nsmall=2),
                                  test="Fisher's test")
            cases.ret <- do.call(rbind, lapply(1:nrow(samples), function(i) {
                contingency.table <- table(variable==samples$variable[i],
                                           covariate==samples$covariate[i])
                data.frame(covariate=name,
                           value=samples$covariate[i], ## covariate level
                           variable=samples$variable[i], ## variable level
                           mean=format(contingency.table[2,2]),
                           var=format(contingency.table[2,2]/sum(contingency.table[2,])*100,digits=3,nsmall=1),
                           p.value=format(fisher.test(contingency.table)$p.value,digits=3,nsmall=2),
                           test="Fisher's test")
            }))
            ret <- rbind(all.ret, cases.ret)
        } else {                
            if (is.factor(variable)) 
                fit <- lm(covariate ~ variable)
            else
                fit <- lm(variable ~ covariate)
            f <- summary(fit)$fstatistic
            p.value <- 1 - pf(f["value"],f["numdf"],f["dendf"])
            
            if (!is.factor(variable) && !is.factor(covariate))
                ret <- data.frame(covariate=name,
                                  value="",
                                  variable="",
                                  mean=format(mean(covariate,na.rm=T)),
                                  var=format(var(covariate,na.rm=T)),
                                  p.value=format(p.value,digits=3,nsmall=2),
                                  test="OLS")
            else {
                if (is.factor(variable)) {
                    var.numeric <- covariate
                    var.factor <- variable
                } else {
                    var.numeric <- variable
                    var.factor <- covariate
                }
                
                all.ret <- data.frame(covariate=name,
                                      value="",
                                      variable="",
                                      mean=format(mean(var.numeric,na.rm=T)),
                                      var=format(sd(var.numeric,na.rm=T)),
                                      p.value=format(p.value,digits=3,nsmall=2),
                                      test="ANOVA")
                cases.ret <- do.call(rbind, lapply(1:nrow(samples), function(i) {
                    if (is.factor(variable)) {
                        var.factor <- variable == samples$variable[i]
                    }
                    else {
                        var.factor <- covariate == samples$covariate[i]
                    }
                    p.value <- t.test(var.numeric[which(var.factor)],
                                      var.numeric[which(!var.factor)])$p.value
                    data.frame(covariate=name,
                               value=samples$covariate[i], ## == "" if covariate is not categorical, otherwise covariate level
                               variable=samples$variable[i],  ## same as above for variable
                               mean=format(mean(var.numeric[which(var.factor)], na.rm=T)),
                               var=format(sd(var.numeric[which(var.factor)], na.rm=T)),
                               p.value=format(p.value,digits=3,nsmall=2),
                               test="t-test")
                }))
                ret <- rbind(all.ret, cases.ret)
            }
        }
        rownames(ret) <- NULL
        ret
    }

    do.call(rbind, lapply(1:ncol(ewas.object$covariates), function(i) {
        summarize.covariate(colnames(ewas.object$covariates)[i],
                            ewas.object$covariates[,i],
                            ewas.object$variable)
    }))
}




