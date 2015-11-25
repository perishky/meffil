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
                           value=samples$covariate[i],
                           variable=samples$variable[i],
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
                        level <- samples$variable[i]
                    }
                    else {
                        var.factor <- covariate == samples$covariate[i]
                        level <- samples$covariate[i]
                    }
                    p.value <- t.test(var.numeric[which(var.factor)],
                                      var.numeric[which(!var.factor)])$p.value
                    data.frame(covariate=name,
                               value=samples$covariate[i],
                               variable=samples$variable[i],                           
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




