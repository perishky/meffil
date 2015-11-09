#' @export
meffil.ewas.sample.characteristics <- function(ewas.object) {
    summarize.variable <- function(name, variable) {
        if (is.factor(variable)) {
            n <- sapply(levels(variable), function(level) sum(variable == level,na.rm=T))
            data.frame(variable=name,
                       value=names(n),
                       mean=n,
                       var=n/length(variable)*100,row.names=NULL)                   
        }
        else {
            data.frame(variable=name,
                       value=NA,
                       mean=mean(variable,na.rm=T),
                       var=sd(variable,na.rm=T),
                       row.names=NULL)
        }
    }

    rbind(summarize.variable("independent variable", ewas.object$variable),
          do.call(rbind, lapply(1:ncol(ewas.object$covariates), function(i) {
              summarize.variable(colnames(ewas.object$covariates)[i],
                                 ewas.object$covariates[,i])
          })))
}

#' @export
meffil.ewas.covariate.associations <- function(ewas.object) {
    summarize.covariate <- function(name,covariate, variable) {
        samples <- expand.grid(variable=if (!is.factor(variable)) NA else levels(variable),
                               covariate=if (!is.factor(covariate)) NA else levels(covariate),
                               stringsAsFactors=F)
        samples$idx <- lapply(1:nrow(samples), function(i)
                              which((is.na(samples$variable[i])
                                     | variable == samples$variable[i]) &
                                    (is.na(samples$covariate[i])
                                     | covariate == samples$covariate[i])))
    
        if (is.factor(variable) && is.factor(covariate)) {
            contingency.table <- table(variable,covariate)
            all.ret <- data.frame(covariate=name,
                                  value=NA,
                                  variable=NA,
                                  mean=NA,
                                  var=NA,
                                  p.value=fisher.test(contingency.table,
                                      simulate.p.value=T)$p.value,
                                  test="Fisher's test")
            cases.ret <- do.call(rbind, lapply(1:nrow(samples), function(i) {
                contingency.table <- table(variable==samples$variable[i],
                                           covariate==samples$covariate[i])
                data.frame(covariate=name,
                           value=samples$covariate[i],
                           variable=samples$variable[i],
                           mean=contingency.table[2,2],
                           var=contingency.table[2,2]/sum(contingency.table[2,])*100,
                           p.value=fisher.test(contingency.table)$p.value,
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
                                  value=NA,
                                  variable=NA,
                                  mean=mean(covariate,na.rm=T),
                                  var=var(covariate,na.rm=T),
                                  p.value=p.value,
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
                                      value=NA,
                                      variable=NA,
                                      mean=mean(var.numeric,na.rm=T),
                                      var=sd(var.numeric,na.rm=T),
                                      p.value=p.value,
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
                               mean=mean(var.numeric[which(var.factor)], na.rm=T),
                               var=sd(var.numeric[which(var.factor)], na.rm=T),
                               p.value=p.value,
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




