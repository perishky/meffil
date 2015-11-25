## A slightly stripped down version of the function
## DoISVA from the ISVA R package.
## https://cran.r-project.org/web/packages/isva/


## require(fastICA)

DoISVA <- function(data.m,pheno.v,cf.m=NULL,factor.log,pvthCF=0.01,th=0.05,ncomp=NULL,
                   verbose=F){
### BEGIN NEW
    stopifnot(is.matrix(data.m)) 
### END NEW
    
    ## Main ISVA function
    isva.o <- isvaFn(data.m,pheno.v,ncomp,verbose=verbose);
    
    if(is.null(cf.m)==FALSE){
        
        ## study pattern of correlation of ISVA components to POI and CFs
        tmp.m <- cbind(pheno.v,cf.m);
        treatfactor <- c(FALSE,factor.log);
                pv.m <- matrix(nrow=ncol(isva.o$isv),ncol=1+ncol(cf.m));
        colnames(pv.m) <- c("POI",colnames(cf.m)); ## POI:phenotype of interest
        for(c in 1:ncol(tmp.m)){

### BEGIN NEW
            if (length(unique(na.omit(tmp.m[,c]))) <= 1) { 
                warning(paste("cf.m[,", c-1, "] is being skipped, doesn't vary.", sep="")) 
                pv.m[,c] <- 1
                next 
            } 
### END NEW
            
            if(treatfactor[c]==FALSE){
                for(sv in 1:ncol(isva.o$isv)){
                    lm.o <- lm(isva.o$isv[,sv] ~ as.numeric(tmp.m[,c]));
                    pv.m[sv,c] <- summary(lm.o)$coeff[2,4];
                }
            }
            else {
                for(sv in 1:ncol(isva.o$isv)){
                    lm.o <- lm(isva.o$isv[,sv] ~ as.factor(tmp.m[,c]));
                    pv.m[sv,c] <- pf(summary(lm.o)$fstat[1],summary(lm.o)$fstat[2],summary(lm.o)$fstat[3],lower.tail=FALSE);
                }
            }
        }

        ## selection of ISVs
        msg("Selecting ISVs", verbose=verbose)
        selisv.idx <- vector();
        for(sv in 1:nrow(pv.m)){
            
            ncf <- length(which(pv.m[sv,2:ncol(pv.m)]< pvthCF)) ## pvth=0.01
                          minpv <- min(pv.m[sv,2:ncol(pv.m)]);
            phpv <- pv.m[sv,1];
            if(ncf > 0){
                if(minpv < phpv){
                    selisv.idx <- c(selisv.idx,sv);
                }
            }
        }
        if (length(selisv.idx) == 0) {
            warning("No ISVs selected because none correlated with the given confounders. Rerun ISVA with cf.m=NULL option")
### BEGIN NEW
            selisv.idx <- 1:ncol(isva.o$isv) 
            pv.m <- NULL
            ## stop
### END NEW
        }        
    }
    else { ### confounder matrix not given, so select all ISVs
        selisv.idx <- 1:ncol(isva.o$isv);
        pv.m <- NULL;
    }
    ## print("Running final multivariate regressions with selected ISVs");
    selisv.m <- matrix(isva.o$isv[,selisv.idx],ncol=length(selisv.idx));
    ## mod <- model.matrix( ~ pheno.v + selisv.m);
    ## modNULL <- model.matrix( ~ selisv.m);
    
    ## df1 <- dim(mod)[2]
    ## df0 <- dim(modNULL)[2]
    ## pv.v <- rep(0, nrow(data.m));
    ## Id <- diag(ncol(data.m))
    ## resid <- data.m %*% (Id - mod %*% solve(t(mod) %*% mod) %*%
    ##                      t(mod))
    ## rss1 <- rowSums(resid * resid)
    ## rm(resid)
    ## residNULL <- data.m %*% (Id - modNULL %*% solve(t(modNULL) %*% modNULL) %*%
    ##                          t(modNULL))
    ## rssNULL <- rowSums(residNULL * residNULL)
    ## rm(residNULL)
    ## fstats <- ((rssNULL - rss1)/(df1 - df0))/(rss1/(ncol(data.m) - df1))
    ## pv.v <- 1 - pf(fstats, df1 = (df1 - df0), df2 = (ncol(data.m) - df1))    
    ## pv.s <- sort(pv.v,decreasing=FALSE,index.return=TRUE);
    ## qv.v <- qvalue(pv.s$x)$qvalue;
    ## ntop <- length(which(qv.v < th));
    ## print(paste("Number of DEGs after ISV adjustment = ",ntop,sep=""));
    ## if(ntop>0){
    ##   pred.idx <- pv.s$ix[1:ntop];
    ## ### find t-stats of significant ones
    ##   lm.o <- lm( t(data.m[pred.idx,]) ~ pheno.v + selisv.m );
    ##   tstats.v <- unlist(lapply(summary(lm.o),function(x){ x$coeff[2,3];}));
    ##   lm.m <- cbind(tstats.v,pv.s$x[1:ntop],qv.v[1:ntop]);
    ##   colnames(lm.m) <- c("t-stat","P-value","q-value");
    ## }
    ## else {
    ##   pred.idx <- NULL;
    ##   lm.m <- NULL;
    ## }    
    ## return(list(spv=pv.s$x,qv=qv.v,rk=pv.s$ix,ndeg=ntop,deg=pred.idx,lm=lm.m,isv=selisv.m,nsv=length(selisv.idx),pvCF=pv.m,selisv=selisv.idx));
    return(list(isv=selisv.m,nsv=length(selisv.idx),pvCF=pv.m,selisv=selisv.idx));
}


EstDimRMT <- function(data.m,plot=TRUE,verbose=F){
    ## standardise matrix
    M <- apply(data.m,2,function(X){ (X - mean(X))/sqrt(var(X))});
    
    sigma2 <- var(as.vector(M));
    Q <- nrow(data.m)/ncol(data.m);
    ns <- ncol(data.m);
    lambdaMAX <- sigma2*(1+1/Q + 2*sqrt(1/Q));
    lambdaMIN <- sigma2*(1+1/Q - 2*sqrt(1/Q));
    delta <- lambdaMAX - lambdaMIN;#  print(delta);
    
    roundN <- 3;
    step <- round(delta/ns,roundN);
    while(step==0){
        roundN <- roundN+1;
        step <- round(delta/ns,roundN);
    }
    
    
    lambda.v <- seq(lambdaMIN,lambdaMAX,by=step);
    dens.v <- vector();
    ii <- 1;
    for(i in lambda.v){
        dens.v[ii] <- (Q/(2*pi*sigma2))*sqrt( (lambdaMAX-i)*(i-lambdaMIN) )/i;
        ii <- ii+1;
    }
    ## theoretical density
    thdens.o <- list(min=lambdaMIN,max=lambdaMAX,step=step,lambda=lambda.v,dens=dens.v);
    C <- 1/nrow(M) * t(M) %*% M;
    eigen.o <- eigen(C,symmetric=TRUE);
    ## empirical density
    estdens.o <- density(eigen.o$values,from=min(eigen.o$values),to=max(eigen.o$values),cut=0);
    intdim <- length(which(eigen.o$values > thdens.o$max));
    ## evalues.v <- eigen.o$values;
    ## plot 
    ## if(plot){
    ##  minx <- min(min(thdens.o$lambda),min(evalues.v));
    ##  maxx <- max(max(thdens.o$lambda),max(evalues.v));
    ##  miny <- min(min(thdens.o$dens),min(estdens.o$y));
    ##  maxy <- max(max(thdens.o$dens),max(estdens.o$y));
    ##  pdf("RMTplot.pdf",width=4,height=4);
    ##  plot(thdens.o$lambda,thdens.o$dens,xlim=c(0.5,maxx),ylim=c(miny,maxy),type="b",col="green",xlab="Folded Eigenvalues",ylab="density",lwd=1.25);
    ##  i <- min(which(estdens.o$x > min(evalues.v)));
    ##  f <- max(which(estdens.o$x < max(evalues.v)));
    ##  points(x=estdens.o$x[i:f],y=estdens.o$y[i:f],type="b",col="red",cex=0.5);
    ##  for(i in 1:intdim){
    ##   abline(v=evalues.v[i],col="red",lwd=2);
    ##  }
    ##  dev.off();
    ## }
 
    return(list(cor=C,dim=intdim,estdens=estdens.o,thdens=thdens.o,evals=eigen.o$values));
}

isvaFn <- function(data.m,pheno.v,ncomp=NULL,verbose=F){
    lm.o <- lm(t(data.m) ~ pheno.v);
    res.m <- t(lm.o$res);
    ## model <- model.matrix(~1+pheno.v); ## 'model' never used!
    if(is.null(ncomp)){
        rmt.o <-  EstDimRMT(res.m, verbose=verbose)
        ncomp <- rmt.o$dim;
        msg(paste("Number of candidate ISVs = ",ncomp,sep=""), verbose=verbose);
    }
    else {
        msg("no need to estimate dimensionality", verbose=verbose);
    }

    ## perform ICA on residual matrix
    fICA.o <- fastICA(res.m,n.comp=ncomp);
    
    ## now construct ISV
    tmp.m <- t(fICA.o$A);
    isv.m <- tmp.m;
    sd <- 1/sqrt(ncol(data.m)-3);
    for(k in 1:ncol(tmp.m)){
        cor.v <- as.vector(cor(t(data.m),tmp.m[, k]))
        z.v <- 0.5*log((1+cor.v)/(1-cor.v));
        pv.v <- 2*pnorm(abs(z.v),0,sd,lower.tail=FALSE)
        tmp.s <- sort(pv.v,decreasing=FALSE,index.return=TRUE);
        ##qv.o <- qvalue(pv.v);
        ##nsig <- length(which(qv.o$qvalues<0.05));
### BEGIN NEW
        qv.o <- p.adjust(pv.v, "fdr")
        nsig <- length(which(qv.o < 0.05))
### END NEW
        if( nsig < 500 ){
            nsig <- 500;
        }
        red.m <- data.m[tmp.s$ix[1:nsig],];
        fICA.o <- fastICA(red.m,n.comp=ncomp);
        cor.v <- abs(cor(tmp.m[,k],t(fICA.o$A)));
        kmax <- which.max(cor.v);
        isv.m[,k] <- t(fICA.o$A)[,kmax]; 
        msg(paste("Built ISV ",k,sep=""), verbose=verbose); 
    }
    return(list(n.isv=ncomp,isv=isv.m));
}

