## A slightly modified version of the function
## DoISVA from the ISVA R package.
## https://cran.r-project.org/web/packages/isva/


## require(fastICA)

isva <- function(data.m,mod,pvthCF=0.01,th=0.05,ncomp=NULL,
                   verbose=F){

    stopifnot(is.matrix(data.m))
    isva.o <- isvaFn(data.m,mod,ncomp,verbose=verbose);
    selisv.idx <- 1:ncol(isva.o$isv);
    selisv.m <- matrix(isva.o$isv[,selisv.idx],ncol=length(selisv.idx));
    return(list(isv=selisv.m,nsv=length(selisv.idx),selisv=selisv.idx));
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
    return(list(cor=C,dim=intdim,estdens=estdens.o,thdens=thdens.o,evals=eigen.o$values));
}

isvaFn <- function(data.m,mod,ncomp=NULL,verbose=F){
## isvaFn <- function(data.m,pheno.v,ncomp=NULL,verbose=F){
    ##lm.o <- lm(t(data.m) ~ pheno.v);
    ##res.m <- t(lm.o$res)

    ## BEGIN NEW
    id <- diag(ncol(data.m))
    res.m <- data.m %*% (id - mod %*% solve(t(mod) %*% mod) %*% t(mod))
    ## END NEW
   
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

