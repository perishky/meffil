#' svaPrep
#'
#' prepare data for sva and/or isva analysis
#' 
#' @param variable Independent variable vector.
#' @param covariates Covariates data frame to include in regression model,
#' one row per sample, one column per covariate (Default: NULL).
#' @param beta.sva Methylation levels matrix, one row per CpG site, one column per sample.
#' @param featureset Name from \code{\link{meffil.list.featuresets}()}  (Default: NA).
#' one row per sample, one column per covariate (Default: NULL).
#' @param most.variable Apply (Independent) Surrogate Variable Analysis to the 
#' given most variable CpG sites (Default: 50000).
svaPrep <- function(variable,covariates,beta.sva,featureset,most.variable){
	autosomal.sites <- meffil.get.autosomal.sites(featureset)
	autosomal.sites <- intersect(autosomal.sites, rownames(beta.sva))
	if (length(autosomal.sites) < most.variable) {
		warning("Probes from the sex chromosomes will be used to calculate surrogate variables.")
	} else {
		beta.sva <- beta.sva[autosomal.sites,]
	}
	var.idx <- order(rowVars(beta.sva, na.rm=T), decreasing=T)[1:most.variable]
	beta.sva <- impute.matrix(beta.sva[var.idx,,drop=F])
	
	if (!is.null(covariates)) {
		cov.frame <- model.frame(~., data.frame(covariates, stringsAsFactors=F), na.action=na.pass)
		mod0 <- model.matrix(~., cov.frame)
	}
	else
		mod0 <- matrix(1, ncol=1, nrow=length(variable))
	mod <- cbind(mod0, variable)
	res <- list("mod" = mod, "beta.sva"=beta.sva)
	return(res)
}

#' svaIsva 
#' 
#' Constructs function calls to surrogate variable analysis or independent surrogate variable analysis
#' 
#' @param beta.sva Methylation levels matrix, one row per CpG site, one column per sample.
#' @param covariates Covariates data frame to include in regression model,
#' one row per sample, one column per covariate (Default: NULL).
#' @param random.seed Value with which to seed the pseudo random number generator for reproducible results
#' @param mod model matrix form \code{\link{svaPrep}}
#' @param n.sv Number of surrogate variables 
#' @param mode Do sva or isva analysis, accepts values "sv" or "isv"
#' @param verbose Set verbosity of isva function
svaIsva <- function(beta.sva,covariates,random.seed,mod,n.sv,mode="sv",verbose){
	set.seed(random.seed)
	if(mode=="sv"){
		sva.ret <- sva(beta.sva, mod=mod, mod0=mod[,1,drop=F], n.sv=n.sv)
	}
	else if (mode=="isv"){
		sva.ret <- isva(beta.sva, mod, ncomp=n.sv, verbose=verbose)
	}
	if (!is.null(covariates))
		svaRes <- data.frame(covariates, sva.ret[[mode]], stringsAsFactors=F)
	else
		svaRes <- as.data.frame(sva.ret[[mode]])
	return(svaRes)
}

#' outliers
#' 
#' selects outliers based on the value of "outlier.iqr.factor"
#' @param beta Methylation levels matrix, one row per CpG site, one column per sample.
#' @param outlier.iqr.factor For each CpG site, prior to fitting regression models,
#' set methylation levels less than
#' \code{Q1 - outlier.iqr.factor * IQR} or more than
#' \code{Q3 + outlier.iqr.factor * IQR} to NA. Here IQR is the inter-quartile
#' range of the methylation levels at the CpG site, i.e. Q3-Q1.
#' Set to NA to skip this step (Default: NA)
outliers <- function(beta,outlier.iqr.factor) {
	q <- rowQuantiles(beta, probs = c(0.25, 0.75), na.rm = T)
	iqr <- q[,2] - q[,1]
	too.hi <- which(beta > q[,2] + outlier.iqr.factor * iqr, arr.ind=T)
	too.lo <- which(beta < q[,1] - outlier.iqr.factor * iqr, arr.ind=T)
	if (nrow(too.hi) > 0) beta[too.hi] <- NA
	if (nrow(too.lo) > 0) beta[too.lo] <- NA
	res <- list("too.lo"=too.lo,"too.hi"=too.hi)
	return(res)
}


#' Epigenome-wide association study
#'
#' Test association with each CpG site.
#'
#' @param beta Methylation levels matrix, one row per CpG site, one column per sample.
#' @param variable Independent variable vector. - NB only handles discrete 2 catagory variables
#' @param covariates Covariates data frame to include in regression model,
#' one row per sample, one column per covariate (Default: NULL).
#' @param batch Batch vector to be included as a random effect (Default: NULL).
#' @param weights Non-negative observation weights.
#' Can be a numeric matrix of individual weights of same dimension as \code{beta},
#' or a numeric vector of weights with length \code{ncol(beta)},
#' or a numeric vector of weights with length \code{nrow(beta)}. 
#' @param cell.counts Proportion of cell counts for one cell type in cases
#' where the samples are mainly composed of two cell types (e.g. saliva) (Default: NULL).
#' @param isva Apply Independent Surrogate Variable Analysis (ISVA) to the
#' methylation levels and include the resulting variables as covariates in a
#' regression model (Default: TRUE).  
#' @param sva Apply Surrogate Variable Analysis (SVA) to the
#' methylation levels and covariates and include
#' the resulting variables as covariates in a regression model (Default: TRUE).
#' @param n.sv Number of surrogate variables to calculate (Default: NULL).
#' @param winsorize.pct Apply all regression models to methylation levels
#' winsorized to the given level. Set to NA to avoid winsorizing (Default: 0.05).
#' @param robust Test associations with the 'robust' option when \code{\link{limma::eBayes}}
#' is called (Default: TRUE).
#' @param rlm Test assocaitions with the 'robust' option when \code{\link{limma:lmFit}}
#' is called (Default: FALSE).
#' @param outlier.iqr.factor For each CpG site, prior to fitting regression models,
#' set methylation levels less than
#' \code{Q1 - outlier.iqr.factor * IQR} or more than
#' \code{Q3 + outlier.iqr.factor * IQR} to NA. Here IQR is the inter-quartile
#' range of the methylation levels at the CpG site, i.e. Q3-Q1.
#' Set to NA to skip this step (Default: NA).
#' @param most.variable Apply (Independent) Surrogate Variable Analysis to the 
#' given most variable CpG sites (Default: 50000).
#' @param featureset Name from \code{\link{meffil.list.featuresets}()}  (Default: NA).
#' @param random.seed Value with which to seed the pseudo random number generator for reproducible results
#' @param verbose Set to TRUE if status updates to be printed (Default: FALSE).
#'
#' @export
meffil.ewas <- function(beta, variable,
						covariates=NULL, batch=NULL, weights=NULL,
						cell.counts=NULL,
						isva=T, sva=T, ## cate?
						n.sv=NULL,
						isva0=F,isva1=F, ## deprecated
						winsorize.pct=0.05,
						robust=TRUE,
						rlm=FALSE,
						outlier.iqr.factor=NA, ## typical value = 3
						most.variable=min(nrow(beta), 50000),
						featureset=NA,
						random.seed=20161123,
						verbose=F) {
	## depreciated Args
	if (isva0 || isva1)
		stop("isva0 and isva1 are deprecated and superceded by isva and sva")
	# warn and set isva, sva values rather than kill?
	
	if (is.na(featureset))
		featureset <- guess.featureset(rownames(beta))
	features <- meffil.get.features(featureset)
	
	## data Validation / type checks
	stopifnot(length(rownames(beta)) > 0 && all(rownames(beta) %in% features$name))
	stopifnot(ncol(beta) == length(variable))
	stopifnot(is.null(covariates) || is.data.frame(covariates) && nrow(covariates) == ncol(beta))
	stopifnot(is.null(batch) || length(batch) == ncol(beta))
	stopifnot(is.null(weights)
			|| is.numeric(weights) && (is.matrix(weights) && nrow(weights) == nrow(beta) && ncol(weights) == ncol(beta)
										|| is.vector(weights) && length(weights) == nrow(beta)
										|| is.vector(weights) && length(weights) == ncol(beta)))
	stopifnot(most.variable > 1 && most.variable <= nrow(beta))
	stopifnot(!is.numeric(winsorize.pct) || winsorize.pct > 0 && winsorize.pct < 0.5)
	
	## cache var and covars
	original.variable <- variable
	original.covariates <- covariates
	if (is.character(variable))
		variable <- as.factor(variable)
	
	stopifnot(!is.factor(variable) || is.ordered(variable) || length(levels(variable)) == 2)
	
	msg("Simplifying any categorical variables.", verbose=verbose)
	variable <- simplify.variable(variable)
	if (!is.null(covariates))
		covariates <- do.call(cbind, lapply(covariates, simplify.variable))
	
	## Getting samples without missing values
	sample.idx <- which(!is.na(variable))
	if (!is.null(covariates))
		sample.idx <- intersect(sample.idx, which(apply(!is.na(covariates), 1, all)))
	
	msg("Removing", ncol(beta) - length(sample.idx), "missing case(s).", verbose=verbose)
	
	if (is.matrix(weights))
		weights <- weights[,sample.idx]
	if (is.vector(weights) && length(weights) == ncol(beta))
		weights <- weights[sample.idx]
	
	beta <- beta[,sample.idx]
	variable <- variable[sample.idx]
	
	if (!is.null(covariates))
		covariates <- covariates[sample.idx,,drop=F]
	
	if (!is.null(batch))
		batch <- batch[sample.idx]
	
	if (!is.null(cell.counts))
		cell.counts <- cell.counts[sample.idx]
	
	if (!is.null(covariates)) {
		pos.var.idx <- which(apply(covariates, 2, var, na.rm=T) > 0)
		msg("Removing", ncol(covariates) - length(pos.var.idx), "covariates with no variance.",
			verbose=verbose)
		covariates <- covariates[,pos.var.idx, drop=F]
	}
	
	## ewas setup
	### Declace list to hold different covariate sets ewas results
	covariate.sets <- list(none=NULL) 
	if (!is.null(covariates))
		covariate.sets$all <- covariates
	
	### winsorisation
	if (is.numeric(winsorize.pct))  {
		msg(winsorize.pct, "- winsorizing the beta matrix.", verbose=verbose)
		beta <- winsorize(beta, pct=winsorize.pct)
	}
	
	### outliers
	too.hi <- too.lo <- NULL
	if (is.numeric(outlier.iqr.factor)) {
		outliersRes <- outliers(beta,outlier.iqr.factor)
		too.hi <- outliersRes$too.hi
		too.lo <- outliersRes$too.lo
	}
	
	## sva/isva
	if (isva || sva) {
		beta.sva <- beta
		res <- svaPrep(variable,covariates,beta.sva,featureset,most.variable)
		mod <- res$mod
		beta.sva <- res$beta.sva
		
		if (isva) {
			msg("ISVA.", verbose=verbose)
			covariate.sets$isva <- svaIsva(beta.sva,covariates,random.seed,mod,n.sv,mode="isv",verbose)
			cat("\n")
		}
		
		if (sva) {
			msg("SVA.", verbose=verbose)
			covariate.sets$sva <- svaIsva(beta.sva,covariates,random.seed,mod,n.sv,mode="sv",verbose)
			cat("\n")
		}
	}
	
	## run EWASs
	analyses <- sapply(names(covariate.sets), function(name) {
		msg("EWAS for covariate set", name, verbose=verbose)
		covariates <- covariate.sets[[name]]
		ewas(variable,
			 beta=beta,
			 covariates=covariates,
			 batch=batch,
			 weights=weights,
			 cell.counts=cell.counts,
			 winsorize.pct=winsorize.pct, 
			 robust=robust, 
			 rlm=rlm)
	}, simplify=F)
	
	## extract results
	p.values <- sapply(analyses, function(analysis) analysis$table$p.value)
	coefficients <- sapply(analyses, function(analysis) analysis$table$coefficient)
	rownames(p.values) <- rownames(coefficients) <- rownames(analyses[[1]]$table)
	
	for (name in names(analyses)) {
		idx <- match(rownames(analyses[[name]]$table), features$name)
		analyses[[name]]$table$chromosome <- features$chromosome[idx]
		analyses[[name]]$table$position <- features$position[idx]
	}
	
	## Return results
	list(class="ewas",
		 version=packageVersion("meffil"),
		 samples=sample.idx,
		 variable=original.variable[sample.idx],
		 covariates=original.covariates[sample.idx,,drop=F],
		 winsorize.pct=winsorize.pct,
		 robust=robust,
		 rlm=rlm,
		 outlier.iqr.factor=outlier.iqr.factor,
		 most.variable=most.variable,
		 p.value=p.values,
		 coefficient=coefficients,
		 analyses=analyses,
		 random.seed=random.seed,
		 too.hi=too.hi,
		 too.lo=too.lo)
}

#' is.ewas.object
#' 
#' @param object a 'meffil' object (not an actual S3/4 R object)
is.ewas.object <- function(object)
	is.list(object) && "class" %in% names(object) && object$class == "ewas"

#' ewas
#' 
#' Test associations between \code{variable} and each row of \code{beta}
#' while adjusting for \code{covariates} (fixed effects) and \code{batch} (random effect).
#' NB only handles discrete 2 catagory variables
#' If \code{cell.counts} is not \code{NULL}, then it is assumed that
#' the methylation data is derived from samples with two cell types.
#' \code{cell.counts} should then be a vector of numbers
#' between 0 and 1 of length equal to \code{variable} corresponding
#' to the proportions of cell of a selected cell type in each sample.
#' The regression model is then modified in order to identify
#' associations specifically in the selected cell type (PMID: 24000956).
ewas <- function(variable, beta, covariates=NULL, batch=NULL, weights=NULL, cell.counts=NULL, winsorize.pct=0.05,
				 robust=TRUE, rlm=FALSE, verbose=F) {
	stopifnot(all(!is.na(variable)))
	stopifnot(length(variable) == ncol(beta))
	stopifnot(is.null(covariates) || nrow(covariates) == ncol(beta))
	stopifnot(is.null(batch) || length(batch) == ncol(beta))
	stopifnot(is.null(cell.counts)
			  || length(cell.counts) == ncol(beta)
			  && all(cell.counts >= 0 & cell.counts <= 1))
  
	method <- "ls"
	if (rlm) method <- "robust"
  
	if (is.null(covariates))
		design <- data.frame(intercept=1, variable=variable)
	else
		design <- data.frame(intercept=1, variable=variable, covariates)
	rownames(design) <- colnames(beta)

	if (!is.null(cell.counts)) {
		## Irizarray method: Measuring cell-type specific differential methylation
		##	  in human brain tissue
		## Mi is methylation level for sample i
		## Xi is the value of the variable of interest for sample i
		## pi is the proportion of the target cell type for sample i
		## thus: Mi ~ target cell methylation * pi + other cell methylation * (1-pi)
		## and target cell methylation = base target methylation + effect of Xi value
		## and other cell methylation = base other methylation + effect of Xi value
		## thus:
		## Mi = (T + U Xi)pi + (O + P Xi)(1-pi) + e
		##	= C + A Xi pi + B Xi (1-pi) + e		
		design <- design[,-which(colnames(design) == "intercept")]
		
		design <- cbind(design * cell.counts,
						typeB=design * (1-cell.counts))
	}

	msg("Linear regression with limma::lmFit", verbose=verbose)
	fit <- NULL
	batch.cor <- NULL
	if (!is.null(batch)) {
		msg("Adjusting for batch effect", verbose=verbose)
		corfit <- duplicateCorrelation(beta, design, block=batch, ndups=1)
		batch.cor <- corfit$consensus
		msg("Linear regression with batch as random effect", verbose=verbose)
		tryCatch(fit <- lmFit(beta, design, method=method, block=batch, cor=batch.cor, weights=weights),
				 error=function(e) {
					 print(e)
					 msg("lmFit failed with random effect batch variable, omitting", verbose=verbose)
				 })
	}
	if (is.null(fit)) {
		msg("Linear regression with only fixed effects", verbose=verbose)
		batch <- NULL
		fit <- lmFit(beta, design, method=method, weights=weights)
	}
	
	msg("Empirical Bayes", verbose=verbose)
	if (is.numeric(winsorize.pct) && robust) {
	  fit.ebayes <- eBayes(fit, robust=T, winsor.tail.p=c(winsorize.pct, winsorize.pct))
	} else {
	  fit.ebayes <- eBayes(fit, robust=robust)
	}

	alpha <- 0.975
	std.error <- (sqrt(fit.ebayes$s2.post) * fit.ebayes$stdev.unscaled[,"variable"])
	margin.error <- (std.error * qt(alpha, df=fit.ebayes$df.total))
	n <- rowSums(!is.na(beta))

	list(design=design,
		 batch=batch,
		 batch.cor=batch.cor,
		 cell.counts=cell.counts,
		 table=data.frame(p.value=fit.ebayes$p.value[,"variable"],
			 fdr=p.adjust(fit.ebayes$p.value[,"variable"], "fdr"),
			 p.holm=p.adjust(fit.ebayes$p.value[,"variable"], "holm"),
			 t.statistic=fit.ebayes$t[,"variable"],
			 coefficient=fit.ebayes$coefficient[,"variable"],
			 coefficient.ci.high=fit.ebayes$coefficient[,"variable"] + margin.error,
			 coefficient.ci.low=fit.ebayes$coefficient[,"variable"] - margin.error,
			 coefficient.se=std.error,
			 n=n))
}


#' Epigenome-wide association study (lm)
#'
#' Test association with each CpG site.
#'
#' @param beta Methylation levels matrix, one row per CpG site, one column per sample.
#' @param variable Independent variable vector.
#' @param covariates Covariates data frame to include in regression model,
#' one row per sample, one column per covariate (Default: NULL).
#' @param weights numeric vector of weights to provide to the lm()/glm() function
#' @param family to perform generalised linear modeling with r's glm() function
#' @param isva Apply Independent Surrogate Variable Analysis (ISVA) to the
#' methylation levels and include the resulting variables as covariates in a
#' regression model (Default: TRUE).  
#' @param sva Apply Surrogate Variable Analysis (SVA) to the
#' methylation levels and covariates and include
#' the resulting variables as covariates in a regression model (Default: TRUE).
#' @param n.sv Number of surrogate variables to calculate (Default: NULL).
#' @param winsorize.pct Apply all regression models to methylation levels
#' winsorized to the given level. Set to NA to avoid winsorizing (Default: 0.05).
#' @param outlier.iqr.factor For each CpG site, prior to fitting regression models,
#' set methylation levels less than
#' \code{Q1 - outlier.iqr.factor * IQR} or more than
#' \code{Q3 + outlier.iqr.factor * IQR} to NA. Here IQR is the inter-quartile
#' range of the methylation levels at the CpG site, i.e. Q3-Q1.
#' Set to NA to skip this step (Default: NA).
#' @param most.variable Apply (Independent) Surrogate Variable Analysis to the 
#' given most variable CpG sites (Default: 50000).
#' @param featureset Name from \code{\link{meffil.list.featuresets}()}  (Default: NA).
#' @param random.seed Value with which to seed the pseudo random number generator for reproducible results
#' @param verbose Set to TRUE if status updates to be printed (Default: FALSE).
#'
#' @export
meffil.ewas.lm <- function(beta, variable,
							covariates=NULL, 
							#batch=NULL, 
							weights=NULL,
							family=NULL,
							#cell.counts=NULL,
							isva=T, sva=T, ## cate?
							n.sv=NULL,
							isva0=F,isva1=F, ## deprecated
							winsorize.pct=0.05,
							outlier.iqr.factor=NA, ## typical value = 3
							most.variable=min(nrow(beta), 50000),
							featureset=NA,
							random.seed=20161123,
							verbose=F) {
	## depreciated Args
	if (isva0 || isva1)
		stop("isva0 and isva1 are deprecated and superceded by isva and sva")
	# warn and set isva, sva values rather than kill?
	
	if (is.na(featureset))
		featureset <- guess.featureset(rownames(beta))
	features <- meffil.get.features(featureset)
	
	phenotype.data <- cbind(variable,covariates)
	
	varName <- colnames(variable)
	variable <- unlist(variable,use.names = F)
	
	## data Validation / type checks
	stopifnot(length(rownames(beta)) > 0 && all(rownames(beta) %in% features$name))
	stopifnot(ncol(beta) == length(variable))
	stopifnot(is.null(covariates) || is.data.frame(covariates) && nrow(covariates) == ncol(beta))
	#stopifnot(is.null(batch) || length(batch) == ncol(beta))
	stopifnot(is.null(weights)|| is.numeric(weights))
	stopifnot(most.variable > 1 && most.variable <= nrow(beta))
	stopifnot(!is.numeric(winsorize.pct) || winsorize.pct > 0 && winsorize.pct < 0.5)
	
	## cache var and covars
	original.variable <- variable
	original.covariates <- covariates

	msg("Simplifying any categorical variables.", verbose=verbose)
	variable <- simplify.variable(variable)
	if (!is.null(covariates))
		covariates <- do.call(cbind, lapply(covariates, simplify.variable))
	
	## Getting samples without missing values
	# sample.idx = all samples without missing value?s
	sample.idx <- which(!is.na(variable))
	
	##!! removing missing values from NB phenotype.data !!##
	##!! get rid ot phenotype.data and merge var/covar to data frame and pass to ewas.lm !!##
	phenotype.data[sample.idx,,drop=F]
	
	if (!is.null(covariates))
		sample.idx <- intersect(sample.idx, which(apply(!is.na(covariates), 1, all)))
	
	msg("Removing", ncol(beta) - length(sample.idx), "missing case(s).", verbose=verbose)
	
	if (is.vector(weights) && length(weights) == ncol(beta))
		weights <- weights[sample.idx]
	
	beta <- beta[,sample.idx]
	variable <- variable[sample.idx]
	
	if (!is.null(covariates))
		covariates <- covariates[sample.idx,,drop=F]
	
	# if (!is.null(batch))
	# 	batch <- batch[sample.idx]
	# 
	# if (!is.null(cell.counts))
	# 	cell.counts <- cell.counts[sample.idx]
	
	if (!is.null(covariates)) {
		pos.var.idx <- which(apply(covariates, 2, var, na.rm=T) > 0)
		msg("Removing", ncol(covariates) - length(pos.var.idx), "covariates with no variance.",
			verbose=verbose)
		covariates <- covariates[,pos.var.idx, drop=F]
	}
	
	## ewas setup
	### Declace list to hold different covariate sets ewas results
	covariate.sets <- list(none=NULL) 
	if (!is.null(covariates))
		covariate.sets$all <- covariates
	
	### winsorisation
	if (is.numeric(winsorize.pct))  {
		msg(winsorize.pct, "- winsorizing the beta matrix.", verbose=verbose)
		beta <- winsorize(beta, pct=winsorize.pct)
	}
	
	### outliers
	too.hi <- too.lo <- NULL
	if (is.numeric(outlier.iqr.factor)) {
		outliersRes <- outliers(beta,outlier.iqr.factor)
		too.hi <- outliersRes$too.hi
		too.lo <- outliersRes$too.lo
	}
	
	## sva/isva
	if (isva || sva) {
		beta.sva <- beta
		res <- svaPrep(variable,covariates,beta.sva,featureset,most.variable)
		mod <- res$mod
		beta.sva <- res$beta.sva
		
		if (isva) {
			msg("ISVA.", verbose=verbose)
			covariate.sets$isva <- svaIsva(beta.sva,covariates,random.seed,mod,n.sv,mode="isv",verbose)
			cat("\n")
		}
		
		if (sva) {
			msg("SVA.", verbose=verbose)
			covariate.sets$sva <- svaIsva(beta.sva,covariates,random.seed,mod,n.sv,mode="sv",verbose)
			cat("\n")
		}
	}
	
	## run EWASs
	analyses <- sapply(names(covariate.sets), function(name) {
		msg("EWAS for covariate set", name, verbose=verbose)
		covariates <- covariate.sets[[name]]
		if(!is.null(covariate.sets[[name]])){
			phenotype.data <- cbind(phenotype.data,covariate.sets[[name]]) ##!!
			#phenotype.data <- cbind(variable,covariate.sets[[name]]) ##!! ??
			#colnames(phenotype.data) <- c(varName,colnames(covariates))
		}
		ewas.lm(variable=varName,
				beta=beta,
				covariates=colnames(covariates),
				#batch=batch,
				weights=weights,
				#cell.counts=cell.counts,
				winsorize.pct=winsorize.pct,
				phenotype.data = phenotype.data,
				family=family
		)
	}, simplify=F)
	
	## extract results
	p.values <- sapply(analyses, function(analysis) analysis$table$p.value)
	coefficients <- sapply(analyses, function(analysis) analysis$table$coefficient)
	rownames(p.values) <- rownames(coefficients) <- rownames(analyses[[1]]$table)
	
	for (name in names(analyses)) {
		idx <- match(rownames(analyses[[name]]$table), features$name)
		analyses[[name]]$table$chromosome <- features$chromosome[idx]
		analyses[[name]]$table$position <- features$position[idx]
	}
	
	list(class="ewas",
		version=packageVersion("meffil"),
		samples=sample.idx,
		variable=original.variable[sample.idx],
		covariates=original.covariates[sample.idx,,drop=F],
		winsorize.pct=winsorize.pct,
		#robust=robust,
		#rlm=rlm,
		outlier.iqr.factor=outlier.iqr.factor,
		most.variable=most.variable,
		p.value=p.values,
		coefficient=coefficients,
		analyses=analyses,
		random.seed=random.seed,
		too.hi=too.hi,
		too.lo=too.lo)
}

#' ewas.lm
#' 
#' performs EWAS using r's lm()/glm() functions
#' 
#' Test associations between \code{variable} and each row of \code{beta}
#' while adjusting for \code{covariates} 
ewas.lm <- function(variable, beta, covariates=NULL,
					weights=NULL, winsorize.pct=0.05,
					verbose=F,phenotype.data=NULL,family=NULL) {
	## Data Validation
	stopifnot(is.character(variable))
	#stopifnot(is.character(covariates))
	stopifnot(is.null(phenotype.data) || nrow(phenotype.data) == ncol(beta))
	n <- rowSums(!is.na(beta))
	
	## Set design
	#(different if no covars)
	design <- NULL
	designFx <- NULL
	if(is.null(covariates)) {
		design <- data.frame(intercept=1, variable=phenotype.data[,variable])
		designFx <- function(variable,covariates) {
			paste0(variable,"~meth.matrix[cpg,]")
		}
	}
	else {
		design <- data.frame(intercept=1, variable=phenotype.data[,variable], covariates=phenotype.data[,covariates])
		designFx <- function(variable,covariates) {
		paste0(variable,"~meth.matrix[cpg,]+",paste0(covariates,collapse="+"))
		}
	}
	rownames(design) <- colnames(beta)
	
	## Set modeling type - fit glm if family provided
	modFx <- NULL
	if(is.null(family)) {
		modFx <- function(FORM,data,family,weights){
			lm(FORM,data,weights)
		}
	}
	else {
		modFx <- function(FORM,data,family,weights){
		glm(FORM,data,family = family,weights=weights)
		}
	}
	
	## Model building and parsing function
	model <- function(cpg, meth.matrix, variable, covariates, phenotype.data,family,weights) {
		FORM=formula(designFx(variable,covariates))
		mod= modFx(FORM,data=phenotype.data,family = family,weights=weights)
		res = summary(mod)$coef[2,c("Estimate","Std. Error","Pr(>|t|)")]
	}
	
	## Run models (in parallel)
	ewas_res <- mclapply(setNames(seq_len(nrow(beta)),dimnames(beta)[[1]]),
						model,
						meth.matrix=beta,
						variable=variable,
						covariates=covariates,
						phenotype.data=phenotype.data,
						family=family,
						weights=weights)
	
	ewas_res <- data.frame(t(sapply(ewas_res,c)))
	
	colnames(ewas_res)<-c("effect.size","std.err","P.value")
	
	## Return results
	list(design=design,
		#batch=batch,
		#batch.cor=batch.cor,
		#cell.counts=cell.counts,
		table=data.frame(p.value=ewas_res$P.value,
						coefficient=ewas_res$effect.size,
						coefficient.se=ewas_res$std.err,
						n=n))
}

