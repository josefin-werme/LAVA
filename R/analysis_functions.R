

#' Perform univariate and bivariate analysis, with univariate filtering
#' 
#' Will only perform bivariate test for phenotypes that pass the univariate significance threshold.
#' If more than two phenotypes are specified, the last will be treated as the target phenotype, and the bivariate test will be performed between this phenotype and all others
#' 
#' @param locus Locus object created using the the \code{\link{process.locus}} function. Contains all the relevant parameters and processed sum-stats for the phenotypes of interest
#' @param phenos Optional argument specifying subset of phenotypes to analyse, and/or their order. If NULL, all phenotypes in the locus object will be analysed (in the order of the locus object). NOTE: the last phenotype will be treated as the phenotype of interest
#' @param univ.thresh P-value threshold for the univariate test, used to determine whether phenotypes exhibit sufficient local heritability for a bivariate test
#' @param adap.thresh The thresholds at which to increase the number of permutations for the p-value generation. 
#' Default number of permutations is 1e+4, but will be increased to 1e+5, and 1e+6 as p-values fall below the respective thresholds.
#' If set to NULL, the maximum number of permutations is capped at the default (Note: this significantly speeds up the analysis, but results in poor accuracy for low p-values)
#' @param p.values Set to F to suppress p-values
#' @param CIs Set to F to suppress 95\% confidence intervals
#' @param param.lim The +- threshold at which estimated parameters are considered to be too far out of bounds (and will be set to NA)
#' @param return.unanalysed If true, the bivariate results will contain all phenotypes, with NAs for those that weren't analysed
#' Otherwise, only the analysed phenotypes will be returned (default).
#' 
#' @return List containing the results from the univariate and bivariate tests (see ?run.univ() and ?run.bivar() for more info)
#' 
#' @export
#' 
run.univ.bivar = function(locus, phenos=NULL, univ.thresh=.05, adap.thresh=c(1e-4, 1e-6), p.values=T, CIs=T, param.lim=1.25, return.unanalysed=F) {
	if (is.null(phenos)) { phenos = locus$phenos } else { if (any(! phenos %in% locus$phenos)) { stop(paste("Invalid phenotype ID provided:", paste(phenos[! phenos %in% locus$phenos]))) } }
	P = length(phenos)
	
	# univariate analysis
	univ = run.univ(locus, phenos)
	
	# bivariate analysis
	if (!return.unanalysed) bivar=NULL else bivar = data.frame(phen1 = phenos[-P], phen2 = phenos[P], rho=NA, rho.lower=NA, rho.upper=NA, r2=NA, r2.lower=NA, r2.upper=NA, p=NA)
	if (any(univ$p[-P] < univ.thresh) & univ$p[P] < univ.thresh) { 
		bivar = run.bivar(locus, phenos = phenos[which(univ$p < univ.thresh)], adap.thresh=adap.thresh, p.values=p.values, CIs=CIs, param.lim=param.lim)
		if (return.unanalysed) { bivar = merge(data.frame(phen1 = phenos[-P]), bivar, all=T, sort=F); bivar = bivar[match(phenos[-P], bivar$phen1),]; bivar$phen2=phenos[P] }
	}
	return(list(univ=univ, bivar=bivar))
}


# Univariate p-values
univariate.test = function(locus, phenos=NULL) {
	if (is.null(phenos)) { phenos = locus$phenos } else { if (any(! phenos %in% locus$phenos)) { stop(paste("Invalid phenotype ID provided:", paste(phenos[! phenos %in% locus$phenos]))) } }
	P = length(phenos)
	
	p = rep(NA, P)
	for (i in 1:P) {
		stat = sum(as.matrix(locus$delta[,phenos])[,i]^2) / as.matrix(locus$sigma[phenos,phenos])[i,i]
		p[i] = ifelse(locus$binary[phenos][i], pchisq(stat, locus$K, lower.tail=F), pf(stat/locus$K, locus$K, locus$N[phenos][i] - locus$K-1, lower.tail=F))
	}
	return(p)
}


#' Run univariate analysis
#' 
#' Performs univariate test to determine the presence of local genetic signal (i.e. the local heritability)
#' 
#' @param locus Locus object created using the the \code{\link{process.locus}} function. Contains all the relevant parameters and processed sum-stats for the phenotypes of interest
#' @param phenos Optional argument specifying subset of phenotypes to analyse. If NULL, all phenotypes will be analysed
#' @param var Set to T to return variance estimate
#' 
#' @return Data frame with the columns:
#' \itemize{
#'     \item phen - analysed phenotypes
#'     \item var - local genetic variance
#'     \item h2.obs - observed heritability
#'     \item h2.latent - population heritability (only relevant for binary phenotypes; requires population prevalence to from input info file)
#'     \item p - p-values from the univariate test (F-test for continuous, Chi-sq for binary)
#' }
#' @export
run.univ = function(locus, phenos=NULL, var=F) {
	if (is.null(phenos)) { phenos = locus$phenos } else { if (any(! phenos %in% locus$phenos)) { stop(paste("Invalid phenotype ID provided:", paste(phenos[! phenos %in% locus$phenos]))) } }
	P = length(phenos)
	
	univ = data.frame(phen = phenos)
	if (var) { univ$var = signif(diag(as.matrix(locus$omega[phenos,phenos])), 6) }
	univ$h2.obs = signif(locus$h2.obs[phenos], 6)
	if (any(locus$binary[phenos])) { univ$h2.latent = signif(locus$h2.latent[phenos],6) }
	univ$p = signif(univariate.test(locus, phenos), 6)
	return(univ)
}


#' Bivariate local genetic correlation analysis
#' 
#' Performs bivariate local genetic correlation analysis between two phenotypes for a single locus.
#' If more than one phenotype is specified, the last phenotype will be treated as the phenotype of interest, 
#' and separate bivariate tests will be performed between this phenotype and all others.
#' Note: this test is symmetric, which phenotype is considered predictor/outcome doesn't matter.
#' 
#' @param locus Locus object created using the the \code{\link{process.locus}} function. Contains all the relevant parameters and processed sum-stats for the phenotypes of interest
#' @param phenos Optional argument specifying subset of phenotypes to analyse, and/or their order. If NULL, all phenotypes in the locus object will be analysed (in the order of the locus object). NOTE: the last phenotype will be treated as the phenotype of interest
#' @param adap.thresh The thresholds at which to increase the number of permutations for the p-value generation. 
#' Default number of permutations is 1e+4, but will be increased to 1e+5, and 1e+6 as p-values fall below the respective thresholds.
#' If set to NULL, the maximum number of permutations is capped at the default (Note: this significantly speeds up the analysis, but results in poor accuracy for low p-values)
#' @param p.values Set to F to suppress p-values
#' @param CIs Set to F to suppress 95\% confidence intervals
#' @param param.lim The +- threshold at which estimated parameters are considered to be too far out of bounds (and will be set to NA)
#' 
#' @return Data frame with the columns:
#' \itemize{
#'     \item phen1 / phen2 - analysed phenotypes
#'     \item rho - standardised coefficient for the local genetic correlation
#'     \item rho.lower / rho.upper - 95\% confidence intervals for rho
#'     \item r2 - proportion of variance in genetic signal for phen1 explained by phen2 (and vice versa)
#'     \item r2.lower / r2.upper - 95\% confidence intervals for the r2
#'     \item p - permutation p-values for the local genetic correlation
#' }
#' @export
run.bivar = function(locus, phenos=NULL, adap.thresh=c(1e-4, 1e-6), p.values=T, CIs=T, param.lim=1.25) {
	if (is.null(phenos)) { phenos = locus$phenos } else { if (any(! phenos %in% locus$phenos)) { stop(paste("Invalid phenotype ID provided:", paste(phenos[! phenos %in% locus$phenos]))) } } ## 20-09-24: added error if faulty phenotype IDs are provided
	P = length(phenos); Px = P-1; Y = phenos[P]
	if (P < 2) { stop("Less than 2 phenotypes provided, cannot perform bivariate analysis") }
	
	bivar = list(); params = c("gamma.std","r2"); ci.params = c("rho.lower","rho.upper","r2.lower","r2.upper")
	bivar = data.frame(matrix(NA, Px, length(params)+7)); colnames(bivar) = c("phen1","phen2",params,ci.params,"p"); bivar$phen2 = phenos[P]
	
	for (i in 1:Px) {
		bivar$phen1[i] = phenos[i]
		
		# estimate params
		mom = estimate.moments(delta=locus$delta[,c(phenos[i],Y)], sigma=locus$sigma[c(phenos[i],Y),c(phenos[i],Y)], K=locus$K)
		for (p in params) { bivar[[p]][i] = signif(mom[[p]], 6) } # store params
		
		# confidence intervals
		if (CIs) {
			ci = ci.bivariate(K = locus$K, omega = mom$omega, sigma = locus$sigma[c(phenos[i],Y),c(phenos[i],Y)])
			bivar$rho.lower[i] = ci$ci.rho.low; bivar$rho.upper[i] = ci$ci.rho.high
			bivar$r2.lower[i] = ci$ci.r2.low; bivar$r2.upper[i] = ci$ci.r2.high
		}
		# p-values
		if (p.values) { bivar$p[i] = signif(integral.p(bivariate.integral, K=locus$K, omega=mom$omega, sigma=locus$sigma[c(phenos[i],Y),c(phenos[i],Y)], adap.thresh=adap.thresh), 6) }	## 20/09/18: updated integral.p()
	}
	# filter any estimates that are too far out of bounds
	bivar = filter.params(data = bivar, params = c(params, ci.params, "p"), param.lim = param.lim)	# first param must be gamma/rho
																									# params just needs to list all that will be set to NA if rho / gamma is too far out of bounds
	# cap out of bounds values (CI's are already capped)
	for (p in params) { bivar[[p]] = cap(bivar[[p]], lim=c(ifelse(p=="r2", 0, -1), 1)) }	# capping rhos at -1/1, and r2s at 0/1
	colnames(bivar)[which(colnames(bivar)=="gamma.std")] = "rho"
	
	return.vars = c("phen1","phen2","rho","rho.lower","rho.upper","r2","r2.lower","r2.upper","p")
	return(bivar[,return.vars])
}





#' Multiple local genetic regression analysis
#' 
#' @param locus Locus object created using the the \code{\link{process.locus}} function. Contains all the relevant parameters and processed sum-stats for the phenotypes of interest
#' @param phenos Optional argument specifying subset of phenotypes to analyse, and/or their order. If NULL, all phenotypes in the locus object will be analysed (in the order of the locus object). NOTE: the last phenotype will be treated as the outcome
#' @param adap.thresh The thresholds at which to increase the number of permutations for the p-value generation. 
#' Default number of permutations is 1e+4, but will be increased to 1e+5, and 1e+6 as p-values fall below the respective thresholds.
#' If set to NULL, the maximum number of permutations is capped at the default (Note: this significantly speeds up the analysis, but results in poor accuracy for low p-values)
#' @param p.values Set to F to suppress p-values
#' @param CIs Set to F to suppress 95\% confidence intervals
#' @param param.lim The +- threshold at which estimated parameters are considered to be too far out of bounds (and will be set to NA)
#' 
#' #' @return Data frame with the columns:
#' \itemize{
#'     \item pretictors / outcome - analysed phenotypes
#'     \item gamma - standardised multiple regression coefficient
#'     \item gamma.lower / gamma.upper - 95\% confidence intervals for gamma
#'     \item r2 - proportion of variance in genetic signal for the outcome explained by all predictors simultaneously
#'     \item r2.lower / r2.upper - 95\% confidence intervals for the r2
#'     \item p - permutation p-values for the gammas
#' }
#' @export
run.multireg = function(locus, phenos=NULL, adap.thresh=c(1e-4, 1e-6), only.full.model=F, p.values=T, CIs=T, param.lim=1.5, suppress.message=F) {
	if (is.null(phenos)) { phenos = locus$phenos } else { if (any(! phenos %in% locus$phenos)) { stop(paste("Invalid phenotype ID provided:", paste(phenos[! phenos %in% locus$phenos]))) } } ## 20-09-24: added error if faulty phenotype IDs are provided
	P = length(phenos); Px = P-1; Y = phenos[P]
	if (P < 3) { stop(paste0("Only ",P," phenotypes provided for conditional analysis; Need at least 3")) }
	cond.idx = 1:Px
	
	if (!suppress.message) print(paste0("~ Running multiple regression for outcome '",Y,"', with predictors '",paste(phenos[1:Px],collapse="', '"),"'"))
	
	cond = list()
	models = list(); for (i in 2:length(cond.idx)) { models[[i-1]] = combn(phenos[1:Px], i) }		# get all unique predictor models
	if (only.full.model) { m=list(); m[[1]] = models[[length(models)]]; models = m }
	
	params = c("gamma.std","r2")
	ci.params = c("gamma.lower","gamma.upper","r2.lower","r2.upper")
	
	for (j in 1:length(models)) {
		cond[[j]] = list()
		for (k in 1:ncol(models[[j]])) {
			mod.idx = models[[j]][,k]		# get index of current model (actually by name, rather than 'index')
			Pm = length(mod.idx)			# no. predictors in curr model
			
			# check invertibility of omega.x
			eig = eigen(locus$omega[mod.idx,mod.idx]); if (any(eig$values / (sum(eig$values)/Pm) < 1e-4)) { stop("omega.X not invertible") }; rm(eig)
			
			# mom
			mom = estimate.moments(delta = locus$delta[,c(mod.idx,Y)], sigma = locus$sigma[c(mod.idx,Y),c(mod.idx,Y)], K = locus$K)
			
			# store mom output
			cond[[j]][[k]] = data.frame(matrix(NA, nrow = Pm, ncol = length(params)+7)); colnames(cond[[j]][[k]]) = c("predictors", "outcome", params, ci.params, "p")
			cond[[j]][[k]]$predictors = mod.idx; cond[[j]][[k]]$outcome = Y
			cond[[j]][[k]]$gamma.std = signif(mom$gamma.std, 6); cond[[j]][[k]]$r2 = rep(signif(mom$r2, 6), Pm)
			
			# 95% confidence intervals
			if (CIs) {
				ci = ci.multivariate(K=locus$K, omega=mom$omega, sigma = locus$sigma[c(mod.idx,Y),c(mod.idx,Y)])	# TODO: why subsetting sigma here but not omega?
				cond[[j]][[k]]$gamma.lower = ci$gamma$ci.low; cond[[j]][[k]]$gamma.upper = ci$gamma$ci.high
				cond[[j]][[k]]$r2.lower = ci$r2$ci.low; cond[[j]][[k]]$r2.upper = ci$r2$ci.high
			}
			# p-values
			if (p.values) {
				pvals = try(integral.p(multivariate.integral, K = locus$K, omega = mom$omega, sigma = locus$sigma[c(mod.idx,Y),c(mod.idx,Y)], adap.thresh=adap.thresh), silent=T) 	# TODO: why subsetting sigma here but not omega?
				cond[[j]][[k]]$p = signif(pvals, 6)
			}
			# filter estimates that are too far out of bounds 
			cond[[j]][[k]] = filter.params(data = cond[[j]][[k]], params = c(params, ci.params, "p"), param.lim = param.lim, multreg=T) # first param must be gamma/rho
			
			# cap r2 to 0 / 1 (CI's are already capped)
			cond[[j]][[k]]$r2 = cap(cond[[j]][[k]]$r2, lim=c(0,1))
			
			cond[[j]][[k]] = cond[[j]][[k]][,c("predictors","outcome","gamma.std","gamma.lower","gamma.upper","r2","r2.lower","r2.upper","p")]
			colnames(cond[[j]][[k]])[which(colnames(cond[[j]][[k]])=="gamma.std")] = "gamma"
		}
	}
	return(cond)
}





#' Partial correlation analysis
#' 
#' Will perform the partial correlation between the first two phenotypes (phen1, phen2) conditioned on the rest (Z). 
#' Phenotype order is based on that within the locus object by default, but can be changed by passing a phenotype vector to the 'phenos' argument.
#' 
#' @param locus Locus object created using the the \code{\link{process.locus}} function. Contains all the relevant parameters and processed sum-stats for the phenotypes of interest
#' @param phenos Optional argument specifying subset of phenotypes to analyse, and/or their order. If NULL, all phenotypes in the locus object will be analysed (in the order of the locus object). NOTE: the first two are p1 & p2 (the phenotypes of interest), the remaining are Z
#' @param adapt.thresh The thresholds at which to increase the number of permutations for the p-value generation. 
#' Default number of permutations is 1e+4, but will be increased to 1e+5, and 1e+6 as p-values fall below the respective thresholds.
#' If set to NULL, the maximum number of permutations is capped at the default (Note: this significantly speeds up the analysis, but results in poor accuracy for low p-values)
#' @param p.values Set to F to suppress p-values
#' @param CIs Set to F to suppress 95\% confidence intervals
#' @param max.r2 Max r2 threshold for the regression of phen1 ~ Z and phen2 ~ Z. If any of these r2's are too high, the partial correlation becomes unstable, and analysis is therefore aborted.
#' 
#' #' @return Data frame with the columns:
#' \itemize{
#'     \item phen1 / phen2 - analysed phenotypes
#'     \item rho - standardised coefficient for the local genetic correlation
#'     \item rho.lower / rho.upper - 95\% confidence intervals for rho
#'     \item r2 - proportion of variance in genetic signal for phen1 explained by phen2 (and vice versa)
#'     \item r2.lower / r2.upper - 95\% confidence intervals for the r2
#'     \item p - permutation p-values for the local genetic correlation
#' }
#' @export
run.partial.cor = function(locus, phenos=NULL, adap.thresh=c(1e-4, 1e-6), p.values=T, CIs=T, max.r2=.95, param.lim=1.25) {
	if (is.null(phenos)) { phenos = locus$phenos } else { if (any(! phenos %in% locus$phenos)) { print(paste("Error: invalid phenotype ID provided:", paste(phenos[! phenos %in% locus$phenos]))); stop() } }
	P = length(phenos); if (P < 3) { stop(paste0("Only ",P," phenotypes provided for partial correlation; need at least 3")) }
	x = 1; y = 2; z = 3:P
	
	print(paste0("~ Running partial correlation for '",phenos[x],"' and '",phenos[y],"', conditioned on '",paste(phenos[z],collapse="' + '"),"'"))
	
	# check invertibility of omega.z
	eig = eigen(locus$omega[phenos,phenos][z,z]); if (any(eig$values / (sum(eig$values)/length(z)) < 1e-4)) { stop("omega.Z not invertible") }; rm(eig)
	
	# get r2 of p1 (x) and p2 (y) on Z
	r2 = data.frame(x = NA, y = NA)
	for (i in 1:2) {
		if (P == 3) {
			r2[i] = run.bivar(locus, phenos[c(3:P,i)], p.values=F, CIs=F)$r2
		} else {
			r2[i] = run.multireg(locus, phenos[c(3:P,i)], only.full.model=T, p.values=F, CIs=F, suppress.message=T)[[1]][[1]]$r2[1]
		}
	}
	out = data.frame(phen1=phenos[x], phen2=phenos[y], z=paste(phenos[z],collapse=";"), r2.phen1_z = r2$x, r2.phen2_z = r2$y, pcor = NA, ci.lower = NA, ci.upper=NA, p = NA)
	out$pcor = partial.cor(locus$omega[phenos,phenos], x, y, z)
	
	if (CIs) {
		ci = ci.pcor(K=locus$K, xy.index=c(x,y), z.index=z, omega=locus$omega[phenos,phenos], sigma=locus$sigma[phenos,phenos])
		out$ci.lower = ci[2]; out$ci.upper = ci[3]
	}
	
	# if r2 is < max.r2, proceed with pvalues
	if (all(out[c("r2.phen1_z","r2.phen2_z")] < max.r2) & p.values) {
		out$p = integral.p(pcov.integral, K=locus$K, omega=locus$omega[phenos,phenos], sigma=locus$sigma[phenos,phenos], adap.thresh=adap.thresh)
	}
	
	# filter any values that are too far out of bounds
	params = c("pcor","ci.lower","ci.upper")
	out = filter.params(data = out, params = params, param.lim = param.lim) # params just needs to list all params that will be set to NA if rho / gamma is too far out of bounds; first param must be gamma/rho
	out$pcor = cap(out$pcor) # cap the pcor at -1/1 (CIs are capped already)
	
	out[,! colnames(out) %in% c("phen1","phen2","z")] = as.data.frame(lapply(out[,! colnames(out) %in% c("phen1","phen2","z")], signif, 6))
	return(out)
}


# remove estimates that are excessively out of bounds
filter.params = function(data, params, param.lim=1.25, multreg=F) {	
	out.of.bounds = abs(data[params[1]]) > abs(param.lim)			# assumes first param is gamma/rho
	
	if (any(out.of.bounds)) {
		# get phenotype colnames
		if (!multreg) { phen.cols = c("phen1","phen2") } else { phen.cols = c("predictors","outcome") }
		# print warning
		message("Warning: Estimates too far out of bounds (+-",param.lim,") for phenotype(s) ",
				paste0(paste0(data[phen.cols[1]][out.of.bounds]," ~ ",data[phen.cols[2]][out.of.bounds]," (",signif(data[params[1]][out.of.bounds],4),")"),collapse=", "),
				" in locus ",locus$id,". Values will be set to NA. To change this threshold, modify the 'param.lim' argument")  # actually dont need to pass locus
		# set to NA
		if (!multreg) { data[out.of.bounds, params] = NA } else { data[,params] = NA }
	}
	return(data)
}

# function for capping values
cap = function(values, lim = c(-1,1)) {
	for (i in 1:length(values)) {
		if (is.na(values[i])) next
		if (values[i] > max(lim)) {
			values[i] = max(lim)
		} else if (values[i] < min(lim)) {
			values[i] = min(lim)
		}
	}
	return(values)
}
