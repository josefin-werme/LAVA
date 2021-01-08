

#' Perform both univariate and bivariate tests, with filtering based on univariate signal
#' 
#' Will only perform bivariate test for all phenotypes that pass the univariate significance threshold. By default, the bivariate test will be performed for all combinations of phenotypes in the locus, 
#' but this can be modified using the 'phenos' and 'target' arguments (see below)
#' 
#' @param locus Locus object created using the the \code{\link{process.locus}} function. Contains all the relevant parameters and processed sum-stats for the phenotypes of interest
#' @param phenos Subset of phenotypes to analyse, and/or their order. If NULL, all phenotypes in the locus object will be analysed (in the order listed within the locus object)
#' @param target Target phenotype of interest. If NULL, bivariate correlations between all pairs of phenotypes will be computed; 
#' Otherwise, this will only be done between the target phenotype and all other phenotypes.
#' @param univ.thresh P-value threshold for the univariate test, used to determine whether phenotypes exhibit sufficient local heritability for a bivariate test
#' @param adap.thresh The thresholds at which to increase the number of iterations for the p-value generation. 
#' Default number of iterations is 1e+4, but will be increased to 1e+5, and 1e+6 as p-values fall below the respective thresholds.
#' If set to NULL, the maximum number of iterations is capped at the default (Note: this significantly speeds up the analysis, but results in poor accuracy for low p-values)
#' @param p.values Set to F to suppress p-values
#' @param CIs Set to F to suppress 95\% confidence intervals
#' @param param.lim The +- threshold at which estimated parameters are considered to be too far out of bounds. If the estimated parameter exceeds this threshold, it is considered unreliable and will be set to NA. 
#' 
#' @return List containing the results from the univariate and bivariate tests (see ?run.univ() and ?run.bivar() for more info)
#' 
#' @export
#' 
run.univ.bivar = function(locus, phenos=NULL, target=NULL, univ.thresh=.05, adap.thresh=c(1e-4, 1e-6), p.values=T, CIs=T, param.lim=1.25) {
	if (is.null(phenos)) { phenos = locus$phenos } else { if (any(! phenos %in% locus$phenos)) { print(paste0("Error: Invalid phenotype ID(s) provided: '",paste0(phenos[! phenos %in% locus$phenos], collapse="', '"),"'")); return(NA) }}
	if (!is.null(target)) {
		if (length(target) > 1) { print(paste0("Error: More than one target phenotype specified")); return(NA) }
		if (! target %in% locus$phenos) { print(paste0("Error: Invalid target phenotype specified: '", target,"'")); return(NA) }
		if (! target %in% phenos) { phenos = c(phenos,target) }		# append target to phenos if not already present
	}
	P = length(phenos)
	
	# univariate analysis
	univ = run.univ(locus, phenos)
	
	# bivariate analysis
	bivar = NA
	if (sum(univ$p < univ.thresh) > 1) {
		if (!is.null(target)) { if (subset(univ, phen==target)$p > univ.thresh) { break() } }	# if target is specified, only proceed if target is sig
		bivar = run.bivar(locus, phenos = subset(univ, p < univ.thresh)$phen, target=target, adap.thresh=adap.thresh, p.values=p.values, CIs=CIs, param.lim=param.lim)
	}
	return(list(univ=univ, bivar=bivar))
}


# Univariate p-values
univariate.test = function(locus, phenos=NULL) {
	if (is.null(phenos)) { phenos = locus$phenos } else { if (any(! phenos %in% locus$phenos)) { print(paste0("Error: Invalid phenotype ID(s) provided: '",paste0(phenos[! phenos %in% locus$phenos], collapse="', '"),"'")); return(NA) }}
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
#' Tests the univariate local genetic signal (i.e. the local heritability) for all phenotypes.
#' 
#' @param locus Locus object created using the the \code{\link{process.locus}} function. Contains all the relevant parameters and processed sum-stats for the phenotypes of interest
#' @param phenos Subset of phenotypes to analyse, and/or their order. If NULL, all phenotypes in the locus object will be analysed (in the order listed within the locus object)
#' @param var Set to T to return variance estimate
#' 
#' @return Data frame with the columns:
#' \itemize{
#'     \item phen - analysed phenotypes
#'     \item var - local genetic variance
#'     \item h2.obs - observed local heritability
#'     \item h2.latent - estimated local population heritability (only relevant for binary phenotypes; requires population prevalence to be specified in input info file)
#'     \item p - p-values from the univariate test (F-test for continuous, Chi-sq for binary)
#' }
#' @export
run.univ = function(locus, phenos=NULL, var=F) {
	if (is.null(phenos)) { phenos = locus$phenos } else { if (any(! phenos %in% locus$phenos)) { print(paste0("Error: Invalid phenotype ID(s) provided: '",paste0(phenos[! phenos %in% locus$phenos], collapse="', '"),"'")); return(NA) } }
	P = length(phenos)
	
	univ = data.frame(phen = phenos)
	if (var) { univ$var = signif(diag(as.matrix(locus$omega[phenos,phenos])), 6) }
	univ$h2.obs = signif(locus$h2.obs[phenos], 6)
	if (any(locus$binary[phenos]) & all(!is.na(locus$h2.latent))) { univ$h2.latent = signif(locus$h2.latent[phenos],6) }
	univ$p = signif(univariate.test(locus, phenos), 6)
	return(univ)
}




#' Bivariate local genetic correlation analysis
#' 
#' Performs bivariate local genetic correlation analysis between two phenotypes.
#' By default, the bivariate test will be performed for all combinations of phenotypes in the locus, 
#' but this can be modified using the 'phenos' and 'target' arguments (see below)
#' 
#' @param locus Locus object created using the the \code{\link{process.locus}} function. Contains all the relevant parameters and processed sum-stats for the phenotypes of interest
#' @param phenos Subset of phenotypes to analyse, and/or their order. If NULL, all phenotypes in the locus object will be analysed (in the order listed within the locus object)
#' @param target Target phenotype of interest. If NULL, bivariate correlations between all pairs of phenotypes will be computed; 
#' Otherwise, this will only be done between the target phenotype and all other phenotypes.
#' @param adap.thresh The thresholds at which to increase the number of iterations for the p-value generation. 
#' Default number of iterations is 1e+4, but will be increased to 1e+5, and 1e+6 as p-values fall below the respective thresholds.
#' If set to NULL, the maximum number of iterations is capped at the default (Note: this significantly speeds up the analysis, but results in poor accuracy for low p-values)
#' @param p.values Set to F to suppress p-values
#' @param CIs Set to F to suppress 95\% confidence intervals
#' @param param.lim The +- threshold at which estimated parameters are considered to be too far out of bounds. If the estimated parameter exceeds this threshold, it is considered unreliable and will be set to NA. 
#' 
#' @return Data frame with the columns:
#' \itemize{
#'     \item phen1 / phen2 - analysed phenotypes
#'     \item rho - standardised coefficient for the local genetic correlation
#'     \item rho.lower / rho.upper - 95\% confidence intervals for rho
#'     \item r2 - proportion of variance in genetic signal for phen1 explained by phen2 (and vice versa)
#'     \item r2.lower / r2.upper - 95\% confidence intervals for the r2
#'     \item p - simulation p-values for the local genetic correlation
#' }
#' @export
run.bivar = function(locus, phenos=NULL, target=NULL, adap.thresh=c(1e-4, 1e-6), p.values=T, CIs=T, param.lim=1.25) {
	if (is.null(phenos)) { phenos = locus$phenos } else { if (any(! phenos %in% locus$phenos)) { print(paste0("Error: Invalid phenotype ID(s) provided: '",paste0(phenos[! phenos %in% locus$phenos], collapse="', '"),"'")); return(NA) }}
	if (!is.null(target)) { 
		if (length(target) > 1) { print(paste0("Error: More than one target phenotype specified")); return(NA) }
		if (! target %in% locus$phenos) { print(paste0("Error: Invalid target phenotype specified: '", target,"'")); return(NA) }
		if (! target %in% phenos) { phenos = c(phenos,target) }		# append target to phenos if not already present
	}
	P = length(phenos)
	if (P < 2) { print("Less than 2 phenotypes provided, cannot perform bivariate analysis"); return(NA) }
	
	if (is.null(target)) {
		pairs = t(combn(phenos,2))	# all unique phenotype pairs
	} else {
		pairs = cbind(phenos[!phenos%in%target], target)
	}
	bivar = list(); params = c("gamma.std","r2"); ci.params = c("rho.lower","rho.upper","r2.lower","r2.upper")
	bivar = data.frame(matrix(NA, nrow(pairs), length(params)+7)); colnames(bivar) = c("phen1","phen2",params,ci.params,"p"); bivar[,c("phen1","phen2")] = pairs
	
	for (i in 1:nrow(pairs)) {
		mom = estimate.moments(delta = locus$delta[,pairs[i,]], sigma = locus$sigma[pairs[i,],pairs[i,]], K = locus$K)  # estimate params
		for (p in params) { bivar[[p]][i] = signif(mom[[p]], 6) } # store params
		
		# confidence intervals
		if (CIs) {
			ci = ci.bivariate(K = locus$K, omega = locus$omega[pairs[i,],pairs[i,]], sigma = locus$sigma[pairs[i,],pairs[i,]])
			bivar$rho.lower[i] = ci$ci.rho.low; bivar$rho.upper[i] = ci$ci.rho.high
			bivar$r2.lower[i] = ci$ci.r2.low; bivar$r2.upper[i] = ci$ci.r2.high
		}
		# p-values
		if (p.values) { bivar$p[i] = signif(integral.p(bivariate.integral, K = locus$K, omega = locus$omega[pairs[i,],pairs[i,]], sigma = locus$sigma[pairs[i,],pairs[i,]], adap.thresh=adap.thresh), 6) }
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





#' Local genetic multiple regression analysis
#' 
#' Will perform a local genetic multiple regression analysis, which models the genetic signal for a single outcome phenotype using two or more predictor phenotypes.
#' Here, the genetic correlations between all predictors will be accounted for, and their genetic relation with the outcome will be conditioned on one another.
#' 
#' @param locus Locus object created using the the \code{\link{process.locus}} function. Contains all the relevant parameters and processed sum-stats for the phenotypes of interest
#' @param phenos Subset of phenotypes to analyse, and/or their order. If NULL, all phenotypes in the locus object will be analysed (in the order listed within the locus object). NOTE: the last phenotype will be treated as the outcome
#' @param adap.thresh The thresholds at which to increase the number of iterations for the p-value generation. 
#' Default number of iterations is 1e+4, but will be increased to 1e+5, and 1e+6 as p-values fall below the respective thresholds.
#' If set to NULL, the maximum number of iterations is capped at the default (Note: this significantly speeds up the analysis, but results in poor accuracy for low p-values)
#' @param p.values Set to F to suppress p-values
#' @param CIs Set to F to suppress 95\% confidence intervals
#' @param param.lim The +- threshold at which estimated parameters are considered to be too far out of bounds. If the estimated parameter exceeds this threshold, it is considered unreliable and will be set to NA. 
#' 
#' @return Data frame with the columns:
#' \itemize{
#'     \item pretictors / outcome - analysed phenotypes
#'     \item gamma - standardised multiple regression coefficient
#'     \item gamma.lower / gamma.upper - 95\% confidence intervals for gamma
#'     \item r2 - proportion of variance in genetic signal for the outcome phenotype explained by all predictor phenotypes simultaneously
#'     \item r2.lower / r2.upper - 95\% confidence intervals for the r2
#'     \item p - simulation p-values for the gammas
#' }
#' @export
run.multireg = function(locus, phenos=NULL, adap.thresh=c(1e-4, 1e-6), only.full.model=F, p.values=T, CIs=T, param.lim=1.5, suppress.message=F) {
	if (is.null(phenos)) { phenos = locus$phenos } else { if (any(! phenos %in% locus$phenos)) { print(paste0("Error: Invalid phenotype ID(s) provided: '",paste0(phenos[! phenos %in% locus$phenos], collapse="', '"),"'")); return(NA) } }
	P = length(phenos); Px = P-1; Y = phenos[P]
	if (P < 3) { print(paste0("Only ",P," phenotypes provided for conditional analysis; Need at least 3")); return(NA) }
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
				ci = ci.multivariate(K = locus$K, omega = locus$omega[c(mod.idx,Y),c(mod.idx,Y)], sigma = locus$sigma[c(mod.idx,Y),c(mod.idx,Y)])
				cond[[j]][[k]]$gamma.lower = ci$gamma$ci.low; cond[[j]][[k]]$gamma.upper = ci$gamma$ci.high
				cond[[j]][[k]]$r2.lower = ci$r2$ci.low; cond[[j]][[k]]$r2.upper = ci$r2$ci.high
			}
			# p-values
			if (p.values) {
				pvals = try(integral.p(multivariate.integral, K = locus$K, omega = locus$omega[c(mod.idx,Y),c(mod.idx,Y)], sigma = locus$sigma[c(mod.idx,Y),c(mod.idx,Y)], adap.thresh=adap.thresh), silent=T)
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




#' Local partial genetic correlation analysis
#' 
#' Will perform a local partial genetic correlation between the first two phenotypes (phen1, phen2) conditioned on the rest (Z). 
#' Phenotype order is based on that within the locus object by default, but can be changed by passing a phenotype vector with the desired order to the 'phenos' argument.
#' 
#' @param locus Locus object created using the the \code{\link{process.locus}} function. Contains all the relevant parameters and processed sum-stats for the phenotypes of interest
#' @param phenos Subset of phenotypes to analyse, and/or their order. If NULL, all phenotypes in the locus object will be analysed (in the order listed within the locus object). NOTE: the first two will be treated as p1 & p2 (the target phenotypes), the remaining as Z (those that are conditioned on)
#' @param adapt.thresh The thresholds at which to increase the number of iterations for the p-value generation. 
#' Default number of iterations is 1e+4, but will be increased to 1e+5, and 1e+6 as p-values fall below the respective thresholds.
#' If set to NULL, the maximum number of iterations is capped at the default (Note: this significantly speeds up the analysis, but results in poor accuracy for low p-values)
#' @param p.values Set to F to suppress p-values
#' @param CIs Set to F to suppress 95\% confidence intervals
#' @param max.r2 Max r2 threshold for the regression of phen1~Z and phen2~Z. If any of these r2's are too high, the partial correlation becomes unstable, and analysis is therefore aborted.
#' @param param.lim The +- threshold at which estimated parameters are considered to be too far out of bounds. If the estimated parameter exceeds this threshold, it is considered unreliable and will be set to NA. 
#' 
#' @return Data frame with the columns:
#' \itemize{
#'     \item phen1 / phen2 - target phenotypes
#'     \item z - phenotype(s) which the genetic correlation between the target phenotypes were conditioned on
#'     \item r2.phen1_z / r2.phen2_z - the proportion of genetic signal in target phenotypes explained by z. Note: if either of these exceed the threshold specified by the max.r2 argument, the analysis will be aborted (to prevent unstable estimates)
#'     \item pcor - the partial correlation between phen1 and phen2 conditioned on z
#'     \item ci.lower / ci.upper - 95\% confidence intervals for the partial genetic correlation
#'     \item p - simulation p-values for the partial genetic correlation
#' }
#' @export
run.partial.cor = function(locus, phenos=NULL, adap.thresh=c(1e-4, 1e-6), p.values=T, CIs=T, max.r2=.95, param.lim=1.25) {
	if (is.null(phenos)) { phenos = locus$phenos } else { if (any(! phenos %in% locus$phenos)) { print(paste0("Error: Invalid phenotype ID(s) provided: '",paste0(phenos[! phenos %in% locus$phenos], collapse="', '"),"'")); return(NA) } }
	P = length(phenos); if (P < 3) { print(paste0("Only ",P," phenotypes provided for partial correlation; need at least 3")); return(NA) }
	x = 1; y = 2; z = 3:P
	
	print(paste0("~ Running partial correlation for '",phenos[x],"' and '",phenos[y],"', conditioned on '",paste(phenos[z],collapse="' + '"),"'"))
	
	# check invertibility of omega.z
	eig = eigen(locus$omega[phenos,phenos][z,z]); if (any(eig$values / (sum(eig$values)/length(z)) < 1e-4)) { stop("omega.Z not invertible") }; rm(eig)
	
	# get r2 of p1 (x) and p2 (y) on Z
	r2 = data.frame(x = NA, y = NA)
	for (i in 1:2) {
		if (P == 3) {
			r2[i] = run.bivar(locus, phenos[c(z,i)], p.values=F, CIs=F)$r2
		} else {
			r2[i] = run.multireg(locus, phenos[c(z,i)], only.full.model=T, p.values=F, CIs=F, suppress.message=T)[[1]][[1]]$r2[1]
		}
	}
	out = data.frame(phen1=phenos[x], phen2=phenos[y], z=paste(phenos[z],collapse=";"), r2.phen1_z = r2$x, r2.phen2_z = r2$y, pcor = NA, ci.lower = NA, ci.upper=NA, p = NA)
	out$pcor = partial.cor(locus$omega[phenos,phenos], x, y, z)
	
	if (CIs) {
		ci = ci.pcor(K=locus$K, xy.index=c(x,y), z.index=z, omega=locus$omega[phenos,phenos], sigma=locus$sigma[phenos,phenos])
		out$ci.lower = ci[2]; out$ci.upper = ci[3]
	}
	
	# if r2 is < max.r2, proceed with p-values
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
				" in locus ",locus$id,". Values will be set to NA. To change this threshold, modify the 'param.lim' argument")  # dont need to pass locus since its an env
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
