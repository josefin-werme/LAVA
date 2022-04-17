#' Process locus info
#' 
#' Processes the summary statistics for all phenotypes within the specified locus and computes all parameters necessary for analysis.
#' 
#' @param loc Locus info for a single locus, obtained using the \code{\link{read.loci}} function. Expects a locus ID ('LOC') together with locus coordinates ('CHR', 'START', 'STOP') and/or a ';' separated SNP list ('SNPS')
#' @param input Input object created with the \code{\link{process.input}} function, containing relevant summary statistics and related info (e.g. sample overlap, case/control ratio)
#' @param phenos Subset of phenotypes from the input object to process. If NULL (default), all phenotypes will be processed
#' @param min.K Minimum number of PCs required to process locus (cannot be less than two). If this criterion is not met, the function will fail and the locus cannot be analysed.
#' @param prune.thresh PC pruning threshold governing the maximum number of PCs to retain.
#' @param max.prop.K Upper bound on retained number of PCs as proportion of lowest sample size of input data.
#' @param drop.failed Determines if failed phenotypes are removed from the output object (default) or retained.
#' PCs are selected as such that the cumulative proportion of variance explained is at least that of the threshold (set to 99 percent by default).
#' 
#' @return Returns an environment containing general locus info, the processed locus related sumstats, and parameters required for analysis. If the function fails (e.g. due to too few SNPs), it will return NULL. 
#' If processing fails for specific phenotypes, only the successful phenotypes will be returned (unless the drop.failed argument is set to true).
#' 
#' \itemize{
#'     \item id - locus ID
#'     \item chr / start / stop - locus coordinates
#'     \item snps - list of locus SNPs
#'     \item n.snps - number of SNPs within locus
#'     \item K - number of PCs retained
#'     \item delta - PC projected joint SNP effects
#'     \item sigma - sampling covariance matrix
#'     \item omega - genetic covariance matrix
#'     \item omega.cor - genetic correlation matrix
#'     \item N - vector of average N across locus SNPs for each phenotype
#'     \item phenos - phenotype IDs
#'     \item binary - boolean vector indicating whether phenotypes are binary
#'     \item h2.obs - observed local heritability
#'     \item h2.latent - estimated local population heritability (only relevant for binary phenotypes; requires population prevalence to be specified in the input info file)
#'     \item failed - boolean vector indicating whether phenotypes failed during processing (only present if drop.failed=F)
#' }
#' 
#' @export


process.locus = function(locus, input, phenos=NULL, min.K=2, prune.thresh=99, max.prop.K=0.75, drop.failed=T) {
	if (is.data.frame(locus) && nrow(locus)!=1) { print("Error: Locus info provided for incorrect number of loci. Please provide only a single locus at a time"); loc=NULL; return(NULL) }  # modified cdl 18/3
	if (!(all(c("LOC","CHR","START","STOP") %in% names(locus)) | all(c("LOC","SNPS") %in% names(locus)))) { print("Error: Locus info data frame is missing some or all of the required headers ('LOC' + 'CHR','START','STOP' and/or 'SNPS')"); loc=NULL; return(NULL) }
	
	min.K = max(min.K, 2) # just to make sure it isn't 1, because that leads to some matrix multiplication dimension issues
	
	# define locus environment & add locus info
	loc = new.env(parent=globalenv())
	loc$id = locus$LOC; loc$chr = locus$CHR; loc$start = locus$START; loc$stop = locus$STOP; loc$snps = locus$SNPS
	
	# add phenotype info (modified + added eQTL check (cdl 18/3))
	if (is.null(phenos)) phenos = as.character(input$info$phenotype) 
	if (any(! phenos %in% as.character(input$info$phenotype))) {print(paste0("Error: Invalid phenotype ID(s) provided: '",paste0(phenos[! phenos %in% as.character(input$info$phenotype)]), collapse="', '","'")); loc=NULL; return(NULL) }
	if ("eqtl" %in% names(input$info) && any(input$info$eqtl[match(phenos, input$info$phenotype)]) && class(locus) != "gene") {print("Error: use function process.eqtl.locus when analyzing eQTL input"); loc=NULL; return(NULL)} #still works on eQTL input if eqtl phenotype not actually analysed
	loc$phenos = phenos; loc$P = length(loc$phenos)
 	loc$binary = input$info$binary[match(loc$phenos, input$info$phenotype)]; names(loc$binary) = loc$phenos

	# get locus SNPs
	if (!is.null(loc$snps)) {
		# if available, use SNP list
		if (!is.list(loc$snps)) {  # added cdl 18/3
			loc$snps = tolower(unlist(strsplit(loc$snps, ';')))
		}
	} else {
		# if not, use bim file coordinates
		loc$snps = tolower(input$ref$bim$snp.name[input$ref$bim$chromosome == loc$chr & input$ref$bim$position >= loc$start & input$ref$bim$position <= loc$stop])
	}
	# subset to SNPs that passed input processing (i.e. exists across data sets + were aligned)
	loc$snps = unique(intersect(input$analysis.snps, loc$snps))	# taking unique in case there are duplicates in loc$snps; using to intersect() to ensure order is same as reference data set
	
	# check that no. SNPs > min.K
	loc$n.snps = length(loc$snps)
	if (loc$n.snps < min.K) { print(paste("Fewer than",min.K,"SNPs in locus",loc$id)); loc=NULL; return(NULL) }
	
	# read in genotype data (only locus SNPs)
	X = read.plink.custom(input$ref.prefix, df.bim=input$ref, select.snps=loc$snps)
	X = as(X$genotypes, "numeric")
	X = scale(X)			# standardise
	X[is.na(X)] = 0			# mean imputation for missing data
	X = scale(X)
	N.ref = nrow(X)
	
	# remove non variant SNPs
	non.var = apply(X, 2, function(x) all(is.na(x)))
	if(any(non.var)) {
		X = X[,!non.var]
		loc$snps = loc$snps[!non.var]
	}
	
	# check that order of SNPs match
	if (!all(tolower(colnames(X))==loc$snps)) { stop(paste0("Program Error: Mismatching SNP order between reference data and analysis SNPs in locus", loc$id,". Please contact developer.")) } # this should never be triggered, but just in case
	
	# prune redundant PCs
	svd = try(svd(X), silent=T)					# try svd
	if (class(svd)=="try-error") { svd = try(svd(X), silent=T) } 	# if fails, try again (some randomness causing error occasionally)
	if (class(svd)=="try-error") {					# if it fails again, do eig
		eig = eigen(cor(X))
		lambda = eig$values
		Q = eig$vectors
	} else {
		lambda = svd$d * svd$d / (N.ref-1)
		Q = svd$v
	}
	cum.perc = cumsum(lambda / sum(lambda) * 100)
	keep = 1:min(which(cum.perc >= prune.thresh))
	
	# Check remaining nr of PCs
	K.max = length(keep) # K renamed to K.max after updating the fit.logistic() function; K is defined within while loop
	if (K.max < min.K) { print(paste("Error: Fewer than",min.K,"PCs in locus",loc$id)); loc=NULL; return(NULL) } # K can drop below min.K during for-loop below, so need to check here; this also ensures the while is guaranteed to terminate
	
	# Subset sum-stats and get locus N
	loc.sum = list(); loc$N = rep(NA, loc$P); names(loc$N) = loc$phenos
	for (i in loc$phenos) {
		# subset sumstats to locus SNPs
		loc.sum[[i]] = input$sum.stats[[i]][input$sum.stats[[i]]$SNP %in% loc$snps,]
		if (!all(loc.sum[[i]]$SNP==loc$snps)) { stop(paste0("Program Error: Mismatching SNP order between sum-stats and reference data for locus", loc$id,". Please contact developer.")) }	# this should never be triggered, but just in case
		# get N
		loc$N[i] = mean(loc.sum[[i]]$N, na.rm=T)	# get mean locus N (for sumstats i)
		if (is.na(loc$N[i])) { loc$N[i] = mean(input$sum.stats[[i]]$N, na.rm=T) }	# if all are NA, set to mean N across all SNPs in sumstats
		loc.sum[[i]]$N[is.na(loc.sum[[i]]$N)] = loc$N[i]	# use mean imputation any for missing per SNP N 
	}

	# cap the number of PCs at a proportion of the (lowest) GWAS input sample size (cdl 16/3)
	if (!is.null(max.prop.K) && !is.na(max.prop.K)) {
		max.K = floor(max.prop.K * min(loc$N))
		if (max.K < min.K) max.K = min.K #min.K setting takes precedence if in conflict
		if (length(keep) > max.K) {
			keep = keep[1:max.K]
			Q = Q[,keep]
			lambda = lambda[keep]
		}
	}

	# Check remaining nr of PCs (moved down, cdl 16/3)
	K.max = length(keep) # K renamed to K.max after updating the fit.logistic() function; K is defined within while loop
	if (K.max < min.K) { print(paste("Error: Fewer than",min.K,"PCs in locus",loc$id)); loc=NULL; return(NULL) } # K can drop below min.K during for-loop below, so need to check here; this also ensures the while is guaranteed to terminate

	# Define G/R for binary phenotypes (used by fit.logistic())
	if (any(loc$binary)) { R = Q[,keep] %*% diag(1/sqrt(lambda[keep])); G = X %*% R } # NOTE: these are further subsetted within the fit.logistic() function if any PCs are dropped, so there is no need to subset again in the while loop
	
	loc$sigma = loc$h2.obs = loc$h2.latent = rep(NA, loc$P); names(loc$sigma) = names(loc$h2.obs) = names(loc$h2.latent) = loc$phenos
	dropped = c()	# any PCs that might be dropped by fit.logistic() due to instability
	
	while (T) {
		keep = which(!(1:K.max %in% dropped))
		loc$K = length(keep)
		if (loc$K < min.K) { print(paste("Error: Fewer than",min.K,"PCs in locus",loc$id)); loc=NULL; return(NULL) } # K can drop below min.K during for-loop below, so need to check here; this also ensures the while is guaranteed to terminate
		
		loc$delta = matrix(NA, loc$K, loc$P); colnames(loc$delta) = loc$phenos # this needs to be defined in here so it updates to right size if further PCs dropped
		rerun = F # will be set to true in case additional PCs dropped during 1:P loop by fit.logistic
		
		for (i in loc$phenos) {
			if (loc$binary[i]) {
				fit = fit.logistic(G, X, R, loc$N[i], subset(input$info,phenotype==i)$prop_cases, loc.sum[[i]]$STAT, loc.sum[[i]]$N, phen.id=i, loc.id=loc$id, dropped=dropped)
				if (is.na(fit[1])) { warning(paste0("Multiple logistic regression model for phenotype '",i,"' in locus ",loc$id," failed to converge (inversion error)")); next() }
				
				# the fit$dropped vector contains IDs of all dropped PCs, including by other phenotypes (or by the same phenotype earlier), so checking here if it got longer
				# if so, this will break out of the current for-loop and re-process all phenotypes with the new dropped vector
				if (length(fit$dropped) > length(dropped)) {
					dropped = fit$dropped
					rerun = T; break
				}
				loc$delta[,i] = fit$beta
				loc$sigma[i] = fit$var
				loc$h2.obs[i] = fit$h2.obs
				
				# compute h2 if pop prevalence info is provided
				if (!is.null(input$info$prevalence)) {
					loc$h2.latent[i] = fit$h2.obs * (subset(input$info,phenotype==i)$prevalence * (1-subset(input$info,phenotype==i)$prevalence) / dnorm(qnorm(subset(input$info,phenotype==i)$prevalence))^2)
				}
			} else {
				r = loc.sum[[i]]$STAT / sqrt(loc.sum[[i]]$STAT^2 + loc.sum[[i]]$N - 2)	# for continuous phenos, convert Z to r
				alpha = Q[,keep] %*% diag(1/lambda[keep]) %*% t(Q[,keep]) %*% r 
				
				loc$delta[,i] = diag(c(sqrt(lambda[keep]))) %*% t(Q[,keep]) %*% alpha	# using keep to filter out dropped PCs (since this one is updated every time a PC is dropped)
				eta = t(r) %*% alpha; eta = (loc$N[i]-1) / (loc$N[i]-loc$K-1) * (1-eta)
				loc$sigma[i] = eta/(loc$N[i]-1)
				
				# get h2
				loc$h2.obs[i] = 1 - (1 - t(loc$delta[,i]) %*% loc$delta[,i]) * ((loc$N[i]-1) / (loc$N[i]-loc$K-1))
			}
		}
		if (!rerun) break # end while loop
	}
	if (loc$P > 1) { loc$sigma = diag(loc$sigma) }; if (!is.null(input$sample.overlap)) { loc$sigma = sqrt(loc$sigma) %*% as.matrix(input$sample.overlap[loc$phenos,loc$phenos]) %*% sqrt(loc$sigma) } # if P > 1 is just because the diag() doesn't work for single phenotype
	loc$sigma = as.matrix(loc$sigma); dimnames(loc$sigma) = rep(list(loc$phenos),2)
	
	# cap any negative h2's at 0
	loc$h2.obs[loc$h2.obs<0] = 0
	loc$h2.latent[loc$h2.latent<0] = 0

	# remove h2's if K/N ratio too high (cdl 22/3)
	thresh.ratio = 0.1
	if (any(loc$K / loc$N > thresh.ratio)) {
		loc$h2.obs[loc$K/loc$N > thresh.ratio] = NA
		loc$h2.latent[loc$K/loc$N > thresh.ratio] = NA
	}

	# get full omega
	loc$omega = t(loc$delta)%*%loc$delta / loc$K - loc$sigma
	loc$omega.cor = suppressWarnings(cov2cor(loc$omega))


	# check if any phenos have negative sigma or omega;
	neg.var = diag(loc$sigma) < 0 | diag(loc$omega) < 0;
	failed = neg.var | is.na(neg.var)	#those that failed due to N < K or wsw.inversion will be NA in the neg.var variable

	if (drop.failed) { # cdl 16/3
		if (all(neg.var, na.rm=T)) { print(paste0("Error: Negative variance estimate for all phenotypes in locus ",loc$id,". This locus cannot be analysed")); loc=NULL; return(NULL) }	# print error if all had negative variance estimate
		if (any(neg.var, na.rm=T)) { print(paste0("Warning: Negative variance estimate for phenotype(s) '",paste(loc$phenos[which(neg.var)],collapse="', '"),"' in locus ",loc$id,"; Dropping these as they cannot be analysed")) }

		# remove all phenotypes that failed (either due to negative variance, N < K, or wsw.inversion problem)
		if (any(failed)) {
			if (all(failed)) { print(paste0("Error: Processing of all phenotypes in locus ",loc$id," failed (see preceeding warning messages for details)")); loc=NULL; return(NULL) }

			for (var in c("N","binary","phenos","h2.obs","h2.latent")) { loc[[var]] = loc[[var]][!failed] }	 # vectors

			loc$delta = as.matrix(loc$delta[,!failed]); colnames(loc$delta) = loc$phenos	# delta

			for (var in c("sigma","omega","omega.cor")) {   # symmetric matrices
				loc[[var]] = as.matrix(loc[[var]][!failed,!failed])
				dimnames(loc[[var]]) = rep(list(loc$phenos),2)	# need to add phenotype IDs again due to as.matrix()
			}
		}
	} else {
		loc$failed = failed
	}
	return(loc)
}



#' Re-process locus to meta-analyse of selected phenotypes
#'
#' Will combine all elements of the requested phenotypes using standard inverse variance weighting, allowing them to be analysed as a single phenotype via the  multivariate analysis functions.
#' Note that the univariate test cannot currently be applied to meta-analysed phenotypes, so please do that beforehand on each phenotype individually.
#'
#' @param locus Locus object defined using the \code{\link{process.locus}} function.
#' @param meta.phenos Phenotypes you want to meta-analyse
#'
#' #'
#' @return This function returns an object just like that \code{\link{process.locus}} function, containing general locus info, the relevant processed sumstats, and info about the input phenotypes.
#'
#' \itemize{
#'     \item id - locus ID
#'     \item chr/start/stop - locus coordinates
#'     \item snps - list of locus SNPs
#'     \item N.snps - number of SNPs
#'     \item K - number of PCs
#'     \item delta - PC projected joint SNP effects for each phenotype
#'     \item sigma - sampling covariance matrix
#'     \item omega - genetic covariance matrix
#'     \item omega.cor - genetic correlation matrix
#'     \item N - vector of average N across locus SNPs for each phenotype
#'     \item phenos - phenotype IDs
#'     \item binary - boolean vector indicating whether phentoypes are binary
#' }
#'
#' @export
meta.analyse.locus = function(locus, meta.phenos) {
	if (any(! meta.phenos %in% locus$phenos)) { stop(paste0("Phenotype(s) not present in locus object: '", paste0(meta.phenos[! meta.phenos %in% locus$phenos],collapse="', '"),"'"))}
	if (length(meta.phenos) < 2) { stop("Specify at least two phenotypes") }
	
	non.meta = locus$phenos[!locus$phenos %in% meta.phenos]
	
	# create new 'meta' locus
	meta = as.environment(as.list(locus, all.names=T))
	meta$phenos = c(paste0("meta.",paste(meta.phenos,collapse="-")), non.meta)
	
	# get params
	w = 1/diag(locus$sigma[meta.phenos, meta.phenos])
	w = as.matrix(w / sum(w))  # this is w.prime in derivs
	
	# create new deltas
	delta.meta = locus$delta[,meta.phenos] %*% w
	meta$delta = cbind(delta.meta, locus$delta[,non.meta])
	colnames(meta$delta) = meta$phenos
	
	# create omega & sigma
	for (v in c("omega","sigma")) {
		v.meta = t(w) %*% locus[[v]][meta.phenos, meta.phenos] %*% w
		v.meta.R = t(w) %*% locus[[v]][non.meta, meta.phenos]
		
		meta[[v]] = rbind(cbind(v.meta, v.meta.R),
						  cbind(t(v.meta.R), locus[[v]][non.meta, non.meta]))
		
		colnames(meta[[v]]) = rownames(meta[[v]]) = meta$phenos
	}
	meta$omega.cor = suppressWarnings(cov2cor(meta$omega))
	
	# set remaining variables
	# for now these are all set to NA, which means they cannot be processed to the univariate analysis function
	for (v in c("binary","N","h2.obs","h2.latent")) {
		meta[[v]] = c(NA, locus[[v]][non.meta])
		names(meta[[v]]) = meta$phenos
	}
	return(meta)
}




### sum-stats read-in and processing functions ###
format.pvalues = function(input, i, min.pval) {
	input$sum.stats[[i]]$P = as.numeric(input$sum.stats[[i]]$P) # if p < numerical limits in R, they will be read in as char
	input$sum.stats[[i]] = input$sum.stats[[i]][!is.na(input$sum.stats[[i]]$P), ]		# remove NA pvalues
	input$sum.stats[[i]]$P[input$sum.stats[[i]]$P < min.pval] = min.pval				# format zero p-values
}

filter.n = function(input, i) {
	input$sum.stats[[i]] = input$sum.stats[[i]][!is.na(input$sum.stats[[i]]$N),]		# remove missing N
	input$sum.stats[[i]] = subset(input$sum.stats[[i]], N > 0)							# remove negative N
}

format.z = function(input, i, min.pval) {
	input$sum.stats[[i]]$STAT = as.numeric(input$sum.stats[[i]]$STAT)
	inf.z = which(abs(input$sum.stats[[i]]$STAT) > abs(qnorm(min.pval/2))) 
	input$sum.stats[[i]]$STAT[inf.z] = abs(qnorm(min.pval/2)) * sign(input$sum.stats[[i]]$STAT[inf.z]) # remove inf z scores
}

process.sumstats = function(input) {
	min.pval = 1e-300
	
	# check if input files exist
	check.files.exist(input$info$filename)
	
	# possible / required headers
	headers = list(); headers$N = c("N","NMISS","N_analyzed"); headers$STAT = c("Z","T","STAT","Zscore"); headers$SNP = c("SNP","ID","SNPID_UKB","SNPID","MarkerName","RSID","RSID_UKB"); headers$B = c("B","BETA"); headers$A1=c("A1","ALT"); headers$A2=c("A2","REF")	 # header variations
	
	# read in sumstats
	print("...Reading in sumstats")
	input$sum.stats = list()
	for (i in 1:input$P) {
		input$sum.stats[[i]] = data.table::fread(input$info$filename[i],data.table=F)				# read in data
		
		for (h in names(headers)) { 																# check for required headers and rename
			if (sum(colnames(input$sum.stats[[i]]) %in% headers[[h]])==0 & h!="STAT" & h!="B") { 	# skip check here for B/Z class headers
				stop(paste0("No valid ",h," header in sumstats file for phenotype: '",input$info$phenotype[i],"'"))
			} else {
				# check if there are more than one valid headers ()
				if (sum(colnames(input$sum.stats[[i]]) %in% headers[[h]]) > 1) { # if so, remove all but the first valid column and print warning
					col.remove = colnames(input$sum.stats[[i]])[colnames(input$sum.stats[[i]]) %in% headers[[h]]][-1]
					input$sum.stats[[i]] = input$sum.stats[[i]][,!colnames(input$sum.stats[[i]]) %in% col.remove]
					print(paste0("Warning: More than one valid ",h," header for phenotype: '",input$info$phenotype[i],"'. Only retainig the first ('",colnames(input$sum.stats[[i]])[colnames(input$sum.stats[[i]]) %in% headers[[h]]],"')."))
				}
				colnames(input$sum.stats[[i]])[colnames(input$sum.stats[[i]]) %in% headers[[h]]] = h		# rename to standard header format
			}
		}
		input$sum.stats[[i]]$SNP = tolower(input$sum.stats[[i]]$SNP)								# set SNP IDs to lower case
		if ("P" %in% colnames(input$sum.stats[[i]])) { format.pvalues(input, i, min.pval) }			# format p-values
		filter.n(input,i)																			# remove missing or negative N
		
		# convert effect sizes to Z if there are none already
		if (! "STAT" %in% colnames(input$sum.stats[[i]])) {
			effect.size = c("B","OR","logOdds")
			param = effect.size[effect.size %in% colnames(input$sum.stats[[i]])][1]	# get relevant effect size parameter (takes first if multiple, e.g. both OR/logOdds)
			# check that both param + P exist
			if (is.na(param) | ! "P" %in% colnames(input$sum.stats[[i]])) {
				stop(paste0("Lack of valid statistics provided (e.g. Z, T, or BETA/OR/logOdds + P) for phenotype: '",input$info$phenotype[i],"'."))
			} else {
				# if so get the sign
				if (param=="OR") { sign = ifelse(input$sum.stats[[i]][[param]] > 1, 1, -1) } else { sign = sign(input$sum.stats[[i]][[param]]) }
			}
			if (all(sign > 1)) { print("Warning: The signs of betas/logOdds are positive; are you sure you did not provide ORs? (if so, please use the appropriate header)") }
			
			# get Z
			Z = -qnorm(input$sum.stats[[i]]$P/2)
			input$sum.stats[[i]]$STAT = Z * sign
		}
		format.z(input, i, min.pval) # format Z stats (truncate out of range ones)
		
		input$sum.stats[[i]] = input$sum.stats[[i]][c("SNP","A1","A2","STAT","N")] 		# retain relevant columns
	}
	names(input$sum.stats) = input$info$phenotype
	
	# reference data bim/afreq file
	print("...Reading in SNP info from reference data")
	input$ref = read.bim.custom(input$ref.prefix, as.env=T)
	input$ref$bim$snp.name = tolower(input$ref$bim$snp.name)		# setting ref SNPs tolower
	
	# get common SNPs
	print("...Extracting common SNPs")
	harmonize.snps(input)

	# align SNPs
	print("...Aligning effect alleles to reference data set")
	align(input)
	
	return(input)
}

harmonize.snps = function(input, check.index=NULL) {
	#determine common snps
	input$analysis.snps = intersect(input$ref$bim$snp.name, input$sum.stats[[1]]$SNP)	# all SNPs will be ordered according to bim file
	if (input$P>1) {
		if (is.null(check.index)) check.index = 2:input$P
		for (i in check.index) { input$analysis.snps = intersect(input$analysis.snps, input$sum.stats[[i]]$SNP) }
	}
	if (length(input$analysis.snps) < 3) { stop("Less than 3 SNPs shared across data sets; make sure you have matching SNP ID formats across sumstats / reference data sets")
	} else { print(paste("...",length(input$analysis.snps),"SNPs shared across data sets")) }

	if (!all(input$ref$bim$snp.name[input$ref$bim$snp.name %in% input$analysis.snps] == input$analysis.snps)) { stop("Program Error: SNPs not ordered according to reference after subsetting. Please contact developer.") }

	# subset sumstats to common snps
	for (i in 1:input$P) {
		input$sum.stats[[i]] = input$sum.stats[[i]][match(input$analysis.snps, input$sum.stats[[i]]$SNP),]
		if(!all(input$sum.stats[[i]]$SNP==input$analysis.snps)) { stop("Program Error: sum-stats SNPs do not match the reference after subsetting. Please contact developer.") }
	}
}

process.sample.overlap = function(sample.overlap.file, phenos) {
	sample.overlap = as.matrix(read.table(sample.overlap.file, check.names=F))	# check.names=F prevents the X in front of numeric phenotype IDs
	idx = match(phenos, colnames(sample.overlap))								# get relevant indices
	if (any(is.na(idx))) { stop(paste0("Phenotype(s) not listed in sample overlap file: '", paste0(phenos[is.na(idx)], collapse="', '"),"'")) }
	return(cov2cor(sample.overlap[idx,idx]))
}

get.input.info = function(input.info.file, sample.overlap.file, ref.prefix, phenos=NULL, input.dir=NULL) {
	check.files.exist(c(input.info.file, sample.overlap.file, paste0(ref.prefix, c(".bim",".bed",".fam"))))
	
	input = new.env(parent=globalenv())
	input$info = read.table(input.info.file, header=T, stringsAsFactors=F) # read in input info file
	if (!all(c("phenotype","cases","controls","filename") %in% colnames(input$info))) { stop("Please provide an input.info file with headers: phenotype, cases, controls, filename") }  # header check on input file
	if (is.null(phenos)) { phenos = input$info$phenotype }
	if (!all(phenos %in% input$info$phenotype)) { stop(paste0("Phenotype(s) not listed in input info file: '", paste0(phenos[!phenos %in% input$info$phenotype], collapse="', '"),"'")) }
	if (!is.null(input.dir)) input$info$filename = paste0(input.dir, "/", input$info$filename) # added (cdl 17/3)
	
	input$info = input$info[match(phenos, input$info$phenotype),]			# match to phenotypes of interest
	input$info$N = input$info$cases + input$info$controls					# total N column
	input$info$prop_cases = input$info$cases / input$info$N 				# get proportion cases
	input$info$binary = !is.na(input$info$prop_cases) & (input$info$prop_cases != 1)	# infer binary phenotypes from prop_cases==1 or NA
	input$P = length(input$info$phenotype)
	input$ref.prefix = ref.prefix
	if (input$P > 1 & !is.null(sample.overlap.file)) { input$sample.overlap = process.sample.overlap(sample.overlap.file, phenos) } else { input$sample.overlap=NULL } # setting sample overlap to null if only one phenotype
	
	return(input)
}



#' Process input
#' 
#' Processes summary statistics and related information (e.g. sample overlap, case/control ratios).
#' Ensures alignment of SNPs across data sets
#' 
#' @param input.info.file Name of info file containing 
#' phenotype IDs ('phenotype'), N cases ('cases'), N controls ('controls'), sumstats file ('filename').
#' For continuous phenotypes, the number of controls should be set to 0, while cases can just  be set to 1 (this is only used for computing the case/control ratio, which should be 1 for continuous phenotypes).
#'
#' @param sample.overlap.file Name of file with sample overlap information.
#' Can be set to NULL if there is no overlap
#' 
#' @param ref.prefix Prefix of reference genotype data in plink format (*.bim, *.bed, *.fam)
#' 
#' @param phenos A vector of phenotype IDs can be provided if only a subset of phenotypes are desired
#' (if NULL, all phenotypes in the input info file will be processed). 
#' This can be convenient if a subset from a larger number of phenotypes are analysed, as only a single input.info / sample overlap file needs to be created.
#'
#' @param input.dir Directory containing the files specified in the info file.
#'
#' @return An object containing processed input data and related info
#' \itemize{
#'     \item info - the processed input info file. Columns added during processing: 
#'     \itemize{
#'        \item N = cases + controls
#'        \item prop_cases = cases / N
#'        \item binary = !is.na(prop_cases) & (prop_cases != 1)
#'     }
#'     \item P - number of phenotypes
#'     \item sample.overlap - sample overlap matrix
#'     \item sum.stats - processed summary statistics (SNP aligned effect sizes, subsetted to common SNPs across data sets, effect sizes converted to Z, etc)
#'     \item ref.prefix - genotype reference data prefix
#'     \item analysis.snps - subset of SNPs that are shared across all data sets and were not removed during alignment
#'     \item unalignable.snps - SNPs removed during alignment (e.g. for being strand ambiguous)
#'     \item ref - environment containing the genotype reference data bim file (ref$bim)
#' }
#' 
#' @export
process.input = function(input.info.file, sample.overlap.file, ref.prefix, phenos=NULL, input.dir=NULL) {
	print("...Processing input info")
	input = get.input.info(input.info.file, sample.overlap.file, ref.prefix, phenos, input.dir)
	if (all(input$info$binary)) { 
		print(paste0("...** All phenotypes treated as BINARY ('",paste0(input$info$phenotype,collapse="', '"),"')"))
	} else if (all(!input$info$binary)) {
		print(paste0("...** All phenotypes treated as CONTINUOUS ('",paste0(input$info$phenotype,collapse="', '"),"')"))
	} else {
		#print(paste0("...** Treating '",paste0(subset(input$info, binary==T)$phenotype,collapse="', '"),"' as BINARY and '",paste0(subset(input$info, binary==F)$phenotype,collapse=", '"),"' as CONTINUOUS"))
		print(paste0("...** BINARY: '",paste0(subset(input$info, binary==T)$phenotype,collapse="', '"),"'"))
		print(paste0("...** CONTINUOUS: '",paste0(subset(input$info, binary==F)$phenotype,collapse=", '"),"'"))
	}
	input = process.sumstats(input)
	return(input)
}


#' Read in locus file
#' 
#' This function reads in the locus file, defining the locus boundaries using either coordinates (headers: 'CHR', 'START', 'STOP') and/or ';' separated SNP lists (header: 'SNPS').
#' A locus ID column is also required (header: 'LOC').\cr
#' \cr
#' If both coordinates and SNP list columns are provided, only the SNP lists will be used for subsetting the reference data (this can be convenient if the SNP coordinates are based on a different genome build version than the reference).\cr
#' 
#' @param loc.file Name of locus file
#' 
#' @export
read.loci = function(loc.file) {
	check.files.exist(loc.file)
	loci = data.table::fread(loc.file, data.table=F)
	if (!(all(c("LOC","CHR","START","STOP") %in% colnames(loci)) | all(c("LOC","SNPS") %in% colnames(loci)))) { stop("Locus file does not contain required headers ('LOC' + 'CHR','START','STOP' and/or 'SNPS')") }
	return(loci)
}


check.files.exist = function(infiles) {
	if (!all(file.exists(infiles))) { stop(paste0("Missing input file(s): '",paste0(infiles[!file.exists(infiles)],collapse=", '"),"'"),". Please ensure correct file names and paths have been provided.") }
}

check.folders.exist = function(folders) {
	if (!all(file.exists(folders))) { stop(paste0("Missing folder(s): '",paste0(folders[!file.exists(folders)],collapse=", '"),"'"),". Please ensure correct paths have been provided.") }
}
