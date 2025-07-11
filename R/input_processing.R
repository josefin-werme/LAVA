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


process.locus = function(locus, input, phenos=NULL, min.K=2, prune.thresh=99, max.prop.K=0.75, drop.failed=T, max.block.size=3000, nref.correction=T, cap.estimates=T) {
	if (is.data.frame(locus) && nrow(locus)!=1) { print("Error: Locus info provided for incorrect number of loci. Please provide only a single locus at a time"); loc=NULL; return(NULL) }  # modified cdl 18/3
	if (!(all(c("LOC","CHR","START","STOP") %in% names(locus)) | all(c("LOC","SNPS") %in% names(locus)))) { print("Error: Locus info data frame is missing some or all of the required headers ('LOC' + 'CHR','START','STOP' and/or 'SNPS')"); loc=NULL; return(NULL) }
	
	min.K = max(min.K, 2) # just to make sure it isn't 1, because that leads to some matrix multiplication dimension issues
	
	# define locus environment & add locus info
	loc = new.env(parent=globalenv())
	loc$id = locus$LOC; loc$chr = locus$CHR; loc$start = locus$START; loc$stop = locus$STOP; loc$snps = locus$SNPS
	
	# add phenotype info (modified + added eQTL check (cdl 18/3))
	if (is.null(phenos)) { 
		phenos = as.character(input$info$phenotype) 
	} else {
		phenos = as.character(phenos) # to guard against unintentional factors
	}
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
		loc$snps = input$reference$snp.info$SNP[input$reference$snp.info$CHR == loc$chr & input$reference$snp.info$POS >= loc$start & input$reference$snp.info$POS <= loc$stop]
	}
	# subset to SNPs that passed input processing (i.e. exists across data sets + were aligned)
	loc$snps = unique(intersect(input$analysis.snps, loc$snps))	# taking unique in case there are duplicates in loc$snps; using to intersect() to ensure order is same as reference data set
	

	# load LD data
	ld = read.ld(input$reference, loc$snps, require.freq=any(loc$binary))
	if (!all(ld$info$SNP %in% loc$snps)) stop("inconsistency in loaded SNPs")
	
	loc$snps = ld$info$SNP; loc$n.snps = length(loc$snps)
	if (loc$n.snps < min.K) { print(paste("Fewer than",min.K,"SNPs in locus",loc$id)); loc=NULL; return(NULL) }

	# Subset sum-stats and get locus N
	loc.sum = list(); loc$N = rep(NA, loc$P); names(loc$N) = loc$phenos; drop.snps = c()
	for (i in loc$phenos) {
		# subset sumstats to locus SNPs
		loc.sum[[i]] = input$sum.stats[[i]][input$sum.stats[[i]]$SNP %in% loc$snps,]
		if (!all(loc.sum[[i]]$SNP==loc$snps)) { stop(paste0("Program Error: Mismatching SNP order between sum-stats and reference data for locus", loc$id,". Please contact developer.")) }	# this should never be triggered, but just in case
		
		# get N
		loc$N[i] = mean(loc.sum[[i]]$N, na.rm=T)	# get mean locus N (for sumstats i)
		if (is.na(loc$N[i])) { loc$N[i] = mean(input$sum.stats[[i]]$N, na.rm=T) }	# if all are NA, set to mean N across all SNPs in sumstats
		loc.sum[[i]]$N[is.na(loc.sum[[i]]$N)] = loc$N[i]	# use mean imputation any for missing per SNP N 
		
		# compute marginal correlations
		if (loc$binary[i]) {
			loc.sum[[i]]$CORR = process.binary(loc.sum[[i]]$STAT, loc.sum[[i]]$N, ld$info$FREQ, input$info$prop_cases[input$info$phenotype == i])
			miss = is.na(loc.sum[[i]]$CORR)
			if (any(miss)) {
				add = loc$snps[miss]; drop.snps = c(drop.snps, add)
				print(paste0("Warning: Unable to reconstruct marginal SNP statistics for binary phenotype '", loc$phenos[i] , "' in locus ", loc$id, " for SNP(s) ", paste0(add, collapse=", "), "; Dropping these as they cannot be analysed")) 
			}
		} else {
			loc.sum[[i]]$CORR = loc.sum[[i]]$STAT / sqrt(loc.sum[[i]]$STAT^2 + loc.sum[[i]]$N - 2)	# convert Z to r	
		}
	}
	
	if (length(drop.snps) > 0) {
		keep = !(loc$snps %in% drop.snps)
		loc$snps = loc$snps[keep]; loc$n.snps = length(loc$snps)
		if (loc$n.snps < min.K) { print(paste("Fewer than",min.K,"SNPs in locus",loc$id)); loc=NULL; return(NULL) }
		
		ld$info = ld$info[keep,]
		if (ld$mode == "plink") ld$data = ld$data[, keep, drop=F]
		else ld$ld = ld$ld[keep, keep, drop=F]
		
		for (i in loc$phenos) loc.sum[[i]] = loc.sum[[i]][keep,]
	}


	# decompose
	R = decompose.ld(ld, prune.thresh, max.block.size)
	
	# cap the number of PCs at a proportion of the (lowest) GWAS input sample size (cdl 16/3)
	if (!is.null(max.prop.K) && !is.na(max.prop.K)) {
		max.K = floor(max.prop.K * min(loc$N))
		if (max.K < min.K) max.K = min.K #min.K setting takes precedence if in conflict
		if (ncol(R) > max.K) R = R[,1:max.K]
	}
		
	# Check remaining nr of PCs 
	loc$K = ncol(R)
	if (loc$K < min.K) { print(paste("Error: Fewer than",min.K,"PCs in locus",loc$id)); loc=NULL; return(NULL) } # K can drop below min.K during for-loop below, so need to check here; this also ensures the while is guaranteed to terminate

	loc$nref.scale = ifelse(nref.correction, (input$reference$sample.size - loc$K - 1) / (input$reference$sample.size - 1), 1)		

	loc$sigma = loc$h2.obs = loc$h2.latent = rep(NA, loc$P); loc$ascertained.h2 = rep(F, loc$P); names(loc$sigma) = names(loc$h2.obs) = names(loc$h2.latent) = names(loc$ascertained.h2) = loc$phenos
	loc$delta = matrix(NA, loc$K, loc$P); colnames(loc$delta) = loc$phenos 
	for (i in loc$phenos) {
		loc$delta[,i] = t(R) %*% loc.sum[[i]]$CORR 

		dtd = sum(loc$delta[,i]^2)
		loc$sigma[i] = (1 - dtd) / (loc$N[i]-loc$K-1)
		loc$h2.obs[i] = (dtd - loc$sigma[i] * loc$K) * loc$nref.scale 
		
		if (loc$binary[i]) {
			case.proportion = input$info$prop_cases[input$info$phenotype == i]
			prevalence = if (!is.null(input$info$prevalence)) input$info$prevalence[input$info$phenotype == i]
			if (!is.null(prevalence) && is.finite(prevalence)) loc$ascertained.h2[i] = TRUE # Lee et al. 2011, formula 23
			else prevalence = case.proportion # fall back to Lee et al. 2011, formula 17
			
			loc$h2.latent[i] = loc$h2.obs[i] / dnorm(qnorm(prevalence))^2 * (prevalence * (1 - prevalence))^2 / (case.proportion * (1 - case.proportion))
		}
	}

	if (loc$P > 1) { loc$sigma = diag(loc$sigma) }; if (!is.null(input$sample.overlap)) { loc$sigma = sqrt(loc$sigma) %*% as.matrix(input$sample.overlap[loc$phenos,loc$phenos]) %*% sqrt(loc$sigma) } # if P > 1 is just because the diag() doesn't work for single phenotype
	loc$sigma = as.matrix(loc$sigma); dimnames(loc$sigma) = rep(list(loc$phenos),2)
	
	# cap any negative h2's at 0
	if (cap.estimates) {for (var in grep("^h2", names(loc), value=T)) loc[[var]][loc[[var]] < 0] = 0}
	
	# remove h2's if K/N ratio too high (cdl 22/3)
	thresh.ratio = 0.1
	if (any(loc$K / loc$N > thresh.ratio)) {for (var in grep("^h2", names(loc), value=T)) loc[[var]][loc$K/loc$N > thresh.ratio] = NA}

	# get full omega
	loc$omega = (t(loc$delta)%*%loc$delta / loc$K - loc$sigma) * loc$nref.scale 
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

			for (var in c("N","binary","phenos","h2.obs","h2.latent", "ascertained.h2")) { loc[[var]] = loc[[var]][!failed] }	 # vectors

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


decompose.ld = function(ld, prune.thresh, max.block.size) {
	if (ld$mode == "plink") {
		X = scale(ld$data); X[is.na(X)] = 0; X = scale(X) 	
		N.ref = nrow(X)
		
		svd = try(svd(X), silent=T)					# try svd
		if (class(svd) == "try-error") { svd = try(svd(X), silent=T) } 	# if fails, try again (some randomness causing error occasionally)
		if (class(svd) != "try-error") {					
			lambda = svd$d * svd$d / (N.ref-1)
			Q = svd$v
		} else { # if it fails again, fall back to direct decomposition
			ld$ld = cor(X)
			ld$mode = "ld"
			return(decompose.ld(ld, prune.thresh, max.block.size))
		}
	} else {
		no.snps = ncol(ld$ld)
		if (no.snps > max.block.size) {
			no.blocks = ceiling(no.snps / max.block.size)
			index = sort(rep(1:no.blocks, length.out = no.snps))
			
			R.base = NULL
			for (i in 1:no.blocks) {
				curr = which(index == i); M = ld$ld[curr,curr] 
				R.curr = decompose.ld(list(ld = M, mode = "ld"), prune.thresh, max.block.size=Inf)
				add = matrix(0, nrow=no.snps, ncol=ncol(R.curr)); add[curr,] = R.curr
				R.base = cbind(R.base, add)
			}
			
			M = t(R.base) %*% ld$ld %*% R.base
			R.block = decompose.ld(list(ld = M, mode = "ld"), prune.thresh, 2*max.block.size)
			return(R.base %*% R.block)
		} else {
			decomp = eigen(ld$ld)
			lambda = decomp$values
			Q = decomp$vectors
		}
	}

	cum.perc = cumsum(lambda); cum.perc = cum.perc / max(cum.perc) * 100
	K = min(which(cum.perc >= prune.thresh))
	return(Q[,1:K] %*% diag(1/sqrt(lambda[1:K])))
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
format.pvalues = function(pvals, min.pval) {
	pvals = as.numeric(pvals) # if p < numerical limits in R, they will be read in as char
	pvals[!is.na(pvals) & pvals < min.pval] = min.pval				# format zero p-values
	return(pvals)
}

format.z = function(zstats, min.pval) {
	zstats = as.numeric(zstats)
	inf.z = which(abs(zstats) > abs(qnorm(min.pval/2)))
	zstats[inf.z] = abs(qnorm(min.pval/2)) * sign(zstats[inf.z]) # remove inf z scores
	return(zstats)
}


read.sumstats.file = function(filename, pheno.label, min.pval=1e-300, N.override=NULL) {
	# possible / required headers
	headers = list(); headers$STAT = c("Z","T","STAT","Zscore"); headers$SNP = c("SNP","ID","SNPID_UKB","SNPID","MarkerName","RSID","RSID_UKB"); headers$B = c("B","BETA"); headers$A1=c("A1","ALT"); headers$A2=c("A2","REF")	 # header variations
	if (is.null(N.override)) headers$N = c("N","NMISS","N_analyzed");

	sum.stats = data.table::fread(filename, data.table=F)				# read in data

	for (h in names(headers)) { 																# check for required headers and rename
		if (sum(colnames(sum.stats) %in% headers[[h]])==0 & h!="STAT" & h!="B") { 	# skip check here for B/Z class headers
			stop(paste0("No valid ",h," header in sumstats file for phenotype: '", pheno.label, "'"))
		} else {
			# check if there are more than one valid headers ()
			if (sum(colnames(sum.stats) %in% headers[[h]]) > 1) { # if so, remove all but the first valid column and print warning
				col.remove = colnames(sum.stats)[colnames(sum.stats) %in% headers[[h]]][-1]
				sum.stats = sum.stats[,!colnames(sum.stats) %in% col.remove]
				print(paste0("Warning: More than one valid ",h," header for phenotype: '", pheno.label, "'. Only retainig the first ('",colnames(sum.stats)[colnames(sum.stats) %in% headers[[h]]],"')."))
			}
			colnames(sum.stats)[colnames(sum.stats) %in% headers[[h]]] = h		# rename to standard header format
		}
	}
	sum.stats$SNP = tolower(sum.stats$SNP)								# set SNP IDs to lower case

	# convert effect sizes to Z if there are none already
	if (! "STAT" %in% colnames(sum.stats)) {
		effect.size = c("B","OR","logOdds")
		param = effect.size[effect.size %in% colnames(sum.stats)][1]	# get relevant effect size parameter (takes first if multiple, e.g. both OR/logOdds)
		# check that both param + P exist
		if (is.na(param) | ! "P" %in% colnames(sum.stats)) {
			stop(paste0("Lack of valid statistics provided (e.g. Z, T, or BETA/OR/logOdds + P) for phenotype: '", pheno.label, "'."))
		} else {
			# if so get the sign
			if (param=="OR") { sign = ifelse(sum.stats[[param]] > 1, 1, -1) } else { sign = sign(sum.stats[[param]]) }
		}
		if (all(sign > 1)) { print("Warning: The signs of betas/logOdds are positive; are you sure you did not provide ORs? (if so, please use the appropriate header)") }

		# get Z
		Z = -qnorm(format.pvalues(sum.stats$P, min.pval)/2)
		sum.stats$STAT = Z * sign
	}
	sum.stats$STAT = format.z(sum.stats$STAT, min.pval) # format Z stats (truncate out of range ones)
	if (!is.null(N.override)) sum.stats$N = N.override

	keep.col = c((if ("GENE" %in% names(sum.stats)) "GENE"),"SNP","A1","A2","STAT","N")
	sum.stats = sum.stats[keep.col] 		# retain relevant columns
	sum.stats = sum.stats[!apply(is.na(sum.stats), 1, any) & sum.stats$N > 0,] # filter missing values and negative sample sizes
  return(sum.stats)
}



process.sumstats = function(input) {
	# check if input files exist
	check.files.exist(input$info$filename)

	# read in sumstats
	print("...Reading in sumstats")
	input$sum.stats = list()
	for (i in 1:input$P) input$sum.stats[[i]] = read.sumstats.file(input$info$filename[i], input$info$phenotype[i])
	names(input$sum.stats) = input$info$phenotype
	
	# reference data bim/afreq file
	print("...Reading in SNP info from reference data")
	input$reference = load.reference(input$reference)
	
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
	input$analysis.snps = intersect(input$reference$snp.info$SNP, input$sum.stats[[1]]$SNP)	# all SNPs will be ordered according to bim file
	if (input$P>1) {
		if (is.null(check.index)) check.index = 2:input$P
		for (i in check.index) { input$analysis.snps = intersect(input$analysis.snps, input$sum.stats[[i]]$SNP) }
	}
	if (length(input$analysis.snps) < 3) { stop("Less than 3 SNPs shared across data sets; make sure you have matching SNP ID formats across sumstats / reference data sets")
	} else { print(paste("...",length(input$analysis.snps),"SNPs shared across data sets")) }

	if (!all(input$reference$snp.info$SNP[input$reference$snp.info$SNP %in% input$analysis.snps] == input$analysis.snps)) { stop("Program Error: SNPs not ordered according to reference after subsetting. Please contact developer.") }

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
	check.files.exist(c(input.info.file, sample.overlap.file))
	prefix.data = check.reference(ref.prefix)
	
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
	input$reference = prefix.data
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
#' @param ref.prefix Prefix of reference genotype data in PLINK format (*.bim, *.bed, *.fam) or LAVA .bcor/.info LD files
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
#'     \item ref.prefix - prefix of reference data
#'     \item analysis.snps - subset of SNPs that are shared across all data sets and were not removed during alignment
#'     \item unalignable.snps - SNPs removed during alignment (e.g. for being strand ambiguous)
#'     \item reference - environment containing the reference data SNP information
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




### Reference data ###


check.reference = function(prefix) {
	prefix.data = list(prefix = prefix)
	
	plink.files = paste0(prefix, c(".bim",".bed",".fam"))
	if (!any(file.exists(plink.files))) {
		prefix.data$mode = "ld"
		if (grepl(".*_chr.+$", prefix)) {
			chromosome = gsub(".*_chr(.+)$", "\\1", prefix)
			if (!(chromosome %in% 1:23) || is.na(suppressWarnings(as.numeric(chromosome)))) stop(paste0("invalid chromosome '", chromosome, "' specified in prefix '", prefix, "'"))
			
			ld.files = paste0(prefix, c(".info", ".bcor"))
			if (any(file.exists(ld.files))) {
				check.files.exist(ld.files)
				prefix.data$prefix = gsub("(.*)_chr[0-9]+$", "\\1", prefix)
				prefix.data$chromosomes = as.numeric(chromosome)
			}
		} else {
			prefix.data$prefix = prefix =  gsub("(.*)_chr$", "\\1", prefix)
			
			file.list = data.frame(file = Sys.glob(paste0(prefix, "_chr*.*")))
			if (nrow(file.list) > 0) {file.list$type = substring(file.list$file, nchar(file.list$file) - 3); file.list = file.list[file.list$type %in% c("bcor", "info"),]}
			if (nrow(file.list) > 0) {
				spec = gsub(paste0(prefix, "_chr"), "", file.list$file, fixed=T)
				file.list$chromosome = substring(spec, 0, nchar(spec)-5)	
				file.list = file.list[file.list$chromosome %in% 1:23,]
				if (nrow(file.list) > 0) {
					file.list$chromosome = suppressWarnings(as.numeric(file.list$chromosome))
					file.list = file.list[!is.na(file.list$chromosome),]
					file.list = file.list[order(file.list$chromosome, file.list$type),]
				}
			}
			
			if (nrow(file.list) > 0) {
				counts = aggregate(file.list$chromosome, list(chromosome=file.list$chromosome), length)
				if (any(counts[,2] > 2)) stop(paste0("too many LD reference files with prefix '", prefix, "' for chromosome(s) ", paste(counts$chromosome[counts[,2] > 2], collapse=", ")))  
				check.files.exist(unlist(lapply(paste0(prefix, "_chr", unique(file.list$chromosome)), function(b) {paste0(b, c(".bcor", ".info"))})))
				
				missing = which(!(1:23 %in% file.list$chromosome))
				if (length(missing) > 0 && missing[1] != 23) warning(paste0("missing LD reference files with prefix '", prefix, "' for chromosome(s) ", paste(missing, collapse=", ")))
				
				prefix.data$chromosomes = unique(file.list$chromosome)					
			}
		}
		if (length(prefix.data$chromosomes) == 0) stop("no PLINK genotype data or LD reference files with prefix '", prefix, "'")
	} else {
		check.files.exist(plink.files)
		prefix.data$mode = "plink"
	}
	return(prefix.data)		
}


check.parameter = function(name, value, min.value=NULL, max.value=NULL) {
	if (length(value) != 1 || !is.numeric(value) || !is.finite(value)) stop(paste0("value for parameter '", name, "' is not a number"))
	if (!is.null(min.value) && value < min.value) stop(paste0("value for parameter '", name, "' should be equal to or greater than ", min.value))
	if (!is.null(max.value) && value > max.value) stop(paste0("value for parameter '", name, "' should be equal to or smaller than ", max.value))
	return(value)
}



load.reference = function(info) {
	if (is.character(info) && length(info) == 1) info = check.reference(info)
	if (is.null(info$mode) || !(info$mode %in% c("plink", "ld"))) stop("invalid reference data specified")
	
	output = as.environment(info)
	if (info$mode == "ld") {
		if (length(info$chromosomes) == 0) stop("invalid reference data specified")
		output$chromosome.snps = lapply(1:23, function(i) {character(0)})
		for (chr in sort(info$chromosomes)) {
			files = paste0(info$prefix, "_chr", chr, c(".info", ".bcor")); check.files.exist(files)
			snp.info = data.table::fread(files[1], data.table=F, showProgress=F); snp.info$SNP = tolower(snp.info$SNP)
			output$snp.info = rbind(output$snp.info, snp.info)
			output$chromosome.snps[[chr]] = snp.info$SNP
		}
		output$sample.size = max(output$snp.info$NOBS)
	} else {
		files = paste0(info$prefix, c(".bed", ".bim", ".fam")); check.files.exist(files)
		output$sample.size = nrow(data.table::fread(files[3], data.table=F, showProgress=F))
		output$snp.info = data.table::fread(files[2], data.table=F, showProgress=T)[,c(2,1,4,5,6)]
		names(output$snp.info) = c("SNP", "CHR", "POS", "A1", "A2")			
		output$snp.info$SNP = tolower(output$snp.info$SNP)		
	}

	return(output)
}


read.ld = function(reference, snp.ids, require.freq=F) {
	if (is.null(reference$mode) || !(reference$mode %in% c("plink", "ld"))) stop("invalid reference data specified")
	
	use = which(reference$snp.info$SNP %in% snp.ids)
	if (length(use) == 0) stop("none of specified SNP IDs are present in reference data")
	
	chromosome = unique(reference$snp.info$CHR[use])
	if (length(chromosome) > 1) stop("SNPs in input span multiple chromosomes")
	
	output = list(mode = reference$mode, info=reference$snp.info[use,])
	if (reference$mode == "plink") {
		parameters = list(maf=0, mac=1, missing=0.05)
		results = load_plink(reference$prefix, nrow(reference$snp.info), reference$sample.size, as.integer(use - 1), parameters)
		if (length(results$include) > 0) {
			use = use[results$include + 1]
			output$info = reference$snp.info[use,]
			output$data = results$data
			colnames(output$data) = output$info$SNP
			
			if (require.freq) output$info$FREQ = apply(results$data, 2, sum, na.rm=T) / (nrow(results$data) - apply(is.na(results$data), 2, sum)) / 2	
		} else stop("no SNPs remaining after filtering")
	} else {
		chr.snps = reference$chromosome.snps[[chromosome]]; use = which(chr.snps %in% snp.ids)
		if (length(use) != nrow(output$info) || !all(chr.snps[use] == output$info$SNP)) stop("inconsistent SNP IDs when loading LD")
		output$ld = load_ld(reference$prefix, chromosome, length(chr.snps), as.integer(use - 1))
		rownames(output$ld) = colnames(output$ld) = chr.snps[use]
	}
	
	return(output)
}
















