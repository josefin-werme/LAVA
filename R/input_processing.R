
#' Process locus info
#' 
#' Processes the summary statistics for all phenotypes within the specified locus and computes all parameters necessary for analysis.
#' 
#' @param loc Locus info for a single locus, obtained using the \code{\link{read.loci}} function. Expects a locus ID ('LOC') together with locus coordinates ('CHR', 'START', 'STOP') and/or a ';' separated SNP list ('SNPS')
#' @param input Input object created with the \code{\link{process.input}} function, containing relevant summary statistics and related info (e.g. sample overlap, case/control ratio)
#' @param min.K Minimum number of PCs required to process locus (cannot be less than two). If this criterion is not met, the function will fail and the locus cannot be analysed.
#' @param prune.thresh PC pruning threshold governing the maximum number of PCs to retain.
#' PCs are selected as such that the cumulative proportion of variance explained is at least that of the threshold (set to 99 percent by default).
#' 
#' @return Returns an environment containing general locus info, the processed sumstats, and parameters required for analysis. If the function fails (e.g. due to too few SNPs), it will return NULL. 
#' If processing fails for specific phenotypes, only the successful phenotypes will be returned.
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
#'     \item binary - boolean vector indicating whether phentoypes are binary
#'     \item h2.obs - observed local heritability
#'     \item h2.latent - estimated local population heritability (only relevant for binary phenotypes; requires population prevalence to be specified in input info file)
#' }
#' 
#' @export
process.locus = function(locus, input, min.K=2, prune.thresh=99) { 
	if (nrow(locus)!=1) { print("Error: Locus info provided incorrect number of loci. Please process only a signle locus at a time"); loc=NULL; return(NULL) }
	if (!(all(c("LOC","CHR","START","STOP") %in% colnames(locus)) | all(c("LOC","SNPS") %in% colnames(locus)))) { print("Locus info data frame is missing the required headers ('LOC' + 'CHR','START','STOP' and/or 'SNPS')"); loc=NULL; return(NULL) }
	
	# define locus environment & add locus info
	loc = new.env(parent=globalenv())
	loc$id = locus$LOC; loc$chr = locus$CHR; loc$start = locus$START; loc$stop = locus$STOP; loc$snps = locus$SNPS
	
	# add phenotype info
	loc$binary = input$info$binary
	loc$phenos = as.character(input$info$phenotype)
	
	min.K = max(min.K, 2) # just to make sure it isn't 1, because that leads to some matrix multiplication dimension issues
	
	# get locus SNPs
	if (!is.null(loc$snps)) {
		# if available, use SNP list
		loc$snps = tolower(unlist(strsplit(loc$snps, ';')))
	} else {
		# if not, use bim file coordinates
		loc$snps = tolower(input$ref$bim$snp.name[input$ref$bim$chromosome == loc$chr & input$ref$bim$position >= loc$start & input$ref$bim$position <= loc$stop])
	}
	# subset to SNPs that passed input processing (i.e. exists across data sets + were aligned)
	loc$snps = unique(intersect(input$analysis.snps, loc$snps))	# taking unique in case there are duplicates in loc$snps; using to intersect() to ensure order is same as reference data set
	
	# check that no SNPs > min.K
	loc$n.snps = length(loc$snps)
	if (loc$n.snps < min.K) { print(paste("Fewer than",min.K,"SNPs in locus",loc$id)); loc=NULL; return(NULL) }
	
	# read in genotype data (only locus SNPs)
	X = read.plink.custom(input$ref.prefix, df.bim=input$ref, select.snps=loc$snps)
	X = as(X$genotypes, "numeric")
	X = scale(X)			# standardise
	X[is.na(X)] = 0			# mean imputation for missing data
	X = scale(X)
	N.ref = nrow(X)
	
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
	# Retain PCs
	cum.perc = cumsum(lambda / sum(lambda) * 100)
	keep = 1:min(which(cum.perc >= prune.thresh))

	# Subset sum-stats and get locus N
	loc.sum = list(); loc$N = rep(NA, input$P)
	for (i in 1:input$P) {
		# subset sumstats to locus SNPs
		loc.sum[[i]] = input$sum.stats[[i]][input$sum.stats[[i]]$SNP %in% loc$snps,]
		if (!all(loc.sum[[i]]$SNP==loc$snps)) { stop(paste0("Program Error: Mismatching SNP order between sum-stats and reference data for locus", loc$id,". Please contact developer.")) }	# this should never be triggered, but just in case
		# get N
		loc$N[i] = mean(loc.sum[[i]]$N, na.rm=T)							# get mean locus N (for sumstats i)
		if (is.na(loc$N[i])) { loc$N[i] = mean(input$sum.stats[[i]]$N, na.rm=T) }			# if all are NA, set to mean N across all SNPs in sumstats
		loc.sum[[i]]$N[is.na(loc.sum[[i]]$N)] = loc$N[i]						# use mean imputation for missing N within locus sumstats
	}
	
	# Check remaining nr of PCs
	K.max = length(keep) # K renamed to K.max after updating the fit.logistic() function; K is defined within while loop
	if (K.max < min.K) { print(paste("Error: Fewer than",min.K,"PCs in locus",loc$id)); loc=NULL; return(NULL) } # K can drop below min.K during for-loop below, so need to check here; this also ensures the while is guaranteed to terminate
	
	# Define G/R for binary phenotypes (used by fit.logistic())
	if (any(loc$binary)) { R = Q[,keep] %*% diag(1/sqrt(lambda[keep])); G = X %*% R } # NOTE: these are further subsetted within the fit.logistic() function if any PCs are dropped, so there is no need to subset again in the while loop
	
	loc$sigma = loc$h2.obs = loc$h2.latent = rep(NA, input$P)
	dropped = c()	# any PCs that might be dropped by fit.logistic() due to instability
	
	while (T) {
		keep = which(!(1:K.max %in% dropped))
		loc$K = length(keep)
		if (loc$K < min.K) { print(paste("Error: Fewer than",min.K,"PCs in locus",loc$id)); loc=NULL; return(NULL) } # K can drop below min.K during for-loop below, so need to check here; this also ensures the while is guaranteed to terminate
		
		loc$delta = matrix(NA, loc$K, input$P); # this needs to be defined in here so it updates to right size if further PCs dropped
		rerun = F # will be set to true in case additional PCs dropped during 1:P loop by fit.logistic
		
		for (i in 1:input$P) {
			# check if N < K+50
			if (loc$N[i] < loc$K + 50) { warning(paste0("N too small for phenotype '",loc$phenos[i],"' in locus ",loc$id)); next() }
			
			if (input$info$binary[i]) {
				fit = fit.logistic(G, X, R, loc$N[i], input$info$prop_cases[i], loc.sum[[i]]$STAT, loc.sum[[i]]$N, phen.id=input$info$phenotype[i], loc.id=loc$id, dropped=dropped)
				if (is.na(fit[1])) { warning(paste0("Multiple logistic regression model for phenotype '",loc$phenos[i],"' in locus ",loc$id," failed to converge (inversion error)")); next() }
				
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
					loc$h2.latent[i] = fit$h2.obs * (input$info$prevalence[i] * (1-input$info$prevalence[i]) / dnorm(qnorm(input$info$prevalence[i]))^2)
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
	if (input$P > 1) { loc$sigma = diag(loc$sigma) }; if (!is.null(input$sample.overlap)) { loc$sigma = sqrt(loc$sigma) %*% as.matrix(input$sample.overlap) %*% sqrt(loc$sigma) } # if P > 1 is just because the diag() doesn't work for single phenotype
	loc$sigma = as.matrix(loc$sigma)  # as.matrix necessary so that the univ test still works with single phenotype
	
	# cap any h2's at 0
	loc$h2.obs[loc$h2.obs<0] = 0
	loc$h2.latent[loc$h2.latent<0] = 0
	
	# get full omega
	loc$omega = t(loc$delta)%*%loc$delta / loc$K - loc$sigma
	loc$omega.cor = suppressWarnings(cov2cor(loc$omega))
	
	# check if any phenos have negative sigma or omega
	neg.var = diag(loc$sigma) < 0 | diag(loc$omega) < 0
	if (all(neg.var, na.rm=T)) { print(paste0("Error: Negative variance estimate for all phenotypes in locus ",loc$id)); loc=NULL; return(NULL) }	# print error if all had negative variance estimate
	if (any(neg.var, na.rm=T)) { warning(paste0("Negative variance estimate for phenotype(s) '",paste(loc$phenos[which(neg.var)],collapse="', '"),"' in locus ",loc$id,"; Dropping these as they cannot be analysed")) }
	
	# remove all phenotypes that failed (either due to negative variance, N < K, or wsw.inversion problem)
	failed = neg.var | is.na(neg.var)	# those that failed due to N < K or wsw.inversion will be NA in the neg.var variable
	if (any(failed)) {
		if (all(failed)) { print(paste0("Error: Processing of all phenotypes in locus ",loc$id," failed (see preceeding warning messages for details)")); loc=NULL; return(NULL) }
		
		loc$delta = as.matrix(loc$delta[,!failed])	# delta
		for (var in c("sigma","omega","omega.cor")) { loc[[var]] = as.matrix(loc[[var]][!failed,!failed]) } 		# symmetric matrices
		for (var in c("N","binary","phenos","h2.obs","h2.latent")) { loc[[var]] = loc[[var]][!failed] }				# vectors
	}
	# Add phenotype IDs to everything
	colnames(loc$delta) = names(loc$N) = names(loc$binary) = names(loc$h2.obs) = names(loc$h2.latent) = loc$phenos
	for (var in c("sigma","omega","omega.cor")) { dimnames(loc[[var]]) = rep(list(loc$phenos),2) }
	
	return(loc)
}





### Sum-stats read-in ###

format.pvalues = function(input, i) {
	input$sum.stats[[i]] = input$sum.stats[[i]][!is.na(input$sum.stats[[i]]$P),]		# remove NA pvalues
	input$sum.stats[[i]]$P[input$sum.stats[[i]]$P<1e-300] = 1e-300						# format zero p-values
}

filter.n = function(input, i) {
	input$sum.stats[[i]] = input$sum.stats[[i]][!is.na(input$sum.stats[[i]]$N),]		# remove missing N
	input$sum.stats[[i]] = subset(input$sum.stats[[i]], N > 0)							# remove negative N
}

process.sumstats = function(input) {
	# check if input files exist
	check.files.exist(input$info$filename)
	
	# possible / required headers
	headers = list(); headers$N = c("N","NMISS","N_analyzed","OBS_CT"); headers$STAT = c("Z","T","STAT","Zscore"); headers$SNP = c("SNP","ID","SNPID_UKB","SNPID","MarkerName","RSID","RSID_UKB"); headers$B = c("B","BETA"); headers$A1=c("A1","ALT"); headers$A2=c("A2","REF")	 # header variations
	
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
					warning(paste0("More than one valid ",h," header for phenotype: '",input$info$phenotype[i],"'. Only retainig the first ('",colnames(input$sum.stats[[i]])[colnames(input$sum.stats[[i]]) %in% headers[[h]]],"')."))
				}
				colnames(input$sum.stats[[i]])[colnames(input$sum.stats[[i]]) %in% headers[[h]]] = h		# rename to standard header format
			}
		}
		input$sum.stats[[i]]$SNP = tolower(input$sum.stats[[i]]$SNP)								# set SNP IDs to lower case
		if ("P" %in% colnames(input$sum.stats[[i]])) { format.pvalues(input,i) }					# format p-values
		filter.n(input,i)																			# remove missing or negative N
		
		# convert effect sizes to Z if there are none already
		if (! "STAT" %in% colnames(input$sum.stats[[i]])) {	# if no Z / T
			effect.size = c("B","OR","logOdds")
			param = effect.size[effect.size %in% colnames(input$sum.stats[[i]])][1]	# get relevant effect size parameter (takes first if multiple, e.g. both OR/logOdds)
			# check that both param + P exist
			if (is.na(param) | ! "P" %in% colnames(input$sum.stats[[i]])) {
				stop(paste0("Lack of valid statistics provided (e.g. Z, T, or BETA/OR/logOdds + P) for phenotype: '",input$info$phenotype[i],"'."))
			} else {
				# if so get the sign
				if (param=="OR") { sign = ifelse(input$sum.stats[[i]][[param]] > 1, 1, -1) } else { sign = sign(input$sum.stats[[i]][[param]]) }
			}
			if (all(sign > 1)) { warning("The signs of betas/logOdds are positive; are you sure you did not provide ORs? (if so, please use the appropriate header)") }
			
			# get Z
			Z = -qnorm(input$sum.stats[[i]]$P/2)
			input$sum.stats[[i]]$STAT = Z * sign
		}
		input$sum.stats[[i]] = input$sum.stats[[i]][c("SNP","A1","A2","STAT","N")] 		# retain relevant columns
	}
	names(input$sum.stats) = input$info$phenotype
	
	# reference data bim/afreq file
	print("...Reading in SNP info from reference data")
	input$ref = read.bim.custom(input$ref.prefix, as.env=T)
	input$ref$bim$snp.name = tolower(input$ref$bim$snp.name)		# setting ref SNPs tolower
	
	# get common SNPs
	print("...Extracting common SNPs")
	input$analysis.snps = intersect(input$ref$bim$snp.name, input$sum.stats[[1]]$SNP)	# all SNPs will be ordered according to bim file
	if (input$P>1) { for (i in 2:input$P) { input$analysis.snps = intersect(input$analysis.snps, input$sum.stats[[i]]$SNP) }}
	if (length(input$analysis.snps) < 3) { stop("Less than 3 SNPs shared across data sets; make sure you have matching SNP ID formats across sumstats / reference data sets")
	} else { print(paste("~",length(input$analysis.snps),"SNPs shared across data sets")) }
	
	if (!all(input$ref$bim$snp.name[input$ref$bim$snp.name %in% input$analysis.snps] == input$analysis.snps)) { stop("Program Error: SNPs not ordered according to reference after subsetting. Please contact developer.") }
	
	# subset sumstats to common snps
	for (i in 1:input$P) {
		input$sum.stats[[i]] = input$sum.stats[[i]][match(input$analysis.snps, input$sum.stats[[i]]$SNP),]
		if(!all(input$sum.stats[[i]]$SNP==input$analysis.snps)) { stop("Program Error: sum-stats SNPs do not match the reference after subsetting. Please contact developer.") }
	}
	
	# align SNPs
	print("...Aligning effect alleles to reference data set")
	align(input)
	
	return(input)
}

process.sample.overlap = function(sample.overlap.file, phenos) {
	sample.overlap = as.matrix(read.table(sample.overlap.file, check.names=F))	# check.names=F prevents the X in front of numeric phenotype IDs
	idx = match(phenos, colnames(sample.overlap))								# get relevant indices
	if (any(is.na(idx))) { stop(paste0("Phenotype(s) not listed in sample overlap file: '", paste0(phenos[is.na(idx)], collapse="', '"),"'")) }
	return(cov2cor(sample.overlap[idx,idx]))
}

get.input.info = function(input.info.file, sample.overlap.file, ref.prefix, phenos=NULL) {
	check.files.exist(c(input.info.file, sample.overlap.file, paste0(ref.prefix, c(".bim",".bed",".fam"))))
	
	input = new.env(parent=globalenv())
	input$info = read.table(input.info.file, header=T, stringsAsFactors=F) # read in input info file
	if (!all(c("phenotype","cases","controls","filename") %in% colnames(input$info))) { stop("Please provide an input.info file with headers: phenotype, cases, controls, filename") }  # header check on input file
	if (is.null(phenos)) { phenos = input$info$phenotype }
	if (!all(phenos %in% input$info$phenotype)) { stop(paste0("Phenotype(s) not listed in input info file: '", paste0(phenos[!phenos %in% input$info$phenotype], collapse="', '"),"'")) }
	
	input$info = input$info[match(phenos, input$info$phenotype),]			# match to phenotypes of interest
	input$info$N = input$info$cases + input$info$controls					# total N column
	input$info$prop_cases = input$info$cases / input$info$N 				# get proportion cases
	input$info$binary = input$info$prop_cases!=1							# infer binary phenotypes from prop_cases
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
#' @param sample.overlap.file Name if file with sample overlap information.
#' Can be set to NULL if there is no overlap
#' 
#' @param ref.prefix Prefix of reference genotype data in plink format (*.bim, *.bed, *.fam)
#' 
#' @param phenos A vector of phenotype IDs can be provided if only a subset of phenotypes are desired
#' (if NULL, all phenotypes in the input info file will be processed). 
#' This can be convenient if a subset from a larger number of phenotypes are analysed, as only a single input.info / sample overlap file needs to be created.
#' 
#' @return An object containing processed input data and related info
#' \itemize{
#'     \item info - the processed input info file. Columns added during processing: 
#'     \itemize{
#'        \item N = cases + controls
#'        \item prop_cases = cases / N
#'        \item binary = prop_cases != 1
#'     }
#'     \item P - number of phenotypes
#'     \item sample.overlap - sample overlap matrix
#'     \item sum.stats - processed summary statistics (SNP aligned effect sizes, subsetted to common SNPs across data sets, effect sizes converted to Z, etc)
#'     \item ref.prefix - genotype reference data prefix
#'     \item analysis.snps - subset of SNPs shared across all data sets
#'     \item unalignable.snps - SNPs removed during alignment
#'     \item ref - environment containing the genotype reference data bim file (ref$bim)
#' }
#' 
#' @export
process.input = function(input.info.file, sample.overlap.file, ref.prefix, phenos=NULL) {
	input = get.input.info(input.info.file, sample.overlap.file, ref.prefix, phenos)
	input = process.sumstats(input)
	return(input)
}


#' Read in locus file
#' 
#' This function reads in the locus file, defining the locus boundaries using either coordinates (headers: 'CHR', 'START', 'STOP') and/or ';' separated SNP lists (header: 'SNPS').
#' A locus ID column is also required (header: 'LOC').\cr
#' \cr
#' If both coordinates and SNP list columns are provided, only the SNP lists will be used for subsetting the reference data (this can be convenient if the SNP coordinates are based on a different GRChX version than the reference).\cr
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
