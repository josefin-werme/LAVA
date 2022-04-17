#### Functions used to align the effect alleles of all sum-stats to the reference data set ###

# allele pair codes
allele.pairs = list( 
	AA=NA, CC=NA, GG=NA, TT=NA,
	AC=1, CA=-1, TG=1, GT=-1,
	AG=2, GA=-2, TC=2, CT=-2,
	AT=11, TA=11, CG=12, GC=12
)

alleles = c("A", "C", "G", "T") 
pair.index = data.frame(a1=rep(alleles, times=4), a2=rep(alleles, each=4), stringsAsFactors=F) 
pair.index$code = unlist(allele.pairs[paste0(pair.index$a1, pair.index$a2)]) 
pair.index$code[pair.index$code > 10] = NA		# setting strand ambiguous alleles to NA

# take input a1 and a2 (string) vectors, output allele pair code
# this uses the alleles and pair.index variables from encompassing scope
# note that this automatically deals with non-ACTG allele codes as well, since the match() calls turn those into NA's
map.alleles = function(a1, a2) {
	# convert to numeric indices
	a1 = match(toupper(a1), alleles)
	a2 = match(toupper(a2), alleles)
	
	return(pair.index$code[a1 + (a2-1)*4]) 
	# this works because of the way the pair.index is set up, as essentially a square matrix that has been flattened into a vector (putting 2nd column below the 1st, etc.)
	# numeric a1 denotes the row index, a2 the column index
}

# takes two allele pair code vectors, returns result of comparison
# return values:
# value of 1 (or T) for matched pairs and same direction
# value of -1 for matched pairs, but different direction
# value of NA for mismatched pairs (including one or both pairs having invalid allele pair code)
com.pair = function(pair1, pair2) {
	res = rep(NA, length(pair1))
	res[abs(pair1) == abs(pair2)] = -1
	res[pair1 == pair2] = 1
	return(res)
}

align = function(input, check.index=NULL) {
	bim.snp.idx = input$ref$bim$snp.name %in% input$analysis.snps 										# index of analysis SNPs in orig bim file (these are ordered acc to the bim file)
	ref.map = map.alleles(input$ref$bim$allele.1[bim.snp.idx], input$ref$bim$allele.2[bim.snp.idx])		# allele map for reference data
  remove = NULL;
	if (is.null(check.index)) { #added check.index code (cdl 17/3)
		input$unalignable.snps = NULL
		check.index = 1:length(input$sum.stats)
	}
	
	for (i in check.index) {
		sum.map = map.alleles(input$sum.stats[[i]]$A1, input$sum.stats[[i]]$A2)				# allele map for sumstats
		input$sum.stats[[i]]$STAT = input$sum.stats[[i]]$STAT * com.pair(ref.map, sum.map)	# flip sign where relevant
		if (any(is.na(input$sum.stats[[i]]$STAT))) { remove = c(remove, which(is.na(input$sum.stats[[i]]$STAT))) }
	}
	# remove any alleles that couldn't be matched/flipped
	if (length(remove)>0) {
		input$unalignable.snps = unique(c(input$unalignable.snps, input$analysis.snps[unique(remove)]))	# store which SNPs couldnt be aligned (modified cdl 17/3)
		input$analysis.snps = input$analysis.snps[-unique(remove)]		# remove unalignable SNPs from reference data SNP list
		for (i in 1:length(input$sum.stats)) {
			input$sum.stats[[i]] = input$sum.stats[[i]][-unique(remove),]	# remove unalignable SNPs from sumstats
			if (!all(input$analysis.snps == input$sum.stats[[i]]$SNP)) { print("Program Error: SNPs do not match after alignment. Please contact developer"); stop() }
		}
		print(paste0("...Removing ",length(input$unalignable.snps)," SNPs which could not be aligned, ", length(input$analysis.snps)," remaining"))
	}
}