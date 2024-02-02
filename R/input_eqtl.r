#' Load eQTL input from preprocessed input files
#'
#' Modifies the GWAS input object, appending eQTL (or sQTL) summary statistics and gene annotation to an existing GWAS input from preprocessed GTEx v8 input files (see \href{https://ctg.cncr.nl/software/lava}{LAVA website}).
#'
#' @param gwas.input GWAS input object generated using \code{\link{process.input}} (object is modified by function).
#'
#' @param eqtl.folder Folder containing the preprocessed eQTL input files.
#'
#' @param file.prefix Prefix of primary input files.
#'
#' @param tissue If set, initializes the input object for analysis of the specified tissue. Equivalent to calling
#' the \code{\link{set.tissue}} function after running process.eqtl.input().
#'
#' @param chromosomes Specifies which chromosomes to load the annotation data for (loading all chromosomes may be time-consuming).
#' Can be set to 'all' to load all chromosomes, 'auto' to load all autosomal chromosomes, or a vector of chromosome numbers or single chromosome number (the X chromosome can be specified as 'X' or 23).
#'
#' @param is.sqtl Set to true if input consists of sQTL rather than eQTL summary statistics.
#'
#' @return No return value, the gwas.input argument is modified in place. An 'eqtl' or 'sqtl' phenotype is added into the existing gwas.input structure, and the following
#' fields are added to gwas.input. Note that if the 'tissue' argument is not set, \code{\link{set.tissue}} must be called on 'gwas.input' to fully initialize the input object for analysis.
#' \itemize{
#'     \item eqtl.files - list of all the file names and folders used
#'     \item eqtl.genes - list containing the names and mapped SNPs for all genes in the eQTL data set
#'     \item eqtl.tissues - data.frame with all the available tissues
#'     \item current.tissue - name of the currently loaded tissue (NULL if not set yet)
#'     \item current.chromosomes - list of currently loaded chromosomes
#'     \item current.genes - vector of genes available for the currently loaded tissue; for sQTL data, these are gene name + splicing ID combinations rather than individual genes
#'     \item current.stats - summary statistics for the currently loaded tissue
#' }
#'
#' @export
process.eqtl.input = function(gwas.input, eqtl.folder, prefix=ifelse(is.sqtl, "gtex_sqtl_v8", "gtex_v8"), tissue=NULL, chromosomes="all", is.sqtl=F, impute.grex=F) {
	print("...Processing eQTL input info (preprocessed input)")

	if (any(gwas.input$info$phenotype == "eqtl")) {print("Error: GWAS input already contains a phenotype 'eqtl'"); return(invisible(NULL))}
	if (any(gwas.input$info$phenotype == "sqtl")) {print("Error: GWAS input already contains a phenotype 'sqtl'"); return(invisible(NULL))}

	gwas.input$eqtl.preprocessed = T

	# check main files
	check.folders.exist(c(eqtl.folder, paste0(eqtl.folder, "/annot")))
   main.files = paste0(eqtl.folder, "/", prefix, ".", c("info", "snps.gz"))
	check.files.exist(main.files)
	gwas.input$eqtl.files = list(
		tissue.file = main.files[1],
		snp.file = main.files[2],
		eqtl.dir = eqtl.folder,
		annot.prefix = paste0(eqtl.folder, "/annot/", prefix)
	)

	# load info file and check
	info = read.table(main.files[1], header=T, stringsAsFactors=F)
	cols = c("tissue", "sampsize", "filename")
	if (!all(cols %in% names(info))) {print(paste0("Error: file '", main.files[1], "' missing columns: ", paste0(cols[!(cols %in% names(info))], collapse=", "))); return(invisible(NULL))}
    
	info$filename = paste0(eqtl.folder, "/", info$filename)
	has.tissue = file.exists(info$filename)
	if (!all(has.tissue)) {
		if (all(!has.tissue)) {print(paste0("Error: all summary statistics files specified in '", main.files[1], "' are missing")); return(invisible(NULL))}
		msg = ifelse(
			mean(has.tissue) > 0.5,
      paste0("missing tissues: ", paste0(info$tissue[!has.tissue], collapse=", ")),
			paste0("available tissues: ", paste0(info$tissue[has.tissue], collapse=", "))
		)
		print(paste0("Warning: summary statistics files not available for all tissues; ", msg))
		info = info[has.tissue,]
	}

 	# load SNP info file and check
	print("...Reading in SNP information")
   sumstats = data.table::fread(main.files[2],data.table=F)
	cols = c("SNP", "A1", "A2")
	if (!all(cols %in% names(sumstats))) {print(paste0("Error: file '", main.files[2], "' missing columns: ", paste0(cols[!(cols %in% names(sumstats))], collapse=", "))); return(invisible(NULL))}
	sumstats$STAT = 1; sumstats$N = NA

	# load gene annotation file into gwas.input (inserted as gwas.input$eqtl.genes)
	print("...Reading in gene annotation file")
	set.chromosomes(gwas.input, chromosomes);
	
	# insert eQTL as dummy phenotype into
   print("...Merging eQTL data into existing GWAS input")
	type.name = ifelse(is.sqtl, "sqtl", "eqtl")
	#gwas.input$P = gwas.input$P + 1; add.index = gwas.input$P
	gwas.input$info[1,] = rep(NA, ncol(gwas.input$info))
	gwas.input$info = gwas.input$info[1,]
	gwas.input$info$phenotype = type.name
	gwas.input$info$binary = FALSE
	gwas.input$info$var.type = type.name

	if (!is.null(gwas.input$sample.overlap)) {
		gwas.input$sample.overlap = cbind(rbind(gwas.input$sample.overlap, 0), 0)
		gwas.input$sample.overlap[add.index,add.index] = 1
	}
	
	# TODO: CHECK what does add.index do in harmonize SNPs????
	# it just selectively checks only the specified indices, so indeed not necessary here
	gwas.input$sum.stats[[type.name]] = sumstats
   #harmonize.snps(gwas.input, add.index) # subset all data sets to shared SNPs
   harmonize.snps(gwas.input) # subset all data sets to shared SNPs
   
#	align(gwas.input, add.index) # this multiplies the dummy constant 1 STAT column by -1 in case of misalignment, storing this in separate DIR column for later use
	align(gwas.input) # this multiplies the dummy constant 1 STAT column by -1 in case of misalignment, storing this in separate DIR column for later use
   gwas.input$sum.stats[[type.name]]$DIR = gwas.input$sum.stats[[type.name]]$STAT
   gwas.input$sum.stats[[type.name]]$STAT = NA
   # TODO: check why these used add.idx instead of type.name (which is used in gwas.input$sumstats above)

	gwas.input$eqtl.tissues = info
	gwas.input$eqtl.type = type.name
	gwas.input$current.tissue = NULL
	gwas.input$current.genes = NULL
	gwas.input$current.stats = NULL

	if (!is.null(tissue)) set.tissue(gwas.input, tissue)
}


#' Process eQTL input from text input files
#'
#' Modifies the GWAS input object, appending eQTL (or sQTL) summary statistics and gene annotation to an existing GWAS input from plain-text input files. Input files should follow the same format
#' as summary statistics input files for other phenotypes, with an additional column GENE to indicate which gene the association results on each line belong to. Input files must be provided separately
#' per chromosome.
#'
#' @param gwas.input GWAS input object generated using \code{\link{process.input}} (object is modified by function).
#'
#' @param eqtl.file.spec File name specifier for eQTL summary statistics input files, with a placeholder [CHR] to indicate the chromosome number in the file name; eg. "my_folder/eqtl.chr[CHR].stats"
#'
#' @param chromosomes Specifies which chromosomes to load the data for (loading all chromosomes may be time-consuming).
#' Can be set to 'all' to load all chromosomes, 'auto' to load all autosomal chromosomes, or a vector of chromosome numbers or single chromosome number (the X chromosome can be specified as 'X' or 23).
#'
#' @param sample.size If set, used as sample size value for each SNP in the input files instead of N column.
#'
#' @return No return value, the gwas.input argument is modified in place. An 'eqtl' phenotype is added into the existing gwas.input structure, and the following
#' fields are added to gwas.input.
#' \itemize{
#'     \item current.chromosomes - list of currently loaded chromosomes
#'     \item current.genes - vector of available genes
#'     \item current.stats - summary statistics for all loaded genes
#' }
#'
#' @export
process.eqtl.text.input = function(gwas.input, eqtl.file.spec, chromosomes="all", sample.size=NULL) {
	print("...Processing eQTL input info (text input)")

	if (any(gwas.input$info$phenotype == "eqtl")) {print("Error: GWAS input already contains a phenotype 'eqtl'"); return(invisible(NULL))}
	if (any(gwas.input$info$phenotype == "sqtl")) {print("Error: GWAS input already contains a phenotype 'sqtl'"); return(invisible(NULL))}

	gwas.input$eqtl.preprocessed = F

	if (!grepl("[CHR]", eqtl.file.spec, fixed=T)) {print("Error: argument 'eqtl.file.spec' does not contain [CHR] placeholder token"); return(invisible(NULL))}
	gwas.input$eqtl.file.specifier = eqtl.file.spec

	if (!is.null(sample.size) && !(is.numeric(sample.size) && is.finite(sample.size) && sample.size > 0)) {print("Error: argument 'sample.size' is negative or not a number"); return(invisible(NULL))}
	gwas.input$eqtl.N = sample.size

	# insert eQTL as dummy phenotype into
	type.name = "eqtl"
	gwas.input$P = gwas.input$P + 1; add.index = gwas.input$P
	gwas.input$info = rbind(gwas.input$info, NA)
	gwas.input$info$phenotype[add.index] = type.name
	gwas.input$info$binary[add.index] = FALSE
	gwas.input$info$var.type = ""
	gwas.input$info$var.type[1:gwas.input$P == add.index] = type.name

	if (!is.null(gwas.input$sample.overlap)) {
		gwas.input$sample.overlap = cbind(rbind(gwas.input$sample.overlap, 0), 0)
		gwas.input$sample.overlap[add.index,add.index] = 1
	}

	gwas.input$eqtl.type = type.name

	# create dummy sum.stats data.frame
	gwas.input$sum.stats[[type.name]] = data.frame(SNP=gwas.input$sum.stats[[1]]$SNP, STAT=NA, N=NA, stringsAsFactors=F)

	set.chromosomes(gwas.input, chromosomes);
}



#' Initialize specific tissue for eQTL input object
#'
#' Modifies the eQTL input object prepared by the process.eqtl.input() function, loading and initializing summary statistics for the specified tissue. Removes data loaded
#' for previously loaded tissue (if any).
#'
#' @param gwas.input GWAS input object preprocessed using \code{\link{process.eqtl.input}} (object is modified by function).
#'
#' @param tissue Tissue to load summary statistics for.
#'
#' @return No return value, the gwas.input argument is modified in place. The following fields are changed
#' \itemize{
#'     \item current.tissue - name of the currently loaded tissue
#'     \item current.genes - vector of genes available for the currently loaded tissue
#'     \item current.stats - summary statistics for the currently loaded tissue
#' }
#'
#' @export
set.tissue = function(gwas.input, tissue) {
	print(paste0("...Loading tissue '", tissue, "'"))
	if (!is.null(gwas.input$eqtl.preprocessed) && !gwas.input$eqtl.preprocessed) {print("Error: eQTL data loaded into gwas.input does not support changing tissues"); return(invisible(NULL))}
	if (!all(c("eqtl.genes", "eqtl.type") %in% names(gwas.input)) || !any(gwas.input$info$var.type == gwas.input$eqtl.type)) {print("Error: gwas.input is not configured for eQTL analysis"); return(invisible(NULL))}
	if (!tissue %in% gwas.input$eqtl.tissues$tissue) {print(paste0("Error: tissue '", tissue, "' is not available in input")); return(invisible(NULL))}
	if (!is.null(gwas.input$current.tissue) && gwas.input$current.tissue == tissue) {print(paste0("Tissue '", tissue, "' already loaded")); return(invisible(NULL))}

	#load summary stats and check
	print("...Reading in summary statistics")
	current = gwas.input$eqtl.tissues[gwas.input$eqtl.tissues == tissue,]
	stats = strsplit(scan(current$filename[1], what="", sep="\n", quiet=TRUE), " ")

	is.sqtl = (gwas.input$eqtl.type == "sqtl")
	if (!all(sapply(stats, length) == ifelse(is.sqtl, 3, 2))) {print(paste0("Error: file '", current$filename[1], "' is improperly formatted")); return(invisible(NULL))}

	make.label = function(s) {ifelse(is.sqtl, paste0(s[1:2], collapse="__"), s[1])}
	genes = list(
		gene.name=sapply(stats, head, 1),
		label=sapply(stats, make.label),
		zstat=sapply(stats, tail, 1)
	)

	#insert into gwas.input object
	gwas.input$sum.stats[[gwas.input$eqtl.type]]$N = current$sampsize[1]
	gwas.input$current.tissue = tissue
	gwas.input$current.genes = genes$label[genes$gene.name %in% gwas.input$eqtl.genes$name] #lists available genes only
	gwas.input$current.stats = genes
}

#' Load gene annotation for eQTL input object
#'
#' Load the eQTL gene annotation for the specified chromosomes into the gwas.input object. Previously loaded gene annotation is replaced
#'
#' @param gwas.input GWAS input object loaded using \code{\link{process.eqtl.input}} or \code{\link{process.eqtl.text.input}} (object is modified by function).
#'
#' @param chromosomes Specifies which chromosomes to load the annotation data for (loading all chromosomes may be time-consuming).
#' Can be set to 'all' to load all chromosomes, 'auto' to load all autosomal chromosomes, or a vector of chromosome numbers or single chromosome number (the X chromosome can be specified as 'X' or 23).
#'
#' @return No return value, the gwas.input argument is modified in place. The following fields are changed
#' \itemize{
#'     \item eqtl.genes - list containing the names and mapped SNPs for all genes in the eQTL data set
#'     \item current.chromosomes - list of currently loaded chromosomes
#'     \item current.genes - vector of genes available for the currently loaded tissue
#' }
#'
#' @export
set.chromosomes = function(gwas.input, chromosomes) {
	if (is.null(gwas.input$eqtl.preprocessed)) {print("Error: gwas.input is not configured for eQTL analysis"); return(invisible(NULL))}

  if (is.null(chromosomes) || chromosomes[1] == "all") {
		chromosomes = c(1:22, "X")
	} else if (chromosomes[1] == "auto") {
		chromosomes = 1:22
	} else {
		chr.all = c(1:22, "X")
		chromosomes[chromosomes %in% c(23, "x")] = "X"
		if (any(!(chromosomes %in% chr.all))) {print(paste0("Error: chromosome codes not recognized: ", paste0(chromosomes[!(chromosomes %in% chr.all)], collapse=", "))); return(invisible(NULL))}
		chromosomes = chr.all[chr.all %in% chromosomes]
	}
	chromosomes = as.character(chromosomes)

	if (gwas.input$eqtl.preprocessed) {
		file.names = paste0(gwas.input$eqtl.files$annot.prefix, ".chr", chromosomes, ".genes.annot.gz")
		check.files.exist(file.names)

		genes = NULL
		for (i in seq_along(chromosomes)) {
			annot = strsplit(scan(file.names[i], what="", sep="\n", quiet=TRUE), " ")
			if (!all(sapply(annot, length) == 2)) {print(paste0("Error: file '", file.names[i], "' is improperly formatted")); return(invisible(NULL))}
			genes = rbind(
	  		genes,
				data.frame(
					name=sapply(annot, head, 1),
					chr=chromosomes[i],
					snps=I(lapply(annot, tail, 1)), #keep this as list, replaced by vector entries later
					processed=rep(F,length(annot)),
					stringsAsFactors=F
				)
			)
		}

    if (!is.null(gwas.input$current.tissue)) gwas.input$current.genes = gwas.input$current.stats$label[gwas.input$current.stats$gene.name %in% genes$name] #reset list of available genes

	} else {
		file.names = sapply(chromosomes, function(chr) {gsub("\\[CHR\\]", chr, gwas.input$eqtl.file.specifier)})
		check.files.exist(file.names)
		genes = NULL

		for (i in seq_along(chromosomes)) {
			chr = chromosomes[i]
			gwas.input$current.stats = list()

			sum.stats = read.sumstats.file(file.names[i], "eqtl", N.override=gwas.input$eqtl.N)
			sum.stats = sum.stats[sum.stats$SNP %in% gwas.input$analysis.snps,] #subset to SNPs available for other phenotypes

			# align to reference data, remove ambiguous / non-matched SNPs
			bim.index = match(sum.stats$SNP, gwas.input$ref$bim$snp.name)
			bim.map = map.alleles(gwas.input$ref$bim$allele.1[bim.index], gwas.input$ref$bim$allele.2[bim.index])
			alignment = com.pair(bim.map, map.alleles(sum.stats$A1, sum.stats$A2))
			sum.stats$STAT = sum.stats$STAT * alignment
			sum.stats = sum.stats[!is.na(sum.stats$STAT),]
			gwas.input$current.stats[[chr]] = sum.stats

			# add to gene index
			genes = rbind(genes, data.frame(GENE=unique(sum.stats$GENE), CHR=chr, stringsAsFactors=F))
		}

		gwas.input$current.genes = genes[,1]
		duplicates = duplicated(gwas.input$current.genes)
    if (any(duplicates)) {
			duplicates = paste(unique(gwas.input$current.genes[duplicates]), collapse=", ")
			print(paste0("Error: following gene IDs encountered for multiple chromosomes: ", duplicates))
			return(invisible(NULL))
		}
	}

	gwas.input$eqtl.genes = genes
	gwas.input$current.chromosomes = chromosomes
}




#' Process eQTL locus info
#'
#' Same as \code{\link{process.locus}}, but for use with eQTL data. Processes the summary statistics for the eQTL data and
#' all phenotypes within the specified genes and computes all parameters necessary for analysis.
#'
#' @param gene Gene to analyse, from the available genes in eqtl.input$current.genes. Can either be a gene name, or a numeric index.
#' @param eqtl.input Input object created with the \code{\link{process.input}} function and modified with \code{\link{process.eqtl.input}}, containing
#' relevant summary statistics and related info (e.g. sample overlap, case/control ratio)
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
process.eqtl.locus = function(gene, eqtl.input, phenos=NULL, min.K=2, prune.thresh=99, max.prop.K=0.75, drop.failed=T, impute.grex=T) {
	if (!all(c("eqtl.genes", "eqtl.type") %in% names(eqtl.input)) || !any(eqtl.input$info$var.type == eqtl.input$eqtl.type)) {print("Error: eqtl.input is not configured for eQTL analysis"); return(invisible(NULL))}
	if (eqtl.input$eqtl.preprocessed && is.null(eqtl.input$current.tissue)) {print("Error: no tissue initialized; run set.tissue function to load tissue data"); return(NULL)}

	if (is.character(gene)) {
		if (!(any(eqtl.input$current.genes == gene))) {print(paste0("Error: gene '", gene, "' not available", if (eqtl.input$eqtl.preprocessed) paste0(" for tissue '", tissue, "'") )); return(NULL)}
	} else {
		gene = suppressWarnings(as.integer(gene))
		if (is.na(gene) || gene <= 0 || gene > length(eqtl.input$current.genes)) {print("Error: specified gene index is non-numeric or out of bounds"); return(NULL)}
		gene = eqtl.input$current.genes[gene]
	}

	if (eqtl.input$eqtl.preprocessed) {
		gene.index.curr = which(eqtl.input$current.stats$label == gene) # offset of input gene in tissue sumstats
		gene.index.main = which(eqtl.input$eqtl.genes$name == eqtl.input$current.stats$gene.name[gene.index.curr]) # offset of input gene in main annotation

		# convert SNP-gene mapping to indices
		if (!eqtl.input$eqtl.genes$processed[gene.index.main]) {
			eqtl.input$eqtl.genes$snps[[gene.index.main]] = match(strsplit(eqtl.input$eqtl.genes$snps[[gene.index.main]], ";")[[1]], eqtl.input$sum.stats[[eqtl.input$eqtl.type]]$SNP)
			eqtl.input$eqtl.genes$processed[gene.index.main] = TRUE
		}
		snp.index = eqtl.input$eqtl.genes$snps[[gene.index.main]] #offsets of SNPs in the gene in sum.stats data.frame (includes NAs)

		# load summary statistics and define locus
		zstat = as.numeric(strsplit(eqtl.input$current.stats$zstat[[gene.index.curr]], ";")[[1]])
		if (length(zstat) != length(snp.index)) {print(paste0("Error: summary statistics for gene '", gene, "' are not consistent with gene annotation")); return(NULL)}
		zstat = zstat * eqtl.input$sum.stats[[eqtl.input$eqtl.type]]$DIR[snp.index]
		valid = !is.na(zstat) # NA's in snp.index will be tagged as invalid here as well
		eqtl.input$sum.stats[[eqtl.input$eqtl.type]]$STAT[snp.index[valid]] = zstat[valid]

		locus = list(LOC=gene, CHR=eqtl.input$eqtl.genes$chr[gene.index.main], SNPS=eqtl.input$sum.stats[[eqtl.input$eqtl.type]]$SNP[snp.index[valid]])
	} else {
  	chr = eqtl.input$eqtl.genes$CHR[eqtl.input$eqtl.genes$GENE == gene]
		curr.stats = eqtl.input$current.stats[[chr]][eqtl.input$current.stats[[chr]]$GENE == gene,]
		stat.index = match(curr.stats$SNP, eqtl.input$sum.stats[[eqtl.input$eqtl.type]]$SNP)
		eqtl.input$sum.stats[[eqtl.input$eqtl.type]]$STAT[stat.index] = curr.stats$STAT
		eqtl.input$sum.stats[[eqtl.input$eqtl.type]]$N[stat.index] = curr.stats$N

  	locus = list(LOC=gene, CHR=chr, SNPS=curr.stats$SNP)
	}

	class(locus) = "gene"

	if (!is.null(phenos)) phenos = unique(c(phenos, eqtl.input$eqtl.type))
	return(process.locus(locus, eqtl.input, phenos, min.K, prune.thresh, max.prop.K, drop.failed, impute.grex))
}


