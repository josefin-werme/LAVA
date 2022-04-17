# These functions are been adapted from the read.bim() and read.plink() functions of the snpStats package
# adaptations were made to prevent the entire bim file from being read in each time a locus is processed

read.bim.custom = function(bim, as.env=T) {
  if (!grepl("\\.bim$", bim)) bim = paste0(bim, ".bim")
  df.bim = data.table::fread(bim, data.table=F)
  if (ncol(df.bim) == 6) {
    names(df.bim) <- c("chromosome", "snp.name", "cM", "position", "allele.1", "allele.2")
    df.bim$snp.name = tolower(df.bim$snp.name) # added cdl 16/3
  } else {
    warning("non-standard .bim file")
  }
  rownames(df.bim) = as.character(df.bim[, 2])

  if (as.env) {
    env = new.env(parent=globalenv()) 
    env$bim = df.bim
    return(env)
  } else {
    return(df.bim)
  }
}

# original code has been commented out, and modifications are indicated with comments
# if calling this with df.bim=NULL as default, it operates almost the same as the original read.plink() function from the snpStats package
# main difference is that it now does the rownames(df.bim) bit before potentially subsetting
#matching on select.snps is case-insensitive, and SNP names from read.bim.custom are returned lower case
read.plink.custom = function (bed, bim, fam, df.bim=NULL, na.strings = c("0","-9"), sep = ".", select.subjects = NULL, select.snps = NULL) {
    lb <- nchar(bed)
    ext <- substr(bed, lb - 3, lb)
    if (ext == ".bed") {
        stub <- substr(bed, 1, lb - 4)
    } else {
        stub <- bed
        bed <- paste(bed, ".bed", sep = "")
    }
    if (missing(bim)) bim <- paste(stub, ".bim", sep = "")
    if (missing(fam)) fam <- paste(stub, ".fam", sep = "")
    
    if (is.null(df.bim)) df.bim = read.bim.custom(bim, as.env=F) # added 
    if (is.environment(df.bim)) { # added
      snps = df.bim$bim$snp.name  # subsetting based on SNP name rather than rownames
    } else {
      snps = rownames(df.bim) # orig code
    }
    
#   df.bim <- read.table(bim, comment.char = "", as.is = TRUE, na.strings = na.strings)
#   snps <- as.character(df.bim[, 2])
    if (!is.null(select.snps)) {
        if (is.character(select.snps)) {
            select.snps <- match(tolower(select.snps), snps) # added tolower (cdl 16/3)
            if (any(is.na(select.snps))) stop("unrecognised snp selected")
        } else if (is.numeric(select.snps)) {
            tot.snps = ifelse(is.environment(df.bim), nrow(df.bim$bim), nrow(df.bim)) # added 
            if (min(select.snps) < 1 || max(select.snps) > tot.snps) stop("selected snp subscript out of range") # modification of line below           
#            if (min(select.snps) < 1 || max(select.snps) > nrow(df.bim)) stop("selected snp subscript out of range") 
            select.snps <- as.integer(select.snps)
        }
        else stop("unrecognised snp selection")
        if (any(duplicated(select.snps))) stop("snp selection contains duplicates")
        select.snps <- sort(select.snps)
        if (is.environment(df.bim)) { # added if-statement; original was just the line in the else clause
          df.bim <- df.bim$bim[select.snps, ]                  
        } else {
          df.bim <- df.bim[select.snps, ]
        }
        snps <- snps[select.snps]
    } else { # added this else clause
      if (is.environment(df.bim)) df.bim = df.bim$env # pulling the data frame out of the env if needed
    }
#    if (ncol(df.bim) == 6) { # moved these into the read.bim function
#        names(df.bim) <- c("chromosome", "snp.name", "cM", "position", "allele.1", "allele.2")
#    } else {
#        warning("non-standard .bim file")
#    }
#    rownames(df.bim) <- snps
    
    if (any(duplicated(snps))) stop("duplicated snp name(s)")
    df.fam <- read.table(fam, comment.char = "", as.is = TRUE, na.strings = na.strings)
    if (!is.null(select.subjects)) {
        if (is.numeric(select.subjects)) {
            if (min(select.snps) < 1 || max(select.snps) > nrow(df.fam))
                stop("selected subject subscript out of range")
            select.subjects <- as.integer(select.subjects)
        } else stop("unrecognised subject selection")
        if (any(duplicated(select.subjects))) stop("subject selection contains duplicates")
        select.subjects <- sort(select.subjects)
#        df.fam <- df.bim[select.subjects, ] 
        df.fam <- df.fam[select.subjects, ] # changed this from line above, clearly that's a bug
    }
    ped <- as.character(df.fam[, 1])
    mem <- as.character(df.fam[, 2])
    if (any(duplicated(ped))) {
        if (any(duplicated(mem))) {
            id <- paste(ped, mem, sep = sep)
            if (any(duplicated(id))) stop("couldn't create unique subject identifiers")
        }
        else id <- mem
    }
    else id <- ped
    names(df.fam) <- c("pedigree", "member", "father", "mother", "sex", "affected")
    rownames(df.fam) <- id
    invisible(getNamespace("snpStats")) # added
    gt <- .Call("readbed", bed, id, snps, select.subjects, select.snps, PACKAGE = "snpStats")
    list(genotypes = gt, fam = df.fam, map = df.bim)
}