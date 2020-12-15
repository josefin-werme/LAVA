---
title: "LAVA TUTORIAL"
author: "Josefin Werme, CTG Lab, VU Amsterdam"
date: "25th Sept, 2020"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{tutorial}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

This tutorial shows you how to read in and analyse data with LAVA (**L**ocal **A**nalysis of [co]**V**ariant **A**ssociation): A tool developed for the analysis of the local genetic correlation (*r~g~*) between traits.

LAVA can analyse the local *r~g~* between two or more phenotypes, analyse both binary and continuous phenotypes, and account for known or estimated sample overlap. In addition, it can test the univariate local genetic signal for all phenotypes of interest in order to filter out unassociated loci.

The tutorial contains example input data that we will use to demonstrate the relevant read-in and analysis functions. You can also inspect the data in the 'lava/vignettes/data' folder.

------------------------------------------------------------------------

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

```{r, results='hide', message=F}
# install.packages('devtools', repos='http://cran.rstudio.com/')
# library(devtools); devtools::install("../../lava")
library(lava)
```

## Input format

#### As input, LAVA needs the following data:

-   **Reference genotype data** in plink format (.bim, .bed, .fam), used for the estimation of LD

-   **Input info file**, containing the columns:

    -   *'*phenotype': phenotype IDs

    -   'cases': number of cases

        -   only relevant for binary phenotypes; set to 1, or the total sample size for continuous phenotypes
        -   the number of cases / controls are used to calculate the case/control ratio for binary phenotypes, which is required in order to properly reconstruct the joint genetic effects for binary phenotypes

    -   'controls': number of controls

        -   set to 0 for continuous phenotypes

    -   'filename': path to summary statistics for all phenotypes

-   **Summary statistics** for all phenotypes of interest. To accomodate common summary statistics formats, a number of different column names are possible for the same column type (but only ONE should be provided; e.g do not enter multiple SNP ID columns):

    -   *'*SNP' / 'ID' / 'SNPID\_UKB'/ 'SNPID' / 'MarkerName' / 'RSID': SNP IDs

    -   *'*A1': effect allele

    -   *'*A2': reference allele

    -   *'*N' / 'NMISS' / 'N\_analyzed': number of samples

    -   *'*Z' / 'T' / 'STAT' / 'Zscore': if provided, no p-values or coefficients are needed; otherwise, please provide:

        -   'B' / 'BETA' / 'OR' / 'logOdds': effect size coefficients.

            -   please provide beta coefficients for continuous phenotypes, and odds ratios / log odds for binary phenotypes

        -   'P': p-values

-   **Locus definition file**: File that defines loci either based on genomic coordinates or a list of SNPs. Required headers:

    -   'LOC': locus ID

    -   'CHR' + 'START' + 'STOP': coordinates

    -   'SNPS': list of SNPS (optional)

        -   note: if a SNP list is provided, coordinates will be ignored, and loci will be subsetted based on SNPs instead. This can be convenient if the locus definition file is based on a different GRChX version than the reference data.

            -   if no SNP list is provided, the GRChX versions of reference and the locus file MUST match! (versions for the summary statistics don't matter, since those are only subsetted based on SNP IDs)

-   **Sample overlap file** (if relevant)

    -   [**TODO**: add info on how to obtain this]

## Process input

Run the script below to process the sample input files. If you wish, you can examine the input object to get an idea of what your data looks like after processing. You can also check out the original data in the 'lava/vingettes/data' folder.

```{r, results='hide'}
### Input file paths
ref.fname = "data/g1000_test"                       # reference genotype data
input.info.file = "data/input.info.txt"             # input info file
sample.overlap.file = "data/overlap.all.phenos.dat" # sample overlap file (can be set to NULL if there is no overlap)

### Read in summary statistics and related info
# if only a subset of phenotypes are needed, set the phenos argument to a character vector indicating which ones
input = process.input(input.info.file, sample.overlap.file, ref.fname, 
                      phenos=c("depression","neuro","bmi"))

# inspect the processed input data
ls(input)				# this is actually an environment; hence ls() rather than str()
ls(input$sum.stats)   # processed summary statistics
input$info            # processed input info file, with additional variables N, prop_cases, binary computed by process.input()
input$sample.overlap  # sample overlap file
head(input$ref$bim)   # bim file from reference data
  
# find out more about this function
?process.input()
```

## Read in locus info and create a locus object

Before analysing the genetic correlation at a locus, we need to convert the marginal SNP effect sizes from GWAS to the corresponding joint effects in order to account for the LD in the locus. Only after we've done this, we can proceed with the analysis.

```{r}
### Read in locus info file
loci = read.loci("data/test.loci")  # read in locus file
head(loci)                          # inspect the locus file
?read.loci()                        # read more about the function and possible data formats

## Create a locus object for the first locus
locus = process.locus(loci[1,], input)
```

If we inspect the locus object, we can see that it contains a lot of info about the locus, such as:

-   Locus coordinates ('chr','start','stop')

-   The number of SNPs ('K.block') / PCs ('K') within the locus

-   The (estimated) PC projected joint SNP effects, $\delta$ ('delta'; see Methods [**paper ref**])

-   The sampling covariance $\sigma^2$ ('sigma')

-   The genetic covariance matrix $\Omega$ ('omega'), with corresponding correlation matrix $\Omega^*$ ('omega.cor'):$\Omega = t(\delta)'\delta / K-\sigma^2$ (Methods, [**paper ref**])

More information about the return variables can be found in the function documentation

```{r}
str(locus)        # inspect locus object
?process.locus()  # check out the function documentation
```

## Perform the univariate test

If we determine the amount of local genetic signal for all phenotypes of interest, we can filter out unassociated loci and dramatically speed up the local r~*g*~ analysis. This is done with the univariate test.

After obtaining the locus object, we can pass this object directly to the run.univ() function in order to perform the univariate test for all phenotypes

-   NOTE: a subset of phenotypes may be analysed by passing an ID vector to the 'pheno' argument

```{r}
run.univ(locus)
```

Here, there was little evidence of any signal for either phenotypes, and there will be little point in testing the local r~*g*~.

If we check another locus however:

```{r}
locus = process.locus(loci[3,], input)
run.univ(locus)
```

we see that there is plenty of signal for all phenotypes, and we can proceed with the bivariate test

## Bivariate local r~*g*~ analysis

To proceed with the bivariane rg analysis, we simply pass the locus object the run.bivar() function

**TODO: MAKE BIVAR FUNC THAT JUST ANALYSES ALL COMBINATIONS OF PHENOS!!**

```{r}
run.bivar(locus)

```

And if we look at the description of this function, we see that this too can take a phenotype vector, should only a subset of phenotypes be analysed. If multiple phenotypes are entered simultaneously, the last phenotype will be treated as the phenotype of interest, and separate bivariate tests will be performed between this phenotype and all others. (Note: this test is symmetric, which phenotype is considered the 'outcome' here doesn't matter).

```{r}
?run.bivar()
```

If we want to analyse the bivariate local r~*g*~ between a subset phenotypes **[...]**

## Obtaining sample overlap from LDSC
