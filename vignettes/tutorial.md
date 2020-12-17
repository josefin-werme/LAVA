LAVA TUTORIAL
================
Josefin Werme, CTG Lab, VU Amsterdam
2020-12-17

This tutorial shows you how to read in and analyse data with LAVA
(**L**ocal **A**nalysis of \[co\]**V**ariant **A**ssociation): A tool
developed for local genetic correlation (*r<sub>g</sub>*) analysis.

LAVA can analyse the local r<sub>*g*</sub> between two or more
phenotypes, analyse both binary and continuous phenotypes, and account
for known or estimated sample overlap. It can also test the univariate
local genetic signal for all phenotypes of interest, which may be used
to filter out unassociated loci.

The tutorial will show you how to install and run LAVA using some
example input data. If you wish, you can inspect the data in the
‘vignettes/data’ folder.

------------------------------------------------------------------------

## Installing LAVA

To read in the genotype data, LAVA uses functions from the snpStats
packge. As this package needs to be installed using BiocManager, its
easiest to do so before installing LAVA

``` r
install.packages("BiocManager")
BiocManager::install("snpStats")
```

You will also need devtools for the installation

``` r
install.packages("devtools")
```

LAVA can then be installed directly from github

``` r
devtools::install_github("https://github.com/josefin-werme/lava.git")   # adapt as necessary
```

Or by downloading the code and installing it from a local directory

``` r
devtools::install("~/Programs/R/lava")      # adapt as necessary
```

``` r
library(lava)
```

## Input format

*NOTE: Example input files can be found in the ‘vignettes/data’
directory*

#### As input, LAVA needs the following data:

-   **Reference genotype data** in plink format (.bim, .bed, .fam), used
    for the estimation of LD

    -   e.g [1000 genomes](https://www.internationalgenome.org/data/)

-   **Input info file**, used for convenient processing of multiple
    phenotypes. Requires the columns:

    -   ***‘phenotype’***: phenotype IDs

    -   ***‘cases’***: number of cases

        -   set to 1 or the total sample size for continuous phenotypes

    -   ***‘controls’:*** number of controls

        -   set to 0 for continuous phenotypes (important!)

    -   ***‘prevalence’*** (optional): the population prevalence of
        binary phenotypes

        -   this is only relevant if you want an estimate of the local
            population h2 for binary phenotypes. Estimates of the
            observed local sample h2 are still provided

    -   ***‘filename’***: paths and file names to the relevant summary
        statistics

-   **Summary statistics** for all phenotypes of interest. To
    accommodate common summary statistics formats, a number of different
    column names are possible for the same column type (but only ONE
    should be provided; i.e do not provide a file that contains multiple
    valid SNP ID columns):

    -   ***‘SNP’ / ‘ID’ / ‘SNPID\_UKB’/ ‘SNPID’ / ‘MarkerName’ /
        ‘RSID’***: SNP IDs

    -   ***‘A1’ / ‘ALT’***: effect allele

    -   ***‘A2’ / ‘REF’***: reference allele

    -   ***‘N’ / ‘NMISS’ / ‘N\_analyzed’***: number of samples

    -   ***‘Z’ / ‘T’ / ‘STAT’ / ‘Zscore’***: if provided, no p-values or
        coefficients are needed; otherwise, please provide:

        -   ***‘B’ / ‘BETA’ / ‘OR’ / ‘logOdds’***: effect size
            coefficients

            -   please provide beta coefficients for continuous
                phenotypes, and odds ratios / log odds for binary
                phenotypes

        -   ***‘P’***: p-values

-   **Locus definition file**: File that defines loci either based on
    genomic coordinates or a list of SNPs. Required headers:

    -   ***‘LOC’***: locus ID

    -   ***‘CHR’, ‘START’, ‘STOP’***: coordinates

    -   ***‘SNPS’***: list of SNPS (optional)

        -   note: if a SNP list is provided, coordinates will be ignored
            and the reference data will be subsetted based on SNP IDs
            instead. This can be convenient if the locus definition file
            is based on a different GRChX version than the reference
            data. If no SNP list is provided, the versions of reference
            and the locus file MUST match!

-   **Sample overlap file** (optional)

    -   This can be obtained using the results from cross-trait LDSC
        (check out the ‘vignettes/sample\_overlap.Rmd’ file for a walk
        through on how to do this)

## Process input

Run the script below to process the sample input files. If you wish, you
can examine the input object to get an idea of what your data looks like
after processing. You can also check out the original data in the
‘vingettes/data’ folder.

``` r
### Read in summarystatistics and related info
input = process.input(input.info.file="data/input.info.txt",           # input info file
                      sample.overlap.file="data/sample.overlap.txt",   # sample overlap file (can be set to NULL if there is no overlap)
                      ref.prefix="data/g1000_test",                    # reference genotype data prefix
                      phenos=c("depression","neuro","bmi"))       # subset of phenotypes listed in the input info file that we want to process

# inspect the processed input data
ls(input)                     # this is actually an environment; hence ls() rather than str()
ls(input$sum.stats)   # processed summary statistics
input$info            # processed input info file, with additional variables N, prop_cases, binary computed by process.input()
input$sample.overlap  # sample overlap file
head(input$ref$bim)   # bim file from reference data
  
# read more about this function
?process.input()
```

## Read in locus definitions and process a locus

Before analysing the genetic correlation at a locus, we need to convert
the marginal SNP effect sizes from GWAS to the corresponding joint
effects in order to account for the LD in the locus. Only after we’ve
done this, we can proceed with the analysis.

``` r
### Read in locus info file
loci = read.loci("data/test.loci")
head(loci)                          # inspect the locus file
#>   LOC CHR     START      STOP
#> 1 100   1 113418038 114664387
#> 2 230   2  26894103  28819510
#> 3 266   2  57952946  59251996
#> 4 374   2 191051955 193033982
#> 5 464   3  47588462  50387742
#> 6 950   6  25684630  26396200
?read.loci()                        # read more about the function and possible data formats

### Create a locus object for the first locus to prepare it for analysis
locus = process.locus(loci[1,], input)
```

If we inspect the locus object, we can see that it contains a lot of
info about the locus, such as:

-   Locus coordinates (‘chr’,‘start’,‘stop’)

-   The number of SNPs (‘n.snps’) / PCs (‘K’) within the locus

-   The (estimated) PC projected joint SNP effects, *δ* (‘delta’)

-   The sampling covariance *σ*<sup>2</sup> (‘sigma’)

-   The genetic covariance matrix
    *Ω* = *t*(*δ*)′*δ*/*K* − *σ*<sup>2</sup> (‘omega’) with
    corresponding correlation matrix *Ω*<sup>\*</sup> (‘omega.cor’)

``` r
ls(locus)                               # inspect locus
#>  [1] "binary"    "chr"       "delta"     "h2.latent" "h2.obs"    "id"       
#>  [7] "K"         "N"         "n.snps"    "omega"     "omega.cor" "phenos"   
#> [13] "sigma"     "snps"      "start"     "stop"
c(locus$chr, locus$start, locus$stop)   # locus coordinates
#> [1]         1 113418038 114664387
str(locus$snps)                         # locus snps
#>  chr [1:2231] "rs2360008" "rs1237670" "rs1235629" "rs1216539" "rs61819971" ...
locus$n.snps                            # N snps
#> [1] 2231
locus$omega                             # genetic covariance matrix
#>               depression         neuro           bmi
#> depression  1.391995e-06 -2.074891e-07 -4.335870e-07
#> neuro      -2.074891e-07  4.477136e-07 -2.674627e-07
#> bmi        -4.335870e-07 -2.674627e-07  7.443673e-07
locus$omega.cor                         # standardised genetic covariance matrix
#>            depression      neuro        bmi
#> depression  1.0000000 -0.2628309 -0.4259550
#> neuro      -0.2628309  1.0000000 -0.4633077
#> bmi        -0.4259550 -0.4633077  1.0000000
locus$phenos                            # locus phenotypes
#> [1] "depression" "neuro"      "bmi"

?process.locus()  # find our more details in the function documentation
```

## Perform the univariate test

If we determine the amount of local genetic signal for all phenotypes of
interest, we can filter out unassociated loci and dramatically speed up
the local r<sub>*g*</sub> analysis. This is done with the univariate
test.

After obtaining the locus object, we pass it directly to the run.univ()
function in order to test the local heritability for all phenotypes (or
a subset)

``` r
# for all phenotypes in the locus
run.univ(locus)
#>         phen      h2.obs h2.latent          p
#> 1 depression 8.05125e-05        NA 0.04240950
#> 2      neuro 1.16406e-04        NA 0.03154340
#> 3        bmi 1.93535e-04        NA 0.00146622

# or just a subset
run.univ(locus, phenos=c("depression","bmi"))
#>         phen      h2.obs h2.latent          p
#> 1 depression 8.05125e-05        NA 0.04240950
#> 2        bmi 1.93535e-04        NA 0.00146622
```

As there was little evidence of any genetic signal for either
phenotypes, and there will be no point in testing the local
r<sub>*g*</sub> (and if we did, the estimates would be unreliable anyway
as we would be dividing some covariance with close to zero variance)

If we check another locus however:

``` r
locus = process.locus(loci[3,], input)
run.univ(locus)
#>         phen      h2.obs h2.latent           p
#> 1 depression 0.000260236        NA 4.13399e-07
#> 2      neuro 0.000573255        NA 1.25993e-14
#> 3        bmi 0.000975242        NA 6.60429e-32
```

we see that there is plenty of signal for all phenotypes, and we can
proceed with the bivariate test!

## Perform a bivariate local r<sub>*g*</sub> analysis

To perform a bivariane r<sub>*g*</sub> analysis, testing the local
r<sub>*g*</sub> between pairs of phenotypes, we simply pass the locus
object the run.bivar() function

``` r
run.bivar(locus)
#>        phen1 phen2       rho rho.lower rho.upper         r2 r2.lower r2.upper
#> 1 depression   bmi 0.0723039  -0.18450   0.32492 0.00522785        0  0.10816
#> 2      neuro   bmi 0.1069490  -0.09658   0.31608 0.01143810        0  0.09990
#>          p
#> 1 0.568418
#> 2 0.301154
```

When multiple phenotypes are entered simultaneously, the last phenotype
will be treated as the phenotype of interest, and separate bivariate
tests will be performed between this phenotype and all others. (Note
that this test is symmetric, the order of the phenotypes doesn’t
matter!).

You can subset and change the order of the phenotypes using the ‘phenos’
argument:

``` r
run.bivar(locus, phenos=c("neuro","bmi","depression"))
#>   phen1      phen2       rho rho.lower rho.upper         r2 r2.lower r2.upper
#> 1 neuro depression 0.5278610   0.27708   0.78421 0.27863700  0.07677  0.61499
#> 2   bmi depression 0.0723039  -0.17989   0.32521 0.00522785  0.00000  0.10824
#>             p
#> 1 0.000481136
#> 2 0.566349000
```

To filter automatically based on the univariate signal, you can also use
the run.univ.bivar() function to perform both tests in one go, with the
condition that the bivariate test is only performed for the phenotypes
that reach the univariate significance:

``` r
# with the default p-value threshold of .05
run.univ.bivar(locus)
#> $univ
#>         phen      h2.obs h2.latent           p
#> 1 depression 0.000260236        NA 4.13399e-07
#> 2      neuro 0.000573255        NA 1.25993e-14
#> 3        bmi 0.000975242        NA 6.60429e-32
#> 
#> $bivar
#>        phen1 phen2       rho rho.lower rho.upper         r2 r2.lower r2.upper
#> 1 depression   bmi 0.0723039  -0.17841   0.32609 0.00522785        0  0.10878
#> 2      neuro   bmi 0.1069490  -0.09581   0.31501 0.01143810        0  0.09923
#>          p
#> 1 0.566060
#> 2 0.301103

# or with a custom p-value threshold
run.univ.bivar(locus, univ.thresh = 1e-8)
#> $univ
#>         phen      h2.obs h2.latent           p
#> 1 depression 0.000260236        NA 4.13399e-07
#> 2      neuro 0.000573255        NA 1.25993e-14
#> 3        bmi 0.000975242        NA 6.60429e-32
#> 
#> $bivar
#>   phen1 phen2      rho rho.lower rho.upper        r2 r2.lower r2.upper        p
#> 1 neuro   bmi 0.106949  -0.09628   0.31176 0.0114381        0  0.09719 0.303243
```

As you can see, this function will only proceed with the bivariate
analysis for all the phenotypes that reach univariate significance at
the specified threshold.

If you would still like to print the bivariate output for unanalysed
phentoypes, you can set the ‘return.unanalysed’ argument to TRUE.

``` r
run.univ.bivar(locus, univ.thresh=1e-8, return.unanalysed=T)
#> $univ
#>         phen      h2.obs h2.latent           p
#> 1 depression 0.000260236        NA 4.13399e-07
#> 2      neuro 0.000573255        NA 1.25993e-14
#> 3        bmi 0.000975242        NA 6.60429e-32
#> 
#> $bivar
#>        phen1 phen2      rho rho.lower rho.upper        r2 r2.lower r2.upper
#> 2 depression   bmi       NA        NA        NA        NA       NA       NA
#> 1      neuro   bmi 0.106949  -0.09986   0.30855 0.0114381        0  0.09542
#>          p
#> 2       NA
#> 1 0.299925
```

## Multivariate approaches

For loci in which you find significant correlations between more than
two phenotypes, you might want to apply one of two possible multivariate
approaches in order to compute conditional local genetic correlations
within the locus

### Multiple regression

With the multiple regression option, you can model the genetic signal
for an outcome phenotype using several other predictor phenotypes
simultaneously. When doing so, any potential local r<sub>*g*</sub>
between the predictor phenotypes will be accounted for, and their
genetic relation with the outcome will be conditional on all the other
predictors.

For this analysis we will select a locus and a phenotype subset which
has a bit more signal:

``` r
input = process.input(input.info.file="data/input.info.txt",
                      sample.overlap.file="data/sample.overlap.txt",
                      ref.prefix="data/g1000_test",
                      phenos=c("asthma","rheuma","diabetes","hypothyroidism"))
#> [1] "...Reading in sumstats"
#> [1] "...Reading in SNP info from reference data"
#> [1] "...Extracting common SNPs"
#> [1] "~ 90342 SNPs shared across data sets"
#> [1] "...Aligning effect alleles to reference data set"
#> [1] "...Removing 12980 SNPs which could not be aligned, 77362 remaining"

locus = process.locus(loci[8,], input)
```

Say that we are interested in modelling hypothyroidism as the outcome.
Lets first check that there is actually some local genetic signal for
all pentoypes, and that all predictors have a r<sub>*g*</sub> with the
outcome

``` r
?run.univ.bivar()
run.univ.bivar(locus)
#> $univ
#>             phen      h2.obs h2.latent           p
#> 1         asthma 0.000130603        NA 3.16573e-04
#> 2         rheuma 0.007015280        NA 1.59844e-62
#> 3       diabetes 0.000045775        NA 8.25651e-02
#> 4 hypothyroidism 0.000607314        NA 3.07248e-22
#> 
#> $bivar
#>    phen1          phen2        rho rho.lower rho.upper          r2 r2.lower
#> 1 asthma hypothyroidism  0.0452930  -0.33905   0.43915 0.002051460        0
#> 2 rheuma hypothyroidism -0.0199968  -0.20411   0.16027 0.000399874        0
#>   r2.upper        p
#> 1  0.20891 0.799961
#> 2  0.04558 0.827897
```

This suggests \[…\]; We can then run the multiple regression for all
predictors. Here, the last phenotype will be treated as the outcome.

``` r
run.multireg(locus)
```

Here you see that, by default, this function does not only return the
full model with all predictors, but also all intermediate models. The
reason for this is that when running a multiple regression, it is
important to be aware of how each predictor affects the fit of all the
other predictors in order to be able to interpret the results.

For example, if I analyse three predictors of which two are very
strongly associated with the outcome but also strongly correlated with
each other, and the third is only moderately associated with the outcome
but uncorrelated with any of the other predictors, the conditional
associations between the strong, but collinear predictors and the
outcome may appear weaker than the moderate, but not collinear
predictors.

If for some reason you only want to return the full model, you can speed
up the analysis by setting the ‘only.full.model’ argument to TRUE

``` r
run.multireg(locus, only.full.model=T)
```

### Partial correlation

In contrast with the multiple regression, the partial correlation allows
you to compute the local r<sub>*g*</sub> between two phenotypes of
interest, conditioned on one or more other phenotypes.

``` r
run.partial.cor()
```

------------------------------------------------------------------------

## Example analysis script for running bivariate r<sub>*g*</sub> analysis across all genomic loci

If you are interested in e.g. the local bivariate genetic correlations
between a large amount of phenotypes, we advice using a cluster
computer. The example analysis script below shows how you might set up
an R script that can be called from the command line

#### Bash script

``` bash
Rscript "g1000_test" "test.loci" "input.info.file" "sample.overlap.file" "depression;bmi" "depression.bmi"
```

#### R script

``` r
# command line arguments, specifying input/output file names and phenotype subset
arg = commandArgs(T); ref.prefix = arg[1]; loc.file = arg[2]; info.file = arg[3]; sample.overlap.file = arg[4]; phenos = unlist(strsplit(arg[5],";")); out.fname = arg[6]

library(data.table); library(lava); print(paste("Running LAVA version",packageVersion("lava")))

### Read in data
loci = read.loci(loc.file); n.loc = nrow(loci)
input = process.input(info.file, sample.overlap.file, ref.prefix, phenos)

print(paste("Starting LAVA analysis for",n.loc,"loci"))
progress = ceiling(quantile(1:n.loc, seq(.05,1,.05)))   # (if you want to print the progress)

### Analyse
u=b=list()
for (i in 1:n.loc) {
        if (i %in% progress) print(paste("..",names(progress[which(progress==i)])))     # (printing progress)
        locus = process.locus(loci[i,], input)                                          # process locus
        
        # It is possible that the locus cannot be defined for various reasons (e.g. too few SNPs, negative variance), so the !is.null(locus) check is necessary before calling the analysis functions.
        # Note that it is also possible that individual phenotypes fail, in which case only the remaining phenotypes will be retured from the process.locus() function. 
        # If you require all phenotypes to be analysed, also include the all(phenos %in% locus$phenos) condition in the if statement as well
        if (!is.null(locus) & all(phenos %in% locus$phenos)) {
                # extract some general locus info for the output
                loc.info = data.frame(locus = locus$id, chr = locus$chr, start = locus$start, stop = locus$stop, n.snps = locus$N.snps, n.pcs = locus$K)
                
                # run the univariate and bivariate tests
                loc.out = run.univ.bivar(locus, univ.thresh=1e-4)
                u[[i]] = cbind(loc.info, loc.out$univ)
                b[[i]] = cbind(loc.info, loc.out$bivar)
        }
}

# save the output
write.table(do.call(rbind,u), paste0(out.fname,".univ.lava"), row.names=F,quote=F,col.names=T)
write.table(do.call(rbind,b), paste0(out.fname,".bivar.lava"), row.names=F,quote=F,col.names=T)

print(paste0("Done! Analysis output written to ",out.fname,".*.lava"))
```
