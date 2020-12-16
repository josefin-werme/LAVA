LAVA TUTORIAL
================
Josefin Werme, CTG Lab, VU Amsterdam
18th Dec, 2020

This tutorial shows you how to read in and analyse data with LAVA
(**L**ocal **A**nalysis of \[co\]**V**ariant **A**ssociation): A tool
developed for local genetic correlation (*r<sub>g</sub>*) analysis
**\[paper ref\]**.

LAVA can analyse the local *r<sub>g</sub>* between two or more
phenotypes, analyse both binary and continuous phenotypes, and account
for known or estimated sample overlap. It can also test the univariate
local genetic signal for all phenotypes of interest, which may be used
to filter out unassociated loci.

The tutorial will show you how to install and run LAVA using some
example input data. If you wish, you can inspect the data in the
‘vignettes/data’ folder.

-----

## Installing LAVA

To read in the genotype data, LAVA uses functions from the snpStats
packge. As this package needs to be installed using BiocManager, its
easiest to do so before installing LAVA

``` r
install.packages("BiocManager") #, repos='http://cran.rstudio.com/')
library("BiocManager")
BiocManager::install("snpStats")
```

You will also need devtools for the installation

``` r
install.packages("devtools")
library("devtools")
```

LAVA can then be installed directly from github

``` r
devtools::install_github("https://github.com/josefin-werme/lava.git")
```

..or by downloading the code and installing it from a local directory

``` r
lava.path = "~/Programs/R/lava"       # edit as necessary
devtools::install(lava.path)
```

``` r
library(lava)
```

## Input format

*NOTE: Example input files can be found in the ‘vignettes/data’
directory*

#### As input, LAVA needs the following data:

  - **Reference genotype data** in plink format (.bim, .bed, .fam), used
    for the estimation of LD
    
      - ``` 
        e.g. [1000 genomes](https://www.internationalgenome.org/data/)
        ```

  - **Input info file**, used for convenient processing of multiple
    phenotypes. Requires the columns:
    
      - *‘phenotype’*: phenotype IDs
    
      - *‘cases’*: number of cases
        
          - set to 1 (or the total sample size) for continuous
            phenotypes
    
      - *‘controls’*: number of controls
        
          - set to 0 for continuous phenotypes
              - *Info about the number of cases and controls is required
                in order to properly reconstruct the joint genetic
                effects for binary phenotypes (see **\[paper ref\]**)*
    
      - *‘prevalence’* (optional)
        
          - this is only relevant if you want an estimate of the local
            population h2 for binary phenotypes. Estimates of the sample
            h2 are still provided
    
      - *‘filename’*: filenames and paths to the relevant summary
        statistics

  - **Summary statistics** for all phenotypes of interest. To
    accommodate common summary statistics formats, a number of different
    column names are possible for the same column type (but only ONE
    should be provided; i.e do not retain multiple SNP ID columns):
    
      - ‘SNP’ / ‘ID’ / ‘SNPID\_UKB’/ ‘SNPID’ / ‘MarkerName’ / ‘RSID’:
        SNP IDs
    
      - ‘A1’ / ‘ALT’: effect allele
    
      - ‘A2’ / ‘REF’: reference allele
    
      - ‘N’ / ‘NMISS’ / ‘N\_analyzed’: number of samples
    
      - ‘Z’ / ‘T’ / ‘STAT’ / ‘Zscore’: if provided, no p-values or
        coefficients are needed; otherwise, please provide:
        
          - ‘B’ / ‘BETA’ / ‘OR’ / ‘logOdds’: effect size coefficients
            
              - please provide beta coefficients for continuous
                phenotypes, and odds ratios / log odds for binary
                phenotypes
        
          - ‘P’: p-values

  - **Locus definition file**: File that defines loci either based on
    genomic coordinates or a list of SNPs. Required headers:
    
      - ‘LOC’: locus ID
    
      - ‘CHR’, ‘START’, ‘STOP’: coordinates
    
      - ‘SNPS’: list of SNPS (optional)
        
          - note: if a SNP list is provided, coordinates will be ignored
            and the reference data will be subsetted based on SNP IDs
            instead. This can be convenient if the locus definition file
            is based on a different GRChX version than the reference
            data. If no SNP list is provided, the GRChX versions of
            reference and the locus file MUST match\!

  - **Sample overlap file** (optional)
    
      - This can be obtained using cross-trait LDSC (check out the
        ‘vignettes/get\_sample\_overlap.Rmd’ file for a walk through
        on how to do this)

## Process input

Run the script below to process the sample input files. If you wish, you
can examine the input object to get an idea of what your data looks like
after processing. You can also check out the original data in the
‘lava/vingettes/data’ folder.

``` r
### Navigate to the vignette directory
# setwd("~/../lava/vignettes")    # adapt as necessary
# setwd("~/Documents/Github/")

### Read in summary statistics and related info
input = process.input(input.info.file="data/input.info.txt",           # input info file
                      sample.overlap.file="data/sample.overlap.dat",   # sample overlap file (can be set to NULL if there is no overlap)
                      ref.prefix="data/g1000_test",                    # reference genotype data prefix
                      phenos=c("depression","neuro","bmi"))            # subset of phenotypes listed in the input info file that we want to process

# inspect the processed input data
ls(input)                     # this is actually an environment; hence ls() rather than str()
ls(input$sum.stats)   # processed summary statistics
input$info            # processed input info file, with additional variables N, prop_cases, binary computed by process.input()
input$sample.overlap  # sample overlap file
head(input$ref$bim)   # bim file from reference data
  
# read more about this function
?process.input()
```

## Read in locus definitions and proces a locus

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

  - Locus coordinates (‘chr’,‘start’,‘stop’)

  - The number of SNPs (‘K.block’) / PCs (‘K’) within the locus

  - The (estimated) PC projected joint SNP effects, \(\delta\) (‘delta’)

  - The sampling covariance \(\sigma^2\) (‘sigma’)

  - The genetic covariance matrix \(\Omega\) (‘omega’), with
    corresponding correlation matrix \(\Omega^*\)
    (‘omega.cor’):\(\Omega = t(\delta)'\delta / K-\sigma^2\)

More information about the return variables can be found in the function
documentation

``` r

ls(locus)                               # inspect locus
#>  [1] "binary"    "chr"       "delta"     "h2.latent" "h2.obs"    "id"       
#>  [7] "K"         "N"         "N.snps"    "omega"     "omega.cor" "phenos"   
#> [13] "sigma"     "snps"      "start"     "stop"
c(locus$chr, locus$start, locus$stop)   # locus coordinates
#> [1]         1 113418038 114664387
str(locus$snps)                         # locus snps
#>  chr [1:2231] "rs2360008" "rs1237670" "rs1235629" "rs1216539" "rs61819971" ...
locus$N.snps                            # N snps
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
locus$sigma
#>              depression        neuro          bmi
#> depression 8.919049e-06 1.323291e-06 1.975528e-07
#> neuro      1.323291e-06 2.645345e-06 4.440069e-08
#> bmi        1.975528e-07 4.440069e-08 2.646752e-06
locus$phenos
#>      [,1]        
#> [1,] "depression"
#> [2,] "neuro"     
#> [3,] "bmi"
locus$N
#> depression      neuro        bmi 
#>   500199.0   377979.6   377749.4

?process.locus()  # find our more details in the function documentation
```

## Perform the univariate test

If we determine the amount of local genetic signal for all phenotypes of
interest, we can filter out unassociated loci and dramatically speed up
the local r<sub>*g*</sub> analysis. This is done with the univariate
test.

After obtaining the locus object, we can pass this object directly to
the run.univ() function in order to perform the univariate test

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

Here, there was little evidence of any signal for either phenotypes, and
there will be little point in testing the local r<sub>*g*</sub>.

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
proceed with the bivariate test

## Bivariate local r<sub>*g*</sub> analysis

To proceed with the bivariane rg analysis, we simply pass the locus
object the run.bivar() function

``` r
run.bivar(locus)
#>        phen1 phen2       rho rho.lower rho.upper         r2 r2.lower r2.upper
#> 1 depression   bmi 0.0723039  -0.18273   0.32590 0.00522785        0  0.10916
#> 2      neuro   bmi 0.1069490  -0.09621   0.31253 0.01143810        0  0.09767
#>          p
#> 1 0.569780
#> 2 0.302601
```

And if we look at the description of this function, we see that this too
can take a phenotype vector, should only a subset of phenotypes be
analysed. If multiple phenotypes are entered simultaneously, the last
phenotype will be treated as the phenotype of interest, and separate
bivariate tests will be performed between this phenotype and all others.
(Note: this test is symmetric, which phenotype is considered the
‘outcome’ here doesn’t matter).

``` r
?run.bivar()
```

To analyse the bivariate test automatically conditioned on univariate
signal, you can also use the run.univ.bivar() function to perform both
tests in one go

``` r
# with the default p-value threshold of XXXX
run.univ.bivar(locus)
#> $univ
#>         phen      h2.obs h2.latent           p
#> 1 depression 0.000260236        NA 4.13399e-07
#> 2      neuro 0.000573255        NA 1.25993e-14
#> 3        bmi 0.000975242        NA 6.60429e-32
#> 
#> $bivar
#>        phen1 phen2       rho rho.lower rho.upper         r2 r2.lower r2.upper
#> 1 depression   bmi 0.0723039  -0.18525   0.32309 0.00522785        0  0.10669
#> 2      neuro   bmi 0.1069490  -0.09788   0.31424 0.01143810        0  0.09875
#>          p
#> 1 0.567976
#> 2 0.299595

# or with a custom p-value threshold
run.univ.bivar(locus)
#> $univ
#>         phen      h2.obs h2.latent           p
#> 1 depression 0.000260236        NA 4.13399e-07
#> 2      neuro 0.000573255        NA 1.25993e-14
#> 3        bmi 0.000975242        NA 6.60429e-32
#> 
#> $bivar
#>        phen1 phen2       rho rho.lower rho.upper         r2 r2.lower r2.upper
#> 1 depression   bmi 0.0723039  -0.18265   0.32806 0.00522785        0  0.10889
#> 2      neuro   bmi 0.1069490  -0.10151   0.31200 0.01143810        0  0.09735
#>          p
#> 1 0.567994
#> 2 0.301485
```

as you can see, this function only proceeds with the bivariate analysis
for all the phenotypes that reach univariate significance at the
specified threshold.

## Multivariate approaches

For loci in which you find significant correlations between more than
two phenotypes, you might want to apply one of two possible multivariate
approaches, which allow you to compute conditional local genetic
correlations

### Multiple regression

The multiple regression allows you to model the genetic signal for an
outcome phenotype using several other predictor phenotypes
simultaneously. When doing so, any potential local r<sub>*g*</sub>
between the predictor phenotypes will be accounted for, i.e., their
genetic association with the outcome will be conditional on all the
other predictors

### Partial correlation

-----

## Example analysis script for bivariate r<sub>*g*</sub> analysis across all loci

If you are interested in e.g. the local bivariate genetic correlations
between a large amount of phenotypes, we advice using e.g. a cluster
computer. The example analysis script below shows how you might set up
an R script that can be called from the command line

``` r
# command line arguments, specifying input/output file names and phenotype subset
arg = commandArgs(T); ref.dat = arg[1]; loc.file = arg[2]; info.file = arg[3]; phenos = unlist(strsplit(arg[4],";")); sample.overlap.file = arg[5]; out.fname = arg[6]
library(data.table); library(lava); print(paste("Running LAVA version",packageVersion("lava")))

### Read in data
loci = read.loci(loc.file); n.loc = nrow(loci)
input = process.input(info.file, sample.overlap.file, ref.dat, phenos)

print(paste("Starting LAVA analysis for",n.loc,"loci"))
progress = ceiling(quantile(1:n.loc, seq(.05,1,.05)))   # (if you want to print the progress)

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
