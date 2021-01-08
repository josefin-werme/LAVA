LAVA TUTORIAL
================
Josefin Werme, CTG Lab, VU Amsterdam
2021-01-08

This tutorial shows you how to read in and analyse data with LAVA
(**L**ocal **A**nalysis of \[co\]**V**ariant **A**ssociation): A tool
developed for local genetic correlation (*r*<sub>g</sub>) analysis.

LAVA can analyse the standard bivariate local *r*<sub>g</sub> between
two phenotypes (binary as well as continuous), and account for known or
estimated sample overlap. It can also test the univariate local genetic
signal for all phenotypes of interest (i.e. the local *h*<sup>2</sup>),
which may be used to filter out non-associated loci. In addition, it can
model the local genetic relations between multiple phenotypes
simultaneously using two possible conditional models: partial
correlation and multiple regression (for more details, see the
[preprint](https://www.biorxiv.org/content/10.1101/2020.12.31.424652v1)).

The tutorial will show you how to install and run LAVA using some
example input data. If you wish, you can inspect the data in the
‘vignettes/data’ folder.

------------------------------------------------------------------------

## Installing LAVA

To read in the genotype data, LAVA uses functions from the snpStats
packge. As this package needs to be installed using BiocManager, you
should do so before installing LAVA

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
devtools::install_github("https://github.com/josefin-werme/lava.git")
```

Or by downloading the code (Download Zip under the Code button on
GitHub), unzipping and installing from disk

``` r
devtools::install("~/Programs/R/lava")      # specify local path here
```

``` r
library(lava)
#> Loading LAVA (v0.0.6)
```

## Input format

*NOTE: Example input files can be found in the ‘vignettes/data’
directory*

#### As input, LAVA needs the following data:

-   **Reference genotype data** in plink format (.bim, .bed, .fam), used
    for the estimation of LD

    -   e.g. [1000 genomes](https://www.internationalgenome.org/data/)
        (pre-processed input files can be found
        [here](https://ctg.cncr.nl/software/magma))

-   **Input info file**, used for convenient processing of multiple
    phenotypes. Requires the columns:

    -   ***phenotype***: phenotype IDs

    -   ***cases***: number of cases

        -   set to 1 for continuous phenotypes

    -   ***controls:*** number of controls

        -   set to 0 for continuous phenotypes (important!); the case /
            control ratio is used to distinguish binary and continuous
            phenotypes, and is required for appropriate processing of
            binary phenotypes.

    -   ***prevalence*** (optional): the population prevalence of binary
        phenotypes

        -   this is only relevant if you want an estimate of the local
            population *h*<sup>2</sup> for binary phenotypes. Estimates
            of the observed local sample *h*<sup>2</sup> are still
            provided

    -   ***filename***: paths and file names to the relevant summary
        statistics

-   **Summary statistics** for all phenotypes of interest. To
    accommodate common summary statistics formats, a number of different
    column names are possible for the same column type (but **only ONE**
    header should be provided for each category; i.e do not provide a
    file that contains multiple valid SNP ID columns):

    -   ***SNP / ID / SNPID\_UKB/ SNPID / MarkerName / RSID***: SNP IDs

    -   ***A1 / ALT***: effect allele

    -   ***A2 / REF***: reference allele

    -   ***N / NMISS / N\_analyzed***: number of samples

    -   ***Z / T / STAT / Zscore***: if provided, no p-values or
        coefficients are needed; otherwise, please provide both:

        -   ***B / BETA / OR / logOdds***: effect size coefficients

        -   ***P***: p-values

-   **Locus definition file**: File that defines loci either based on
    genomic coordinates or a list of SNPs. Required headers:

    -   ***LOC***: locus ID

    -   ***CHR, START, STOP***: coordinates (i.e. basepair positions)

    -   ***SNPS***: list of SNPS. ‘;’ separated, no white-space
        (optional)

        -   note: if a SNP list is provided, coordinates will be ignored
            and the reference data will be subsetted based on SNP IDs
            instead. This can be convenient if the locus definition file
            is based on a different GRChX version than the reference
            data. If no SNP list is provided, the versions of reference
            and the locus file MUST match!

-   **Sample overlap file** (optional)

    -   This can be obtained using the results from cross-trait LDSC
        (see the ‘vignettes/sample\_overlap.Rmd’ file for a walk through
        on how to do this)

## Process input

Run the script below to process the sample input files. If you wish, you
can examine the input object to get an idea of what your data looks like
after processing. You can also check out the original data in the
‘vignettes/data’ folder.

``` r
### Read in summary statistics and related info
input = process.input(input.info.file="vignettes/data/input.info.txt",           # input info file
                      sample.overlap.file="vignettes/data/sample.overlap.txt",   # sample overlap file (can be set to NULL if there is no overlap)
                      ref.prefix="vignettes/data/g1000_test",                    # reference genotype data prefix
                      phenos=c("depression","neuro","bmi"))       # subset of phenotypes listed in the input info file that we want to process

# inspect the processed input data
ls(input)                     # this is actually an environment; hence ls() rather than str()
ls(input$sum.stats)   # processed summary statistics
head(input$sum.stats$bmi)
input$info            # processed input info file, with additional variables N, prop_cases, binary computed by process.input()
input$sample.overlap  # sample overlap file
head(input$ref$bim)   # bim file from reference data
  
# read more about this function
?process.input()
```

As you can see in the input$info data frame, additional columns have
been added to indicate whether phenotypes are binary, and what the the
proportion of cases are for binary phenotypes (computed from the cases /
controls columns). *(NOTE: the N column here is not used for the
analysis, that currently needs to be provided via the summary
statistics).*

## Prepare a locus for analysis

Before analysing the genetic correlation at a locus, we need to convert
the marginal GWAS SNP effects within the locus to their corresponding
joint effects (in order to account for the LD between SNPs), and
estimate certain key parameter that are used for the analysis.

``` r
### Read in locus info file
loci = read.loci("vignettes/data/test.loci")
head(loci)  # inspect the locus file
#>   LOC CHR     START      STOP
#> 1 100   1 113418038 114664387
#> 2 230   2  26894103  28819510
#> 3 266   2  57952946  59251996
#> 4 374   2 191051955 193033982
#> 5 464   3  47588462  50387742
#> 6 950   6  25684630  26396200
?read.loci  # read more about the function and possible data formats

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
c(locus$chr, locus$start, locus$stop)   # locus coordinates
str(locus$snps)                         # locus snps
locus$n.snps                            # N snps
locus$omega                             # genetic covariance matrix
locus$omega.cor                         # standardised genetic covariance matrix
locus$phenos                            # locus phenotypes

?process.locus  # find our more details in the function documentation
```

## Perform the univariate test

To determine the amount of local genetic signal for all phenotypes of
interest and filter out non-associated loci, we use the univariate test.

After obtaining the locus object, we pass it directly to the run.univ()
function in order to test the heritability within that locus for all
phenotypes (or a subset)

``` r
# for all phenotypes in the locus
run.univ(locus)
#>         phen      h2.obs          p
#> 1 depression 8.05125e-05 0.04240950
#> 2      neuro 1.16406e-04 0.03154340
#> 3        bmi 1.93535e-04 0.00146622

# or just a subset
run.univ(locus, phenos=c("depression","bmi"))
#>         phen      h2.obs          p
#> 1 depression 8.05125e-05 0.04240950
#> 2        bmi 1.93535e-04 0.00146622
```

As there was little evidence of any genetic signal for either
phenotypes, there will be no point in testing the local *r*<sub>g</sub>
(if we did, the estimates would be unreliable and uninterpretable)..

If we check another locus however:

``` r
locus = process.locus(loci[3,], input)
run.univ(locus)
#>         phen      h2.obs           p
#> 1 depression 0.000260236 4.13399e-07
#> 2      neuro 0.000573255 1.25993e-14
#> 3        bmi 0.000975242 6.60429e-32
```

we see that there is plenty of signal for all phenotypes, and we can
proceed with the bivariate test!

## Analyse the bivariate local *r*<sub>g</sub> between pairs of phenotypes

To perform the bivariate test and obtain the local *r*<sub>g</sub>
between pairs of phenotypes, we simply pass the locus object the
run.bivar() function

``` r
run.bivar(locus)
#>        phen1 phen2       rho rho.lower rho.upper         r2 r2.lower r2.upper
#> 1 depression neuro 0.5278610   0.27394   0.78031 0.27863700  0.07504  0.60888
#> 2 depression   bmi 0.0723039  -0.18270   0.32772 0.00522785  0.00000  0.10956
#> 3      neuro   bmi 0.1069490  -0.09711   0.31301 0.01143810  0.00000  0.09820
#>             p
#> 1 0.000502283
#> 2 0.567621000
#> 3 0.300836000
```

When multiple phenotypes are entered simultaneously, bivariate
*r*<sub>g</sub>‘s will be computed between all unique phenotype pairs by
default. If you want to analyse only a subset of phenotypes, you can use
the ’phenos’ argument

``` r
run.bivar(locus, phenos=c("neuro","depression"))
#>   phen1      phen2      rho rho.lower rho.upper       r2 r2.lower r2.upper
#> 1 neuro depression 0.527861   0.27603   0.77648 0.278637  0.07619  0.60293
#>             p
#> 1 0.000485853
```

Additionally, should you only be interested in the *r*<sub>g</sub>‘s
between one target phenotype and all others, you can use the ’target’
argument to prevent all pairwise tests from being computed and focus
only on those relevant to the target

``` r
run.bivar(locus, target="bmi")
#>        phen1 phen2       rho rho.lower rho.upper         r2 r2.lower r2.upper
#> 1 depression   bmi 0.0723039  -0.18811   0.32174 0.00522785        0  0.10588
#> 2      neuro   bmi 0.1069490  -0.09725   0.31273 0.01143810        0  0.09795
#>          p
#> 1 0.565684
#> 2 0.301970
# ?run.bivar  # check the function description for more options
```

*\[ Note that as we use a simulation procedure to obtain the p-values
and confidence intervals here, the results might differ slightly between
different runs \]*

To **filter automatically based on the univariate signal**, you can also
use the run.univ.bivar() function. This will perform both tests at once,
with the condition that the bivariate test is only performed for the
phenotypes that reach the desired univariate significance threshold

``` r
# with the default p-value threshold of .05
run.univ.bivar(locus)
#> $univ
#>         phen      h2.obs           p
#> 1 depression 0.000260236 4.13399e-07
#> 2      neuro 0.000573255 1.25993e-14
#> 3        bmi 0.000975242 6.60429e-32
#> 
#> $bivar
#>        phen1 phen2       rho rho.lower rho.upper         r2 r2.lower r2.upper
#> 1 depression neuro 0.5278610   0.27516   0.77893 0.27863700  0.07571  0.60673
#> 2 depression   bmi 0.0723039  -0.18733   0.32510 0.00522785  0.00000  0.10755
#> 3      neuro   bmi 0.1069490  -0.09721   0.31441 0.01143810  0.00000  0.09885
#>             p
#> 1 0.000504686
#> 2 0.565774000
#> 3 0.303979000

# or with a custom p-value threshold
run.univ.bivar(locus, univ.thresh = 1e-8)
#> $univ
#>         phen      h2.obs           p
#> 1 depression 0.000260236 4.13399e-07
#> 2      neuro 0.000573255 1.25993e-14
#> 3        bmi 0.000975242 6.60429e-32
#> 
#> $bivar
#>   phen1 phen2      rho rho.lower rho.upper        r2 r2.lower r2.upper        p
#> 1 neuro   bmi 0.106949  -0.09058   0.31889 0.0114381        0  0.10169 0.301688
```

*\[ You can also use the ‘phenos’ and ‘target’ arguments here to subset
and specify target phenotypes \]*

## Analyse conditional genetic relations between several phenotypes

For loci in which you find significant correlations between more than
two phenotypes, you might want to apply one of two possible multivariate
approaches in order to compute conditional local genetic relations
within the locus

### Multiple regression

With the multiple regression option, you can model the genetic signal
for an outcome phenotype using several other predictor phenotypes
simultaneously. When doing so, any potential local *r*<sub>g</sub>
between the predictor phenotypes will be accounted for, and their
relation with the outcome will be conditional on all the other
predictors.

For this analysis we will select a locus and a phenotype subset which
has a bit more signal:

``` r
input = process.input(input.info.file="vignettes/data/input.info.txt",
                      sample.overlap.file="vignettes/data/sample.overlap.txt",
                      ref.prefix="vignettes/data/g1000_test",
                      phenos=c("asthma","rheuma","diabetes","hypothyroidism"))
#> [1] "...Reading in sumstats"
#> [1] "...Reading in SNP info from reference data"
#> [1] "...Extracting common SNPs"
#> [1] "~ 90342 SNPs shared across data sets"
#> [1] "...Aligning effect alleles to reference data set"
#> [1] "...Removing 12980 SNPs which could not be aligned, 77362 remaining"

locus = process.locus(loci[19,], input)
```

Say that we are interested in modelling hypothyroidism as the outcome.
Lets first check that there is actually some local genetic signal for
all phenotypes, and that all predictors are genetically correlated with
the outcome within this locus

``` r
run.univ.bivar(locus, target="hypothyroidism") 
#> $univ
#>             phen      h2.obs            p
#> 1         asthma 0.001778430 4.66411e-148
#> 2         rheuma 0.054520700  0.00000e+00
#> 3       diabetes 0.000661558  1.54348e-50
#> 4 hypothyroidism 0.001717760  6.33045e-84
#> 
#> $bivar
#>      phen1          phen2      rho rho.lower rho.upper       r2 r2.lower
#> 1   asthma hypothyroidism 0.847161   0.77659   0.91018 0.717682  0.60309
#> 2   rheuma hypothyroidism 0.525462   0.43744   0.60737 0.276111  0.19136
#> 3 diabetes hypothyroidism 0.859903   0.77335   0.93920 0.739434  0.59807
#>   r2.upper           p
#> 1  0.82843 4.09981e-45
#> 2  0.36890 1.29213e-24
#> 3  0.88210 1.89228e-25
```

They are! We can then run the multiple regression for all predictors.

With this function, the last phenotype will be treated as the outcome by
default (order can be changed with the ‘phenos’ argument)

``` r
run.multireg(locus, adap.thresh=NULL) # We can set the adap.thresh argument to NULL for the sake of speeding up the analysis in this example, though note that this migh lead to some loss of accuracy for very low p-values, and we do not recommended doing this for a proper analysis (see function manual for more details)
#> [1] "~ Running multiple regression for outcome 'hypothyroidism', with predictors 'asthma', 'rheuma', 'diabetes'"
#> [[1]]
#> [[1]][[1]]
#>   predictors        outcome     gamma gamma.lower gamma.upper       r2 r2.lower
#> 1     asthma hypothyroidism 0.8184750     0.70381     0.94126 0.719291   0.6136
#> 2     rheuma hypothyroidism 0.0493099    -0.10675     0.19171 0.719291   0.6136
#>   r2.upper           p
#> 1  0.83339 6.33993e-24
#> 2  0.83339 5.09336e-01
#> 
#> [[1]][[2]]
#>   predictors        outcome    gamma gamma.lower gamma.upper       r2 r2.lower
#> 1     asthma hypothyroidism 0.405003    -0.07631     0.76953 0.778212  0.68283
#> 2   diabetes hypothyroidism 0.505998     0.11995     0.98241 0.778212  0.68283
#>   r2.upper         p
#> 1  0.89727 0.1039200
#> 2  0.89727 0.0376151
#> 
#> [[1]][[3]]
#>   predictors        outcome     gamma gamma.lower gamma.upper       r2 r2.lower
#> 1     rheuma hypothyroidism -0.195568    -0.50284     0.02426 0.757829  0.60818
#> 2   diabetes hypothyroidism  1.000800     0.80113     1.27018 0.757829  0.60818
#>   r2.upper           p
#> 1  0.93001 1.12548e-01
#> 2  0.93001 2.03808e-12
#> 
#> 
#> [[2]]
#> [[2]][[1]]
#>   predictors        outcome     gamma gamma.lower gamma.upper       r2 r2.lower
#> 1     asthma hypothyroidism  0.372955    -0.31407     0.76345 0.790052  0.69434
#> 2     rheuma hypothyroidism -0.158500    -0.50786     0.04452 0.790052  0.69434
#> 3   diabetes hypothyroidism  0.648194     0.15998     1.51102 0.790052  0.69434
#>   r2.upper         p
#> 1  0.93588 0.1889980
#> 2  0.93588 0.1876560
#> 3  0.93588 0.0559055
```

By default, this function does not only return the full model with all
predictors, but also all intermediate models. The reason for this is
that it is important to be aware of how each predictor affects the fit
of the other predictors in order to be able to interpret the results
correctly.

For example, in this case we can see that both asthma and diabetes are
quite significant in the two-predictor models with rheuma, but not in
the full model, or the two-predictor model with asthma and diabetes
together. This suggests that there is some collinearity between these
predictors, which explains why neither appear very significant in the
full model (despite both showing strong bivariate correlations with the
hypothyroidism).

We can also look directly at the correlations between all predictors,
which indeed confirms our suspicion

``` r
run.bivar(locus, phenos=c("asthma","rheuma","diabetes"))
#>    phen1    phen2      rho rho.lower rho.upper       r2 r2.lower r2.upper
#> 1 asthma   rheuma 0.581756   0.51601   0.64475 0.338440  0.26627  0.41570
#> 2 asthma diabetes 0.873835   0.79919   0.94282 0.763587  0.63870  0.88892
#> 3 rheuma diabetes 0.720453   0.63190   0.80713 0.519053  0.39930  0.65146
#>             p
#> 1 1.24045e-46
#> 2 2.68922e-31
#> 3 2.85072e-29
```

Clues of this collinearity are also evident in the multivariate model
*r*<sup>2</sup>’s, which show that the proportion of genetic signal for
hypothyroidism that is explained by the predictors (jointly) is quite
high in all cases, but changes little regardless of whether either
asthma or diabetes (or both) are in the model, and is generally similar
to their individual *r*<sup>2</sup>’s with hypothyroidism (which we
determined with the bivariate previously).

Had we not looked at these intermediate models beforehand (and also
didn’t already know the bivariate correlations between the phenotypes),
this pattern would have not been obvious, and we might have made the
erroneous conclusion that neither are relevant to hypothyroidism.

For this reason, we strongly recommend examining the bivariate
correlations and computing all intermediate conditional models when
running this type of analysis. If for some reason you only want to
return the full model, however, you do so by setting the
‘only.full.model’ argument to TRUE

``` r
run.multireg(locus, only.full.model=T)
#> [1] "~ Running multiple regression for outcome 'hypothyroidism', with predictors 'asthma', 'rheuma', 'diabetes'"
#> [[1]]
#> [[1]][[1]]
#>   predictors        outcome     gamma gamma.lower gamma.upper       r2 r2.lower
#> 1     asthma hypothyroidism  0.372955    -0.31258     0.75702 0.790052   0.6914
#> 2     rheuma hypothyroidism -0.158500    -0.49925     0.04538 0.790052   0.6914
#> 3   diabetes hypothyroidism  0.648194     0.15569     1.52003 0.790052   0.6914
#>   r2.upper         p
#> 1  0.93569 0.1909840
#> 2  0.93569 0.1877720
#> 3  0.93569 0.0554024
```

### Partial correlation

In contrast with the multiple regression, the partial correlation allows
you to compute the local r<sub>*g*</sub> between two phenotypes of
interest, conditioned on one or more other phenotypes.

For this example, we will use the same locus and set of phenotypes as
before. Lets say we are now interested in the partial correlation
between hypothyroidism and diabetes. Recall that there was also some
collinearity between asthma and diabetes. We can therefore expect that
if we compute the partial correlation between hypothyroidism and
diabetes conditioned on asthma, that asthma might account for a notable
proportion of the r<sub>*g*</sub> between hypothyroidism and asthma
(which we determined to be .86, from the bivariate test earlier)

``` r
run.partial.cor(locus, phenos=c("hypothyroidism","diabetes","asthma"))
#> [1] "~ Running partial correlation for 'hypothyroidism' and 'diabetes', conditioned on 'asthma'"
#>            phen1    phen2      z r2.phen1_z r2.phen2_z     pcor ci.lower
#> 1 hypothyroidism diabetes asthma   0.717682   0.763587 0.463037  0.11456
#>   ci.upper         p
#> 1   0.7674 0.0143628
```

Indeed, here you see that the partial correlation has been almost
halved, and is no longer significant. The columns r2.phen1\_z and
r2.phen2\_z confirm that asthma explains a substantial proportion of
genetic variance for both hypothyroidism and diabetes.

If we instead condition on rheuma

``` r
run.partial.cor(locus, phenos=c("hypothyroidism","diabetes","rheuma"))
#> [1] "~ Running partial correlation for 'hypothyroidism' and 'diabetes', conditioned on 'rheuma'"
#>            phen1    phen2      z r2.phen1_z r2.phen2_z     pcor ci.lower
#> 1 hypothyroidism diabetes rheuma   0.276111   0.519053 0.815756  0.68126
#>   ci.upper           p
#> 1  0.94583 5.33825e-13
```

The partial correlation (.82) is now only slightly lower than the
bivariate correlation (.86).

As seen previously, this phenotype was far less strongly genetically
correlated with hypothyroidism (at .53), and as indicated by the r2’s
here, it is also explains less of the genetic variance in hypothyroidism
and diabetes compared to asthma.

We can also condition on both asthma and rheuma (or any number of
phenotypes) in one go:

``` r
run.partial.cor(locus, phenos=c("hypothyroidism","diabetes","rheuma","asthma"))
#> [1] "~ Running partial correlation for 'hypothyroidism' and 'diabetes', conditioned on 'rheuma' + 'asthma'"
#>            phen1    phen2             z r2.phen1_z r2.phen2_z     pcor ci.lower
#> 1 hypothyroidism diabetes rheuma;asthma   0.719291   0.831584 0.502074  0.12937
#>   ci.upper         p
#> 1   0.8692 0.0146603
```

------------------------------------------------------------------------

## Example analysis script for running bivariate *r*<sub>g</sub> analysis across all genomic loci

If you are interested in analysing the bivariate local *r*<sub>g</sub>’s
between a large amount of phenotypes, we advice using a cluster
computer. The example analysis script below shows how you might set up
an R script that can be called from the command line, analysing the
local *r*<sub>g</sub>’s across all loci in the locus file.

#### Bash

This is how you may call the R script from the command line

``` bash
Rscript "g1000_test" "test.loci" "input.info.file" "sample.overlap.file" "depression;bmi" "depression.bmi"
```

#### R script

``` r
# command line arguments, specifying input/output file names and phenotype subset
arg = commandArgs(T); ref.prefix = arg[1]; loc.file = arg[2]; info.file = arg[3]; sample.overlap.file = arg[4]; phenos = unlist(strsplit(arg[5],";")); out.fname = arg[6]

### Load package
library(lava)

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
        
        # It is possible that the locus cannot be defined for various reasons (e.g. too few SNPs), so the !is.null(locus) check is necessary before calling the analysis functions.
        if (!is.null(locus)) {
                # extract some general locus info for the output
                loc.info = data.frame(locus = locus$id, chr = locus$chr, start = locus$start, stop = locus$stop, n.snps = locus$n.snps, n.pcs = locus$K)
                
                # run the univariate and bivariate tests
                loc.out = run.univ.bivar(locus, univ.thresh=1e-4)
                u[[i]] = cbind(loc.info, loc.out$univ)
                if(!is.null(loc.out$bivar)) b[[i]] = cbind(loc.info, loc.out$bivar)
        }
}

# save the output
write.table(do.call(rbind,u), paste0(out.fname,".univ.lava"), row.names=F,quote=F,col.names=T)
write.table(do.call(rbind,b), paste0(out.fname,".bivar.lava"), row.names=F,quote=F,col.names=T)

print(paste0("Done! Analysis output written to ",out.fname,".*.lava"))
```
