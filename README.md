LAVA DOCUMENTATION
================
Josefin Werme (<j.werme@vu.nl>), Christiaan de Leeuw
(<c.a.de.leeuw@vu.nl>), CTG Lab, VU Amsterdam
2025-06-27

**NOTE**: as of LAVA version 0.1.5, custom LAVA LD reference files based on a large UK Biobank cohort is
now available for use with LAVA (linked below). We **strongly** recommend that users switch to using the UK Biobank reference instead of
1,000 Genomes for any analysis of European ancestry data, as use of the 1,000 Genomes data has been found ***\[TODO: preprint link\]*** to lead to
possible inflation of type 1 error rates and bias in the local heritability estimates due to its small sample size.

------------------------------------------------------------------------

The documentation and tutorial below shows you how to read in and
analyse data with LAVA (**L**ocal **A**nalysis of \[co\]**V**ariant
**A**ssociation): A tool developed for local genetic correlation
(*r*<sub>g</sub>) analysis. The current version of LAVA is version 0.1.5 (see [here](UPDATES.md) for
update history).

LAVA can analyse the standard bivariate local *r*<sub>g</sub> between
two phenotypes (binary as well as continuous), and account for known or
estimated sample overlap. It can also test the univariate local genetic
signal for all phenotypes of interest (i.e. the local *h*<sup>2</sup>),
which may be used to filter out non-associated loci. In addition, it can
model the local genetic relations between multiple phenotypes
simultaneously using two possible conditional models: partial
correlation and multiple regression (for more details, see the [LAVA
paper](https://www.nature.com/articles/s41588-022-01017-y)).

The tutorial will show you how to install and run LAVA using some
example input data. If you wish, you can inspect the data in the
‘vignettes/data’ folder.

To be kept informed about important LAVA updates, you can join the mailing list [here](mailto:c.a.de.leeuw@vu.nl?subject=lava-mailing-list).

------------------------------------------------------------------------

## Installing LAVA

LAVA can be installed directly from github using the following command;
add the argument ‘build_vignettes=TRUE’ to also install the vignettes
(this may require installing Pandoc first, if not yet available on your
system):

``` r
if (!require("remotes", quietly=T)) install.packages("remotes")
remotes::install_github("josefin-werme/LAVA")
```

After installation, load the LAVA library to use it:

``` r
library(LAVA)
```

## Input format

*NOTE: Example input files can be found in the ‘vignettes/data’
directory*

#### As input, LAVA needs the following data:

- **LD reference data** in custom LAVA LD format or genotype data in
  PLINK format (.bed/.bim/.fam)

  - Preprocessed reference data files can be downloaded
    [here](REFERENCE.md)

- **Input info file**, used for convenient processing of multiple
  phenotypes. Requires the columns:

  - ***phenotype***: phenotype IDs

  - ***cases***: number of cases (set to NA for continuous phenotypes)

  - ***controls:*** number of controls (set to NA for continuous
    phenotypes)

  - ***prevalence*** (optional): the population prevalence of binary
    phenotypes

    - this is only relevant if you want an estimate of the local
      population *h*<sup>2</sup> for binary phenotypes. Estimates of the
      observed local sample *h*<sup>2</sup> are still provided

  - ***filename***: paths and file names to the relevant summary
    statistics

- **Summary statistics** for all phenotypes of interest. To accommodate
  common summary statistics formats, a number of different column names
  are possible for the same column type (but **only ONE** header should
  be provided for each category; i.e do not provide a file that contains
  multiple valid SNP ID columns):

  - ***SNP / ID / SNPID_UKB/ SNPID / MarkerName / RSID / RSID_UKB***:
    SNP IDs

  - ***A1 / ALT***: effect allele

  - ***A2 / REF***: reference allele

  - ***N / NMISS / N_analyzed***: number of samples

  - ***Z / T / STAT / Zscore***: if provided, no p-values or
    coefficients are needed; otherwise, please provide both:

    - ***B / BETA / OR / logOdds***: effect size coefficients

    - ***P***: p-values

  - ***GENE***: (optional) outcome gene, for use with eQTL data

- **Locus definition file**: File that defines loci either based on
  genomic coordinates or a list of SNPs (the locus file that we used in
  the LAVA preprint can be found in the support_data folder; this file
  was obtained via <https://github.com/cadeleeuw/lava-partitioning>
  using the g1000 data phase 3, build GRCh37/hg19). The locus file
  requires the following headers:

  - ***LOC***: locus ID

  - ***CHR, START, STOP***: coordinates (i.e. basepair positions)

  - ***SNPS***: list of SNPS. ‘;’ separated, no white-space (optional)

    - note: if a SNP list is provided, coordinates will be ignored and
      the reference data will be subsetted based on SNP IDs instead.
      This can be convenient if the locus definition file is based on a
      different GRChX version than the reference data. If no SNP list is
      provided, the versions of reference and the locus file MUST match!

- **Sample overlap file** (optional)

  - This can be obtained using the results from cross-trait LDSC (see
    the ‘vignettes/sample_overlap.Rmd’ file for a walk through on how to
    do this)

## Process input

Run the script below to process the sample input files. If you wish, you
can examine the input object to get an idea of what your data looks like
after processing. You can also check out the original data in the
‘vignettes/data’ folder.

``` r
### Read in summary statistics and related info
input = process.input(input.info.file="vignettes/data/input.info.txt",           # input info file
                      sample.overlap.file="vignettes/data/sample.overlap.txt",   # sample overlap file (can be set to NULL if there is no overlap)
                      ref.prefix="vignettes/data/g1000_test",                    # reference data prefix
                      phenos=c("depression","neuro","bmi"))       # subset of phenotypes listed in the input info file that we want to process

# inspect the processed input data
ls(input)                     # this is actually an environment; hence ls() rather than str()
ls(input$sum.stats)   # processed summary statistics
head(input$sum.stats$bmi)
input$info            # processed input info file, with additional variables computed by process.input()
input$sample.overlap  # sample overlap file


head(input$ref$bim)   # bim file from reference data
```

``` r
# read more about this function
?process.input()
```

As you can see in the input\$info data frame, additional columns have
been added to indicate whether phenotypes are binary, and what the the
proportion of cases are for binary phenotypes (computed from the cases /
controls columns) *(NOTE: the N column here is not used for the
analysis, that currently needs to be provided via the summary
statistics).*

When using LAVA binary LD format reference data, the chromosome part
should be omitted from the provided prefix (eg. “lava-ukb-v1.1” instead
of “lava-ukb-v1.1_chr1”) if you are intending to analyse multiple
chromosomes.

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
```

``` r
?read.loci  # read more about the function and possible data formats
```

``` r
### Create a locus object for the first locus to prepare it for analysis
locus = process.locus(loci[1,], input)
```

If we inspect the locus object, we can see that it contains a lot of
info about the locus, such as:

- Locus coordinates (‘chr’,‘start’,‘stop’)

- The number of SNPs (‘n.snps’) / PCs (‘K’) within the locus

- The estimated PC projected joint SNP effects, $\hat{\delta}$ (‘delta’)

- The estimated sampling covariance matrix $\hat{\Sigma}$ (‘sigma’)

- The estimated genetic covariance matrix
  $\hat{\Omega} = \hat{\delta}^T\hat{\delta} / K-\hat{\Sigma}$ (‘omega’)
  with corresponding correlation matrix $\hat{\Omega}^*$ (‘omega.cor’)

``` r
ls(locus)                               # inspect locus
c(locus$chr, locus$start, locus$stop)   # locus coordinates
str(locus$snps)                         # locus snps
locus$n.snps                            # N snps
locus$omega                             # genetic covariance matrix
locus$omega.cor                         # standardised genetic covariance matrix
locus$phenos                            # locus phenotypes
```

``` r
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
#>         phen      h2.obs   h2.latent ascertained          p
#> 1 depression 8.45733e-05 0.000141198       FALSE 0.03655800
#> 2      neuro 1.16406e-04          NA       FALSE 0.03154340
#> 3        bmi 1.93535e-04          NA       FALSE 0.00146622
```

``` r

# or just a subset
run.univ(locus, phenos=c("depression","bmi"))
#>         phen      h2.obs   h2.latent ascertained          p
#> 1 depression 8.45733e-05 0.000141198       FALSE 0.03655800
#> 2        bmi 1.93535e-04          NA       FALSE 0.00146622
```

As there was little evidence of any genetic signal for either
phenotypes, there will be no point in testing the local *r*<sub>g</sub>
(if we did, the estimates would be unreliable and uninterpretable)..

If we check another locus however:

``` r
locus = process.locus(loci[3,], input)
run.univ(locus)
#>         phen      h2.obs   h2.latent ascertained           p
#> 1 depression 0.000263286 0.000439566       FALSE 3.19957e-07
#> 2      neuro 0.000573255          NA       FALSE 1.25993e-14
#> 3        bmi 0.000975242          NA       FALSE 6.60429e-32
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
#> 1 depression neuro 0.5258740   0.26976   0.77509 0.27654300  0.07277  0.60077
#> 2 depression   bmi 0.0692754  -0.18337   0.31952 0.00479909  0.00000  0.10342
#> 3      neuro   bmi 0.1069490  -0.09500   0.30912 0.01143810  0.00000  0.09556
#>             p
#> 1 0.000474846
#> 2 0.581570000
#> 3 0.300358000
```

When multiple phenotypes are entered simultaneously, bivariate
*r*<sub>g</sub>‘s will be computed between all unique phenotype pairs by
default. If you want to analyse only a subset of phenotypes, you can use
the ’phenos’ argument

``` r
run.bivar(locus, phenos=c("neuro","depression"))
#>   phen1      phen2      rho rho.lower rho.upper       r2 r2.lower r2.upper
#> 1 neuro depression 0.525874   0.26693   0.77677 0.276543  0.07125  0.60338
#>             p
#> 1 0.000507518
```

Additionally, should you only be interested in the *r*<sub>g</sub>‘s
between one target phenotype and all others, you can use the ’target’
argument to prevent all pairwise tests from being computed and focus
only on those relevant to the target

``` r
run.bivar(locus, target="bmi")
#>        phen1 phen2       rho rho.lower rho.upper         r2 r2.lower r2.upper
#> 1 depression   bmi 0.0692754  -0.18681   0.32451 0.00479909        0  0.10761
#> 2      neuro   bmi 0.1069490  -0.09720   0.31104 0.01143810        0  0.09677
#>          p
#> 1 0.581177
#> 2 0.300670
```

``` r
?run.bivar  # check the function description for more options
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
#>         phen      h2.obs   h2.latent ascertained           p
#> 1 depression 0.000263286 0.000439566       FALSE 3.19957e-07
#> 2      neuro 0.000573255          NA       FALSE 1.25993e-14
#> 3        bmi 0.000975242          NA       FALSE 6.60429e-32
#> 
#> $bivar
#>        phen1 phen2       rho rho.lower rho.upper         r2 r2.lower r2.upper
#> 1 depression neuro 0.5258740   0.26787   0.77017 0.27654300  0.07176  0.59316
#> 2 depression   bmi 0.0692754  -0.18531   0.32807 0.00479909  0.00000  0.10951
#> 3      neuro   bmi 0.1069490  -0.09730   0.31002 0.01143810  0.00000  0.09611
#>             p
#> 1 0.000511901
#> 2 0.581637000
#> 3 0.301516000
```

``` r

# or with a custom p-value threshold
run.univ.bivar(locus, univ.thresh = 1e-8)
#> $univ
#>         phen      h2.obs   h2.latent ascertained           p
#> 1 depression 0.000263286 0.000439566       FALSE 3.19957e-07
#> 2      neuro 0.000573255          NA       FALSE 1.25993e-14
#> 3        bmi 0.000975242          NA       FALSE 6.60429e-32
#> 
#> $bivar
#>   phen1 phen2      rho rho.lower rho.upper        r2 r2.lower r2.upper        p
#> 1 neuro   bmi 0.106949  -0.09611    0.3137 0.0114381        0  0.09846 0.300953
```

*\[ You can use the ‘phenos’ and ‘target’ arguments with this function
too \]*

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
#> [1] "...Processing input info"
#> [1] "...** All phenotypes treated as BINARY ('asthma', 'rheuma', 'diabetes', 'hypothyroidism')"
#> [1] "...Reading in sumstats"
#> [1] "...Reading in SNP info from reference data"
#> [1] "...Extracting common SNPs"
#> [1] "... 90342 SNPs shared across data sets"
#> [1] "...Aligning effect alleles to reference data set"
#> [1] "...Removing 12980 SNPs which could not be aligned, 77362 remaining"
```

``` r

locus = process.locus(loci[19,], input)
```

Say that we are interested in modelling hypothyroidism as the outcome.
Lets first check that there is actually some local genetic signal for
all phenotypes, and that all predictors are genetically correlated with
the outcome within this locus

``` r
run.univ.bivar(locus, target="hypothyroidism") 
#> $univ
#>             phen      h2.obs  h2.latent ascertained            p
#> 1         asthma 0.001832050 0.00495321       FALSE 1.77873e-143
#> 2         rheuma 0.060403700 0.11280300       FALSE  0.00000e+00
#> 3       diabetes 0.000664802 0.00305030       FALSE  1.60790e-48
#> 4 hypothyroidism 0.001736150 0.00744354       FALSE  8.48011e-84
#> 
#> $bivar
#>      phen1          phen2      rho rho.lower rho.upper       r2 r2.lower
#> 1   asthma hypothyroidism 0.839093   0.76826   0.90496 0.704077  0.59022
#> 2   rheuma hypothyroidism 0.454805   0.36246   0.54278 0.206848  0.13138
#> 3 diabetes hypothyroidism 0.862047   0.77788   0.94112 0.743125  0.60510
#>   r2.upper           p
#> 1  0.81895 1.33436e-42
#> 2  0.29461 7.50381e-19
#> 3  0.88571 3.08935e-25
```

They are! We can then run the multiple regression with all predictors.

``` r
run.multireg(locus, target='hypothyroidism', adap.thresh=NULL) # We can set the adap.thresh argument to NULL for the sake of speeding up the analysis in this example, though note that this migh lead to some loss of accuracy for very low p-values, and we do not recommended doing this for a proper analysis (see function manual for more details)
#> [1] "~ Running multiple regression for outcome 'hypothyroidism', with predictors 'asthma', 'rheuma', 'diabetes'"
#> [[1]]
#> [[1]][[1]]
#>   predictors        outcome     gamma gamma.lower gamma.upper       r2 r2.lower
#> 1     asthma hypothyroidism 0.7950650     0.69949     0.89504 0.711771   0.6039
#> 2     rheuma hypothyroidism 0.0981469    -0.03592     0.22070 0.711771   0.6039
#>   r2.upper           p
#> 1  0.82595 3.47167e-32
#> 2  0.82595 1.41994e-01
#> 
#> [[1]][[2]]
#>   predictors        outcome    gamma gamma.lower gamma.upper       r2 r2.lower
#> 1     asthma hypothyroidism 0.387956    -0.01046     0.69098 0.786513  0.68822
#> 2   diabetes hypothyroidism 0.534753     0.21720     0.92602 0.786513  0.68822
#>   r2.upper         p
#> 1  0.90409 0.0700524
#> 2  0.90409 0.0126032
#> 
#> [[1]][[3]]
#>   predictors        outcome     gamma gamma.lower gamma.upper       r2 r2.lower
#> 1     rheuma hypothyroidism -0.150299    -0.39228     0.03910 0.756684  0.60939
#> 2   diabetes hypothyroidism  0.957072     0.79879     1.15923 0.756684  0.60939
#>   r2.upper           p
#> 1  0.91843 1.46802e-01
#> 2  0.91843 1.28128e-15
#> 
#> 
#> [[2]]
#> [[2]][[1]]
#>   predictors        outcome      gamma gamma.lower gamma.upper       r2
#> 1     asthma hypothyroidism  0.3586490    -0.15148     0.68187 0.792225
#> 2     rheuma hypothyroidism -0.0996357    -0.37206     0.07799 0.792225
#> 3   diabetes hypothyroidism  0.6224720     0.21184     1.24739 0.792225
#>   r2.lower r2.upper        p
#> 1  0.69201  0.93043 0.140059
#> 2  0.69201  0.93043 0.329050
#> 3  0.69201  0.93043 0.028268
```

Here, the ‘target’ argument specifies the outcome phenotype of interest,
and all others will be treated as predictors by default (if you want
only a subset of predictors, use the ‘phenos’ argument).

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

We can also look directly at the bivariate correlations between all
predictors, which indeed confirms this

``` r
run.bivar(locus, phenos=c("asthma","rheuma","diabetes"))
#>    phen1    phen2      rho rho.lower rho.upper       r2 r2.lower r2.upper
#> 1 asthma   rheuma 0.448590   0.37588   0.52025 0.201233  0.14128  0.27066
#> 2 asthma diabetes 0.843636   0.76135   0.92109 0.711721  0.57966  0.84841
#> 3 rheuma diabetes 0.632245   0.53112   0.73093 0.399734  0.28209  0.53427
#>             p
#> 1 1.16328e-27
#> 2 4.62370e-29
#> 3 2.16009e-22
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
run.multireg(locus, target='hypothyroidism', only.full.model=T)
#> [1] "~ Running multiple regression for outcome 'hypothyroidism', with predictors 'asthma', 'rheuma', 'diabetes'"
#> [[1]]
#> [[1]][[1]]
#>   predictors        outcome      gamma gamma.lower gamma.upper       r2
#> 1     asthma hypothyroidism  0.3586490    -0.15599     0.68721 0.792225
#> 2     rheuma hypothyroidism -0.0996357    -0.36418     0.07686 0.792225
#> 3   diabetes hypothyroidism  0.6224720     0.21572     1.24415 0.792225
#>   r2.lower r2.upper         p
#> 1  0.69311  0.92856 0.1528650
#> 2  0.69311  0.92856 0.3310510
#> 3  0.69311  0.92856 0.0292865
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

We use the ‘target’ argument to denote which phenotype pair we are
interested in, and use the ‘phenos’ argument to exclude rheuma

``` r
run.pcor(locus, target=c("hypothyroidism","diabetes"), phenos='asthma')
#> [1] "~ Running partial correlation for 'hypothyroidism' and 'diabetes', conditioned on 'asthma'"
#>            phen1    phen2      z r2.phen1_z r2.phen2_z     pcor ci.lower
#> 1 hypothyroidism diabetes asthma   0.704077   0.711721 0.527801  0.22299
#>   ci.upper          p
#> 1  0.80226 0.00303935
```

Indeed, here you see that the partial correlation has been almost
halved, and is no longer significant. The columns r2.phen1_z and
r2.phen2_z confirm that asthma explains a substantial proportion of
genetic variance for both hypothyroidism and diabetes.

If we instead condition on rheuma

``` r
run.pcor(locus, target=c("hypothyroidism","diabetes"), phenos="rheuma", adap.thresh=NULL)
#> [1] "~ Running partial correlation for 'hypothyroidism' and 'diabetes', conditioned on 'rheuma'"
#>            phen1    phen2      z r2.phen1_z r2.phen2_z     pcor ci.lower
#> 1 hypothyroidism diabetes rheuma   0.206848   0.399734 0.832604  0.71944
#>   ci.upper           p
#> 1  0.94413 9.55212e-17
```

The partial correlation (.82) is now only slightly lower than the
bivariate correlation (.86), and still highly significant.

As seen previously, this phenotype was far less strongly genetically
correlated with hypothyroidism (at .53), and as indicated by the r2’s
here, it is also explains less of the genetic variance in hypothyroidism
and diabetes compared to asthma.

We can also condition on both asthma and rheuma (or any number of
phenotypes) in one go. Note that unless subsetting is done with the
‘phenos’ argument, all phenotypes will be conditioned on by default:

``` r
run.pcor(locus, target=c("hypothyroidism","diabetes"))
#> [1] "~ Running partial correlation for 'hypothyroidism' and 'diabetes', conditioned on 'asthma' + 'rheuma'"
#>            phen1    phen2             z r2.phen1_z r2.phen2_z     pcor ci.lower
#> 1 hypothyroidism diabetes asthma;rheuma   0.711771   0.792362 0.528328  0.18768
#>   ci.upper          p
#> 1  0.84762 0.00639349
```

------------------------------------------------------------------------

## Analysis of eQTL input

LAVA can also be used in conjunction with eQTL input data to estimate
and test local genetic correlations of phenotypes with the gene
expression of specific genes, akin to TWAS, and additional input
processing functions have been provided to facilite this.

To perform the analysis, first run the process.input() function as
described above to load the summary statistics for all non-eQTL input
phenotypes. The resulting `input` object returned by process.input() is
then run through a secondary eQTL input processing function to add the
eQTL summary statistics to it.

### Loading input

The eQTL summary statistics should be provided in general data files,
which should have the same format as the input files for regular
phenotypes, plus an additional **GENE** column indicating which gene the
SNP associations on each line belong to. Input files must be provided
separately for each chromosome, with file names that are identical
except for the chromosome code. Data is loaded as follows:

``` r
process.eqtl.input(input, "data/eqtl_text/whole_blood.chr[CHR].stats", chromosomes=14, sample.size=750)
```

The second argument of process.eqtl.input() specifies the eQTL file
format, with the `[CHR]` a placeholder for the chromosome code; in this
case, LAVA will look for the file
`data/eqtl_text/whole_blood.chr14.stats`. The sample.size argument is
optional, and can be omitted if a sample size column is present in the
input files.

### Performing the analysis

After loading the eQTL data into the `input` object, analysis can be
performed by generating a locus object for the gene to be analysed:

``` r
locus = process.eqtl.locus("ENSG00000005700.14", input)
```

This function works the same way as process.locus(), except that the
first argument should specify a gene rather than a locus. Genes can be
specified either by name (as above) or by numeric index, and the list of
all loaded genes is stored in `input$current.genes`. Subsequent analysis
on the `locus` object is performed in the same way as in a regular LAVA
analysis, and so a simple analysis script would look something like:

``` r
for (i in 1:length(input$current.genes)) {
  locus = process.eqtl.locus(i, input)
  univ = run.univ(locus)
  bivar = run.bivar(locus)
  
  # ...
}
```

------------------------------------------------------------------------

## Example analysis script for running bivariate *r*<sub>g</sub> analysis across all genomic loci

If you are interested in analysing the bivariate local *r*<sub>g</sub>’s
between a large amount of phenotypes, we advice using a cluster
computer. The example analysis script below shows how you might set up
an R script that can be called from the command line, analysing the
local *r*<sub>g</sub>’s across all loci in the locus file.

#### Bash

This is how you may call the R script from the command line (**note**:
file paths may need to be adapted depending on your set-up)

``` bash
Rscript lava_script.R "g1000_test" "test.loci" "input.info.txt" "sample.overlap.txt" "depression;bmi" "depression.bmi"
```

#### R script

``` r
# command line arguments, specifying input/output file names and phenotype subset
arg = commandArgs(T); ref.prefix = arg[1]; loc.file = arg[2]; info.file = arg[3]; sample.overlap.file = arg[4]; phenos = unlist(strsplit(arg[5],";")); out.fname = arg[6]

### Load package
library(LAVA)

### Read in data
loci = read.loci(loc.file); n.loc = nrow(loci)
input = process.input(info.file, sample.overlap.file, ref.prefix, phenos)

### Set univariate pvalue threshold
univ.p.thresh = #[SET]

### Analyse
print(paste("Starting LAVA analysis for",n.loc,"loci"))
progress = ceiling(quantile(1:n.loc, seq(.05,1,.05)))   # (if you want to print the progress)

u=b=list()
for (i in 1:n.loc) {
        if (i %in% progress) print(paste("..",names(progress[which(progress==i)])))     # (printing progress)
        locus = process.locus(loci[i,], input)                                          # process locus
        
        # It is possible that the locus cannot be defined for various reasons (e.g. too few SNPs), so the !is.null(locus) check is necessary before calling the analysis functions.
        if (!is.null(locus)) {
                # extract some general locus info for the output
                loc.info = data.frame(locus = locus$id, chr = locus$chr, start = locus$start, stop = locus$stop, n.snps = locus$n.snps, n.pcs = locus$K)
                
                # run the univariate and bivariate tests
                loc.out = run.univ.bivar(locus, univ.thresh = univ.p.thresh)
                u[[i]] = cbind(loc.info, loc.out$univ)
                if(!is.null(loc.out$bivar)) b[[i]] = cbind(loc.info, loc.out$bivar)
        }
}

# save the output
write.table(do.call(rbind,u), paste0(out.fname,".univ.lava"), row.names=F,quote=F,col.names=T)
write.table(do.call(rbind,b), paste0(out.fname,".bivar.lava"), row.names=F,quote=F,col.names=T)

print(paste0("Done! Analysis output written to ",out.fname,".*.lava"))
```

------------------------------------------------------------------------

For more elaborate scripts intended for whole genome bivariate rg
analyses on a SLURM cluster (which also **parallelise** across loci),
please see the cluster_setup/scripts folder here:
<https://surfdrive.surf.nl/files/index.php/s/rtmNm8YZERGwl7f>
