# Update history

### v0.1.5
- added support for custom LAVA LD reference file format
	- made UK Biobank based LD files available for analysis of European ancestry data
- added (partial) correction for low reference data sample sizes
- updated analysis of binary phenotypes
  - changed univariate analysis model from logistic regression to linear probability
  - changed computation of latent-scale local heritability estimates (using formulas from [Lee et al. 2011](https://pubmed.ncbi.nlm.nih.gov/21376301/))
    - if population prevalence not provided, assumes no ascertainment and uses formula 17
    - if population prevalence is provided, uses formula 23
- removed dependency on snpStats package


### v0.1.3
- added support for eQTL analysis

