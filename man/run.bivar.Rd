% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis_functions.R
\name{run.bivar}
\alias{run.bivar}
\title{Bivariate local genetic correlation analysis}
\usage{
run.bivar(
  locus,
  phenos = NULL,
  target = NULL,
  adap.thresh = c(1e-04, 1e-06),
  p.values = T,
  CIs = T,
  param.lim = 1.25,
  cap.estimates = T
)
}
\arguments{
\item{locus}{Locus object created using the the \code{\link{process.locus}} function. Contains all the relevant parameters and processed sum-stats for the phenotypes of interest}

\item{phenos}{Subset of phenotypes to analyse. If NULL, all phenotypes in the locus object will be analysed}

\item{target}{Target phenotype of interest. If NULL, bivariate correlations between all pairs of phenotypes will be computed; 
Otherwise, only the relations between the target phenotype and the other phenotypes will be tested.}

\item{adap.thresh}{The thresholds at which to increase the number of iterations for the p-value generation. 
Default number of iterations is 1e+4, but will be increased to 1e+5, and 1e+6 as p-values fall below the respective thresholds.
If set to NULL, the maximum number of iterations is capped at the default (Note: this significantly speeds up the analysis, but results in poor accuracy for low p-values)}

\item{p.values}{Set to F to suppress p-values}

\item{CIs}{Set to F to suppress 95\% confidence intervals}

\item{param.lim}{The +- threshold at which estimated parameters are considered to be too far out of bounds. If the estimated parameter exceeds this threshold, it is considered unreliable and will be set to NA.}
}
\value{
Data frame with the columns:
\itemize{
    \item phen1 / phen2 - analysed phenotypes
    \item rho - standardised coefficient for the local genetic correlation
    \item rho.lower / rho.upper - 95\% confidence intervals for rho
    \item r2 - proportion of variance in genetic signal for phen1 explained by phen2 (and vice versa)
    \item r2.lower / r2.upper - 95\% confidence intervals for the r2
    \item p - simulation p-values for the local genetic correlation
}
}
\description{
Performs bivariate local genetic correlation analysis between two phenotypes.
By default, the bivariate test will be performed for all combinations of phenotypes in the locus, 
but this can be modified using the 'phenos' and 'target' arguments (see below)
}
