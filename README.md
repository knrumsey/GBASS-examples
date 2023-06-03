
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GBASS Examples

<!-- badges: start -->
<!-- badges: end -->

This repo contains scripts used to generate figures and tables for the
examples found in the manuscript *Generalized Bayesian MARS: Tools for
Emulating Stochastic Computer Models* (link to paper when published).

- **mcycle:** Script for quantile regression using the well-know
  motorcycle dataset.
- **Piston:** Scripts for performing simulation study comparison of
  GBASS package with several competitors on the stochatic Piston
  function.
- **Fish:** Scripts for quantile regression comparison on the LogoNet
  “Fish capture recapture” agent based model.
- **Epi:** Scripts for quantile-based sensitivity analysis of a
  stochastic SIR model with intervention.
- **qrsvm** A local copy of the qrsvm package which is (at the time of
  writing) inaccesible on CRAN.

To run these examples, the `GBASS` package ([found
here](https://github.com/knrumsey/GBASS)) is needed. The `GBASS` package
can be installed by typing:

``` r
# install.packages("devtools")
devtools::install_github("knrumsey/GBASS")
```
