
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

# Copyright Notice

© 2021. Triad National Security, LLC. All rights reserved. This program
was produced under U.S. Government contract 89233218CNA000001 for Los
Alamos National Laboratory (LANL), which is operated by Triad National
Security, LLC for the U.S. Department of Energy/National Nuclear
Security Administration. All rights in the program are reserved by Triad
National Security, LLC, and the U.S. Department of Energy/National
Nuclear Security Administration. The Government is granted for itself
and others acting on its behalf a nonexclusive, paid-up, irrevocable
worldwide license in this material to reproduce, prepare derivative
works, distribute copies to the public, perform publicly and display
publicly, and to permit others to do so.
