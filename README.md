
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rhosa: Higher-Order Spectral Analysis in R

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/rhosa)](https://CRAN.R-project.org/package=rhosa)
[![R-CMD-check](https://github.com/tabe/rhosa/workflows/R-CMD-check/badge.svg)](https://github.com/tabe/rhosa/actions)
<!-- badges: end -->

This package aims to provide functions to estimate higher-order spectra
or polyspectra of multivariate time series, such as bispectrum and
bicoherence ([Brillinger and Irizarry
1998](#ref-brillinger_investigation_1998)). They are useful for
e.g. detecting nonlinear interaction between stationary time series
driven by periodic signals ([Abe et al. 2024](#ref-abe_detecting_2024)).

## Installation

You can install the released version of rhosa from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("rhosa")
```

Alternatively, the development version from
[GitHub](https://github.com/) with
[remotes](https://cran.r-project.org/package=remotes):

``` r
# install.packages("remotes")
remotes::install_github("tabe/rhosa")
```

## Acknowledgement

The author thanks Alessandro E. P. Villa for his generous support to
this project.

## References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-abe_detecting_2024" class="csl-entry">

Abe, Takeshi, Yoshiyuki Asai, Alessandra Lintas, and Alessandro E. P.
Villa. 2024. “Detection of Quadratic Phase Coupling by Cross-Bicoherence
and Spectral Granger Causality in Bifrequencies Interactions.”
*Scientific Reports* 14 (1): 8521.
<https://doi.org/10.1038/s41598-024-59004-8>.

</div>

<div id="ref-brillinger_investigation_1998" class="csl-entry">

Brillinger, D. R., and R. A. Irizarry. 1998. “An Investigation of the
Second- and Higher-Order Spectra of Music.” *Signal Processing* 65 (2):
161–79. <https://doi.org/10.1016/S0165-1684(97)00217-X>.

</div>

</div>
