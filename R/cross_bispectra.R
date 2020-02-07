## -*- mode: R -*-
##
## Copyright (C) 2020 Takeshi Abe <tabe@fixedpoint.jp>
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <https://www.gnu.org/licenses/>.

.generate_1st_region <- function(n) {
    stopifnot(length(n) == 1, n >= 1)

    ymax <- function(x) ifelse(x <= n / 4, x, n / 2 - x)
    xs <- 0:(n %/% 2)
    do.call(rbind, Map(function(x, u) data.frame(x1 = x, x2 = 0:u), xs, ymax(xs)))
}

#' Estimate cross-bispectrum from time series data.
#'
#' Estimate cross-bispectrum from three real-valued time series data.
#'
#' @param x Given 1st time series, as a data frame or matrix with which columns
#' correspond to sampled stretches.
#' @param y Given 2nd time series, with the same dimension as x.
#' @param z Given 3rd time series, with the same dimension as x (and thus as y).
#' @param dft_given If TRUE, suppose that DFTs is given instead of time series
#' data and skip the fast fourier transform.
#'
#' @return A data frame including the following columns:
#' \describe{
#' \item{f1:}{
#' The first elements of frequency pairs.
#' }
#' \item{f2:}{
#' The second elements of frequency pairs.
#' }
#' \item{value:}{
#' The estimated cross-bispectrum at each frequency pair.
#' }
#' }
#'
#' @references
#' K. S. Lii and K. N. Helland. 1981. Cross-Bispectrum Computation and Variance Estimation. ACM Trans. Math. Softw. 7, 3 (September 1981), 284–294. DOI:https://doi.org/10.1145/355958.35596
#'
#' @export
cross_bispectrum <- function(x, y, z,
                             dft_given = FALSE) {

    ## Make data a matrix
    if (!is.matrix(x))
        x <- as.matrix(x)
    if (!is.matrix(y))
        y <- as.matrix(y)
    if (!is.matrix(z))
        z <- as.matrix(z)

    ## the number of stretch
    L <- ncol(x)
    ## the length of each stretch
    V <- nrow(x)
    if (V == 0)
        stop("row of length 0 given")
    if (!all(dim(x) == dim(y)))
        stop("x's dimension is different from y's")
    if (!all(dim(x) == dim(z)))
        stop("x's dimension is different from z's")

    dft_x <- if (dft_given) x else stats::mvfft(x)
    dft_y <- if (dft_given) y else stats::mvfft(y)
    dft_z <- if (dft_given) z else stats::mvfft(z)

    r1 <- .generate_1st_region(V)

    value <- sapply(seq_len(nrow(r1)), function(i) {
        f1 <- r1$x1[i] + 1
        f2 <- r1$x2[i] + 1
        f3 <- r1$x1[i] + r1$x2[i] + 1
        mean(dft_x[f1,] * dft_y[f2,] * Conj(dft_z[f3,])) / ((2 * pi)^2 * V)
    })

    data.frame(f1 = r1$x1 / V,
               f2 = r1$x2 / V,
               value = value)
}
