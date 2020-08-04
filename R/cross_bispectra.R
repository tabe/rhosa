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

## Return a data frame of frequency pairs in the 1st quadrant.
## The unit of frequency x1 and x2 is cycles.
##
## @param n the number of samples
.generate_1st_quadrant <- function(n) {
    stopifnot(length(n) == 1, n >= 1)

    if (n < 4) {
        data.frame(x1 = integer(), x2 = integer())
    } else {
        ymax <- function(x) n %/% 2 - x
        xs <- seq_len(n %/% 2 - 1)
        do.call(rbind, Map(function(x, u) data.frame(x1 = x, x2 = 1:u), xs, ymax(xs)))
    }
}

## Return a data frame of frequency pairs in the 4th quadrant.
## The x2 is positive, yet represents a negative frequency of the same absolute value.
## Note that it excludes the row of x1 == x2, at which the third frequency is zero.
##
## @inheritParams .generate_1st_quadrant
.generate_4th_quadrant <- function(n) {
    stopifnot(length(n) == 1, n >= 1)

    if (n < 4) {
        data.frame(x1 = integer(), x2 = integer())
    } else {
        xs <- seq_len(n %/% 2 - 1)
        do.call(rbind, Map(function(x, u) expand.grid(x1 = x, x2 = setdiff(xs, x)), xs))
    }
}

## Return a data frame of frequency pairs in the 3rd region of Fig. 1 (a) in Lii and
## Helland (1981). Like .generate_4th_quadrant(), the resulting x2 is positive, yet
## represents a negative frequency of the same absolute value.
##
## @inheritParams .generate_1st_quadrant
.generate_3rd_region <- function(n) {
    data <- .generate_4th_quadrant(n)
    data[data$x1 / 2 >= data$x2,]
}

#' Estimate cross-bispectrum from time series data.
#'
#' Estimate cross-bispectrum from three real-valued time series data.
#'
#' @param x Given 1st time series, as a data frame or matrix with which columns
#' correspond to sampled stretches.
#' @param y Given 2nd time series, with the same dimension as x.
#' @param z Optional 3rd time series, with the same dimension as x (and thus as y).
#' If omitted, \code{y} is used instead.
#' @param dft_given If TRUE, suppose that DFTs is given instead of time series
#' data and skip the fast fourier transform. Default: FALSE.
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
#' K. S. Lii and K. N. Helland. 1981. Cross-Bispectrum Computation and Variance Estimation. ACM Trans. Math. Softw. 7, 3 (September 1981), 284–294. DOI:https://doi.org/10.1145/355958.355961
#'
#' @examples
#' x <- seq_len(1280)
#' v1 <- sapply(x, function(x) {sin(2 * x)}) + rnorm(1280)
#' v2 <- sapply(x, function(x) {sin(3 * x + 1)}) + rnorm(1280)
#' v3 <- sapply(x, function(x) {cos(2 * x) * cos(3 * x + 1)}) + rnorm(1280)
#' m1 <- matrix(v1, nrow = 128)
#' m2 <- matrix(v2, nrow = 128)
#' m3 <- matrix(v3, nrow = 128)
#' cross_bispectrum(m1, m2, m3)
#'
#' d1 <- stats::mvfft(m1)
#' d2 <- stats::mvfft(m2)
#' d3 <- stats::mvfft(m3)
#' cross_bispectrum(d1, d2, d3, dft_given = TRUE)
#'
#' @export
cross_bispectrum <- function(x, y, z = y,
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

    q1 <- .generate_1st_quadrant(V)

    v1 <- vapply(seq_len(nrow(q1)), function(i) {
        f1 <- q1$x1[i] + 1
        f2 <- q1$x2[i] + 1
        f3 <- q1$x1[i] + q1$x2[i] + 1
        mean(dft_x[f1,] * dft_y[f2,] * Conj(dft_z[f3,])) / ((2 * pi)^2 * V)
    }, complex(1))

    if (identical(y, z)) {
        r3 <- .generate_3rd_region(V)

        v3 <- vapply(seq_len(nrow(r3)), function(i) {
            f1 <- r3$x1[i] + 1
            f2 <- r3$x2[i] + 1
            f3 <- r3$x1[i] - r3$x2[i] + 1
            mean(dft_x[f1,] * Conj(dft_y[f2,]) * Conj(dft_z[f3,])) / ((2 * pi)^2 * V)
        }, complex(1))

        data.frame(f1 = c(q1$x1, r3$x1) / V,
                   f2 = c(q1$x2, -r3$x2) / V,
                   value = c(v1, v3))
    } else {
        q4 <- .generate_4th_quadrant(V)

        v4 <- vapply(seq_len(nrow(q4)), function(i) {
            f1 <- q4$x1[i] + 1
            f2 <- q4$x2[i] + 1
            if (q4$x1[i] > q4$x2[i]) { # in the 3rd or 4th region
                f3 <- q4$x1[i] - q4$x2[i] + 1
                mean(dft_x[f1,] * Conj(dft_y[f2,]) * Conj(dft_z[f3,])) / ((2 * pi)^2 * V)
            } else { # in the 5th or 6th region
                f3 <- q4$x2[i] - q4$x1[i] + 1
                mean(dft_x[f1,] * Conj(dft_y[f2,]) * dft_z[f3,]) / ((2 * pi)^2 * V)
            }
        }, complex(1))

        data.frame(f1 = c(q1$x1, q4$x1) / V,
                   f2 = c(q1$x2, -q4$x2) / V,
                   value = c(v1, v4))
    }
}

#' Estimate cross-coherence from time series data.
#'
#' Estimate cross-coherence from three real-valued time series data.
#'
#' @inheritParams cross_bispectrum
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
#' The estimated value of magnitude-squared cross-bicoherence at the respective
#' frequency pair.
#' }
#' }
#'
#' @inherit cross_bispectrum details references
#'
#' @examples
#' x <- seq_len(1280)
#' v1 <- sapply(x, function(x) {sin(2 * x)}) + rnorm(1280)
#' v2 <- sapply(x, function(x) {sin(3 * x + 1)}) + rnorm(1280)
#' v3 <- sapply(x, function(x) {cos(2 * x) * cos(3 * x + 1)}) + rnorm(1280)
#' m1 <- matrix(v1, nrow = 128)
#' m2 <- matrix(v2, nrow = 128)
#' m3 <- matrix(v3, nrow = 128)
#' cross_bicoherence(m1, m2, m3)
#'
#' d1 <- stats::mvfft(m1)
#' d2 <- stats::mvfft(m2)
#' d3 <- stats::mvfft(m3)
#' cross_bicoherence(d1, d2, d3, dft_given = TRUE)
#'
#' @export
cross_bicoherence <- function(x, y, z = y,
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

    q1 <- .generate_1st_quadrant(V)

    v1 <- vapply(seq_len(nrow(q1)), function(i) {
        f1 <- q1$x1[i] + 1
        f2 <- q1$x2[i] + 1
        f3 <- q1$x1[i] + q1$x2[i] + 1
        tp <- dft_x[f1,] * dft_y[f2,] * Conj(dft_z[f3,])
        abs(sum(tp)) / sum(abs(tp))
    }, numeric(1))

    if (identical(y, z)) {
        r3 <- .generate_3rd_region(V)

        v3 <- vapply(seq_len(nrow(r3)), function(i) {
            f1 <- r3$x1[i] + 1
            f2 <- r3$x2[i] + 1
            f3 <- r3$x1[i] - r3$x2[i] + 1
            tp <- dft_x[f1,] * Conj(dft_y[f2,]) * Conj(dft_z[f3,])
            abs(sum(tp)) / sum(abs(tp))
        }, numeric(1))

        data.frame(f1 = c(q1$x1, r3$x1) / V,
                   f2 = c(q1$x2, -r3$x2) / V,
                   value = c(v1, v3))
    } else {
        q4 <- .generate_4th_quadrant(V)

        v4 <- vapply(seq_len(nrow(q4)), function(i) {
            f1 <- q4$x1[i] + 1
            f2 <- q4$x2[i] + 1
            if (q4$x1[i] > q4$x2[i]) { # in the 3rd or 4th region
                f3 <- q4$x1[i] - q4$x2[i] + 1
                tp <- dft_x[f1,] * Conj(dft_y[f2,]) * Conj(dft_z[f3,])
                abs(sum(tp)) / sum(abs(tp))
            } else { # in the 5th or 6th region
                f3 <- q4$x2[i] - q4$x1[i] + 1
                tp <- dft_x[f1,] * Conj(dft_y[f2,]) * dft_z[f3,]
                abs(sum(tp)) / sum(abs(tp))
            }
        }, numeric(1))

        data.frame(f1 = c(q1$x1, q4$x1) / V,
                   f2 = c(q1$x2, -q4$x2) / V,
                   value = c(v1, v4))
    }
}
