## -*- mode: R -*-
##
## Copyright (C) 2019-2020 Takeshi Abe <tabe@fixedpoint.jp>
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

## On unit of frequency
## We suppose that values of frequencies given to our API is in normalized
## frequency i.e. in cycles/sample (or, more precisely, in cycles per sampling
## interval).  So its domain is [0, 1).
## In other words, the normalized Nyquist frequency for real signals is 0.5
## cycles/sample, and the normalized sample rate is 1.

.generate_triangle <- function(n) {
    stopifnot(length(n) == 1, n >= 1)

    ymax <- function(x) ifelse(x <= n / 3, x, n - 2 * x)
    xs <- 0:(n %/% 2)
    do.call(rbind, Map(function(x, u) data.frame(x1 = x, x2 = 0:u), xs, ymax(xs)))
}

## Return the window function of given name
##
## @param name The name of a window function.
## @return the function.
.named_window_function <- function(name) {
    if (identical(name, "hamming"))
        return(.hamming_window)
    if (identical(name, "hann"))
        return(.hann_window)
    if (identical(name, "blackman"))
        return(.blackman_window)
    stop(sprintf("%s: no such function in rhosa", name))
}

## Build a taper function
##
## @param h A window function, or NULL for no tapering.
## @param V The length of the array given to the resulting function.
## @return The taper function.
.taper_function <- function(h, V) {
    if (is.null(h)) {
        identity
    } else {
        f <- .named_window_function(h)
        function(x) sapply(1:V, function(v) f(v / (V + 1)) * x[v])
    }
}

## tapered discrete Fourier transform of the l-th stretch
##
## Note that l is 1-based.
## d <- function(f, l) {
##     sum(Map(function(v) h((v + 1) / (V + 1)) * data[v + 1, l] * exp(2 * pi * -1i * v * f), 0:(V-1)))
## }
## Apply taper to each column if necessary.
.tdft <- function(data, window_function, V) {
    tapered_data <- as.matrix(apply(data, 2, .taper_function(window_function, V)))
    ## tdft is a complex-valued matrix of dimension (V, L).
    tdft <- stats::mvfft(tapered_data)
    function(x, l) {
        v <- x + 1
        stopifnot(1 <= v, v <= V)
        tdft[v, l]
    }
}

.h <- function(k, window_function) {
    if (is.null(window_function)) {
        1
    } else {
        f <- .named_window_function(window_function)
        a <- stats::integrate(function(x) f(x)^k, 0, 1)
        a$value
    }
}

.bispectrum <- function(L, V, h3, d, tr) {
    ## 3rd order periodogram of the l-th stretch
    ## Again, l is 1-based.
    ## See equation (A.5) in [1]'s Appendix A.
    I3 <- function(lambda, mu, l) {
        r <- d(lambda, l) * d(mu, l) * Conj(d(lambda + mu, l))
        r / (((2*pi)^2) * V * h3)
    }

    ## The estimate of the bispectrum
    ## See equation (A.6) in [1]'s Appendix A.
    f3 <- function(lambda, mu) {
        mean(sapply(1:L, function(l) I3(lambda, mu, l)))
    }

    data.frame(f1 = tr$x1 / V,
               f2 = tr$x2 / V,
               value = mapply(f3, tr$x1, tr$x2, SIMPLIFY = TRUE))
}

#' Estimate bispectrum from time series data.
#'
#' Estimate bispectrum from real- or complex-valued time series data.
#'
#' @param data Given time series, as a data frame or matrix with which columns
#' correspond to sampled stretches.
#' @param window_function A window function's name for tapering. Defaults to
#' \code{NULL} ("no tapering").
#'
#' Currently the following window functions are available: Hamming window ("hamming"),
#' Hann window ("hann"), and Blackman window ("blackman").
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
#' The estimated bispectrum at each frequency pair.
#' }
#' }
#'
#' @references
#' Brillinger, D.R. and Irizarry, R.A.
#' "An investigation of the second- and higher-order spectra of music."
#' Signal Processing, Volume 65, Issue 2, 30 March 1998, Pages 161-179.
#'
#' @examples
#' f <- function(x) {
#'     sin(2 * x) + sin(3 * x + 1) + sin(2 * x) * sin(3 * x + 1)
#' }
#' v <- sapply(seq_len(1280), f) + rnorm(1280)
#' m <- matrix(v, nrow = 128)
#' bispectrum(m)
#' bispectrum(m, "hamming")
#'
#' @export
bispectrum <- function(data, window_function = NULL) {

    ## Make data a matrix
    if (!is.matrix(data))
        data <- as.matrix(data)

    ## the number of stretch
    L <- ncol(data)
    ## the length of each stretch
    V <- nrow(data)
    if (V == 0)
        stop("row of length 0 given")

    h3 <- .h(3, window_function)
    d <- .tdft(data, window_function, V)
    tr <- .generate_triangle(V)
    .bispectrum(L, V, h3, d, tr)
}

#' Estimate bicoherence from given time series data.
#'
#' Estimate magnitude-squared bicoherence from given real- or complex-valued
#' time series data.
#'
#' @inheritParams bispectrum
#' @param alpha The alpha level of the hypotesis test. Defaults to 0.05.
#' @param p_adjust_method The correction method for p-values, given to
#' \code{\link[stats]{p.adjust}()}. Defaults to "BH" (Benjamini and Hochberg).
#' No correction if a non-character is given.
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
#' The estimate of magnitude-squared bicoherence at the respective frequency
#' pair.
#' }
#' \item{p_value:}{
#' The (corrected, if requested) p-value for hypothesis testing under null
#' hypothesis that bicoherence is 0.
#' }
#' \item{significance:}{
#' TRUE if the null hypothesis of the above hypothesis test is rejected
#' with given \code{alpha} level.
#' }
#' }
#'
#' @inherit bispectrum details references
#'
#' @examples
#' f <- function(x) {
#'     sin(2 * x) + sin(3 * x + 1) + sin(2 * x) * sin(3 * x + 1)
#' }
#' v <- sapply(seq_len(1280), f) + rnorm(1280)
#' m <- matrix(v, nrow = 128)
#' bicoherence(m)
#' bicoherence(m, "hamming")
#'
#' @export
bicoherence <- function(data,
                        window_function = NULL,
                        alpha = 0.05,
                        p_adjust_method = "BH") {

    ## Make data a matrix
    if (!is.matrix(data))
        data <- as.matrix(data)

    ## the number of stretch
    L <- ncol(data)
    ## the length of each stretch
    V <- nrow(data)
    if (V == 0)
        stop("row of length 0 given")

    h2 <- .h(2, window_function)
    h3 <- .h(3, window_function)
    h6 <- .h(6, window_function)

    d <- .tdft(data, window_function, V)

    ## 2nd order periodogram of the l-th strech
    ## Again, l is 1-based.
    ## The same as equation (5) in [1], but with adjusting taper's effect.
    I2 <- function(f, l) {
        abs(d(f, l))^2 / (2*pi * V * h2)
    }

    ## Estimate the power spectrum.
    ## See [1]'s definition of the estimated power spectrum.
    f2 <- function(f) {
        stopifnot(length(f) == 1)
        mean(sapply(1:L, function(l) I2(f, l)))
    }

    denom <- function(x, y) {
        stopifnot(length(x) == length(y))
        mapply(function(fx, fy) {f2(fx) * f2(fy) * f2(fx + fy)}, x, y)
    }

    tr <- .generate_triangle(V)
    bs <- .bispectrum(L, V, h3, d, tr)
    msbc <- abs(bs$value)^2 / denom(tr$x1, tr$x2)

    ## The mean of approximated distribution under null hypothesis that
    ## bicoherence = 0 is said to be exponentially distributed.
    ## See (A.11) in [1]'s Appendix A.
    m <- h6 * V / ((h3^2) * L * (2*pi))

    ## p-value, corrected if requested
    p_value <- stats::pexp(msbc, rate = 1 / m, lower.tail = FALSE)
    if (is.character(p_adjust_method))
        p_value <- stats::p.adjust(p_value, p_adjust_method)

    data.frame(f1 = bs$f1,
               f2 = bs$f2,
               value = msbc,
               p_value = p_value,
               significance = (p_value < alpha))
}

## References:
## [2] David R. Brillinger, "Asymptotic Properties of Spectral Estimates of Second Order",
##     Biometrika, Volume 56, Issue 2, August 1969, Pages 375–390, https://doi.org/10.1093/biomet/56.2.375
## [3] David R. Brillinger, "Time Series: Data Analysis and Theory (Classics in Applied Mathematics)",
##     Society for Industrial and Applied Mathematics, 2001, ISBN: 978-0-89871-501-9, https://doi.org/10.1137/1.9780898719246
