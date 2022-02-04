## -*- mode: R -*-
##
## Copyright (C) 2022 Takeshi Abe <tabe@fixedpoint.jp>
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

#' Calculate biperiodogram
#'
#' Calculate the biperiodogram of real-valued time series
#'
#' @param x Given time series (or its DFT), as a data frame or matrix with
#' which columns correspond to sampled stretches
#' @param dft_given If TRUE, suppose that DFTs are given instead of time series
#' data and skip the fast fourier transform. Default: \code{FALSE}.
#' @param mc If \code{TRUE}, calculation is done in parallel computation.
#' Defaults to \code{FALSE}.
#' @param mc_cores The number of cores in use for parallel computation, passed
#' \code{\link[parallel:mclapply]{parallel::mcmapply}()} etc. as \code{mc.cores}.
#'
#' @return A list with names
#' \describe{
#' \item{f1:}{
#' The first elements of frequency pairs.
#' }
#' \item{f2:}{
#' The second elements of frequency pairs.
#' }
#' \item{value:}{
#' The biperiodogram as a matrix.
#' Each of its rows is for a frequency pair; its columns correspond to stretches.
#' }
#' }
#'
#' @references
#' Hinich, M.J., 1994. Higher order cumulants and cumulant spectra. Circuits Systems and Signal Process 13, 391–402. doi:10.1007/BF01183737
#'
#' @examples
#' f <- function(x) {
#'     sin(2 * x) + sin(3 * x + 1) + sin(2 * x) * sin(3 * x + 1)
#' }
#' v <- sapply(seq_len(1280), f) + rnorm(1280)
#' m <- matrix(v, nrow = 128)
#' bp <- biperiodogram(m)
#'
#' m2 <- stats::mvfft(m)
#' bp2 <- biperiodogram(m2, dft_given = TRUE)
#'
#' @export
biperiodogram <- function(x,
                          dft_given = FALSE,
                          mc = FALSE,
                          mc_cores = getOption("mc.cores", 2L)) {

    ## Make data a matrix
    if (!is.matrix(x))
        x <- as.matrix(x)

    ## the number of stretch
    L <- ncol(x)
    ## the length of each stretch
    V <- nrow(x)
    if (V == 0)
        stop("row of length 0 given")

    dft_x <- if (dft_given) x else stats::mvfft(x)

    tr <- .generate_triangle(V)

    f <- function(x1, x2) {
        i <- x1 + 1
        j <- x2 + 1
        k <- x1 + x2 + 1
        g <- function(l) {
            dft_x[i, l] * dft_x[j, l] * Conj(dft_x[k, l]) / V
        }
        vapply(1:L, g, complex(1))
    }

    value <- if (mc)
                 parallel::mcmapply(f, tr$x1, tr$x2, SIMPLIFY = TRUE, mc.cores = mc_cores)
             else
                 mapply(f, tr$x1, tr$x2, SIMPLIFY = TRUE)

    list(f1 = tr$x1 / V,
         f2 = tr$x2 / V,
         value = value)
}
