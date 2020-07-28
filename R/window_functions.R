## -*- mode: R -*-
##
## Copyright (C) 2019 Takeshi Abe <tabe@fixedpoint.jp>
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

## Window functions
## They takes an argument, which domain is [0, 1].

#' Hamming window function
#'
#' Calculate
#' \deqn{\frac{25}{46} - \frac{21}{46} \cos(2 \pi x).}{%
#' 25/46 - (21/46) * cos(2 * pi * x).
#' }
#'
#' @param x A real number in [0, 1].
#' @return A real number in [0, 1].
#' @seealso \code{\link{.hann_window}} and \code{\link{.blackman_window}}.
#' @keywords internal
.hamming_window <- function(x) {
    25/46 - (21/46) * cos(2 * pi * x)
}

#' Hann window function
#'
#' Calculate
#' \deqn{\frac{1 - \cos(2 * \pi * x)}{2}.}{%
#' 0.5 * (1 - cos(2 * pi * x)).
#' }
#'
#' @inheritParams .hamming_window
#' @inherit .hamming_window return
#' @seealso \code{\link{.hamming_window}} and \code{\link{.blackman_window}}.
#' @keywords internal
.hann_window <- function(x) {
    ## sin(pi * x)^2
    0.5 * (1 - cos(2 * pi * x))
}

#' Blackman window function
#'
#' Calculate
#' \deqn{\frac{42}{100} - \frac{\cos(2 \pi x)}{2} + \frac{8 \cos(4 \pi x)}{100}.}{%
#' 0.42 - 0.5 * cos(2 * pi * x) + 0.08 * cos(4 * pi * x).
#' }
#'
#' @inheritParams .hamming_window
#' @inherit .hamming_window return
#' @seealso \code{\link{.hamming_window}} and \code{\link{.hann_window}}.
#' @keywords internal
.blackman_window <- function(x) {
    0.42 - 0.5 * cos(2 * pi * x) + 0.08 * cos(4 * pi * x)
}
