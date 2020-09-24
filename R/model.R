## -*- mode: R; fill-column: 80 -*-
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

#' A three-channel model of quadratic phase coupling
#'
#' Simulate observations by a three-channel model of quadratic phase coupling.
#'
#' Given three periodic functions, this function generate a list of three data
#' frames in which each column represents a simulated observation at a channel.
#' The phase is chosen at random from \eqn{[0, 2 \pi]}{[0, 2 * pi]} for each
#' observation and each channel.
#'
#' @param f1 A function of period \eqn{2 \pi}{2 * pi} for the first channel.
#' @param f2 A function of period \eqn{2 \pi}{2 * pi} for the second channel.
#' @param f3 A function of period \eqn{2 \pi}{2 * pi} for the third channel.
#' @param num_samples The number of sampling points in an observation.
#' @param num_observations The number of observations.
#' @param input_freq The scaling factor for the frequencies of input periodic functions.
#' It can be a scalar or a vector of length three.
#' If a scalar is given, the same frequency is used for all of inputs.
#' @param noise_sd The standard deviation of a Gaussian noise perturbing samples.
#' It can be a scalar or a vector of length three.
#' If a scalar is given, the same value is used for all of noises.
#' Giving 0 is possible and specifies no noise.
#'
#' @return A list of three data frames: \code{c1}, \code{c2}, and \code{c3}.
#' Each element has \code{num_observations} columns and \code{num_samples} rows.
#'
#' @examples
#' sawtooth <- function(r) {
#'     x <- r/(2*pi)
#'     x - floor(x) - 0.5
#' }
#' data <- three_channel_model(cos, sin, sawtooth,
#'                             input_freq = c(0.2, 0.3, 0.4),
#'                             noise_sd = 0.9)
#'
#' @export
three_channel_model <- function(f1, f2, f3,
                                num_samples = 256,
                                num_observations = 100,
                                input_freq = c(1.2, 0.7, 0.8),
                                noise_sd = 1) {
    if (length(input_freq) == 1)
        input_freq <- rep_len(input_freq, 3)
    if (length(noise_sd) == 1)
        noise_sd <- rep_len(noise_sd, 3)

    i1 <- function(x, p) {f1(input_freq[1] * (2 * pi * x) + p)}
    i2 <- function(x, p) {f2(input_freq[2] * (2 * pi * x) + p)}
    i3 <- function(x, p) {f3(input_freq[3] * (2 * pi * x) + p)}

    tc <- function(k) {
        set.seed(k)
        ps <- stats::runif(3, min = 0, max = 2 * pi)
        function(x) {
            n <- length(x)
            c1 <- i1(x, ps[1]) + stats::rnorm(n, mean = 0, sd = noise_sd[1])
            c2 <- i2(x, ps[2]) + stats::rnorm(n, mean = 0, sd = noise_sd[2])
            c3 <- i3(x, ps[3]) + stats::rnorm(n, mean = 0, sd = noise_sd[3]) + c1 * c2
            data.frame(c1, c2, c3)
        }
    }

    sample_tc <- function() {
        Map(function(f) {f(seq_len(num_samples))}, Map(tc, seq_len(num_observations)))
    }

    c1_data_frame <- function(y) {
        do.call(cbind, Map(function(k) {y[[k]]$c1}, seq_len(num_observations)))
    }

    c2_data_frame <- function(y) {
        do.call(cbind, Map(function(k) {y[[k]]$c2}, seq_len(num_observations)))
    }

    c3_data_frame <- function(y) {
        do.call(cbind, Map(function(k) {y[[k]]$c3}, seq_len(num_observations)))
    }

    y <- sample_tc()
    list(c1 = c1_data_frame(y),
         c2 = c2_data_frame(y),
         c3 = c3_data_frame(y))
}
