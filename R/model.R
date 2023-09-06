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
#' Given three periodic functions, this function generates a list of three data
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
#' @return A list of six data frames:
#' \code{i1}, \code{i2}, \code{i3}, \code{o1}, \code{o2}, and \code{o3}.
#' Each element has \code{num_observations} columns and \code{num_samples} rows.
#' \code{i1}, \code{i2}, and \code{i3} are observations of input signals;
#' \code{o1}, \code{o2}, and \code{o3} are of output.
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

    g1 <- function(x, p) {f1(input_freq[1] * (2 * pi * x) + p)}
    g2 <- function(x, p) {f2(input_freq[2] * (2 * pi * x) + p)}
    g3 <- function(x, p) {f3(input_freq[3] * (2 * pi * x) + p)}

    tc <- function(k) {
        set.seed(k)
        ps <- stats::runif(3, min = 0, max = 2 * pi)
        function(x) {
            n <- length(x)
            i1 <- g1(x, ps[1])
            i2 <- g2(x, ps[2])
            i3 <- g3(x, ps[3])
            o1 <- i1 + stats::rnorm(n, mean = 0, sd = noise_sd[1])
            o2 <- i2 + stats::rnorm(n, mean = 0, sd = noise_sd[2])
            o3 <- i3 + stats::rnorm(n, mean = 0, sd = noise_sd[3]) + o1 * o2
            data.frame(i1, i2, i3, o1, o2, o3)
        }
    }

    sample_tc <- function() {
        Map(function(f) {f(seq_len(num_samples))}, Map(tc, seq_len(num_observations)))
    }

    get_data_frame <- function(y, name) {
        do.call(cbind, Map(function(k) {y[[k]][[name]]}, seq_len(num_observations)))
    }

    y <- sample_tc()
    list(i1 = get_data_frame(y, "i1"),
         i2 = get_data_frame(y, "i2"),
         i3 = get_data_frame(y, "i3"),
         o1 = get_data_frame(y, "o1"),
         o2 = get_data_frame(y, "o2"),
         o3 = get_data_frame(y, "o3"))
}

#' A test signal of the phase coherence between three oscillators
#'
#' Generate test signals which involve three oscillators described in Kim and Powers (1979).
#'
#' This function produces a list of numeric vectors; its each element represents
#' a test signal in which three oscillators b, c, and d are superimposed.
#' The ratio of the frequency of b (f1) to the Nyquist frequency is 0.220 and
#' the ratio of the frequency of c (f2) to the Nyquist frequency is 0.375, by default.
#' The d's frequency f3 is equal to f1 + f2 unless specified otherwise.
#' Optionally the product of b and c is also added to signals.
#'
#' @param fbfN b's frequency divided by the Nyquist frequency; \code{0.220} by default.
#' @param fcfN c's frequency divided by the Nyquist frequency; \code{0.375} by default.
#' @param fdfN d's frequency divided by the Nyquist frequency; \code{fbfN + fcfN} by default.
#' @param num_points The number of sampling points in a record; 128 by default.
#' @param num_records The number of records; 64 by default.
#' @param noise_sd The standard deviation of a Gaussian noise perturbing samples; 0.1 (-20dB) by default.
#' @param phase_coherence If TRUE (default), the phase coherence in the signal d is on; otherwise off.
#' @param product_term If TRUE, the product of b and c is included in the model; FALSE by default.
#'
#' @return A matrix of \code{num_points} rows x \code{num_records} columns.
#'
#' @examples
#' data <- kim_and_powers_model()
#'
#' @export
kim_and_powers_model <- function(fbfN = 0.220,
                                 fcfN = 0.375,
                                 fdfN = fbfN + fcfN,
                                 num_points = 128,
                                 num_records = 64,
                                 noise_sd = 0.1,
                                 phase_coherence = TRUE,
                                 product_term = FALSE) {
    do.call(cbind, Map(function(i) {
        t <- seq_len(num_points)
        tb <- stats::runif(1, min = -pi, max = pi)
        tc <- stats::runif(1, min = -pi, max = pi)
        td <- ifelse(phase_coherence, tb + tc, stats::runif(1, min = -pi, max = pi))
        b <- cos(pi * fbfN * t + tb)
        c <- cos(pi * fcfN * t + tc)
        d <- 0.5 * cos(pi * fdfN * t + td)
        n <- stats::rnorm(num_points, sd = noise_sd)
        if (product_term)
            b + c + d + b * c + n
        else
            b + c + d + n
    }, seq_len(num_records)))
}
