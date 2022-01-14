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

## Tests for bispectrum()

assert_bispectrum_result <- function(x) {
    expect_s3_class(x, "data.frame")
    expect_named(x, c("f1", "f2", "value"))
}

test_that("bispectrum accepts a vector", {
    v <- rnorm(64)
    bs <- expect_silent(bispectrum(v))
    assert_bispectrum_result(bs)
})

test_that("bispectrum accepts an array", {
    a <- array(rnorm(64))
    bs <- expect_silent(bispectrum(a))
    assert_bispectrum_result(bs)
})

test_that("bispectrum accepts a matrix", {
    v <- rnorm(1024)
    m <- matrix(v, ncol = 16)
    bs <- expect_silent(bispectrum(m))
    assert_bispectrum_result(bs)
})

test_that("bispectrum accepts a data.frame", {
    df <- data.frame(a = runif(64),
                     b = rnorm(64),
                     c = runif(64, min = -2, max = 2),
                     d = rnorm(64, sd = 3))
    bs <- expect_silent(bispectrum(df))
    assert_bispectrum_result(bs)
})

test_that("bispectrum accepts a data.matrix", {
    df <- data.frame(a = runif(64),
                     b = rnorm(64),
                     c = runif(64, min = -2, max = 2),
                     d = rnorm(64, sd = 3))
    dm <- data.matrix(df)
    bs <- expect_silent(bispectrum(dm))
    assert_bispectrum_result(bs)
})

test_that("bispectrum accepts a time-series", {
    v <- rnorm(64)
    bs <- expect_silent(bispectrum(stats::ts(v)))
    assert_bispectrum_result(bs)
})

test_that("bispectrum rejects data of length 0", {
    expect_error(bispectrum(numeric()), "row of length 0 given")
})

test_that("bispectrum of a time-series of length 1", {
    bs <- expect_silent(bispectrum(stats::ts(c(0.777))))
    assert_bispectrum_result(bs)
})

test_that("bispectrum of an alternating series", {
    v <- rep(c(1, -1), 32)
    bs <- expect_silent(bispectrum(v))
    assert_bispectrum_result(bs)
    expect_equal(bs$value, rep(0+0i, 374))
})

test_that("bispectrum of series with the opposite sign", {
    v <- runif(64, min = -1, max = 1)
    bs1 <- expect_silent(bispectrum(v))
    bs2 <- expect_silent(bispectrum(-v))
    assert_bispectrum_result(bs1)
    assert_bispectrum_result(bs2)
    expect_equal(bs1$value, -bs2$value)
})

test_that("bispectrum returns the same result regardless of mc", {
    skip_on_os("windows") # where mc.cores > 1 is not allowed
    v <- runif(64, min = -1, max = 1)
    bs1 <- expect_silent(bispectrum(v))
    bs2 <- expect_silent(bispectrum(v, mc = TRUE, mc_cores = 2))
    assert_bispectrum_result(bs1)
    assert_bispectrum_result(bs2)
    expect_equal(bs1, bs2)
})

## Tests for bicoherence()

assert_bicoherence_result <- function(x) {
    expect_s3_class(x, "data.frame")
    expect_named(x, c("f1", "f2", "value", "p_value", "significance"))
}

test_that("bicoherence accepts a vector", {
    v <- rnorm(64)
    bs <- expect_silent(bicoherence(v))
    assert_bicoherence_result(bs)
})

test_that("bicoherence accepts an array", {
    a <- array(rnorm(64))
    bs <- expect_silent(bicoherence(a))
    assert_bicoherence_result(bs)
})

test_that("bicoherence accepts a matrix", {
    v <- rnorm(1024)
    m <- matrix(v, ncol = 16)
    bs <- expect_silent(bicoherence(m))
    assert_bicoherence_result(bs)
})

test_that("bicoherence accepts a data.frame", {
    df <- data.frame(a = runif(64),
                     b = rnorm(64),
                     c = runif(64, min = -2, max = 2),
                     d = rnorm(64, sd = 3))
    bs <- expect_silent(bicoherence(df))
    assert_bicoherence_result(bs)
})

test_that("bicoherence accepts a data.matrix", {
    df <- data.frame(a = runif(64),
                     b = rnorm(64),
                     c = runif(64, min = -2, max = 2),
                     d = rnorm(64, sd = 3))
    dm <- data.matrix(df)
    bs <- expect_silent(bicoherence(dm))
    assert_bicoherence_result(bs)
})

test_that("bicoherence accepts a time-series", {
    v <- rnorm(64)
    bs <- expect_silent(bicoherence(stats::ts(v)))
    assert_bicoherence_result(bs)
})

test_that("bicoherence rejects data of length 0", {
    expect_error(bicoherence(numeric()), "row of length 0 given")
})

test_that("bicoherence of a time-series of length 1", {
    bs <- expect_silent(bicoherence(stats::ts(c(0.777))))
    assert_bicoherence_result(bs)
})

test_that("bicoherence returns the same result regardless of mc", {
    skip_on_os("windows") # where mc.cores > 1 is not allowed
    v <- runif(64, min = -1, max = 1)
    bc1 <- expect_silent(bicoherence(v))
    bc2 <- expect_silent(bicoherence(v, mc = TRUE, mc_cores = 2))
    assert_bicoherence_result(bc1)
    assert_bicoherence_result(bc2)
    expect_equal(bc1, bc2)
})
