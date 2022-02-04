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

## Tests for biperiodogram()

assert_biperiodogram_result <- function(x) {
    expect_type(x, "list")
    expect_named(x, c("f1", "f2", "value"))
}

expect_equal_biperiodogram <- function(x, y) {
    expect_equal(x$f1, y$f1)
    expect_equal(x$f2, y$f2)
    expect_equal(x$value, y$value)
}

test_that("biperiodogram accepts a vector", {
    v <- rnorm(64)
    bp <- expect_silent(biperiodogram(v))
    assert_biperiodogram_result(bp)
})

test_that("biperiodogram accepts an array", {
    a <- array(rnorm(64))
    bp <- expect_silent(biperiodogram(a))
    assert_biperiodogram_result(bp)
})

test_that("biperiodogram accepts a matrix", {
    v <- rnorm(1024)
    m <- matrix(v, ncol = 16)
    bp <- expect_silent(biperiodogram(m))
    assert_biperiodogram_result(bp)
})

test_that("biperiodogram accepts a data.frame", {
    df <- data.frame(a = runif(64),
                     b = rnorm(64),
                     c = runif(64, min = -2, max = 2),
                     d = rnorm(64, sd = 3))
    bp <- expect_silent(biperiodogram(df))
    assert_biperiodogram_result(bp)
})

test_that("biperiodogram accepts a data.matrix", {
    df <- data.frame(a = runif(64),
                     b = rnorm(64),
                     c = runif(64, min = -2, max = 2),
                     d = rnorm(64, sd = 3))
    dm <- data.matrix(df)
    bp <- expect_silent(biperiodogram(dm))
    assert_biperiodogram_result(bp)
})

test_that("biperiodogram accepts a time-series", {
    v <- rnorm(64)
    bp <- expect_silent(biperiodogram(stats::ts(v)))
    assert_biperiodogram_result(bp)
})

test_that("biperiodogram rejects data of length 0", {
    expect_error(biperiodogram(numeric()), "row of length 0 given")
})

test_that("biperiodogram of a time-series of length 1", {
    bp <- expect_silent(biperiodogram(stats::ts(rnorm(1))))
    assert_biperiodogram_result(bp)
})

test_that("biperiodogram of an alternating series", {
    v <- rep(c(1, -1), 32)
    bp <- expect_silent(biperiodogram(v))
    assert_biperiodogram_result(bp)
    expect_equal(bp$value, rep(0+0i, 374))
})

test_that("biperiodogram of series with the opposite sign", {
    v <- runif(64, min = -1, max = 1)
    bp1 <- expect_silent(biperiodogram(v))
    bp2 <- expect_silent(biperiodogram(-v))
    assert_biperiodogram_result(bp1)
    assert_biperiodogram_result(bp2)
    expect_equal(bp1$value, -bp2$value)
})

test_that("option dft_given does not change the result", {
    v <- rnorm(1024)
    m <- matrix(v, ncol = 16)
    bp1 <- expect_silent(biperiodogram(m))
    m2 <- stats::mvfft(m)
    bp2 <- expect_silent(biperiodogram(m2, dft_given = TRUE))
    assert_biperiodogram_result(bp1)
    assert_biperiodogram_result(bp2)
    expect_equal_biperiodogram(bp1, bp2)
})

test_that("biperiodogram returns the same result regardless of mc", {
    skip_on_os("windows") # where mc.cores > 1 is not allowed
    x <- runif(64, min = -1, max = 1)
    bp1 <- expect_silent(biperiodogram(x))
    bp2 <- expect_silent(biperiodogram(x, mc = TRUE, mc_cores = 2))
    assert_biperiodogram_result(bp1)
    assert_biperiodogram_result(bp2)
    expect_equal_biperiodogram(bp1, bp2)
})
