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

## Tests for cross_bispectrum()

assert_cross_bispectrum_result <- function(x) {
    expect_s3_class(x, "data.frame")
    expect_named(x, c("f1", "f2", "value"))
}

test_that("cross_bispectrum accepts a vector", {
    v1 <- rnorm(64)
    v2 <- rnorm(64)
    v3 <- rnorm(64)
    bs <- expect_silent(cross_bispectrum(v1, v2, v3))
    assert_cross_bispectrum_result(bs)
})

test_that("cross_bispectrum accepts an array", {
    a1 <- array(rnorm(64))
    a2 <- array(rnorm(64))
    a3 <- array(rnorm(64))
    bs <- expect_silent(cross_bispectrum(a1, a2, a3))
    assert_cross_bispectrum_result(bs)
})

test_that("cross_bispectrum accepts a matrix", {
    v1 <- rnorm(1024)
    v2 <- rnorm(1024)
    v3 <- rnorm(1024)
    m1 <- matrix(v1, ncol = 16)
    m2 <- matrix(v2, ncol = 16)
    m3 <- matrix(v3, ncol = 16)
    bs <- expect_silent(cross_bispectrum(m1, m2, m3))
    assert_cross_bispectrum_result(bs)
})

test_that("cross_bispectrum accepts a data.frame", {
    df1 <- data.frame(a = runif(64),
                     b = rnorm(64),
                     c = runif(64, min = -2, max = 2),
                     d = rnorm(64, sd = 3))
    df2 <- data.frame(a = runif(64),
                     b = rnorm(64),
                     c = runif(64, min = -2, max = 2),
                     d = rnorm(64, sd = 3))
    df3 <- data.frame(a = runif(64),
                     b = rnorm(64),
                     c = runif(64, min = -2, max = 2),
                     d = rnorm(64, sd = 3))
    bs <- expect_silent(cross_bispectrum(df1, df2, df3))
    assert_cross_bispectrum_result(bs)
})

test_that("cross_bispectrum accepts a data.matrix", {
    df1 <- data.frame(a = runif(64),
                     b = rnorm(64),
                     c = runif(64, min = -2, max = 2),
                     d = rnorm(64, sd = 3))
    df2 <- data.frame(a = runif(64),
                     b = rnorm(64),
                     c = runif(64, min = -2, max = 2),
                     d = rnorm(64, sd = 3))
    df3 <- data.frame(a = runif(64),
                     b = rnorm(64),
                     c = runif(64, min = -2, max = 2),
                     d = rnorm(64, sd = 3))
    dm1 <- data.matrix(df1)
    dm2 <- data.matrix(df3)
    dm3 <- data.matrix(df3)
    bs <- expect_silent(cross_bispectrum(dm1, dm2, dm3))
    assert_cross_bispectrum_result(bs)
})

test_that("cross_bispectrum accepts a time-series", {
    v1 <- rnorm(64)
    v2 <- rnorm(64)
    v3 <- rnorm(64)
    bs <- expect_silent(cross_bispectrum(stats::ts(v1), stats::ts(v2), stats::ts(v3)))
    assert_cross_bispectrum_result(bs)
})

test_that("cross_bispectrum rejects data of length 0", {
    expect_error(cross_bispectrum(numeric(), numeric(), numeric()), "row of length 0 given")
})

test_that("cross_bispectrum of a time-series of length 1", {
    bs <- expect_silent(cross_bispectrum(stats::ts(c(0.777)), stats::ts(c(0.888)), stats::ts(c(0.999))))
    assert_cross_bispectrum_result(bs)
})

## Tests for cross_bicoherence()

assert_cross_bicoherence_result <- function(x) {
    expect_s3_class(x, "data.frame")
    expect_named(x, c("f1", "f2", "value"))
}

test_that("cross_bicoherence accepts a vector", {
    v1 <- rnorm(64)
    v2 <- rnorm(64)
    v3 <- rnorm(64)
    bs <- expect_silent(cross_bicoherence(v1, v2, v3))
    assert_cross_bicoherence_result(bs)
})

test_that("cross_bicoherence accepts an array", {
    a1 <- array(rnorm(64))
    a2 <- array(rnorm(64))
    a3 <- array(rnorm(64))
    bs <- expect_silent(cross_bicoherence(a1, a2, a3))
    assert_cross_bicoherence_result(bs)
})

test_that("cross_bicoherence accepts a matrix", {
    v1 <- rnorm(1024)
    v2 <- rnorm(1024)
    v3 <- rnorm(1024)
    m1 <- matrix(v1, ncol = 16)
    m2 <- matrix(v2, ncol = 16)
    m3 <- matrix(v3, ncol = 16)
    bs <- expect_silent(cross_bicoherence(m1, m2, m3))
    assert_cross_bicoherence_result(bs)
})

test_that("cross_bicoherence accepts a data.frame", {
    df1 <- data.frame(a = runif(64),
                     b = rnorm(64),
                     c = runif(64, min = -2, max = 2),
                     d = rnorm(64, sd = 3))
    df2 <- data.frame(a = runif(64),
                     b = rnorm(64),
                     c = runif(64, min = -2, max = 2),
                     d = rnorm(64, sd = 3))
    df3 <- data.frame(a = runif(64),
                     b = rnorm(64),
                     c = runif(64, min = -2, max = 2),
                     d = rnorm(64, sd = 3))
    bs <- expect_silent(cross_bicoherence(df1, df2, df3))
    assert_cross_bicoherence_result(bs)
})

test_that("cross_bicoherence accepts a data.matrix", {
    df1 <- data.frame(a = runif(64),
                     b = rnorm(64),
                     c = runif(64, min = -2, max = 2),
                     d = rnorm(64, sd = 3))
    df2 <- data.frame(a = runif(64),
                     b = rnorm(64),
                     c = runif(64, min = -2, max = 2),
                     d = rnorm(64, sd = 3))
    df3 <- data.frame(a = runif(64),
                     b = rnorm(64),
                     c = runif(64, min = -2, max = 2),
                     d = rnorm(64, sd = 3))
    dm1 <- data.matrix(df1)
    dm2 <- data.matrix(df3)
    dm3 <- data.matrix(df3)
    bs <- expect_silent(cross_bicoherence(dm1, dm2, dm3))
    assert_cross_bicoherence_result(bs)
})

test_that("cross_bicoherence accepts a time-series", {
    v1 <- rnorm(64)
    v2 <- rnorm(64)
    v3 <- rnorm(64)
    bs <- expect_silent(cross_bicoherence(stats::ts(v1), stats::ts(v2), stats::ts(v3)))
    assert_cross_bicoherence_result(bs)
})

test_that("cross_bicoherence rejects data of length 0", {
    expect_error(cross_bicoherence(numeric(), numeric(), numeric()), "row of length 0 given")
})

test_that("cross_bicoherence of a time-series of length 1", {
    bs <- expect_silent(cross_bicoherence(stats::ts(c(0.777)), stats::ts(c(0.888)), stats::ts(c(0.999))))
    assert_cross_bicoherence_result(bs)
})
