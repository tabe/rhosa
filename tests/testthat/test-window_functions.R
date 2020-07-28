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

assert_monotonically_increasing <- function(f, from, to) {
    x <- runif(1, min = from, max = to)
    y <- runif(1, min = x, max = to)
    expect_lte(f(x), f(y))
}

assert_monotonically_decreasing <- function(f, from, to) {
    x <- runif(1, min = from, max = to)
    y <- runif(1, min = x, max = to)
    expect_gte(f(x), f(y))
}

test_that("Hamming window function", {
    expect_gt(.hamming_window(0), 0)
    expect_lt(.hamming_window(0), 0.1)
    assert_monotonically_increasing(.hamming_window, 0, 0.5)
    expect_equal(.hamming_window(0.5), 1)
    assert_monotonically_decreasing(.hamming_window, 0.5, 1)
    expect_gt(.hamming_window(1), 0)
    expect_lt(.hamming_window(1), 0.1)
})

test_that("Hann window function", {
    expect_equal(.hann_window(0), 0)
    assert_monotonically_increasing(.hann_window, 0, 0.5)
    expect_equal(.hann_window(0.5), 1)
    assert_monotonically_decreasing(.hann_window, 0.5, 1)
    expect_equal(.hann_window(1), 0)
})

test_that("Blackman window function", {
    expect_equal(.blackman_window(0), 0)
    assert_monotonically_increasing(.blackman_window, 0, 0.5)
    expect_equal(.blackman_window(0.5), 1)
    assert_monotonically_decreasing(.blackman_window, 0.5, 1)
    expect_equal(.blackman_window(1), 0)
})
