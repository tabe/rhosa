## -*- mode: R -*-
##
## Copyright (C) 2023 Takeshi Abe <tabe@fixedpoint.jp>
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

test_that("kim_and_powers_model", {
    data1 <- expect_silent(kim_and_powers_model())
    expect_type(data1, "double")
    expect_equal(nrow(data1), 128)
    expect_equal(ncol(data1), 64)

    data2 <- expect_silent(kim_and_powers_model(phase_coherence = FALSE, product_term = FALSE))
    expect_type(data2, "double")
    expect_equal(nrow(data2), 128)
    expect_equal(ncol(data2), 64)

    data3 <- expect_silent(kim_and_powers_model(phase_coherence = FALSE, product_term = TRUE))
    expect_type(data3, "double")
    expect_equal(nrow(data3), 128)
    expect_equal(ncol(data3), 64)

    data4 <- expect_silent(kim_and_powers_model(phase_coherence = TRUE, product_term = TRUE))
    expect_type(data4, "double")
    expect_equal(nrow(data4), 128)
    expect_equal(ncol(data4), 64)
})
