## -*- mode: R -*-
##
## Copyright (C) 2021 Takeshi Abe <tabe@fixedpoint.jp>
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

.phi <- function(eta1, eta2) {
    ## log(gamma(eta2 + 1)/(-eta1)^(eta2+1))
    log(gamma(eta2 + 1)) - (eta2 + 1) * log(-eta1)
}

.revpsum <- function(x) {
    rev(cumsum(rev(x)))
}

#' Estimate cross-bicoherence's empirical null distribution by a mode matching method
#'
#' Estimate false discovery rate by fitting scaled chi-squared distribution as an empirical null of cross-bicoherence with Schwartzman's mode matching method.
#'
#' @param xbc cross-bicoherence, returned from \code{cross_bicoherence}.
#' @param t_max the upper limit of interval \deqn{S_0}, see the reference.
#' @param d the bin width of the tuning parameter.
#' 
#' @references
#' Schwartzman, Armin. “Empirical Null and False Discovery Rate Inference for Exponential Families.” Annals of Applied Statistics 2, no. 4 (December 2008): 1332–59. https://doi.org/10.1214/08-AOAS184.
#'
#' @export
mode_matching <- function(xbc, t_max = NULL, d = 0.001) {
    n <- nrow(xbc)
    bin_max <- seq(d, 1, by = d)
    K <- length(bin_max)
    bin_center <- bin_max - d/2
    hc <- unlist(Map(function(t) {
        tt <- ifelse(t == 1, t+d, t)
        sum((t-d <= xbc$value) & (xbc$value < tt))
    }, bin_max), use.names = FALSE)

    if (is.null(t_max))
        t_max <- 0.9
    repeat {
        s0 <- bin_max <= t_max
        n0 <- sum(s0)
        hc0 <- hc[s0]
        h <- log(n * d)
        data <- data.frame(lambda = hc0,
                           x1 = bin_center[s0],
                           x2 = log(bin_center[s0]),
                           h = h)

        m <- stats::glm(lambda ~ x1 + x2 + offset(h), family = "poisson", data = data)
        if (!m$converged) {
            if (t_max > 0.2) {
                warning("failed Poisson regression")
                t_max <- t_max - 0.1
                next
            } else {
                stop("failed Poisson regression")
            }
        }
        C <- m$coefficients[1]
        eta1 <- m$coefficients[2]
        if (eta1 >= 0) {
            if (t_max > 0.2) {
                warning("non-negative eta1")
                t_max <- t_max - 0.1
                next
            }
        }
        eta2 <- m$coefficients[3]
        log_p0 <- C + .phi(eta1, eta2)

        y_hat <- exp(bin_center * eta1 + log(bin_center) * eta2 + C + h)
        y_hat2 <- y_hat/2
        hc2 <- hc/2
        fdr <- (.revpsum(y_hat) - y_hat2) / (.revpsum(hc) - hc2) # div-by-zero is OK

        k <- floor(xbc$value/d)
        idx <- ifelse(k == 0, 1, k)

        return(list(t_max = t_max,
                    d = d,
                    C = C,
                    eta1 = eta1,
                    eta2 = eta2,
                    a = -1/(2 * eta1),
                    nu = 2 * (eta2 + 1),
                    log_p0 = log_p0,
                    fdr = fdr[idx]))
    }
}
