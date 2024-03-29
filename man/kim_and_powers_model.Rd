% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model.R
\name{kim_and_powers_model}
\alias{kim_and_powers_model}
\title{A test signal of the phase coherence between three oscillators}
\usage{
kim_and_powers_model(
  fbfN = 0.22,
  fcfN = 0.375,
  fdfN = fbfN + fcfN,
  num_points = 128,
  num_records = 64,
  noise_sd = 0.1,
  phase_coherence = TRUE,
  product_term = FALSE
)
}
\arguments{
\item{fbfN}{b's frequency divided by the Nyquist frequency; \code{0.220} by default.}

\item{fcfN}{c's frequency divided by the Nyquist frequency; \code{0.375} by default.}

\item{fdfN}{d's frequency divided by the Nyquist frequency; \code{fbfN + fcfN} by default.}

\item{num_points}{The number of sampling points in a record; 128 by default.}

\item{num_records}{The number of records; 64 by default.}

\item{noise_sd}{The standard deviation of a Gaussian noise perturbing samples; 0.1 (-20dB) by default.}

\item{phase_coherence}{If TRUE (default), the phase coherence in the signal d is on; otherwise off.}

\item{product_term}{If TRUE, the product of b and c is included in the model; FALSE by default.}
}
\value{
A matrix of \code{num_points} rows x \code{num_records} columns.
}
\description{
Generate test signals which involve three oscillators described in Kim and Powers (1979).
}
\details{
This function produces a list of numeric vectors; its each element represents
a test signal in which three oscillators b, c, and d are superimposed.
The ratio of the frequency of b (f1) to the Nyquist frequency is 0.220 and
the ratio of the frequency of c (f2) to the Nyquist frequency is 0.375, by default.
The d's frequency f3 is equal to f1 + f2 unless specified otherwise.
Optionally the product of b and c is also added to signals.
}
\examples{
data <- kim_and_powers_model()

}
