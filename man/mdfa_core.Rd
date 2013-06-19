\name{mdfa_core}
\alias{mdfa_core}
\title{Core calculation for mdfa}
\usage{
  mdfa_core(L, cutoff, spectral_estimate, Gamma = NULL,
    K = NROW(spectral_estimate) - 1, lambda = 0, Lag = 0,
    expweight = 0, i1 = FALSE, i2 = FALSE,
    weight_constraint = 1, lambda_cross = 0,
    lambda_decay = c(0, 0), lambda_smooth = 0.1,
    lin_expweight = FALSE, shift_constraint = 0,
    grand_mean = TRUE, ...)
}
\description{
  Wrapper around Marc's code that sets defaults and is
  available for dispatch from within the mdfa family of
  functions.
}

