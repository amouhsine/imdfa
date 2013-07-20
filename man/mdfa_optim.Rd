\name{mdfa_optim}
\alias{mdfa_optim}
\title{mdfa set up for optimization}
\usage{
  mdfa_optim(x, L, spectral_estimate, Gamma, criterion,
    crit_args = NULL, ...)
}
\arguments{
  \item{x}{a vector of parameters to optimize on}

  \item{spectral_estimate}{a spectral estimate object}

  \item{Gamma}{}

  \item{criterion}{a string selecting one of the fit
  statistics or a function returning one number}

  \item{crit_args}{a list of further named arguments to
  pass to criteron}

  \item{...}{further arguments to pass to mdfa_core}
}
\description{
  WARNING: This function is very alpha, parameters may
  change. A wrapper for mdfa that takes a vector of
  parameters so that it can be used in optim() calls. The
  order of parameters in x: 1. cutoff 2. lambda 3.
  expweight 4. lambda_cross 5. lambda_decay[[1]] 6.
  lambda_decay[[2]] 7. lambda_smooth
}
\details{
  If the criteron is a function, it should take a mdfa
  object and return a numeric vector of length one.
}

