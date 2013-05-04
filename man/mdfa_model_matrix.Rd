\name{mdfa_model_matrix}
\alias{mdfa_model_matrix}
\title{Create a matrix appropriate for passing to mdfa().}
\usage{
  mdfa_model_matrix(formula, data)
}
\arguments{
  \item{formula}{a formula describing the desired model}

  \item{a}{matrix, dataframe, list or xts object}
}
\description{
  Creates a matrix suitable for use in an imdfa model fit.
  The formula can be structured like this "y ~ y + x" or in
  the more standard "y ~ x" format, as desired.
}

