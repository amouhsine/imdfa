\name{mdfa}
\alias{mdfa}
\title{Fit an i-mdfa model.}
\usage{
  mdfa(x, ...)
}
\arguments{
  \item{data}{a data.frame, list, matrix, vector or xts
  data series}

  \item{formula}{a symbolic description of the model as a
  formula object. If a formula is used, the data object
  must be a data frome or list.}

  \item{spectral_estimate}{The spectral estimate. Usually
  the output of a DFT for each of the variables.  If a
  spectral_estimate is passed into any of these methods,
  data and formula parameters are ignored.}

  \item{...}{further arguments to be passed to mdfa_core}
}
\description{
  Main function to calculate filter using I-MDFA code. If a
  formula is not used, the first column of the data is the
  dependent, and the remaining columns are the predictors.
}
\examples{
# Reproduce Wildi's 2012-06 EURI example
m <- mdfa(gdp_monthly ~ gdp_monthly + production_sa_1 + order_book_sa_2 + stocks_sa_4 + prod_expect_sa_5,
        data = euri)
summary(m)
plot(m)
coef(m)
mp <- predict(m)
plot(mp)
}
\keyword{dfa}
\keyword{imdfa}
\keyword{mdfa}

