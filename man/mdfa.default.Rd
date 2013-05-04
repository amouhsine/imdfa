\name{mdfa.default}
\alias{mdfa.default}
\title{Fit an i-mdfa model.}
\usage{
  \method{mdfa}{default} (data, formula = NULL,
    keep_data = TRUE, spectral_estimate = NULL, d = 0,
    cutoff = pi/12, Gamma = NULL, L = 24, ...)
}
\arguments{
  \item{data}{a data.frame, list, matrix, vector or xts
  data series}

  \item{formula}{a symbolic description of the model as a
  formula object. If a formula is used, the data object
  must be a data frome or list.}

  \item{...}{further arguments to be passed to mdfa_core}
}
\description{
  Main function to calculate filter using I-MDFA code.
}
\examples{
m <- mdfa(y ~ a + b, df)
summary(m)
plot(m)
coef(m)
mp <- predict(m)
plot(mp)
}
\keyword{dfa}
\keyword{imdfa}
\keyword{mdfa}

