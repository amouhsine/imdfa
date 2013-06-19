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

  \item{...}{further arguments to be passed to mdfa_core}
}
\description{
  Main function to calculate filter using I-MDFA code. If a
  formula is not used, the first column of the data is the
  dependent, and the remaining columns are the predictors.
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

