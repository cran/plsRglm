% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coefs.plsRglm.R
\name{coefs.plsRglm}
\alias{coefs.plsRglm}
\title{Coefficients for bootstrap computations of PLSGLR models}
\usage{
coefs.plsRglm(
  dataset,
  ind,
  nt,
  modele,
  family = NULL,
  maxcoefvalues,
  ifbootfail,
  verbose
)
}
\arguments{
\item{dataset}{dataset to resample}

\item{ind}{indices for resampling}

\item{nt}{number of components to use}

\item{modele}{type of modele to use, see \link{plsRglm}}

\item{family}{glm family to use, see \link{plsRglm}}

\item{maxcoefvalues}{maximum values allowed for the estimates of the
coefficients to discard those coming from singular bootstrap samples}

\item{ifbootfail}{value to return if the estimation fails on a bootstrap
sample}

\item{verbose}{should info messages be displayed ?}
}
\value{
estimates on a bootstrap sample or \code{ifbootfail} value if the
bootstrap computation fails.
}
\description{
A function passed to \code{boot} to perform bootstrap.
}
\examples{

data(Cornell)

# (Y,X) bootstrap of a PLSGLR model
# statistic=coefs.plsRglm is the default for (Y,X) bootstrap of a PLSGLR models.
set.seed(250)
modplsglm <- plsRglm(Y~.,data=Cornell,1,modele="pls-glm-family",family=gaussian)
Cornell.bootYX <- bootplsglm(modplsglm, R=250, typeboot="plsmodel", 
statistic=coefs.plsRglm, verbose=FALSE)

}
\seealso{
See also \code{\link{bootplsglm}}.
}
\author{
Frédéric Bertrand\cr
\email{frederic.bertrand@utt.fr}\cr
\url{https://fbertran.github.io/homepage/}
}
\keyword{models}
