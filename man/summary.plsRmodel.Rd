% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.plsRmodel.R
\name{summary.plsRmodel}
\alias{summary.plsRmodel}
\title{Summary method for plsR models}
\usage{
\method{summary}{plsRmodel}(object, ...)
}
\arguments{
\item{object}{an object of the class \code{"plsRmodel"}}

\item{\dots}{further arguments to be passed to or from methods.}
}
\value{
\item{call }{function call of plsRmodel}
}
\description{
This function provides a summary method for the class \code{"plsRmodel"}
}
\examples{

data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
modpls <- plsR(yCornell,XCornell,3,modele="pls")
class(modpls)
summary(modpls)
rm(list=c("XCornell","yCornell","modpls"))

}
\references{
Nicolas Meyer, Myriam Maumy-Bertrand et
Frédéric Bertrand (2010). Comparing the linear and the
logistic PLS regression with qualitative predictors: application to
allelotyping data. \emph{Journal de la Societe Francaise de Statistique},
151(2), pages 1-18.
\url{http://publications-sfds.math.cnrs.fr/index.php/J-SFdS/article/view/47}
}
\seealso{
\code{\link{summary}}
}
\author{
Frédéric Bertrand\cr
\email{frederic.bertrand@utt.fr}\cr
\url{https://fbertran.github.io/homepage/}
}
\keyword{methods}
\keyword{print}
