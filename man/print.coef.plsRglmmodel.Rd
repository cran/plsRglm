% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.coef.plsRglmmodel.R
\name{print.coef.plsRglmmodel}
\alias{print.coef.plsRglmmodel}
\title{Print method for plsRglm models}
\usage{
\method{print}{coef.plsRglmmodel}(x, ...)
}
\arguments{
\item{x}{an object of the class \code{"coef.plsRglmmodel"}}

\item{\dots}{not used}
}
\value{
\code{NULL}
}
\description{
This function provides a print method for the class
\code{"coef.plsRglmmodel"}
}
\examples{

data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
modplsglm <- plsRglm(yCornell,XCornell,3,modele="pls-glm-family",family=gaussian())
class(modplsglm)
print(coef(modplsglm))
rm(list=c("XCornell","yCornell","modplsglm"))

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
\code{\link{print}}
}
\author{
Frédéric Bertrand\cr
\email{frederic.bertrand@utt.fr}\cr
\url{https://fbertran.github.io/homepage/}
}
\keyword{methods}
\keyword{print}
