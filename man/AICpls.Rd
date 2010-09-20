\name{AICpls}
\alias{AICpls}
\alias{AICpls2}
\title{AIC functions for plsR models}
\description{
This function provides AIC computations for an univariate plsR model.
}
\usage{
AICpls(nn, ncomp, residpls)
AICpls2(ncomp, dataYpls, YChapeaupls, residpls)
}
\arguments{
  \item{nn}{Number of rows of the dataset}
  \item{ncomp}{Number of components}
  \item{dataYpls}{Vector of the response}
  \item{YChapeaupls}{Fitted values of a plsR univariate plsR model}
  \item{residpls}{Residuals of a fitted univariate plsR model}
}
\details{
AIC functions for plsR models with univariate response.
}
\value{
  \item{real}{AIC value}
}
\references{
Baibing Li, Julian Morris, Elaine B. Martin,
Model selection for partial least squares regression, 
Chemometrics and Intelligent Laboratory Systems 64 (2002) 79-89.
\url{http://dx.doi.org/10.1016/S0169-7439(02)00051-5}
}
\author{\enc{Frédéric}{Fr\'ed\'eric} Bertrand\cr
\email{frederic.bertrand@math.unistra.fr}\cr
\url{http://www-irma.u-strasbg.fr/~fbertran/}
}
\seealso{\code{\link{loglikpls}} and \code{\link{loglikpls2}} and \code{\link{AIC}} for loglikelihood computations for linear models}
\examples{
data(pine)
ypine <- pine[,11]
Xpine <- pine[,1:10]
(Pinscaled <- as.data.frame(cbind(scale(log(ypine)),scale(as.matrix(Xpine)))))
colnames(Pinscaled)[1] <- "yy"

lm(yy~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10,data=Pinscaled)

modpls <- plsR(log(ypine),Xpine,10)
modpls$Std.Coeffs
lm(yy~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10,data=Pinscaled)

AIC(lm(yy~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10,data=Pinscaled))
print(logLik(lm(yy~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10,data=Pinscaled)))

sum(dnorm(modpls$RepY, modpls$Std.ValsPredictY, sqrt(mean(modpls$residY^2)), log=TRUE))
sum(dnorm(Pinscaled$yy,fitted(lm(yy~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10,data=Pinscaled)),sqrt(mean(residuals(lm(yy~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10,data=Pinscaled))^2)), log=TRUE))
loglikpls2(modpls$RepY,modpls$Std.ValsPredictY,modpls$residY)
loglikpls2(Pinscaled$yy,fitted(lm(yy~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10,data=Pinscaled)),residuals(lm(yy~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10,data=Pinscaled)))
loglikpls(33,10,modpls$residY)
loglikpls(33,10,residuals(lm(yy~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10,data=Pinscaled)))
AICpls(33,10,residuals(lm(yy~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10,data=Pinscaled)))
AICpls(33,10,modpls$residY)
}
\keyword{models}
\keyword{regression}
\keyword{utilities}
