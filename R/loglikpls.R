#' loglikelihood function for plsR models
#' 
#' This function provides loglikelihood computation for an univariate plsR
#' model.
#' 
#' Loglikelihood functions for plsR models with univariate response.
#' 
#' @param residpls Residuals of a fitted univariate plsR model
#' @param weights Weights of observations
#' @return \item{real}{Loglikelihood value}
#' @author Frédéric Bertrand\cr
#' \email{frederic.bertrand@@utt.fr}\cr
#' \url{https://fbertran.github.io/homepage/}
#' @seealso \code{\link{AICpls}} for AIC computation and \code{\link{logLik}}
#' for loglikelihood computations for linear models
#' @references Baibing Li, Julian Morris, Elaine B. Martin, Model selection for
#' partial least squares regression, Chemometrics and Intelligent Laboratory
#' Systems 64 (2002) 79-89, \doi{10.1016/S0169-7439(02)00051-5}.
#' @keywords models regression utilities
#' @examples
#' 
#' data(pine)
#' ypine <- pine[,11]
#' Xpine <- pine[,1:10]
#' (Pinscaled <- as.data.frame(cbind(scale(ypine),scale(as.matrix(Xpine)))))
#' colnames(Pinscaled)[1] <- "yy"
#' 
#' lm(yy~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10,data=Pinscaled)
#' 
#' modpls <- plsR(ypine,Xpine,10)
#' modpls$Std.Coeffs
#' lm(yy~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10,data=Pinscaled)
#' 
#' AIC(lm(yy~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10,data=Pinscaled))
#' print(logLik(lm(yy~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10,data=Pinscaled)))
#' 
#' sum(dnorm(modpls$RepY, modpls$Std.ValsPredictY, sqrt(mean(modpls$residY^2)), log=TRUE))
#' sum(dnorm(Pinscaled$yy,fitted(lm(yy~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10,data=Pinscaled)),
#' sqrt(mean(residuals(lm(yy~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10,data=Pinscaled))^2)), log=TRUE))
#' loglikpls(modpls$residY)
#' loglikpls(residuals(lm(yy~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10,data=Pinscaled)))
#' AICpls(10,residuals(lm(yy~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10,data=Pinscaled)))
#' AICpls(10,modpls$residY)
#' 
#' @export loglikpls
loglikpls <- function(residpls,weights=rep.int(1,length(residpls))) {
if(is.null(weights)){rep.int(1,length(residpls))}
nn=length(residpls)
return(1/2*(sum(log(weights))-nn*(log(2*pi) + 1 - log(nn) + log(sum(weights*residpls^2)))))}
