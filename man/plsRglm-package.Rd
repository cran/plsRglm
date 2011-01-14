\name{plsRglm-package}
\alias{plsRglm-package}
\docType{package}
\title{
Partial least squares Regression for generalized linear models
}
\description{
This packages provides Partial least squares Regression for generalized linear models and kfold crossvalidation of such models using various criteria. It allows for missing data in the eXplanatory variables. Bootstrap confidence intervals constructions are also available.
}
\details{
\tabular{ll}{
Package: \tab plsRglm\cr
Version: \tab 0.6.3\cr
Date: \tab 2011-01-08\cr
Depends: \tab R (>= 2.4.0)\cr
Imports: \tab mvtnorm, boot\cr
Enhances: \tab pls\cr
Suggests: \tab MASS\cr
License: \tab GPL-3\cr
Encoding: \tab latin1\cr
URL: \tab http://www-irma.u-strasbg.fr/~fbertran/\cr
Classification/MSC: \tab 62J12, 62J99\cr
Packaged: \tab 2011-01-08 14:47:17 UTC; Bertrand\cr
Built: \tab R 2.11.1; ; 2011-01-08 14:47:17 UTC; windows\cr
}

Index:
\preformatted{
AICpls                  AIC functions for plsR models
aze                     Microsat Dataset
aze_compl               As aze without missing values
bootpls                 Non-parametric Bootstrap for PLS models
bootplsglm              Non-parametric Bootstrap for PLS generalized
                        linear models
bordeaux                Quality of wine dataset
boxplots.bootpls        Boxplot bootstrap distributions
coefs.plsR              Coefficients for bootstrap computations
coefs.plsRglm           Coefficients for bootstrap computations
confints.bootpls        Bootstrap confidence intervals
CorMat                  Correlation matrix for simulating plsR datasets
Cornell                 Cornell dataset
dicho                   Dichotomization
fowlkes                 Fowlkes dataset
kfolds2Chisq            Computes Predicted Chisquare for kfold cross
                        validated partial least squares regression
                        models.
kfolds2Chisqind         Computes individual Predicted Chisquare for
                        kfold cross validated partial least squares
                        regression models.
kfolds2coeff            Extracts coefficients from kfold cross
                        validated partial least squares regression
                        models
kfolds2CVinfos_glm      Extracts and computes information criteria and
                        fits statistics for kfold cross validated
                        partial least squares glm models
kfolds2CVinfos_v2       Extracts and computes information criteria and
                        fits statistics for kfold cross validated
                        partial least squares glm models
kfolds2Mclassed         Number of missclassified individuals for kfold
                        cross validated partial least squares
                        regression models.
kfolds2Mclassedind      Number of missclassified individuals per group
                        for kfold cross validated partial least squares
                        regression models.
kfolds2Press            Computes PRESS for kfold cross validated
                        partial least squares regression models.
kfolds2Pressind         Computes individual PRESS for kfold cross
                        validated partial least squares regression
                        models.
loglikpls               loglikelihood functions for plsR models
permcoefs.plsR          Coefficients computation for permutation
                        bootstrap
permcoefs.plsRglm       Coefficients computation for permutation
                        bootstrap
pine                    Pine dataset
pine_full               Full pine dataset
pine_sup                Supplementary data for pine dataset
plots.confints.bootpls
                        Plot bootstrap confidence intervals
PLS_glm                 Partial least squares Regression generalized
                        linear models
PLS_glm_formula         Partial least squares Regression generalized
                        linear models
PLS_glm_kfoldcv         Partial least squares regression glm models
                        with kfold cross validation
PLS_glm_kfoldcv_formula
                        Partial least squares regression glm models
                        with kfold cross validation
PLS_glm_wvc             Light version of PLS\_glm for cross validation
                        purposes
PLS_lm                  Partial least squares Regression models with
                        leave one out cross validation
PLS_lm_formula          Partial least squares Regression models with
                        leave one out cross validation
PLS_lm_kfoldcv          Partial least squares regression models with
                        kfold cross validation
PLS_lm_kfoldcv_formula
                        Partial least squares regression models with
                        kfold cross validation
PLS_lm_wvc              Light version of PLS\_lm for cross validation
                        purposes
PLS_v1                  Partial least squares Regression models with
                        leave one out cross validation
PLS_v1_kfoldcv          Partial least squares regression models with
                        kfold cross validation
PLS_v1_wvc              Light version of PLS\_v1 for cross validation
                        purposes
PLS_v2                  Partial least squares Regression generalized
                        linear models
PLS_v2_kfoldcv          Partial least squares regression glm models
                        with kfold cross validation
PLS_v2_wvc              Light version of PLS\_v2 for cross validation
                        purposes
plsR                    Partial least squares Regression models with
                        leave one out cross validation
plsRglm                 Partial least squares Regression generalized
                        linear models
plsRglm-package         Partial least squares Regression for
                        generalized linear models
print.plsRglmmodel      Print method for plsRglm models
print.plsRmodel         Print method for plsR models
print.summary.plsRglmmodel
                        Print method for summaries of plsRglm models
print.summary.plsRmodel
                        Print method for summaries of plsR models
simul_data_complete     Data generating detailed process for
                        multivariate plsR models
simul_data_UniYX        Data generating function for univariate plsR
                        models
simul_data_YX           Data generating function for multivariate plsR
                        models
summary.plsRglmmodel    Summary method for plsRglm models
summary.plsRmodel       Summary method for plsR models
tilt.bootpls            Tilted bootstrap for PLS models
tilt.bootplsglm         Tilted bootstrap for PLS models
XbordeauxNA             Missing data analysis for the quality of wine
                        dataset
XpineNAX21              Missing data analysis for the pine dataset
}
}
\author{
Frederic Bertrand <frederic.bertrand@math.unistra.fr>, Nicolas Meyer <Nicolas.Meyer@nmeyer@unistra.fr>, Myriam Maumy-Bertrand <myriam.maumy-bertrand@math.unistra.fr>.

Maintainer: Frederic Bertrand <frederic.bertrand@math.unistra.fr>
}
\references{
Nicolas Meyer, Myriam Maumy-Bertrand et \enc{Frédéric}{Fr\'ed\'eric} Bertrand (2010). Comparaison de la \enc{régression}{r\'egression} PLS et de la \enc{régression}{r\'egression} logistique PLS : application aux \enc{données}{donn\'ees} \enc{d'allélotypage}{d'all\'elotypage}. \emph{Journal de la Société Française de Statistique}, 151(2), pages 1-18.
\url{http://smf4.emath.fr/Publications/JSFdS/151_2/pdf/sfds_jsfds_151_2_1-18.pdf}
}
\keyword{ package }
\examples{
data(pine)
ypine <- pine[,11]
Xpine <- pine[,1:10]
(Pinscaled <- as.data.frame(cbind(scale(log(ypine)),scale(as.matrix(Xpine)))))
colnames(Pinscaled)[1] <- "yy"

modpls <- plsR(log(ypine),Xpine,10)
modpls$Std.Coeffs
}
