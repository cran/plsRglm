\name{plsRglm-package}
\alias{plsRglm-package}
\docType{package}
\title{
Partial least squares Regression for generalized linear models
}
\description{
This packages provides Partial least squares Regression for generalized linear models and kfold crossvalidation of such models using various criteria. It works either for complete datasets or datasets with missing values.
}
\details{
\tabular{ll}{
Package: \tab plsRglm\cr
Version: \tab 0.1.2\cr
Date: \tab 2009-09-06\cr
Requires: \tab \cr
Enhances: \tab \cr
Suggests: \tab MASS, mvtnorm\cr
Depends: \tab R (>= 2.4.0)\cr
License: \tab GPL-3\cr
Encoding: \tab latin1\cr
URL: \tab http://www-irma.u-strasbg.fr/~fbertran/\cr
Classification/MSC: \tab 62J12, 62J99\cr
Packaged: \tab Tue Sep  8 01:22:14 2009; Administrateur\cr
Built: \tab R 2.8.0; ; 2009-09-08 01:22:15; windows\cr
}

Index:
\preformatted{
AICpls                  ~~function to do ... ~~
AICpls2                 ~~function to do ... ~~
aze                     ~~ data name/kind ... ~~
aze_compl               ~~ data name/kind ... ~~
bordeaux                ~~ data name/kind ... ~~
CorMat                  ~~ data name/kind ... ~~
Cornell                 ~~ data name/kind ... ~~
dicho                   ~~function to do ... ~~
fowlkes                 ~~ data name/kind ... ~~
kfolds2Chisq            ~~function to do ... ~~
kfolds2Chisqind         ~~function to do ... ~~
kfolds2coeff            ~~function to do ... ~~
kfolds2CVinfos_glm      ~~function to do ... ~~
kfolds2CVinfos_lm       ~~function to do ... ~~
kfolds2CVinfos_v1       ~~function to do ... ~~
kfolds2CVinfos_v2       ~~function to do ... ~~
kfolds2Press            ~~function to do ... ~~
kfolds2Pressind         ~~function to do ... ~~
loglikpls               ~~function to do ... ~~
loglikpls2              ~~function to do ... ~~
pine                    ~~ data name/kind ... ~~
pine_full               ~~ data name/kind ... ~~
pine_sup                ~~ data name/kind ... ~~
PLS_glm                 ~~function to do ... ~~
PLS_glm_kfoldcv         ~~function to do ... ~~
PLS_glm_wvc             ~~function to do ... ~~
PLS_lm                  ~~function to do ... ~~
PLS_lm_kfoldcv          ~~function to do ... ~~
PLS_lm_wvc              ~~function to do ... ~~
PLS_v1                  ~~function to do ... ~~
PLS_v1_kfoldcv          ~~function to do ... ~~
PLS_v1_wvc              ~~function to do ... ~~
PLS_v2                  ~~function to do ... ~~
PLS_v2_kfoldcv          ~~function to do ... ~~
PLS_v2_wvc              ~~function to do ... ~~
plsR                    ~~function to do ... ~~
plsRglm                 ~~function to do ... ~~
plsRglmmodel.default    ~~function to do ... ~~
plsRmodel.default       ~~function to do ... ~~
print.plsRglmmodel      ~~function to do ... ~~
print.plsRmodel         ~~function to do ... ~~
print.summary.plsRglmmodel
                        ~~function to do ... ~~
print.summary.plsRmodel
                        ~~function to do ... ~~
simul_data_complete     ~~function to do ... ~~
simul_data_UniYX        ~~function to do ... ~~
simul_data_UniYX_old    ~~function to do ... ~~
simul_data_YX           ~~function to do ... ~~
summary.plsRglmmodel    ~~function to do ... ~~
summary.plsRmodel       ~~function to do ... ~~
XbordeauxNA             ~~ data name/kind ... ~~
XpineNAX21              ~~ data name/kind ... ~~
}
~~ An overview of how to use the package, including the most important functions ~~
}
\author{
\enc{Frédéric}{Fr\'ed\'eric} Bertrand <frederic.bertrand@math.unistra.fr>, Nicolas Meyer <Nicolas.Meyer@medecine.u-strasbg.fr>, Myriam Maumy-Bertrand <myriam.maumy-bertrand@math.unistra.fr>.

Maintainer: \enc{Frédéric}{Fr\'ed\'eric} Bertrand <frederic.bertrand@math.unistra.fr>
}
\references{
Nicolas Meyer, Myriam Maumy-Bertrand et \enc{Frédéric}{Fr\'ed\'eric} Bertrand (2010). Comparaison de la \enc{régression}{r\'egression} PLS et de la \enc{régression}{r\'egression} logistique PLS : application aux \enc{données}{donn\'ees} \enc{d'allélotypage}{d'all\'elotypage}. \emph{Journal de la Société Française de Statistique}, 151(2), pages 1-18.
\url{http://smf4.emath.fr/Publications/JSFdS/151_2/pdf/sfds_jsfds_151_2_1-18.pdf}
}
\examples{
data(pine)
ypine <- pine[,11]
Xpine <- pine[,1:10]
(Pinscaled <- as.data.frame(cbind(scale(log(ypine)),scale(as.matrix(Xpine)))))
colnames(Pinscaled)[1] <- "yy"

modpls <- plsR(log(ypine),Xpine,10)
modpls$Std.Coeffs
}
\keyword{package}
\seealso{\code{\link[pls:mvr]{pls}}}
