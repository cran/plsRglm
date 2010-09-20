\name{plsRglmmodel.default}
\alias{plsRglmmodel.default}
\title{Partial least squares Regression generalized linear models}
\description{
This function implements Partial least squares Regression generalized linear models complete or incomplete datasets.
}
\usage{
plsRglmmodel.default(dataY, dataX, nt = 2, limQ2set = 0.0975, dataPredictY = dataX, modele = "pls", family = NULL, typeVC = "none", EstimXNA = FALSE, scaleX = TRUE, scaleY = NULL, pvals.expli = FALSE, alpha.pvals.expli = 0.05, MClassed = FALSE, tol_Xi = 10^(-12))
}
\arguments{
  \item{dataY}{response (training) dataset}
  \item{dataX}{predictor(s) (training) dataset}
  \item{nt}{number of components to be extracted}
  \item{limQ2set}{limit value for the Q2}
  \item{dataPredictY}{predictor(s) (testing) dataset}
  \item{modele}{name of the PLS glm model to be fitted (\code{"pls"}, \code{"pls-glm-gaussian"}, \code{"pls-glm-logistic"}, \code{"pls-glm-polr"}).}
  \item{family}{for the present moment the family argumlent is ignored and set thanks to the value of modele.}
  \item{typeVC}{type of leave one out crossed validation. For back compatibility purpose.
    \describe{
      \item{\code{none}}{no crossed validation}
      \item{\code{standard}}{no crossed validation}
      \item{\code{missingdata}}{no crossed validation}
      \item{\code{adaptative}}{no crossed validation}
      }
  }
  \item{EstimXNA}{only for \code{modele="pls"}. Set whether the missing X values have to be estimated.}
  \item{scaleX}{scale the predictor(s) : must be set to TRUE for \code{modele="pls"} and should be for glms pls.}
  \item{scaleY}{scale the response : Yes/No. Ignored since non always possible for glm responses.}
  \item{pvals.expli}{should individual p-values be reported to tune model selection ?}
  \item{alpha.pvals.expli}{level of significance for predictors when pvals.expli=TRUE}
  \item{MClassed}{number of missclassed cases, should only be used for binary responses}
  \item{tol_Xi}{minimal value for Norm2(Xi) and \eqn{\mathrm{det}(pp' \times pp)}{det(pp'*pp)} if there is any missing value in the \code{dataX}. It defaults to \eqn{10^{-12}}{10^{-12}}}
}
\details{
There are four different models available :
\describe{
\item{\code{"pls"}}{ordinary pls models}
\item{\code{"pls-glm-gaussian"}}{glm gaussian with identity link pls models}
\item{\code{"pls-glm-logistic"}}{glm binomial with logit link pls models}
\item{\code{"pls-glm-polr"}}{glm polr with logit link pls models}
}
}
\value{
Depends on the model that was used to fit the model.
}
\references{
Nicolas Meyer, Myriam Maumy-Bertrand et \enc{Frédéric}{Fr\'ed\'eric} Bertrand (2010). Comparaison de la \enc{régression}{r\'egression} PLS et de la \enc{régression}{r\'egression} logistique PLS : application aux \enc{données}{donn\'ees} \enc{d'allélotypage}{d'all\'elotypage}. \emph{Journal de la Société Française de Statistique}, 151(2), pages 1-18.
\url{http://smf4.emath.fr/Publications/JSFdS/151_2/pdf/sfds_jsfds_151_2_1-18.pdf}
}
\author{\enc{Frédéric}{Fr\'ed\'eric} Bertrand\cr
\email{frederic.bertrand@math.unistra.fr}\cr
\url{http://www-irma.u-strasbg.fr/~fbertran/}
}
\note{Use \code{plsRglm} instead.}
\seealso{\code{\link{plsRmodel.default}}}
\examples{
data(aze_compl)
Xaze_compl<-aze_compl[,2:34]
yaze_compl<-aze_compl$y
plsRglmmodel.default(yaze_compl,Xaze_compl,nt=10,modele="pls",MClassed=TRUE)$InfCrit
modpls <- plsRglm(yaze_compl,Xaze_compl,nt=10,modele="pls-glm-logistic",MClassed=TRUE,pvals.expli=TRUE)
modpls$InfCrit
rm(list=c("Xaze_compl","yaze_compl","modpls"))


data(bordeaux)
Xbordeaux<-bordeaux[,1:4]
ybordeaux<-factor(bordeaux$Quality,ordered=TRUE)
modpls <- plsRglmmodel.default(ybordeaux,Xbordeaux,10,modele="pls-glm-polr")
modpls$InfCrit

XbordeauxNA<-Xbordeaux
XbordeauxNA[1,1] <- NA
modplsNA <- plsRglmmodel.default(ybordeaux,XbordeauxNA,10,modele="pls-glm-polr")
modplsNA$InfCrit
rm(list=c("Xbordeaux","XbordeauxNA","ybordeaux"))
}
\keyword{models}
\keyword{regression}
