\name{OLD_PLS_v2_vc}
\alias{OLD_PLS_v2_vc}
\title{Former plsRglm function}
\description{
This function in a beta version of the plsRglm with leave one out cross validation.
}
\usage{
OLD_PLS_v2_vc(dataY, dataX, nt = 2, limQ2set = 0.0975, dataPredictY = dataX, modele = "pls", family = NULL, typeVC = "none", EstimXNA = FALSE, scaleX = TRUE, scaleY = NULL, pvals.expli = FALSE, alpha.pvals.expli = 0.05)
}
\arguments{
  \item{dataY}{response (training) dataset}
  \item{dataX}{predictor(s) (training) dataset}
  \item{nt}{number of components to be extracted}
  \item{limQ2set}{limit value for the Q2}
  \item{dataPredictY}{predictor(s) (testing) dataset}
  \item{modele}{name of the PLS model to be fitted (\code{"pls"}, \code{"pls-glm-gaussian"}, \code{"pls-glm-logistic"}, \code{"pls-glm-polr"}).}
  \item{family}{for the present moment the family argumlent is ignored and set thanks to the value of modele.}
  \item{typeVC}{type of crossed validation (only leave-one-out now). Several procedures are available and may be forced.
    \describe{
      \item{\code{none}}{no crossed validation}
      \item{\code{standard}}{as in SIMCA for datasets without missing values and with all values predicted as those with missing values for datasets with any missing values}
      \item{\code{missingdata}}{all values predicted as those with missing values for datasets with any missing values}
      \item{\code{adaptative}}{predict a response value for an x with any missing value as those with missing values and for an x without any missing value as those without missing values.}
      }
  }
  \item{EstimXNA}{only for \code{modele="pls"}. Set whether the missing X values have to be estimated.}
  \item{scaleX}{scale the predictor(s) : must be set to TRUE for \code{modele="pls"} and should be for glms pls.}
  \item{scaleY}{scale the response : Yes/No. Ignored since non always possible for glm responses.}
  \item{pvals.expli}{should individual p-values be reported to tune model selection ?}
  \item{alpha.pvals.expli}{level of significance for predictors when pvals.expli=TRUE}
}
\details{
Cross validation not the same way for PLS than for PLS-glm : pls VC only on last components effect as in former SIMCA versions.\cr

There are four different models available :
\describe{
\item{\code{"pls"}}{ordinary pls models}
\item{\code{"pls-glm-gaussian"}}{glm gaussian with identity link pls models}
\item{\code{"pls-glm-logistic"}}{glm binomial with logit link pls models}
\item{\code{"pls-glm-polr"}}{glm polr with logit link pls models}
}
}
\value{
Depends on the chosen model.
}
\references{
Nicolas Meyer, Myriam Maumy-Bertrand et \enc{Frédéric}{Fr\'ed\'eric} Bertrand (2010). Comparaison de la \enc{régression}{r\'egression} PLS et de la \enc{régression}{r\'egression} logistique PLS : application aux \enc{données}{donn\'ees} \enc{d'allélotypage}{d'all\'elotypage}. \emph{Journal de la Société Française de Statistique}, 151(2), pages 1-18.
\url{http://smf4.emath.fr/Publications/JSFdS/151_2/pdf/sfds_jsfds_151_2_1-18.pdf}
}
\author{\enc{Frédéric}{Fr\'ed\'eric} Bertrand\cr
\email{frederic.bertrand@math.unistra.fr}\cr
\url{http://www-irma.u-strasbg.fr/~fbertran/}
}
\note{Should not be used with incomplete datasets. \code{plsRglm} and \code{PLS_glm_kfoldcv} must be preferred.}
\section{Warning}{Deprecated function.}
\seealso{\code{\link{plsR}}, \code{\link{plsRglm}} and \code{\link{PLS_glm_kfoldcv}}.}
\examples{\donttest{
data(pine)
Xpine<-pine[,1:10]
ypine<-pine[,11]
OLD_PLS_v2_vc(log(ypine),Xpine,4,typeVC="standard",modele="pls")$CVinfos
OLD_PLS_v2_vc(log(ypine),Xpine,4,typeVC="missingdata",modele="pls")$CVinfos

OLD_PLS_v2_vc(log(ypine),Xpine,4,typeVC="standard",modele="pls-glm-gaussian")$InfCrit
OLD_PLS_v2_vc(log(ypine),Xpine,4,typeVC="missingdata",modele="pls-glm-gaussian")$InfCrit


data(aze_compl)
Xaze_compl<-aze_compl[,2:34]
yaze_compl<-aze_compl$y
OLD_PLS_v2_vc(yaze_compl,Xaze_compl,10,modele="pls-glm-logistic",typeVC="none")$InfCrit
OLD_PLS_v2_vc(yaze_compl,Xaze_compl,10,modele="pls-glm-logistic",typeVC="standard")$InfCrit
OLD_PLS_v2_vc(yaze_compl,Xaze_compl,10,modele="pls-glm-logistic",typeVC="none",pvals.expli=TRUE)$valpvalstep


data(bordeaux)
Xbordeaux<-bordeaux[,1:4]
ybordeaux<-factor(bordeaux$Quality,ordered=TRUE)
OLD_PLS_v2_vc(as.numeric(ybordeaux),Xbordeaux,4,modele="pls",typeVC="standard")$CVinfos
OLD_PLS_v2_vc(ybordeaux,Xbordeaux,4,modele="pls-glm-polr",typeVC="none")$InfCrit


# plsR and gaussian plsRglm with missing data
XpineNAX21 <- Xpine
XpineNAX21[1,2] <- NA
cbind((log(ypine)),OLD_PLS_v2_vc(log(ypine),XpineNAX21,4,typeVC="none",modele="pls")$YChapeau,OLD_PLS_v2_vc(log(ypine),XpineNAX21,4,typeVC="none",modele="pls-glm-gaussian")$YChapeau)
cbind((log(ypine)),OLD_PLS_v2_vc(log(ypine),XpineNAX21,4,typeVC="none",modele="pls")$ValsPredictY,OLD_PLS_v2_vc(log(ypine),XpineNAX21,4,typeVC="none",modele="pls-glm-gaussian")$ValsPredictY)
rm("ypine","Xpine","XpineNAX21","yaze_compl","Xaze_compl","ybordeaux","Xbordeaux")
}
}\keyword{internal}
\keyword{models}
\keyword{regression}

