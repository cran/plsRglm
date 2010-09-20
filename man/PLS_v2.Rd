\name{PLS_v2}
\alias{PLS_v2}
\title{Partial least squares Regression generalized linear models}
\description{
This function implements Partial least squares Regression generalized linear models complete or incomplete datasets.
}
\usage{
PLS_v2(dataY, dataX, nt = 2, limQ2set = 0.0975, dataPredictY = dataX, modele = "pls", family = NULL, typeVC = "none", EstimXNA = FALSE, scaleX = TRUE, scaleY = NULL, pvals.expli = FALSE, alpha.pvals.expli = 0.05, MClassed = FALSE, tol_Xi = 10^(-12))
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
  \item{scaleY}{scale the response : Yes/No. Ignored since not always possible for glm responses.}
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
\note{
Use \code{plsRglm} instead.
}
\seealso{\code{\link{PLS_v2_wvc}} and \code{\link{PLS_v2_kfoldcv}}}
\examples{
data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
PLS_v2(yCornell,XCornell,3)$uscores
PLS_v2(yCornell,XCornell,3)$pp
PLS_v2(yCornell,XCornell,3)$Coeffs
PLS_v2(yCornell,XCornell,10)$InfCrit
PLS_v2(yCornell,XCornell,10,modele="pls-glm-gaussian")$InfCrit
rm(list=c("XCornell","yCornell"))


data(pine)
Xpine<-pine[,1:10]
ypine<-pine[,11]
PLS_v2(log(ypine),Xpine,1)$Std.Coeffs
PLS_v2(log(ypine),Xpine,1)$Coeffs
PLS_v2(log(ypine),Xpine,4)$Std.Coeffs
PLS_v2(log(ypine),Xpine,4)$Coeffs
PLS_v2(log(ypine),Xpine,4)$PredictY[1,]
PLS_v2(log(ypine),Xpine,4,dataPredictY=Xpine[1,])$PredictY[1,]

XpineNAX21 <- Xpine
XpineNAX21[1,2] <- NA
str(PLS_v2(log(ypine),XpineNAX21,2))
PLS_v2(log(ypine),XpineNAX21,4)$Std.Coeffs
PLS_v2(log(ypine),XpineNAX21,4)$YChapeau[1,]
PLS_v2(log(ypine),Xpine,4)$YChapeau[1,]
PLS_v2(log(ypine),XpineNAX21,4)$CoeffC
PLS_v2(log(ypine),XpineNAX21,4,EstimXNA=TRUE)$XChapeau
PLS_v2(log(ypine),XpineNAX21,4,EstimXNA=TRUE)$XChapeauNA

# compare pls-glm-gaussian with classic plsR
cbind(PLS_v2(log(ypine),Xpine,4,modele="pls")$Std.Coeffs,PLS_v2(log(ypine),Xpine,4,modele="pls-glm-gaussian")$Std.Coeffs)

# without missing data
cbind(log(ypine),PLS_v2(log(ypine),Xpine,4,modele="pls")$YChapeau,PLS_v2(log(ypine),Xpine,4,modele="pls-glm-gaussian")$YChapeau)
cbind(log(ypine),PLS_v2(log(ypine),XpineNAX21,4,modele="pls")$YChapeau,PLS_v2(log(ypine),XpineNAX21,4,modele="pls-glm-gaussian")$YChapeau)

# with missing data
cbind((log(ypine)),PLS_v2(log(ypine),XpineNAX21,4,modele="pls")$YChapeau,PLS_v2(log(ypine),XpineNAX21,4,modele="pls-glm-gaussian")$YChapeau)
cbind((log(ypine)),PLS_v2(log(ypine),XpineNAX21,4,modele="pls")$ValsPredictY,PLS_v2(log(ypine),XpineNAX21,4,modele="pls-glm-gaussian")$ValsPredictY)
rm(list=c("Xpine","ypine"))


\donttest{
data(fowlkes)
Xfowlkes <- fowlkes[,2:13]
yfowlkes <- fowlkes[,1]
modpls <- PLS_v2(yfowlkes,Xfowlkes,4,modele="pls-glm-logistic",pvals.expli=TRUE)
modpls$pvalstep
rm(list=c("Xfowlkes","yfowlkes","modpls"))
}

data(aze_compl)
Xaze_compl<-aze_compl[,2:34]
yaze_compl<-aze_compl$y
PLS_v2(yaze_compl,Xaze_compl,nt=10,modele="pls",MClassed=TRUE)$InfCrit
modpls <- PLS_v2(yaze_compl,Xaze_compl,nt=10,modele="pls-glm-logistic",MClassed=TRUE,pvals.expli=TRUE)
modpls$InfCrit
modpls$valpvalstep
modpls$Coeffsmodel_vals

plot(PLS_v2(yaze_compl,Xaze_compl,4,modele="pls-glm-logistic")$FinalModel)
PLS_v2(yaze_compl[-c(99,72)],Xaze_compl[-c(99,72),],4,modele="pls-glm-logistic",pvals.expli=TRUE)$pvalstep
plot(PLS_v2(yaze_compl[-c(99,72)],Xaze_compl[-c(99,72),],4,modele="pls-glm-logistic",pvals.expli=TRUE)$FinalModel)
rm(list=c("Xaze_compl","yaze_compl","modpls"))


data(bordeaux)
Xbordeaux<-bordeaux[,1:4]
ybordeaux<-factor(bordeaux$Quality,ordered=TRUE)
modpls <- PLS_v2(ybordeaux,Xbordeaux,10,modele="pls-glm-polr")
modpls$Coeffsmodel_vals
modpls$InfCrit

XbordeauxNA<-Xbordeaux
XbordeauxNA[1,1] <- NA
modplsNA <- PLS_v2(ybordeaux,XbordeauxNA,10,modele="pls-glm-polr")
modplsNA$Coeffsmodel_vals
modplsNA$InfCrit
rm(list=c("Xbordeaux","XbordeauxNA","ybordeaux"))


\donttest{
dimX <- 6
Astar <- 4
dataAstar4 <- t(replicate(250,simul_data_UniYX(dimX,Astar)))
ysimbin1 <- dicho(dataAstar4)[,1]
Xsimbin1 <- dicho(dataAstar4)[,2:(dimX+1)]
modplsglm <- PLS_v2(ysimbin1,Xsimbin1,10,modele="pls-glm-logistic")
modplsglm$computed_nt
modplsglm$InfCrit
rm(list=c("dimX","Astar","dataAstar4","ysimbin1","Xsimbin1","modplsglm"))


dimX <- 24
Astar <- 2
dataAstar2 <- t(replicate(250,simul_data_UniYX(dimX,Astar)))
ysimbin1 <- dicho(dataAstar2)[,1]
Xsimbin1 <- dicho(dataAstar2)[,2:(dimX+1)]
modplsglm <- PLS_v2(ysimbin1,Xsimbin1,10,modele="pls-glm-logistic")
modplsglm$computed_nt
modplsglm$InfCrit
rm(list=c("dimX","Astar","dataAstar2","ysimbin1","Xsimbin1","modplsglm"))


dimX <- 24
Astar <- 3
dataAstar3 <- t(replicate(250,simul_data_UniYX(dimX,Astar)))
ysimbin1 <- dicho(dataAstar3)[,1]
Xsimbin1 <- dicho(dataAstar3)[,2:(dimX+1)]
modplsglm <- PLS_v2(ysimbin1,Xsimbin1,10,modele="pls-glm-logistic")
modplsglm$computed_nt
modplsglm$InfCrit
rm(list=c("dimX","Astar","dataAstar3","ysimbin1","Xsimbin1","modplsglm"))


dimX <- 24
Astar <- 4
dataAstar4 <- t(replicate(250,simul_data_UniYX(dimX,Astar)))
ysimbin1 <- dicho(dataAstar4)[,1]
Xsimbin1 <- dicho(dataAstar4)[,2:(dimX+1)]
modplsglm <- PLS_v2(ysimbin1,Xsimbin1,10,modele="pls-glm-logistic")
modplsglm$computed_nt
modplsglm$InfCrit
rm(list=c("dimX","Astar","dataAstar4","ysimbin1","Xsimbin1","modplsglm"))


dimX <- 24
Astar <- 5
dataAstar5 <- t(replicate(250,simul_data_UniYX(dimX,Astar)))
ysimbin1 <- dicho(dataAstar5)[,1]
Xsimbin1 <- dicho(dataAstar5)[,2:(dimX+1)]
modplsglm <- PLS_v2(ysimbin1,Xsimbin1,10,modele="pls-glm-logistic")
modplsglm$computed_nt
modplsglm$InfCrit
rm(list=c("dimX","Astar","dataAstar5","ysimbin1","Xsimbin1","modplsglm"))


dimX <- 24
Astar <- 6
dataAstar6 <- t(replicate(250,simul_data_UniYX(dimX,Astar)))
ysimbin1 <- dicho(dataAstar6)[,1]
Xsimbin1 <- dicho(dataAstar6)[,2:(dimX+1)]
modplsglm <- PLS_v2(ysimbin1,Xsimbin1,10,modele="pls-glm-logistic")
modplsglm$computed_nt
modplsglm$InfCrit
rm(list=c("dimX","Astar","dataAstar6","ysimbin1","Xsimbin1","modplsglm"))
}
}\keyword{models}
\keyword{regression}
