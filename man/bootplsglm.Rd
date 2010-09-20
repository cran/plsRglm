\name{bootplsglm}
\alias{bootplsglm}
\title{Non-parametric Bootstrap for PLS generalized linear models}
\description{
Provides a wrapper for the bootstrap function \code{boot} from the \code{boot} R package.\cr Implements non-parametric bootstrap for PLS generalized linear models by case resampling.
}
\usage{
bootplsglm(object, typeboot="plsmodel", R=250, statistic=coefs.plsRglm, sim="ordinary", stype="i")
}
\arguments{
  \item{object}{ ~~Explain \code{object} here~~ }
  \item{typeboot}{ ~~Explain \code{typeboot} here~~ }
  \item{R}{ ~~Explain \code{R} here~~ }
  \item{statistic}{ ~~Explain \code{statistic} here~~ }
  \item{sim}{ ~~Explain \code{sim} here~~ }
  \item{stype}{ ~~Explain \code{stype} here~~ }
}
\details{
  ~~ More details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{\enc{Frédéric}{Fr\'ed\'eric} Bertrand\cr
\email{frederic.bertrand@math.unistra.fr}\cr
\url{http://www-irma.u-strasbg.fr/~fbertran/}
}
\note{ ~~some notes~~ 
}
\seealso{\code{\link[boot:boot]{boot}}}
\examples{\donttest{
data(aze_compl)
Xaze_compl<-aze_compl[,2:34]
yaze_compl<-aze_compl$y

dataset <- cbind(y=yaze_compl,Xaze_compl)

library(boot)
# Lazraq-Cléroux PLS bootstrap Classic

aze_compl.boot2 <- boot(data=dataset, statistic=coefs.plsRglm, sim="ordinary", stype="i", R=250, nt=3, modele="pls-glm-logistic")
aze_compl.boot <- bootplsglm(plsRglm(yaze_compl,Xaze_compl,3,modele="pls-glm-logistic"), sim="ordinary", stype="i", R=250)
boot.ci(aze_compl.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=1)
boot.ci(aze_compl.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=2)
boot.ci(aze_compl.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=3)
boot.ci(aze_compl.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=4)
boot.ci(aze_compl.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=5)
boot.ci(aze_compl.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=6)

boxplots.bootpls(aze_compl.boot)
confints.bootpls(aze_compl.boot)
plots.confints.bootpls(confints.bootpls(aze_compl.boot))


plot(aze_compl.boot,index=2)
jack.after.boot(aze_compl.boot, index=2, useJ=TRUE, nt=3)
plot(aze_compl.boot, index=2,jack=TRUE)
# tilt.boot(data=dataset, statistic=coefs.plsRglm, R=c(499, 100, 100), alpha=c(0.025, 0.975), sim="ordinary", stype="i", index=1, nt=3, modele="pls-glm-logistic")
aze_compl.tilt.boot <- tilt.bootplsglm(plsRglm(yaze_compl,Xaze_compl,3, modele="pls-glm-logistic"), statistic=coefs.plsR, R=c(499, 100, 100), alpha=c(0.025, 0.975), sim="ordinary", stype="i", index=1)

# PLS bootstrap balanced

aze_compl.boot <- boot(data=dataset, statistic=coefs.plsRglm, sim="balanced", stype="i", R=250, nt=3, modele="pls-glm-logistic")
aze_compl.boot <- bootplsglm(plsRglm(yaze_compl,Xaze_compl,3,modele="pls-glm-logistic"), sim="balanced", stype="i", R=250)
boot.ci(aze_compl.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=1)
boot.ci(aze_compl.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=2)
boot.ci(aze_compl.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=3)
boot.ci(aze_compl.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=4)
boot.ci(aze_compl.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=5)
boot.ci(aze_compl.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=6)


boxplots.bootpls(aze_compl.boot)
confints.bootpls(aze_compl.boot)
plots.confints.bootpls(confints.bootpls(aze_compl.boot))



plot(aze_compl.boot)
jack.after.boot(aze_compl.boot, index=1, useJ=TRUE, nt=3)
plot(aze_compl.boot,jack=TRUE)
# tilt.boot(data=dataset, statistic=coefs.plsRglm, R=c(499, 250, 250), alpha=c(0.025, 0.975), sim="balanced", stype="i", index=1, nt=3, modele="pls-glm-logistic")
aze_compl.tilt.boot <- tilt.bootplsglm(plsRglm(yaze_compl,Xaze_compl,3, modele="pls-glm-logistic"), statistic=coefs.plsR, R=c(499, 100, 100), alpha=c(0.025, 0.975), sim="balanced", stype="i", index=1)


# PLS permutation bootstrap

aze_compl.boot <- boot(data=dataset, statistic=permcoefs.plsRglm, sim="permutation", stype="i", R=250, nt=3, modele="pls-glm-logistic")
aze_compl.boot <- bootplsglm(plsRglm(yaze_compl,Xaze_compl,3,modele="pls-glm-logistic"), sim="permutation", stype="i", R=250)
boot.ci(aze_compl.boot, conf = c(0.90,0.95), type = c("norm","basic","perc"), index=1)
boot.ci(aze_compl.boot, conf = c(0.90,0.95), type = c("norm","basic","perc"), index=2)
boot.ci(aze_compl.boot, conf = c(0.90,0.95), type = c("norm","basic","perc"), index=3)
boot.ci(aze_compl.boot, conf = c(0.90,0.95), type = c("norm","basic","perc"), index=4)
boot.ci(aze_compl.boot, conf = c(0.90,0.95), type = c("norm","basic","perc"), index=5)
boot.ci(aze_compl.boot, conf = c(0.90,0.95), type = c("norm","basic","perc"), index=6)


boxplots.bootpls(aze_compl.boot)


plot(aze_compl.boot)



data(aze)
Xaze<-aze[,2:34]
yaze<-aze$y

dataset <- cbind(y=yaze,Xaze)

library(boot)
aze.boot2 <- boot(data=dataset, statistic=coefs.plsRglm, sim="ordinary", stype="i", R=250, nt=3, modele="pls-glm-logistic")
aze.boot <- bootplsglm(plsRglm(yaze,Xaze,3,modele="pls-glm-logistic"), sim="ordinary", stype="i", R=250)
boot.ci(aze.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=1)
boot.ci(aze.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=2)
boot.ci(aze.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=3)
boot.ci(aze.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=4)
boot.ci(aze.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=5)
boot.ci(aze.boot, conf = c(0.90,0.95), type = c("norm","basic","perc","bca"), index=6)


boxplots.bootpls(aze_compl.boot)
confints.bootpls(aze_compl.boot)
plots.confints.bootpls(confints.bootpls(aze_compl.boot))



plot(aze.boot)
jack.after.boot(aze.boot, index=1, useJ=TRUE, nt=3)
plot(aze.boot,jack=TRUE)
# tilt.boot(data=dataset, statistic=coefs.plsRglm, R=c(499, 100, 100), alpha=c(0.025, 0.975), sim="ordinary", stype="i", index=1, nt=3, modele="pls-glm-logistic")
aze.tilt.boot <- tilt.bootplsglm(plsRglm(yaze,Xaze,3, modele="pls-glm-logistic"), statistic=coefs.plsR, R=c(499, 100, 100), alpha=c(0.025, 0.975), sim="ordinary", stype="i", index=1)


data(bordeaux)
Xbordeaux<-bordeaux[,1:4]
ybordeaux<-factor(bordeaux$Quality,ordered=TRUE)

dataset <- cbind(y=ybordeaux,Xbordeaux)

# PLS permutation bootstrap

library(boot)
bordeaux.boot2 <- boot(data=dataset, statistic=permcoefs.plsRglm, sim="permutation", stype="i", R=250, nt=3, modele="pls-glm-polr")
bordeaux.boot <- bootplsglm(plsRglm(ybordeaux,Xbordeaux,3,modele="pls-glm-polr"), sim="permutation", stype="i", R=250)
boot.ci(bordeaux.boot, conf = c(0.90,0.95), type = c("norm","basic","perc"), index=1)
boot.ci(bordeaux.boot, conf = c(0.90,0.95), type = c("norm","basic","perc"), index=2)
boot.ci(bordeaux.boot, conf = c(0.90,0.95), type = c("norm","basic","perc"), index=3)
boot.ci(bordeaux.boot, conf = c(0.90,0.95), type = c("norm","basic","perc"), index=4)
boot.ci(bordeaux.boot, conf = c(0.90,0.95), type = c("norm","basic","perc"), index=5)
boot.ci(bordeaux.boot, conf = c(0.90,0.95), type = c("norm","basic","perc"), index=6)

boxplots.bootpls(bordeaux.boot)



plot(bordeaux.boot)
}
}\keyword{models}
