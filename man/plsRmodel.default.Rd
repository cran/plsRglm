\name{plsRmodel.default}
\alias{plsRmodel.default}
\title{Partial least squares Regression models with leave one out cross validation}
\description{
This function implements Partial least squares Regression models with leave one out cross validation for complete or incomplete datasets.
}
\usage{
plsRmodel.default(dataY, dataX, nt = 2, limQ2set = 0.0975, dataPredictY = dataX, modele = "pls", family = NULL, typeVC = "none", EstimXNA = FALSE, scaleX = TRUE, scaleY = NULL, pvals.expli = FALSE, alpha.pvals.expli = 0.05, MClassed = FALSE, tol_Xi = 10^(-12))
}
\arguments{
  \item{dataY}{response (training) dataset}
  \item{dataX}{predictor(s) (training) dataset}
  \item{nt}{number of components to be extracted}
  \item{limQ2set}{limit value for the Q2}
  \item{dataPredictY}{predictor(s) (testing) dataset}
  \item{modele}{name of the PLS model to be fitted, only (\code{"pls"} available for this fonction.}
  \item{family}{for the present moment the family argument is ignored and set thanks to the value of modele.}
  \item{typeVC}{type of leave one out crossed validation. Several procedures are available and may be forced.
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
  \item{MClassed}{number of missclassified cases, should only be used for binary responses}
  \item{tol_Xi}{minimal value for Norm2(Xi) and \eqn{\mathrm{det}(pp' \times pp)}{det(pp'*pp)} if there is any missing value in the \code{dataX}. It defaults to \eqn{10^{-12}}{10^{-12}}}
}
\details{
There are several ways to deal with missing values that leads to different computations of leave one out cross validation criteria.
}
\value{
  \item{nr}{Description of 'comp1'}
  \item{nc}{Description of 'comp2'}
  \item{nt}{Description of 'comp2'}
  \item{ww}{Description of 'comp2'}
  \item{wwnorm}{Description of 'comp2'}
  \item{wwetoile}{Description of 'comp2'}
  \item{tt}{Description of 'comp2'}
  \item{pp}{Description of 'comp2'}
  \item{CoeffC}{Description of 'comp2'}
  \item{uscores}{Description of 'comp2'}
  \item{YChapeau}{Description of 'comp2'}
  \item{residYChapeau}{Description of 'comp2'}
  \item{RepY}{Description of 'comp2'}
  \item{na.miss.Y}{Description of 'comp2'}
  \item{YNA}{Description of 'comp2'}
  \item{residY}{Description of 'comp2'}
  \item{ExpliX}{Description of 'comp2'}
  \item{na.miss.X}{Description of 'comp2'}
  \item{XXNA}{Description of 'comp2'}
  \item{residXX}{Description of 'comp2'}
  \item{PredictY}{Description of 'comp2'}
  \item{press.ind}{Description of 'comp2'}
  \item{press.tot}{Description of 'comp2'}
  \item{familly}{Description of 'comp2'}
  \item{ttPredictY}{Description of 'comp2'}
  \item{typeVC}{Description of 'comp2'}
  \item{dataY}{Description of 'comp2'}
  \item{computed_nt}{Description of 'comp2'}
  \item{CoeffCFull}{Description of 'comp2'}
  \item{CoeffConstante}{Description of 'comp2'}
  \item{Std.Coeffs}{Description of 'comp2'}
  \item{press.ind2}{Description of 'comp2'}
  \item{RSSresidY}{Description of 'comp2'}
  \item{Coeffs}{Description of 'comp2'}
  \item{Yresidus}{Description of 'comp2'}
  \item{RSS}{Description of 'comp2'}
  \item{residusY}{Description of 'comp2'}
  \item{AIC.std}{AIC.std vs number of components (AIC computed for the standardized model}
  \item{AIC}{AIC vs nber of components}
  \item{optional}{If the response is assumed to be binary:\cr
        i.e. \code{MClassed=TRUE}.    
  \describe{
  \item{\code{MissClassed}}{Number of miss classed results}
  \item{\code{Probs}}{Probability predicted by the model. These atre not true probabilities since they may lay outside of [0,1]}
  \item{\code{Probs.trc}}{Probability predicted by the model constrined in [0,1]}
  }
  }
  \item{ttPredictFittedMissingY}{Description of 'comp2'}
  \item{optional}{If cross validation was requested:\cr
        i.e. \code{typeVC="standard"}, \code{typeVC="missingdata"} or \code{typeVC="adaptative"}.
  \describe{
  \item{\code{R2residY}}{Description of 'comp2'}
  \item{\code{R2}}{Description of 'comp2'}
  \item{\code{press.tot2}}{Description of 'comp2'}
  \item{\code{Q2}}{Description of 'comp2'}
  \item{\code{limQ2}}{Description of 'comp2'}
  \item{\code{Q2_2}}{Description of 'comp2'}
  \item{\code{Q2cum}}{Description of 'comp2'}
  \item{\code{Q2cum_2}}{Description of 'comp2'}
  \item{\code{CVinfos}}{Description of 'comp2'}
  }
  }
  \item{optional}{If no cross validation was requested :\cr
        i.e. \code{typeVC="none"}.
  \describe{
  \item{\code{InfCrit}}{Several information criteria of fit statistics}
  }
  }
  \item{Std.ValsPredictY}{Description of 'comp2'}
  \item{ValsPredictY}{Description of 'comp2'}
  \item{Std.XChapeau}{Description of 'comp2'}
  \item{XXwotNA}{Description of 'comp2'}
}
\references{
Nicolas Meyer, Myriam Maumy-Bertrand et \enc{Frédéric}{Fr\'ed\'eric} Bertrand (2010). Comparaison de la \enc{régression}{r\'egression} PLS et de la \enc{régression}{r\'egression} logistique PLS : application aux \enc{données}{donn\'ees} \enc{d'allélotypage}{d'all\'elotypage}. \emph{Journal de la Société Française de Statistique}, 151(2), pages 1-18.
\url{http://smf4.emath.fr/Publications/JSFdS/151_2/pdf/sfds_jsfds_151_2_1-18.pdf}
}
\author{\enc{Frédéric}{Fr\'ed\'eric} Bertrand\cr
\email{frederic.bertrand@math.unistra.fr}\cr
\url{http://www-irma.u-strasbg.fr/~fbertran/}
}
\note{Use \code{plsR} instead.}
\seealso{\code{\link{plsRglmmodel.default}}}
\examples{
data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
plsRmodel.default(yCornell,XCornell,10)$InfCrit
rm(list=c("XCornell","yCornell"))
}
\keyword{models}
\keyword{regression}
