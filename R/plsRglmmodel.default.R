plsRglmmodel.default <- function(dataY,dataX,nt=2,limQ2set=.0975,dataPredictY=dataX,modele="pls",family=NULL,typeVC="none",EstimXNA=FALSE,scaleX=TRUE,scaleY=NULL,pvals.expli=FALSE,alpha.pvals.expli=.05,MClassed=FALSE,tol_Xi=10^(-12))
{
  estmodel <- PLS_glm(dataY=dataY,dataX=dataX,nt=nt,limQ2set=limQ2set,dataPredictY=dataPredictY,modele=modele,family=family,typeVC=typeVC,EstimXNA=EstimXNA,scaleX=scaleX,scaleY=scaleY,pvals.expli=pvals.expli,alpha.pvals.expli=alpha.pvals.expli,MClassed=MClassed,tol_Xi=tol_Xi)
  if (typeVC!="none") {print("Use plsRglm_kfoldcv for applying kfold cross validation to glms");break}
  class(estmodel) <- "plsRglmmodel"
  estmodel$call <- match.call()
  estmodel
}
