#' @rdname cv.plsR
#' @export
cv.plsRmodel.formula <- function(formula,data=NULL,nt=2,limQ2set=.0975,modele="pls", K=5, NK=1, grouplist=NULL, random=TRUE, scaleX=TRUE, scaleY=NULL, keepcoeffs=FALSE, keepfolds=FALSE, keepdataY=TRUE, keepMclassed=FALSE, tol_Xi=10^(-12), weights,subset,contrasts=NULL,verbose=TRUE) 
{
if (!(modele %in% c("pls"))) {stop("Use cv.plsRglm to cross-validate PLSRGLRs")}
if (missing(data)) {data <- environment(formula)}
mf <- match.call(expand.dots = FALSE)
m <- match(c("formula","data","nt","limQ2set","modele","K","NK","grouplist","random","scaleX","scaleY","keepcoeffs","keepfolds","keepdataY","keepMclassed","tol_Xi","weights","subset","contrasts","verbose"), names(mf), 0L)
mf <- mf[c(1L, m)]
mf[[1L]] <- as.name("PLS_lm_kfoldcv_formula")
cvmodel <- eval(mf, parent.frame())
  class(cvmodel) <- "cv.plsRmodel"
  cvmodel$call <- match.call()
  return(cvmodel)
}
