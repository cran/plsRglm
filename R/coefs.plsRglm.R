coefs.plsRglm <- function(dataset,ind,nt,modele){
tempcoefs <- PLS_glm_wvc(dataY =dataset[ind,1], dataX=dataset[ind,-1], nt=nt, modele=modele, keepstd.coeffs=TRUE)$std.coeffs
if(!is.null(tempcoefs)) {return(tempcoefs)} else {return(as.matrix(as.numeric(c(0,rep(NA,ncol(dataset))))))}
}
