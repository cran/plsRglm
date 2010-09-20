coefs.plsR <- function(dataset,ind,nt,modele){
tempcoefs <- PLS_lm_wvc(dataY =dataset[ind,1], dataX=dataset[ind,-1], nt=nt, keepstd.coeffs=TRUE)$std.coeffs
if(!is.null(tempcoefs)) {return(tempcoefs)} else {return(as.matrix(as.numeric(c(0,rep(NA,ncol(dataset))))))}
}
