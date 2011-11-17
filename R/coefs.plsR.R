coefs.plsR <- function(dataset,ind,nt,modele){
tempcoefs <- try(PLS_lm_wvc(dataY =dataset[ind,1], dataX=dataset[ind,-1], nt=nt, keepstd.coeffs=TRUE)$std.coeffs,silent=TRUE)
    if (is.numeric(tempcoefs)) {
        return(tempcoefs)
    }
    else {
        return(as.matrix(as.numeric(rep(NA, ncol(dataset)))))
    }
}
