coefs.plsRglm <- function(dataset, ind, nt, modele, family=NULL) 
{
    tempcoefs <- try(PLS_glm_wvc(dataY = dataset[ind, 1], dataX = dataset[ind, 
        -1], nt = nt, modele = modele, family=family, keepstd.coeffs = TRUE)$std.coeffs, silent=TRUE)
    if (is.numeric(tempcoefs)) {
        return(tempcoefs)
    }
    else {
        return(as.matrix(as.numeric(ifelse(any(class(dataset[,1])=="factor"),rep(NA, ncol(dataset)+nlevels(dataset[,1])-1),rep(NA, ncol(dataset))))))
    }
}
