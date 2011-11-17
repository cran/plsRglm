permcoefs.plsRglm <- function(dataset,ind,nt,modele,family=NULL,method="logistic"){
PLS_glm_wvc(dataY =dataset[,1], dataX=dataset[ind,-1], nt=nt, modele=modele, family=family, keepstd.coeffs=TRUE, method=method)$std.coeffs
}
