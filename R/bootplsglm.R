bootplsglm <- function(object, typeboot="plsmodel", R=250, statistic=coefs.plsRglm, sim="ordinary", stype="i",...){
callplsRglm <- object$call
#dataset <- cbind(y = eval(callplsRglm$dataY),eval(callplsRglm$dataX))
dataset <- cbind(y = object$dataY,object$dataX)
nt <- eval(callplsRglm$nt)
if(!is.null(callplsRglm$modele)){modele <- eval(callplsRglm$modele)} else {modele <- "pls"}
if(!is.null(callplsRglm$family)){family <- eval(callplsRglm$family)} else {family <- NULL}
if(typeboot=="plsmodel"){
temp.bootplsRglm <- if(!(sim=="permutation")){boot(data=dataset, statistic=coefs.plsRglm, sim=sim, stype=stype, R=R, nt=nt, modele=modele, family=family,...)} else {
boot(data=dataset, statistic=permcoefs.plsRglm, sim=sim, stype=stype, R=R, nt=nt, modele=modele, family=family)}
indices.temp.bootplsRglm <- !is.na(temp.bootplsRglm$t[,1])
temp.bootplsRglm$t=temp.bootplsRglm$t[indices.temp.bootplsRglm,]
temp.bootplsRglm$R=sum(indices.temp.bootplsRglm)
temp.bootplsRglm$call$R<-sum(indices.temp.bootplsRglm)
return(temp.bootplsRglm)
}
if(typeboot=="fmodel_np"){
temp.bootplsRglm <- if(!(sim=="permutation")){boot(data=dataset, statistic=coefs.plsRglm, sim=sim, stype=stype, R=R, nt=nt, modele=modele, family=family,...)} else {
boot(data=dataset, statistic=permcoefs.plsRglm, sim=sim, stype=stype, R=R, nt=nt, modele=modele, family=family)}
indices.temp.bootplsRglm <- !is.na(temp.bootplsRglm$t[,1])
temp.bootplsRglm$t=temp.bootplsRglm$t[indices.temp.bootplsRglm,]
temp.bootplsRglm$R=sum(indices.temp.bootplsRglm)
temp.bootplsRglm$call$R<-sum(indices.temp.bootplsRglm)
return(temp.bootplsRglm)
}
if(typeboot=="fmodel_par"){
temp.bootplsRglm <- if(!(sim=="permutation")){boot(data=dataset, statistic=coefs.plsRglm, sim=sim, stype=stype, R=R, nt=nt, modele=modele, family=family,...)} else {
boot(data=dataset, statistic=permcoefs.plsRglm, sim=sim, stype=stype, R=R, nt=nt, modele=modele, family=family)}
indices.temp.bootplsRglm <- !is.na(temp.bootplsRglm$t[,1])
temp.bootplsRglm$t=temp.bootplsRglm$t[indices.temp.bootplsRglm,]
temp.bootplsRglm$R=sum(indices.temp.bootplsRglm)
temp.bootplsRglm$call$R<-sum(indices.temp.bootplsRglm)
return(temp.bootplsRglm)
}
}