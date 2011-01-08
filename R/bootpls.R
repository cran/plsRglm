bootpls <- function(object, typeboot="plsmodel", R=250, statistic=coefs.plsR, sim="ordinary", stype="i",...){
callplsR <- object$call
#dataset <- cbind(y = eval(callplsR$dataY),eval(callplsR$dataX))
dataset <- cbind(y = object$dataY,object$dataX)
nt <- eval(callplsR$nt)
if(typeboot=="plsmodel"){
#return(boot(data=dataset, statistic=if(!(sim=="permutation")){coefs.plsR} else {permcoefs.plsR}, sim=sim, stype=stype, R=R, nt=nt, ...))
temp.bootplsR <- if(!(sim=="permutation")){boot(data=dataset, statistic=coefs.plsR, sim=sim, stype=stype, R=R, nt=nt, ...)} else {
boot(data=dataset, statistic=permcoefs.plsR, sim=sim, stype=stype, R=R, nt=nt)}
indices.temp.bootplsR <- !is.na(temp.bootplsR$t[,1])
temp.bootplsR$t=temp.bootplsR$t[indices.temp.bootplsR,]
temp.bootplsR$R=sum(indices.temp.bootplsR)
temp.bootplsR$call$R<-sum(indices.temp.bootplsR)
return(temp.bootplsR)
}
if(typeboot=="fmodel_np"){
#return(boot(data=dataset, statistic=if(!(sim=="permutation")){coefs.plsR} else {permcoefs.plsR}, sim=sim, stype=stype, R=R, nt=nt, ...))
temp.bootplsR <- if(!(sim=="permutation")){boot(data=dataset, statistic=coefs.plsR, sim=sim, stype=stype, R=R, nt=nt, ...)} else {
boot(data=dataset, statistic=permcoefs.plsR, sim=sim, stype=stype, R=R, nt=nt)}
indices.temp.bootplsR <- !is.na(temp.bootplsR$t[,1])
temp.bootplsR$t=temp.bootplsR$t[indices.temp.bootplsR,]
temp.bootplsR$R=sum(indices.temp.bootplsR)
temp.bootplsR$call$R<-sum(indices.temp.bootplsR)
return(temp.bootplsR)
}
if(typeboot=="fmodel_par"){
#return(boot(data=dataset, statistic=if(!(sim=="permutation")){coefs.plsR} else {permcoefs.plsR}, sim=sim, stype=stype, R=R, nt=nt, ...))
temp.bootplsR <- if(!(sim=="permutation")){boot(data=dataset, statistic=coefs.plsR, sim=sim, stype=stype, R=R, nt=nt, ...)} else {
boot(data=dataset, statistic=permcoefs.plsR, sim=sim, stype=stype, R=R, nt=nt)}
indices.temp.bootplsR <- !is.na(temp.bootplsR$t[,1])
temp.bootplsR$t=temp.bootplsR$t[indices.temp.bootplsR,]
temp.bootplsR$R=sum(indices.temp.bootplsR)
temp.bootplsR$call$R<-sum(indices.temp.bootplsR)
return(temp.bootplsR)
}
}
