bootpls <- function(object, typeboot="plsmodel", R=250, statistic=coefs.plsR, sim="ordinary", stype="i",...){
callplsR <- object$call
dataset <- cbind(y = eval(callplsR$dataY),eval(callplsR$dataX))
nt <- eval(callplsR$nt)
if(typeboot=="plsmodel"){
return(boot(data=dataset, statistic=if(!(sim=="permutation")){coefs.plsR} else {permcoefs.plsR}, sim=sim, stype=stype, R=R, nt=nt, ...))
}
if(typeboot=="fmodel_np"){
return(boot(data=dataset, statistic=if(!(sim=="permutation")){coefs.plsR} else {permcoefs.plsR}, sim=sim, stype=stype, R=R, nt=nt, ...))
}
if(typeboot=="fmodel_par"){
return(boot(data=dataset, statistic=if(!(sim=="permutation")){coefs.plsR} else {permcoefs.plsR}, sim=sim, stype=stype, R=R, nt=nt, ...))
}
}
