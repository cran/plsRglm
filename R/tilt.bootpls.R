tilt.bootpls <- function(object, typeboot="plsmodel", statistic=coefs.plsR, R=c(499, 250, 250), alpha=c(0.025, 0.975), sim="ordinary", stype="i", index=1){
callplsR <- object$call
dataset <- cbind(y = eval(callplsR$dataY),eval(callplsR$dataX))
nt <- eval(callplsR$nt)
if(typeboot=="plsmodel"){
return(tilt.boot(data=dataset, statistic=if(!(sim=="permutation")){coefs.plsR} else {permcoefs.plsR}, sim=sim, stype=stype, R=R, nt=nt))
}
if(typeboot=="fmodel_np"){
return(tilt.boot(data=dataset, statistic=if(!(sim=="permutation")){coefs.plsRnp} else {permcoefs.plsRnp}, sim=sim, stype=stype, R=R, nt=nt))
}
if(typeboot=="fmodel_par"){
return(tilt.boot(data=dataset, statistic=if(!(sim=="permutation")){coefs.plsR} else {permcoefs.plsR}, sim=sim, stype=stype, R=R, nt=nt))
}
}
