tilt.bootplsglm <- function(object, typeboot="plsmodel", statistic=coefs.plsRglm, R=c(499, 250, 250), alpha=c(0.025, 0.975), sim="ordinary", stype="i", index=1){
callplsRglm <- object$call
dataset <- cbind(y = eval(callplsRglm$dataY),eval(callplsRglm$dataX))
nt <- eval(callplsRglm$nt)
if(!is.null(callplsRglm$modele)){modele <- eval(callplsRglm$modele)} else {modele <- "pls"}
if(typeboot=="plsmodel"){
return(tilt.boot(data=dataset, statistic=if(!(sim=="permutation")){coefs.plsRglm} else {permcoefs.plsRglm}, sim=sim, stype=stype, R=R, nt=nt, modele=modele))
}
if(typeboot=="fmodel_np"){
return(tilt.boot(data=dataset, statistic=if(!(sim=="permutation")){coefs.plsRglm} else {permcoefs.plsRglm}, sim=sim, stype=stype, R=R, nt=nt, modele=modele))
}
if(typeboot=="fmodel_par"){
return(tilt.boot(data=dataset, statistic=if(!(sim=="permutation")){coefs.plsRglm} else {permcoefs.plsRglm}, sim=sim, stype=stype, R=R, nt=nt, modele=modele))
}
}
