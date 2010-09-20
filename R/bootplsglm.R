bootplsglm <- function(object, typeboot="plsmodel", R=250, statistic=coefs.plsRglm, sim="ordinary", stype="i"){
callplsRglm <- object$call
dataset <- cbind(y = eval(callplsRglm$dataY),eval(callplsRglm$dataX))
nt <- eval(callplsRglm$nt)
if(!is.null(callplsRglm$modele)){modele <- eval(callplsRglm$modele)} else {modele <- "pls"}
if(typeboot=="plsmodel"){
return(boot(data=dataset, statistic=if(!(sim=="permutation")){coefs.plsRglm} else {permcoefs.plsRglm}, sim=sim, stype=stype, R=R, nt=nt, modele=modele))
}
if(typeboot=="fmodel_np"){
return(boot(data=dataset, statistic=if(!(sim=="permutation")){coefs.plsRglm} else {permcoefs.plsRglm}, sim=sim, stype=stype, R=R, nt=nt, modele=modele))
}
if(typeboot=="fmodel_par"){
return(boot(data=dataset, statistic=if(!(sim=="permutation")){coefs.plsRglm} else {permcoefs.plsRglm}, sim=sim, stype=stype, R=R, nt=nt, modele=modele))
}
}