data(bordeaux)
Xbordeaux<-bordeaux[,1:4]
ybordeaux<-factor(bordeaux$Quality,ordered=TRUE)
modpls <- PLS_glm(ybordeaux,Xbordeaux,10,modele="pls-glm-polr")
modpls$Coeffsmodel_vals
modpls$InfCrit
modpls <- plsRglm(ybordeaux,Xbordeaux,10,modele="pls-glm-polr")
modpls$Coeffsmodel_vals
modpls$InfCrit


XbordeauxNA<-Xbordeaux
XbordeauxNA[1,1] <- NA
modplsNA <- PLS_glm(ybordeaux,XbordeauxNA,10,modele="pls-glm-polr")
modplsNA$Coeffsmodel_vals
modplsNA$InfCrit
modplsNA <- plsRglm(ybordeaux,XbordeauxNA,10,modele="pls-glm-polr")
modplsNA$Coeffsmodel_vals
modplsNA$InfCrit


data(bordeaux)
Xbordeaux<-bordeaux[,1:4]
ybordeaux<-factor(bordeaux$Quality,ordered=TRUE)
dataset <- cbind(y=ybordeaux,Xbordeaux)
options(contrasts = c("contr.treatment", "contr.poly"))
library(boot)
bordeaux.boot3<- bootplsglm(plsRglm(ybordeaux,Xbordeaux,1,modele="pls-glm-polr"), sim="permutation", stype="i", R=250)
#bordeaux.boot3 <- boot(data=dataset, statistic=permcoefs.plsRglm, sim="permutation", stype="i", R=250, nt=1, modele="pls-glm-polr")
boxplots.bootpls(bordeaux.boot3)
boxplots.bootpls(bordeaux.boot3,ranget0=TRUE)

bordeaux.boot2<- bootplsglm(plsRglm(ybordeaux,Xbordeaux,1,modele="pls-glm-polr"), sim="permutation", stype="i", R=250, strata=unclass(ybordeaux))
bordeaux.boot2 <- boot(data=dataset, statistic=permcoefs.plsRglm, sim="permutation", stype="i", R=250, nt=1, modele="pls-glm-polr", strata=unclass(ybordeaux))
boxplots.bootpls(bordeaux.boot2)


bordeaux.boot4<- bootplsglm(plsRglm(ybordeaux,Xbordeaux,1,modele="pls-glm-polr"), sim="balanced", stype="i", R=250)
#bordeaux.boot4 <- boot(data=dataset, statistic=coefs.plsRglm, sim="balanced", stype="i", R=250, nt=1, modele="pls-glm-polr")
#tempboot4 <- bordeaux.boot4
#tempboot4$t0=bordeaux.boot4$t0
#tempboot4$t=bordeaux.boot4$t[!is.na(bordeaux.boot4$t[,1]),]
#tempboot4$R=sum(!is.na(bordeaux.boot4$t[,1]))
#tempboot4$call$R<-sum(!is.na(bordeaux.boot4$t[,1]))
#tempboot4
boxplots.bootpls(bordeaux.boot4)
boot.ci(bordeaux.boot4)
confints.bootpls(bordeaux.boot4)
plots.confints.bootpls(confints.bootpls(bordeaux.boot4))

bordeaux.boot5<- bootplsglm(plsRglm(ybordeaux,Xbordeaux,1,modele="pls-glm-polr"), sim="balanced", stype="i", R=250, strata=unclass(ybordeaux))
#bordeaux.boot5 <- boot(data=dataset, statistic=coefs.plsRglm, sim="balanced", stype="i", R=250, nt=1, modele="pls-glm-polr", strata=unclass(ybordeaux))
#tempboot5 <- bordeaux.boot5
#tempboot5$t0=bordeaux.boot5$t0
#tempboot5$t=bordeaux.boot5$t[!is.na(bordeaux.boot5$t[,1]),]
#tempboot5$R=sum(!is.na(bordeaux.boot5$t[,1]))
#tempboot5$call$R<-sum(!is.na(bordeaux.boot5$t[,1]))
#tempboot5
boxplots.bootpls(bordeaux.boot5)
boot.ci(bordeaux.boot5)
confints.bootpls(bordeaux.boot5)
plots.confints.bootpls(confints.bootpls(bordeaux.boot5))


#install.packages("chemometrics")
if("chemometrics" %in% installed.packages()){
library(chemometrics)
data(hyptis)
hyptis
yhyptis <- factor(hyptis$Group,ordered=TRUE)
Xhyptis <- as.data.frame(hyptis[,c(1:6)])
dataset <- cbind(y=yhyptis,Xhyptis)
options(contrasts = c("contr.treatment", "contr.poly"))
modpls2 <- plsRglm(yhyptis,Xhyptis,6,modele="pls-glm-polr")
modpls2$Coeffsmodel_vals
modpls2$InfCrit
modpls2$Coeffs
modpls2$std.coeffs

table(yhyptis,predict(modpls2$FinalModel,type="class"))


modpls3 <- PLS_glm(yhyptis[-c(1,2,3)],Xhyptis[-c(1,2,3),],3,modele="pls-glm-polr",dataPredictY=Xhyptis[c(1,2,3),])
bbb <- PLS_glm_kfoldcv(yhyptis,Xhyptis,nt=4,K=10,random=TRUE,modele="pls-glm-polr",keepcoeffs=T)
kfolds2CVinfos_glm(bbb,MClassed=TRUE)



hyptis.boot3<- bootplsglm(plsRglm(yhyptis,Xhyptis,3,modele="pls-glm-polr"), sim="permutation", stype="i", R=250)
#library(boot)
#hyptis.boot3 <- boot(data=dataset, statistic=permcoefs.plsRglm, sim="permutation", stype="i", R=250, nt=3, modele="pls-glm-polr")
rownames(hyptis.boot3$t0)<-c("1|2\n","2|3\n","3|4\n","Sabi\nnene","Pin\nene","Cine\nole","Terpi\nnene","Fenc\nhone","Terpi\nnolene")
boxplots.bootpls(hyptis.boot3)
boxplots.bootpls(hyptis.boot3,xaxisticks=FALSE)
boxplots.bootpls(hyptis.boot3,ranget0=TRUE)
boxplots.bootpls(hyptis.boot3,ranget0=TRUE,xaxisticks=FALSE)

hyptis.boot4<- bootplsglm(plsRglm(yhyptis,Xhyptis,3,modele="pls-glm-polr"), sim="balanced", stype="i", R=250)
#hyptis.boot4 <- boot(data=dataset, statistic=coefs.plsRglm, sim="balanced", stype="i", R=500, nt=3, modele="pls-glm-polr")
#tempboot4 <- hyptis.boot4
#tempboot4$t0=hyptis.boot4$t0
#tempboot4$t=hyptis.boot4$t[!is.na(hyptis.boot4$t[,1]),]
#tempboot4$R=sum(!is.na(hyptis.boot4$t[,1]))
#tempboot4$call$R<-sum(!is.na(hyptis.boot4$t[,1]))
#tempboot4
rownames(hyptis.boot4$t0)<-c("1|2\n","2|3\n","3|4\n","Sabi\nnene","Pin\nene","Cine\nole","Terpi\nnene","Fenc\nhone","Terpi\nnolene")
boxplots.bootpls(hyptis.boot4)
boxplots.bootpls(hyptis.boot4,xaxisticks=FALSE)
confints.bootpls(hyptis.boot4)
plots.confints.bootpls(confints.bootpls(hyptis.boot4),legendpos = "bottomleft",xaxisticks=TRUE)
plots.confints.bootpls(confints.bootpls(hyptis.boot4),legendpos = "bottomleft",xaxisticks=FALSE)
