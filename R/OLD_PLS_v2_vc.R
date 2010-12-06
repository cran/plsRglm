OLD_PLS_v2_vc <- function(dataY,dataX,nt=2,limQ2set=.0975,dataPredictY=dataX,modele="pls",
                    family=NULL,typeVC="none",EstimXNA=FALSE,scaleX=TRUE,scaleY=NULL,pvals.expli=FALSE,alpha.pvals.expli=.05) {

#
# dataY : response(s) (training) dataset
# dataX : predictor(s) (training) dataset
# nt : number of components to be extracted
# limQ2set : limit value for the Q2
# dataPredictY : predictor(s) (testing) dataset
# modele : type of PLS model to be fitted ("pls", "pls-glm-gaussian", "pls-glm-logistic", "pls-glm-polr")
# family : for the moment the family argumlent is ignored and set thanks to the value of modele
# typeVC : type of crossed validation (only leave-one-out now). Several procedures are available and may be forced.
#    - "none" : no crossed validation
#    - "standard" : as in SIMCA for datasets without missing values and with all values predicted as those with missing values for datasets with any missing values
#    - "missingdata" : all values predicted as those with missing values for datasets with any missing values
#    - "adaptative" : predict a response value for an x with any missing value as those with missing values and for an x without any missing value as those without missing values.
# EstimaXNA : only for modele="pls". Set whether the missing X values have to be estimated.
# scaleX : scale the predictor(s) : must be set to TRUE for modele="pls" and should be for glms pls.
# scaleY : scale the response : Yes/No. Ignored since non always possible for glm responses.
# pvals.expli : should individual p-values be reported to enable model selection as in Tenenhaus
# alpha.pvals.expli : level of significance for predictors when pvals.expli=TRUE
#
# Cross validation not the same way for PLS than for PLS-glm : pls VC only on last components effect as in SIMCA. PLS Glm -> full cross validation.





##################################################
#                                                #
#    Initialization and formatting the inputs    #
#                                                #
##################################################

cat("____************************************************____\n")
if (!is.data.frame(dataX)) {dataX <- data.frame(dataX)}
if (!(modele %in% c("pls","pls-glm-logistic","pls-glm-gaussian","pls-glm-polr"))) {break}
if (modele=="pls-glm-logistic") {family<-binomial()}
if (modele=="pls-glm-gaussian") {family<-gaussian()}
scaleY <- NULL
if (is.null(scaleY)) {
if (!(modele %in% c("pls"))) {scaleY <- FALSE} else {scaleY <- TRUE}
}
if (scaleY) {RepY <- scale(dataY)}
else {
    RepY <- dataY
    attr(RepY,"scaled:center") <- 0
    attr(RepY,"scaled:scale") <- 1
}
if (scaleX) {ExpliX <- scale(dataX)
    PredictY <- sweep(sweep(dataPredictY, 2, attr(ExpliX,"scaled:center")), 2 ,attr(ExpliX,"scaled:scale"), "/")
}
else {
    ExpliX <- dataX
    attr(ExpliX,"scaled:center") <- rep(0,ncol(dataX))
    attr(ExpliX,"scaled:scale") <- rep(1,ncol(dataX))
    PredictY <- (dataPredictY)
}

if (any(is.na(dataX))) {na.miss.X <- TRUE} else na.miss.X <- FALSE
if (any(is.na(dataY))) {na.miss.Y <- TRUE} else na.miss.Y <- FALSE
if (any(is.na(PredictY))) {na.miss.PredictY <- TRUE} else na.miss.PredictY <- FALSE

XXNA <- !(is.na(ExpliX))
YNA <- !(is.na(RepY))
PredictYNA <- !is.na(PredictY)

ExpliXwotNA <- as.matrix(ExpliX)
ExpliXwotNA[!XXNA] <- 0

XXwotNA <- as.matrix(ExpliX)
XXwotNA[!XXNA] <- 0

dataXwotNA <- as.matrix(dataX)
dataXwotNA[!XXNA] <- 0

YwotNA <- as.matrix(RepY)
YwotNA[!YNA] <- 0

dataYwotNA <- as.matrix(dataY)

dataYwotNA[!YNA] <- 0

PredictYwotNA <- as.matrix(PredictY)
PredictYwotNA [is.na(PredictY)] <- 0

if (modele == "pls-glm-polr") {
dataY <- as.factor(dataY)
YwotNA <- as.factor(YwotNA)}

res <- list(nr=nrow(ExpliX),nc=ncol(ExpliX),ww=NULL,wwnorm=NULL,wwetoile=NULL,tt=NULL,pp=NULL,CoeffC=NULL,uscores=NULL,YChapeau=NULL,residYChapeau=NULL,RepY=RepY,na.miss.Y=na.miss.Y,YNA=YNA,residY=RepY,ExpliX=ExpliX,na.miss.X=na.miss.X,XXNA=XXNA,residXX=ExpliX,PredictY=PredictYwotNA,RSS=rep(NA,nt),RSSresidY=rep(NA,nt),R2=rep(NA,nt),R2residY=rep(NA,nt),press.ind=NULL,press.tot=NULL,Q2cum=rep(NA, nt),family=family,ttPredictY = NULL,typeVC=typeVC) 
res$temppred <- NULL

##############################################
######                PLS               ######
##############################################
if (modele == "pls") {
if (scaleY) {res$YChapeau=rep(attr(RepY,"scaled:center"),nrow(ExpliX))
res$residYChapeau=rep(0,nrow(ExpliX))}
else
{res$YChapeau=rep(mean(RepY),nrow(ExpliX))
res$residYChapeau=rep(mean(RepY),nrow(ExpliX))}
}





################################################
################################################
##                                            ##
##  Beginning of the loop for the components  ##
##                                            ##
################################################
################################################


for (kk in 1:nt) {
XXwotNA <- as.matrix(res$residXX)
XXwotNA[!XXNA] <- 0
YwotNA <- as.matrix(res$residY)
YwotNA[!YNA] <- 0
tempww <- rep(0,res$nc)




##############################################
#                                            #
#     Weight computation for each model      #
#                                            #
##############################################

##############################################
######                PLS               ######
##############################################
if (modele == "pls") {
tempww <- t(XXwotNA)%*%YwotNA/(t(XXNA)%*%YwotNA^2)
}

##############################################
######              PLS-GLM             ######
##############################################
if (modele %in% c("pls-glm-gaussian","pls-glm-logistic")) {
if (!pvals.expli) {
XXwotNA[!XXNA] <- NA
for (jj in 1:(res$nc)) {
    tempww[jj] <- coef(glm(YwotNA~cbind(res$tt,XXwotNA[,jj]),family=family))[kk+1]
}
XXwotNA[!XXNA] <- 0
rm(jj)}
else {
XXwotNA[!XXNA] <- NA
tempvalpvalstep <- rep(0,res$nc)
temppvalstep <- rep(0,res$nc)
for (jj in 1:(res$nc)) {
    tmww <- summary(glm(YwotNA~cbind(res$tt,XXwotNA[,jj]),family=family))$coefficients[kk+1,]
    tempww[jj] <- tmww[1]
    tempvalpvalstep[jj] <- tmww[4]
    temppvalstep[jj] <- (tmww[4] < alpha.pvals.expli)
}
XXwotNA[!XXNA] <- 0
rm(jj)
res$valpvalstep <- cbind(res$valpvalstep,tempvalpvalstep)
res$pvalstep <- cbind(res$pvalstep,temppvalstep)
}
}

##############################################
######           PLS-GLM-POLR           ######
##############################################
if (modele %in% c("pls-glm-polr")) {
YwotNA <- as.factor(YwotNA)
XXwotNA[!XXNA] <- NA
library(MASS)
for (jj in 1:(res$nc)) {
    tempww[jj] <- -1*polr(YwotNA~cbind(res$tt,XXwotNA[,jj]),na.action=na.exclude)$coef[kk]
}
XXwotNA[!XXNA] <- 0
rm(jj)}




##############################################
#                                            #
# Computation of the components (model free) #
#                                            #
##############################################

res$ww <- cbind(res$ww,tempww)

tempwwnorm <- tempww/sqrt(drop(crossprod(tempww)))
res$wwnorm <- cbind(res$wwnorm,tempwwnorm)

temptt <- XXwotNA%*%tempwwnorm/(XXNA%*%(tempwwnorm^2))
res$tt <- cbind(res$tt,temptt)       #tt = scores

temppp <- rep(0,res$nc)
for (jj in 1:(res$nc)) {
     temppp[jj] <- crossprod(temptt,XXwotNA[,jj])/drop(crossprod(XXNA[,jj],temptt^2))
}
res$residXX <- XXwotNA-temptt%*%temppp
res$pp <- cbind(res$pp,temppp)   #Contrib of X




##############################################
#                                            #
#      Computation of the coefficients       #
#      of the model with kk components       #
#                                            #
##############################################

##############################################
######                PLS               ######
##############################################
if (modele == "pls") {
if (kk==1) {
tempCoeffC <- solve(t(res$tt[YNA])%*%res$tt[YNA])%*%t(res$tt[YNA])%*%YwotNA[YNA]
res$CoeffCFull <- matrix(c(tempCoeffC,rep(NA,nt-kk)),ncol=1)
tempCoeffConstante <- 0
} else {
if (!(na.miss.X | na.miss.Y)) {# This disjonction is unnecessary
tempCoeffC <- c(rep(0,kk-1),solve(t(res$tt[YNA,kk])%*%res$tt[YNA,kk])%*%t(res$tt[YNA,kk])%*%YwotNA[YNA])
#tempCoeffC <- solve(t(res$tt[YNA,])%*%res$tt[YNA,])%*%t(res$tt[YNA,])%*%YwotNA[YNA]  
tempCoeffConstante <- 0
res$CoeffCFull <- cbind(res$CoeffCFull,c(tempCoeffC,rep(NA,nt-kk)))
}
else
{
tempCoeffC <- c(rep(0,kk-1),solve(t(res$tt[YNA,kk])%*%res$tt[YNA,kk])%*%t(res$tt[YNA,kk])%*%YwotNA[YNA])  
#tempCoeffC <- solve(t(res$tt[YNA,])%*%res$tt[YNA,])%*%t(res$tt[YNA,])%*%YwotNA[YNA]  
tempCoeffConstante <- 0
res$CoeffCFull <- cbind(res$CoeffCFull,c(tempCoeffC,rep(NA,nt-kk)))
}
}

res$wwetoile <- (res$wwnorm)%*%solve(t(res$pp)%*%res$wwnorm)
res$CoeffC <- diag(res$CoeffCFull)
res$CoeffConstante <- tempCoeffConstante
res$Std.Coeffs <- rbind(tempCoeffConstante,res$wwetoile%*%res$CoeffC)
rownames(res$Std.Coeffs) <- c("Intercept",colnames(ExpliX))
}


##############################################
######              PLS-GLM             ######
##############################################
if (modele %in% c("pls-glm-logistic","pls-glm-gaussian")) {
if (kk==1) {
tempconstglm <- glm(YwotNA~1,family=family)
res$AIC <- AIC(tempconstglm)
res$BIC <- AIC(tempconstglm, k = log(res$nr))
res$Coeffsmodel_vals <- rbind(summary(tempconstglm)$coefficients,matrix(rep(NA,4*nt),ncol=4))
res$ChisqPearson <- crossprod(residuals.glm(tempconstglm,type="pearson"))  
rm(tempconstglm)
tempregglm <- glm(YwotNA~res$tt,family=family)
res$AIC <- cbind(res$AIC,AIC(tempregglm))
res$BIC <- cbind(res$BIC,AIC(tempregglm, k = log(res$nr)))
res$Coeffsmodel_vals <- cbind(res$Coeffsmodel_vals,rbind(summary(tempregglm)$coefficients,matrix(rep(NA,4*(nt-kk)),ncol=4)))
res$ChisqPearson <- c(res$ChisqPearson,crossprod(residuals.glm(tempregglm,type="pearson")))
tempCoeffC <- as.vector(coef(tempregglm))
res$CoeffCFull <- matrix(c(tempCoeffC,rep(NA,nt-kk)),ncol=1)
tempCoeffConstante <- tempCoeffC[1]
res$CoeffConstante <- tempCoeffConstante
tempCoeffC <- tempCoeffC[-1]
} else {
if (!(na.miss.X | na.miss.Y)) {
tempregglm <- glm(YwotNA~res$tt,family=family)
res$AIC <- cbind(res$AIC,AIC(tempregglm))
res$BIC <- cbind(res$BIC,AIC(tempregglm, k = log(res$nr)))
res$Coeffsmodel_vals <- cbind(res$Coeffsmodel_vals,rbind(summary(tempregglm)$coefficients,matrix(rep(NA,4*(nt-kk)),ncol=4)))
res$ChisqPearson <- c(res$ChisqPearson,crossprod(residuals.glm(tempregglm,type="pearson")))
tempCoeffC <- as.vector(coef(tempregglm))  
res$CoeffCFull <- cbind(res$CoeffCFull,c(tempCoeffC,rep(NA,nt-kk)))
tempCoeffConstante <- tempCoeffC[1]
res$CoeffConstante <- cbind(res$CoeffConstante,tempCoeffConstante)
tempCoeffC <- tempCoeffC[-1]
}
else
{
tempregglm <- glm(YwotNA~res$tt,family=family)
res$AIC <- cbind(res$AIC,AIC(tempregglm))
res$BIC <- cbind(res$BIC,AIC(tempregglm, k = log(res$nr)))
res$Coeffsmodel_vals <- cbind(res$Coeffsmodel_vals,rbind(summary(tempregglm)$coefficients,matrix(rep(NA,4*(nt-kk)),ncol=4)))
res$ChisqPearson <- c(res$ChisqPearson,crossprod(residuals.glm(tempregglm,type="pearson")))
tempCoeffC <- as.vector(coef(tempregglm))  
res$CoeffCFull <- cbind(res$CoeffCFull,c(tempCoeffC,rep(NA,nt-kk)))
tempCoeffConstante <- tempCoeffC[1]
res$CoeffConstante <- cbind(res$CoeffConstante,tempCoeffConstante)
tempCoeffC <- tempCoeffC[-1]
}
}

res$wwetoile <- (res$wwnorm)%*%solve(t(res$pp)%*%res$wwnorm)
res$CoeffC <- tempCoeffC
res$Std.Coeffs <- rbind(tempCoeffConstante,res$wwetoile%*%res$CoeffC)
rownames(res$Std.Coeffs) <- c("Intercept",colnames(ExpliX))
}


##############################################
######           PLS-GLM-POLR           ######
##############################################

if (modele %in% c("pls-glm-polr")) {
if (kk==1) {
tempconstpolr <- polr(YwotNA~1,na.action=na.exclude,Hess=TRUE)
res$AIC <- AIC(tempconstpolr)
res$BIC <- AIC(tempconstpolr, k = log(res$nr))
res$Coeffsmodel_vals <- rbind(summary(tempconstpolr)$coefficients,matrix(rep(NA,3*nt),ncol=3))
rm(tempconstpolr)
tempregpolr <- polr(YwotNA~res$tt,na.action=na.exclude,Hess=TRUE)
res$AIC <- cbind(res$AIC,AIC(tempregpolr))
res$BIC <- cbind(res$BIC,AIC(tempregpolr, k = log(res$nr)))
res$Coeffsmodel_vals <- cbind(res$Coeffsmodel_vals,rbind(summary(tempregpolr)$coefficients,matrix(rep(NA,3*(nt-kk)),ncol=3)))
tempCoeffC <- -1*as.vector(tempregpolr$coef)
tempCoeffConstante <- as.vector(tempregpolr$zeta)
res$CoeffCFull <- matrix(c(tempCoeffConstante,tempCoeffC,rep(NA,nt-kk)),ncol=1)
res$CoeffConstante <- tempCoeffConstante
} else {
if (!(na.miss.X | na.miss.Y)) {
tempregpolr <- polr(YwotNA~res$tt,na.action=na.exclude,Hess=TRUE)
res$AIC <- cbind(res$AIC,AIC(tempregpolr))
res$BIC <- cbind(res$BIC,AIC(tempregpolr, k = log(res$nr)))
res$Coeffsmodel_vals <- cbind(res$Coeffsmodel_vals,rbind(summary(tempregpolr)$coefficients,matrix(rep(NA,3*(nt-kk)),ncol=3)))
tempCoeffC <- -1*as.vector(tempregpolr$coef)  
tempCoeffConstante <- as.vector(tempregpolr$zeta)
res$CoeffCFull <- cbind(res$CoeffCFull,c(tempCoeffConstante,tempCoeffC,rep(NA,nt-kk)))
res$CoeffConstante <- cbind(res$CoeffConstante,tempCoeffConstante)
}
else
{
tempregpolr <- polr(YwotNA~res$tt,na.action=na.exclude,Hess=TRUE)
res$AIC <- cbind(res$AIC,AIC(tempregpolr))
res$BIC <- cbind(res$BIC,AIC(tempregpolr, k = log(res$nr)))
res$Coeffsmodel_vals <- cbind(res$Coeffsmodel_vals,rbind(summary(tempregpolr)$coefficients,matrix(rep(NA,3*(nt-kk)),ncol=3)))
tempCoeffC <- -1*as.vector(tempregpolr$coef)  
tempCoeffConstante <- as.vector(tempregpolr$zeta)
res$CoeffCFull <- cbind(res$CoeffCFull,c(tempCoeffConstante,tempCoeffC,rep(NA,nt-kk)))
res$CoeffConstante <- cbind(res$CoeffConstante,tempCoeffConstante)
}
}

res$wwetoile <- (res$wwnorm)%*%solve(t(res$pp)%*%res$wwnorm)
res$CoeffC <- tempCoeffC
res$Std.Coeffs <- as.matrix(rbind(as.matrix(tempCoeffConstante),res$wwetoile%*%res$CoeffC))
rownames(res$Std.Coeffs) <- c(names(tempregpolr$zeta),colnames(ExpliX))
}




##############################################
#                                            #
#       Prediction of the components         #
#     as if missing values (model free)      #
#       For cross-validating the GLM         #
#                                            #
##############################################

res$ttPredictFittedMissingY <- NULL
for (ii in 1:res$nr) {
res$ttPredictFittedMissingY <- rbind(res$ttPredictFittedMissingY,t(solve(t(res$pp[XXNA[ii,],])%*%res$pp[XXNA[ii,],])%*%t(res$pp[XXNA[ii,],])%*%(ExpliXwotNA[ii,])[XXNA[ii,]]))
}
rm(ii)




if (!(na.miss.X | na.miss.Y)) {

##############################################
#                                            #
#             Cross validation               #
#           without missing value            #
#                                            #
##############################################

##############################################
######                PLS               ######
##############################################
if (modele == "pls") {
if (typeVC == "none") {} else {
if (typeVC %in% c("standard","adaptative")) {
if (kk==1) {
cat(paste("____TypeVC____",typeVC,"____\n"))
}
temppred <- rep(0, res$nr)
#temppred2 <- rep(0, res$nr)
for (i in 1:(res$nr)) {     
                tempww.cv <- t(XXwotNA[-i, ])%*%YwotNA[-i]/(t(XXNA[-i, ])%*%YwotNA[-i]^2)
                tempwwnorm.cv <- tempww.cv/sqrt(drop(crossprod(tempww.cv)))
                temptt.cv <- XXwotNA[-i, ]%*%tempwwnorm.cv/(XXNA[-i, ]%*%(tempwwnorm.cv^2))
                tempc.cv <- solve(t(temptt.cv)%*%temptt.cv)%*%t(temptt.cv)%*%(YwotNA[-i])
                tempc.cv <- as.vector(tempc.cv)
                temppred[i] <- tempc.cv * XXwotNA[i, ] %*% tempwwnorm.cv
#                temppred2[i] <- res$YChapeau[i]+attr(res$RepY,"scaled:scale")*temppred[i]  
}
rm(i)
res$press.ind <- cbind(res$press.ind,(YwotNA-temppred)^2) 
res$press.ind2 <- cbind(res$press.ind2,(dataYwotNA-res$YChapeau-attr(res$RepY,"scaled:scale")*temppred)^2)
}
else {
if (typeVC == "missingdata") {
if (kk==1) {
cat(paste("____TypeVC____",typeVC,"____\n"))
}
temppp.cv <- res$pp  
temppred <- rep(0, res$nr)
for (i in 1:(res$nr)) {     
                tempww.cv <- t(XXwotNA[-i, ])%*%YwotNA[-i]/(t(XXNA[-i, ])%*%YwotNA[-i]^2)
                tempwwnorm.cv <- tempww.cv/sqrt(drop(crossprod(tempww.cv)))
                temptt.cv <- XXwotNA[-i, ]%*%tempwwnorm.cv/(XXNA[-i, ]%*%(tempwwnorm.cv^2))
                tempc.cv <- solve(t(temptt.cv)%*%temptt.cv)%*%t(temptt.cv)%*%(YwotNA[-i])    
                for (jj in 1:(res$nc)) {
                     temppp.cv[jj,kk] <- crossprod(temptt.cv,(XXwotNA[-i,])[,jj])/drop(crossprod((XXNA[-i, ])[,jj],temptt.cv^2))
                }
                ttPredictY.cv <- (solve(t(temppp.cv[XXNA[i,],])%*%temppp.cv[XXNA[i,],])%*%t(temppp.cv[XXNA[i,],])%*%(ExpliXwotNA[i,])[XXNA[i,]])[kk]
                temppred[i] <- tempc.cv*ttPredictY.cv   
}
rm(i)
res$press.ind <- cbind(res$press.ind,(YwotNA-temppred)^2)
res$press.ind2 <- cbind(res$press.ind2,(dataYwotNA-res$YChapeau-attr(res$RepY,"scaled:scale")*temppred)^2)
}
else {
if (kk==1) {
cat(paste("____TypeVC____",typeVC,"____inexistant____\n"))
}
}
}
}
res$residYChapeau <- res$tt%*%tempCoeffC
res$RSSresidY[kk] <- crossprod(res$residY-res$residYChapeau)


tempCoeffs <- res$wwetoile%*%res$CoeffC*attr(res$RepY,"scaled:scale")/attr(res$ExpliX,"scaled:scale")
tempConstante <- attr(res$RepY,"scaled:center")-sum(tempCoeffs*attr(res$ExpliX,"scaled:center"))
res$Coeffs <- rbind(tempConstante,tempCoeffs)

res$YChapeau <- attr(res$RepY,"scaled:center")+attr(res$RepY,"scaled:scale")*res$tt%*%res$CoeffC             
res$Yresidus <- dataY-res$YChapeau
res$RSS[kk] <- crossprod(res$Yresidus)
}
##############################################


##############################################
######              PLS-GLM             ######
##############################################
if (modele %in% c("pls-glm-gaussian","pls-glm-logistic")) {
if (typeVC == "none") {} else {
if (typeVC %in% c("standard","adaptative")) {
if (kk==1) {
cat(paste("____TypeVC____",typeVC,"____\n"))
cat("Not available, please use PLS_v2_kfoldcv\n")
}
}
else {
if (typeVC == "missingdata") {
if (kk==1) {
cat(paste("____TypeVC____",typeVC,"____\n"))
}
}
else {
if (kk==1) {
cat(paste("____TypeVC____",typeVC,"____inexistant____\n"))
cat("Not available, please use PLS_v2_kfoldcv\n")
}
}
}
}
res$residYChapeau <- tempregglm$linear.predictors
res$RSSresidY[kk] <- crossprod(res$residY-res$residYChapeau)


tempCoeffs <- res$wwetoile%*%res$CoeffC*attr(res$RepY,"scaled:scale")/attr(res$ExpliX,"scaled:scale")
tempConstante <- attr(res$RepY,"scaled:center")-sum(tempCoeffs*attr(res$ExpliX,"scaled:center"))+attr(res$RepY,"scaled:scale")*res$Std.Coeffs[1]
res$Coeffs <- rbind(tempConstante,tempCoeffs)

res$YChapeau <- tempregglm$fitted.values          
res$Yresidus <- dataY-res$YChapeau
res$RSS[kk] <- crossprod(res$Yresidus)
}
if (modele %in% c("pls-glm-polr")) {
if (typeVC == "none") {} else {
if (typeVC %in% c("standard","adaptative")) {
if (kk==1) {
cat(paste("____TypeVC____",typeVC,"____\n"))
cat("Not available, please use PLS_v2_kfoldcv\n")
}
}
}
}
##############################################
}

else {
if (na.miss.X & !na.miss.Y) {


##############################################
#                                            #
#             Cross validation               #
#           with missing value(s)            #
#                                            #
##############################################


if (kk==1) {
cat("____There are some NAs in X but not in Y____\n")
}

##############################################
######                PLS               ######
##############################################
if (modele == "pls") {
if (typeVC == "none") {} else {
if (typeVC %in% c("standard","missingdata")) {
if (kk==1) {
cat(paste("____TypeVC____",typeVC,"____\n"))
}
temppp.cv <- res$pp  
temppred <- rep(0, res$nr)
for (i in 1:(res$nr)) {     
                tempww.cv <- t(XXwotNA[-i, ])%*%YwotNA[-i]/(t(XXNA[-i, ])%*%YwotNA[-i]^2)
                tempwwnorm.cv <- tempww.cv/sqrt(drop(crossprod(tempww.cv)))
                temptt.cv <- XXwotNA[-i, ]%*%tempwwnorm.cv/(XXNA[-i, ]%*%(tempwwnorm.cv^2))
                tempc.cv <- solve(t(temptt.cv)%*%temptt.cv)%*%t(temptt.cv)%*%(YwotNA[-i])    
                for (jj in 1:(res$nc)) {
                     temppp.cv[jj,kk] <- crossprod(temptt.cv,(XXwotNA[-i,])[,jj])/drop(crossprod((XXNA[-i, ])[,jj],temptt.cv^2))
                }
                ttPredictY.cv <- (solve(t(temppp.cv[XXNA[i,],])%*%temppp.cv[XXNA[i,],])%*%t(temppp.cv[XXNA[i,],])%*%(XXwotNA[i,])[XXNA[i,]])[kk]
                temppred[i] <- tempc.cv*ttPredictY.cv   
}
res$press.ind <- cbind(res$press.ind,(YwotNA-temppred)^2)
res$press.ind2 <- cbind(res$press.ind2,(dataYwotNA-res$YChapeau-attr(res$RepY,"scaled:scale")*temppred)^2)
}
else {
if (typeVC == "adaptative") {
if (kk==1) {
cat(paste("____TypeVC____",typeVC,"____\n"))
}
temppp.cv <- res$pp  
temppred <- rep(0, res$nr)
for (i in 1:(res$nr)) {     
if (all(XXNA[i,])) {
                tempww.cv <- t(XXwotNA[-i, ])%*%YwotNA[-i]/(t(XXNA[-i, ])%*%YwotNA[-i]^2)
                tempwwnorm.cv <- tempww.cv/sqrt(drop(crossprod(tempww.cv)))
                temptt.cv <- XXwotNA[-i, ]%*%tempwwnorm.cv/(XXNA[-i, ]%*%(tempwwnorm.cv^2))
                tempc.cv <- solve(t(temptt.cv)%*%temptt.cv)%*%t(temptt.cv)%*%(YwotNA[-i])
                tempc.cv <- as.vector(tempc.cv)
                temppred[i] <- tempc.cv * XXwotNA[i, ] %*% tempwwnorm.cv
}
else {
                tempww.cv <- t(XXwotNA[-i, ])%*%YwotNA[-i]/(t(XXNA[-i, ])%*%YwotNA[-i]^2)
                tempwwnorm.cv <- tempww.cv/sqrt(drop(crossprod(tempww.cv)))
                temptt.cv <- XXwotNA[-i, ]%*%tempwwnorm.cv/(XXNA[-i, ]%*%(tempwwnorm.cv^2))
                tempc.cv <- solve(t(temptt.cv)%*%temptt.cv)%*%t(temptt.cv)%*%(YwotNA[-i])    
                for (jj in 1:(res$nc)) {
                     temppp.cv[jj,kk] <- crossprod(temptt.cv,(XXwotNA[-i,])[,jj])/drop(crossprod((XXNA[-i, ])[,jj],temptt.cv^2))
                }
                ttPredictY.cv <- (solve(t(temppp.cv[XXNA[i,],])%*%temppp.cv[XXNA[i,],])%*%t(temppp.cv[XXNA[i,],])%*%(XXwotNA[i,])[XXNA[i,]])[kk]
                temppred[i] <- tempc.cv*ttPredictY.cv   
}
}
res$press.ind <- cbind(res$press.ind,(YwotNA-temppred)^2)
res$press.ind2 <- cbind(res$press.ind2,(dataYwotNA-res$YChapeau-attr(res$RepY,"scaled:scale")*temppred)^2)
}
else {
if (kk==1) {
cat(paste("____TypeVC____",typeVC,"____inexistant____\n"))
}
}
}
}
res$residYChapeau <- res$tt%*%tempCoeffC
res$RSSresidY[kk] <- crossprod(res$residY-res$residYChapeau)


tempCoeffs <- res$wwetoile%*%res$CoeffC*attr(res$RepY,"scaled:scale")/attr(res$ExpliX,"scaled:scale")
tempConstante <- attr(res$RepY,"scaled:center")-sum(tempCoeffs*attr(res$ExpliX,"scaled:center"))
res$Coeffs <- rbind(tempConstante,tempCoeffs)

res$YChapeau <- attr(res$RepY,"scaled:center")+attr(res$RepY,"scaled:scale")*res$tt%*%res$CoeffC            
res$Yresidus <- dataY-res$YChapeau
res$RSS[kk] <- crossprod(res$Yresidus)
}
##############################################



##############################################
######              PLS-GLM             ######
##############################################
if (modele %in% c("pls-glm-gaussian","pls-glm-logistic")) {
if (typeVC == "none") {} else {
if (typeVC %in% c("standard","missingdata")) {
if (kk==1) {
cat(paste("____TypeVC____",typeVC,"____inexistant____\n"))
cat("Not available, please use PLS_v2_kfoldcv\n")
}
}
else {
if (typeVC == "adaptative") {
if (kk==1) {
cat(paste("____TypeVC____",typeVC,"____inexistant____\n"))
cat("Not available, please use PLS_v2_kfoldcv\n")
}
}
else {
if (kk==1) {
cat(paste("____TypeVC____",typeVC,"____inexistant____\n"))
cat("Not available, please use PLS_v2_kfoldcv\n")
}
}
}
}
res$residYChapeau <- tempregglm$linear.predictors
res$RSSresidY[kk] <- crossprod(res$residY-res$residYChapeau)


tempCoeffs <- res$wwetoile%*%res$CoeffC*attr(res$RepY,"scaled:scale")/attr(res$ExpliX,"scaled:scale")
tempConstante <- attr(res$RepY,"scaled:center")-sum(tempCoeffs*attr(res$ExpliX,"scaled:center"))
res$Coeffs <- rbind(tempConstante,tempCoeffs)

res$YChapeau <- tempregglm$fitted.values                      
res$Yresidus <- dataY-res$YChapeau
res$RSS[kk] <- crossprod(res$Yresidus)
}
if (modele %in% c("pls-glm-polr")) {
if (typeVC == "none") {} else {
if (typeVC %in% c("standard","adaptative")) {
if (kk==1) {
cat(paste("____TypeVC____",typeVC,"____\n"))
cat("Not available, please use PLS_v2_kfoldcv\n")
}
}
}
}
##############################################
}

else {
if (kk==1) {
cat("____There are some NAs both in X and Y____\n")
}
}
}


##############################################
#                                            #
#      Update and end of loop cleaning       #
#        (Especially useful for PLS)         #
#                                            #
##############################################


##############################################
######                PLS               ######
##############################################
if (modele == "pls") {
res$uscores <- cbind(res$uscores,res$residY/res$CoeffC[kk])
res$residY <- res$residY - res$tt%*%tempCoeffC 
res$residusY <- cbind(res$residusY,res$residY)

rm(tempww)
rm(tempwwnorm)
rm(temptt)
rm(temppp)
rm(tempCoeffC)
rm(tempCoeffs)
rm(tempConstante)
}

##############################################
######              PLS-GLM             ######
##############################################
if (modele %in% c("pls-glm-gaussian","pls-glm-logistic")) {
res$residY <- res$residY 
res$residusY <- cbind(res$residusY,res$residY)

rm(tempww)
rm(tempwwnorm)
rm(temptt)
rm(temppp)
rm(tempCoeffC)
rm(tempCoeffs)
rm(tempConstante)
}

##############################################
######           PLS-GLM-POLR           ######
##############################################
if (modele %in% c("pls-glm-polr")) {
res$residY <- res$residY 
res$residusY <- cbind(res$residusY,res$residY)

rm(tempww)
rm(tempwwnorm)
rm(temptt)
rm(temppp)
rm(tempCoeffC)
}
cat("____Component____",kk,"____\n")
}




##############################################
##############################################
##                                          ##
##    End of the loop on the components     ##
##                                          ##
##############################################
##############################################




##############################################
#                                            #
#           Predicting components            #
#                                            #
##############################################

if (!(na.miss.X | na.miss.Y)) {
cat("____Predicting X without NA neither in X or Y____\n")
res$ttPredictY <- PredictYwotNA%*%res$wwetoile 
colnames(res$ttPredictY) <- paste("tt",1:nt,sep="")
}
else {
if (na.miss.X & !na.miss.Y) {
cat("____Predicting X with NA in X and not in Y____\n")
for (ii in 1:nrow(PredictYwotNA)) {  
      res$ttPredictY <- rbind(res$ttPredictY,t(solve(t(res$pp[PredictYNA[ii,],])%*%res$pp[PredictYNA[ii,],])%*%t(res$pp[PredictYNA[ii,],])%*%(PredictYwotNA[ii,])[PredictYNA[ii,]]))
}
colnames(res$ttPredictY) <- paste("tt",1:nt,sep="")
}
else {
cat("____There are some NAs both in X and Y____\n")
}
}




##############################################
#                                            #
#          Computing RSS, PRESS,             #
#           Chi2, Q2 and Q2cum               #
#                                            #
##############################################

##############################################
######                PLS               ######
##############################################
if (modele == "pls") {
res$RSSresidY<- c(crossprod(RepY),res$RSSresidY)
res$R2residY <- 1-res$RSSresidY[2:(nt+1)]/res$RSSresidY[1]
res$RSS <- c(crossprod(dataY-mean(dataY)),res$RSS)
res$R2 <- 1-res$RSS[2:(nt+1)]/res$RSS[1]
if (typeVC %in% c("standard","missingdata","adaptative")) {
res$press.tot <- colSums(res$press.ind)
res$press.tot2 <- colSums(res$press.ind2)
res$Q2 <- 1-res$press.tot/res$RSSresidY[-(nt+1)]
res$limQ2 <- rep(limQ2set,nt)
res$Q2_2 <- 1-res$press.tot2/res$RSS[-(nt+1)]


for (k in 1:nt) {res$Q2cum[k] <- prod(res$press.tot[1:k])/prod(res$RSSresidY[1:k])}
res$Q2cum <- 1 - res$Q2cum
for (k in 1:nt) {res$Q2cum_2[k] <- prod(res$press.tot2[1:k])/prod(res$RSS[1:k])}
res$Q2cum_2 <- 1 - res$Q2cum_2
res$CVinfos <- cbind(res$Q2cum_2, res$limQ2, res$Q2_2[1:nt], res$press.tot2[1:nt], res$RSS[1:nt], res$R2[1:nt], res$R2residY[1:nt], res$RSSresidY[1:nt], res$press.tot[1:nt], res$Q2[1:nt], res$limQ2, res$Q2cum)
dimnames(res$CVinfos) <- list(1:nt, c("Q2cum_Y", "LimQ2_Y", "Q2_Y", "PRESS_Y", "RSS_Y", "R2_Y", "R2_residY", "RSS_residY", "PRESS_residY", "Q2_residY", "LimQ2", "Q2cum_residY"))
}
res$InfCrit <- cbind(res$RSS, c(0,res$R2), c(0,res$R2residY), res$RSSresidY)
dimnames(res$InfCrit) <- list(paste("Nb_Comp_",0:nt), c("RSS_Y", "R2_Y", "R2_residY", "RSS_residY"))
}


##############################################
######              PLS-GLM             ######
##############################################
if (modele %in% c("pls-glm-gaussian", "pls-glm-logistic")) {
res$RSSresidY<- c(crossprod(RepY-mean(RepY)),res$RSSresidY)
res$R2residY <- 1-res$RSSresidY[2:(nt+1)]/res$RSSresidY[1]
res$RSS <- c(crossprod(dataY-mean(dataY)),res$RSS)
res$R2 <- 1-res$RSS[2:(nt+1)]/res$RSS[1]
if (typeVC %in% c("standard","missingdata","adaptative")) {
}
res$InfCrit <- cbind(t(res$AIC), t(res$BIC), res$ChisqPearson, res$RSS, c(0,res$R2), c(0,res$R2residY), res$RSSresidY)
dimnames(res$InfCrit) <- list(paste("Nb_Comp_",0:nt), c("AIC", "BIC", "Chi2_Pearson_Y", "RSS_Y", "R2_Y", "R2_residY", "RSS_residY"))
}


##############################################
######           PLS-GLM-POLR           ######
##############################################
if (modele == "pls-glm-polr") {
if (typeVC %in% c("standard","missingdata","adaptative")) {
}
res$InfCrit <- cbind(t(res$AIC), t(res$BIC))
dimnames(res$InfCrit) <- list(paste("Nb_Comp_",0:nt), c("AIC", "BIC"))
}




##########################################
#                                        #
#          Predicting responses          #
#                                        #
##########################################


##############################################
######               PLS                ######
##############################################
if (modele == "pls") {
res$YChapeau <- attr(res$RepY,"scaled:center")+attr(res$RepY,"scaled:scale")*res$tt%*%res$CoeffC            
rownames(res$YChapeau) <- rownames(ExpliX)

res$Std.ValsPredictY <- res$ttPredictY%*%res$CoeffC
res$ValsPredictY <- attr(res$RepY,"scaled:center")+attr(res$RepY,"scaled:scale")*res$ttPredictY%*%res$CoeffC

res$Std.XChapeau <- res$tt%*%t(res$pp)
rownames(res$Std.XChapeau) <- rownames(ExpliX)
if (EstimXNA) {
res$XChapeau <- sweep(sweep(res$Std.XChapeau,2,attr(res$ExpliX,"scaled:scale"),FUN="*"),2,attr(res$ExpliX,"scaled:center"),FUN="+")
rownames(res$XChapeau) <- rownames(ExpliX)
colnames(res$XChapeau) <- colnames(ExpliX)

res$XChapeauNA <- sweep(sweep(res$Std.XChapeau,2,attr(res$ExpliX,"scaled:scale"),FUN="*"),2,attr(res$ExpliX,"scaled:center"),FUN="+")*!XXNA
rownames(res$XChapeau) <- rownames(ExpliX)
colnames(res$XChapeau) <- colnames(ExpliX)
}
names(res$CoeffC) <- paste("Coeff_Comp_Reg",1:nt)
rownames(res$Coeffs) <- c("Intercept",colnames(ExpliX))
}


##############################################
######              PLS-GLM             ######
##############################################
if (modele %in% c("pls-glm-gaussian","pls-glm-logistic")) {
res$YChapeau <- as.matrix(tempregglm$fitted.values)            
rownames(res$YChapeau) <- rownames(ExpliX)

res$Std.ValsPredictY <- predict(tempregglm,newdata=as.data.frame(res$ttPredictY))
res$ValsPredictY <- predict(tempregglm,newdata=as.data.frame(res$ttPredictY),type = "response")

res$Std.XChapeau <- res$tt%*%t(res$pp)
rownames(res$Std.XChapeau) <- rownames(ExpliX)
names(res$CoeffC) <- paste("Coeff_Comp_Reg",1:nt)
rownames(res$Coeffs) <- c("Intercept",colnames(ExpliX))
res$FinalModel <- tempregglm
}


##############################################
######           PLS-GLM-POLR           ######
##############################################
if (modele %in% c("pls-glm-polr")) {
res$YChapeau <- tempregpolr$fitted.values
res$YChapeauCat <- predict(tempregpolr,type="class")
rownames(res$YChapeau) <- rownames(ExpliX)


res$Std.XChapeau <- res$tt%*%t(res$pp)
rownames(res$Std.XChapeau) <- rownames(ExpliX)
names(res$CoeffC) <- paste("Coeff_Comp_Reg",1:nt)
res$MissClassedFinal <- sum(!(unclass(res$YChapeauCat)==unclass(res$RepY)))
res$FinalModel <- tempregpolr
}


rownames(res$pp) <- colnames(ExpliX)
colnames(res$pp) <- paste("Comp_",1:nt)
rownames(res$ww) <- colnames(ExpliX)
colnames(res$ww) <- paste("Comp_",1:nt)
rownames(res$wwnorm) <- colnames(ExpliX)
colnames(res$wwnorm) <- paste("Comp_",1:nt)
rownames(res$wwetoile) <- colnames(ExpliX)
colnames(res$wwetoile) <- paste("Coord_Comp_",1:nt)
rownames(res$tt) <- rownames(ExpliX)
colnames(res$tt) <- paste("Comp_",1:nt)
res$XXwotNA <- XXwotNA
cat("****________________________________________________****\n")
cat("\n")
return(res)
}
