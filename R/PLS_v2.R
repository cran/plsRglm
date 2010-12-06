PLS_v2 <- function(dataY,dataX,nt=2,limQ2set=.0975,dataPredictY=dataX,modele="pls",
                    family=NULL,typeVC="none",EstimXNA=FALSE,scaleX=TRUE,scaleY=NULL,pvals.expli=FALSE,alpha.pvals.expli=.05,MClassed=FALSE,tol_Xi=10^(-12)) {  

##################################################
#                                                #
#    Initialization and formatting the inputs    #
#                                                #
##################################################

cat("____************************************************____\n")
if (!is.data.frame(dataX)) {dataX <- data.frame(dataX)}
if (!(modele %in% c("pls","pls-glm-logistic","pls-glm-gaussian","pls-glm-polr","pls-glm-lrm"))) {break}
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
if(is.null(colnames(ExpliX))){colnames(ExpliX)<-paste("X",1:ncol(ExpliX),sep=".")}
if(is.null(rownames(ExpliX))){rownames(ExpliX)<-1:nrow(ExpliX)}

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

res <- list(nr=nrow(ExpliX),nc=ncol(ExpliX),nt=nt,ww=NULL,wwnorm=NULL,wwetoile=NULL,tt=NULL,pp=NULL,CoeffC=NULL,uscores=NULL,YChapeau=NULL,residYChapeau=NULL,RepY=RepY,na.miss.Y=na.miss.Y,YNA=YNA,residY=RepY,ExpliX=ExpliX,na.miss.X=na.miss.X,XXNA=XXNA,residXX=ExpliX,PredictY=PredictYwotNA,RSS=rep(NA,nt),RSSresidY=rep(NA,nt),R2=rep(NA,nt),R2residY=rep(NA,nt),press.ind=NULL,press.tot=NULL,Q2cum=rep(NA, nt),family=family,ttPredictY = NULL,typeVC=typeVC) 
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

res$computed_nt <- 0
break_nt <- FALSE
break_nt_vc <- FALSE

for (kk in 1:nt) {
XXwotNA <- as.matrix(res$residXX)
XXwotNA[!XXNA] <- 0
YwotNA <- as.matrix(res$residY)
YwotNA[!YNA] <- 0
tempww <- rep(0,res$nc)


temptest <- sqrt(colSums(res$residXX^2, na.rm=TRUE))
if(any(temptest<tol_Xi)) {
break_nt <- TRUE
if (is.null(names(which(temptest<tol_Xi)))) {
print(paste("Warning : ",paste(names(which(temptest<tol_Xi)),sep="",collapse=" ")," < 10^{-12}",sep=""))
} else {
print(paste("Warning : ",paste((which(temptest<tol_Xi)),sep="",collapse=" ")," < 10^{-12}",sep=""))
}
print(paste("Warning only ",res$computed_nt," components could thus be extracted",sep=""))
rm(temptest)
break
}

res$computed_nt <- kk

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
tts <- res$tt
for (jj in 1:(res$nc)) {
    tempww[jj] <- -1*polr(YwotNA~cbind(tts,XXwotNA[,jj]),na.action=na.exclude)$coef[kk]
}
XXwotNA[!XXNA] <- 0
rm(jj,tts)}




##############################################
#                                            #
# Computation of the components (model free) #
#                                            #
##############################################

tempwwnorm <- tempww/sqrt(drop(crossprod(tempww)))

temptt <- XXwotNA%*%tempwwnorm/(XXNA%*%(tempwwnorm^2))

temppp <- rep(0,res$nc)
for (jj in 1:(res$nc)) {
     temppp[jj] <- crossprod(temptt,XXwotNA[,jj])/drop(crossprod(XXNA[,jj],temptt^2))
}
res$residXX <- XXwotNA-temptt%*%temppp

if (na.miss.X & !na.miss.Y) {
for (ii in 1:res$nr) {
if(1/kappa(t(cbind(res$pp,temppp)[XXNA[ii,],])%*%cbind(res$pp,temppp)[XXNA[ii,],])<tol_Xi) {
break_nt <- TRUE
res$computed_nt <- kk-1
cat(paste("Warning : determinant of t(cbind(res$pp,temppp)[XXNA[",ii,",],])%*%cbind(res$pp,temppp)[XXNA[",ii,",],] < 10^{-12}\n",sep=""))
cat(paste("Warning only ",res$computed_nt," components could thus be extracted\n",sep=""))
break
}
}
rm(ii)
if(break_nt==TRUE) {break}
}

if (na.miss.X & !na.miss.Y) {
for (ii in 1:nrow(PredictYwotNA)) {
if(1/kappa(t(cbind(res$pp,temppp)[PredictYNA[ii,],])%*%cbind(res$pp,temppp)[PredictYNA[ii,],])<tol_Xi) {
break_nt <- TRUE
res$computed_nt <- kk-1
cat(paste("Warning : determinant of t(cbind(res$pp,temppp)[PredictYNA[",ii,",],])%*%cbind(res$pp,temppp)[PredictYNA[",ii,",],] < 10^{-12}\n",sep=""))
cat(paste("Warning only ",res$computed_nt," components could thus be extracted\n",sep=""))
break
}
}
rm(ii)
if(break_nt==TRUE) {break}
}

res$ww <- cbind(res$ww,tempww)
res$wwnorm <- cbind(res$wwnorm,tempwwnorm)
res$tt <- cbind(res$tt,temptt)       
res$pp <- cbind(res$pp,temppp)   




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
if (!(na.miss.X | na.miss.Y)) {
tempCoeffC <- c(rep(0,kk-1),solve(t(res$tt[YNA,kk])%*%res$tt[YNA,kk])%*%t(res$tt[YNA,kk])%*%YwotNA[YNA])  
tempCoeffConstante <- 0
res$CoeffCFull <- cbind(res$CoeffCFull,c(tempCoeffC,rep(NA,nt-kk)))
}
else
{
tempCoeffC <- c(rep(0,kk-1),solve(t(res$tt[YNA,kk])%*%res$tt[YNA,kk])%*%t(res$tt[YNA,kk])%*%YwotNA[YNA])  
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
if (modele %in% c("pls-glm-logistic")) {
res$MissClassed <- sum(unclass(res$RepY)!=ifelse(predict(tempconstglm,type="response") < 0.5, 0,1))
}
rm(tempconstglm)
tt<-res$tt
tempregglm <- glm(YwotNA~tt,family=family)
rm(tt)
res$AIC <- cbind(res$AIC,AIC(tempregglm))
res$BIC <- cbind(res$BIC,AIC(tempregglm, k = log(res$nr)))
res$Coeffsmodel_vals <- cbind(res$Coeffsmodel_vals,rbind(summary(tempregglm)$coefficients,matrix(rep(NA,4*(nt-kk)),ncol=4)))
res$ChisqPearson <- c(res$ChisqPearson,crossprod(residuals.glm(tempregglm,type="pearson")))
if (modele %in% c("pls-glm-logistic")) {
res$MissClassed <- cbind(res$MissClassed,sum(unclass(res$RepY)!=ifelse(predict(tempregglm,type="response") < 0.5, 0,1)))
}
tempCoeffC <- as.vector(coef(tempregglm))
res$CoeffCFull <- matrix(c(tempCoeffC,rep(NA,nt-kk)),ncol=1)
tempCoeffConstante <- tempCoeffC[1]
res$CoeffConstante <- tempCoeffConstante
tempCoeffC <- tempCoeffC[-1]
} else {
if (!(na.miss.X | na.miss.Y)) {
tt<-res$tt
tempregglm <- glm(YwotNA~tt,family=family)
rm(tt)
res$AIC <- cbind(res$AIC,AIC(tempregglm))
res$BIC <- cbind(res$BIC,AIC(tempregglm, k = log(res$nr)))
res$Coeffsmodel_vals <- cbind(res$Coeffsmodel_vals,rbind(summary(tempregglm)$coefficients,matrix(rep(NA,4*(nt-kk)),ncol=4)))
res$ChisqPearson <- c(res$ChisqPearson,crossprod(residuals.glm(tempregglm,type="pearson")))
if (modele %in% c("pls-glm-logistic")) {
res$MissClassed <- cbind(res$MissClassed,sum(unclass(res$RepY)!=ifelse(predict(tempregglm,type="response") < 0.5, 0,1)))
}
tempCoeffC <- as.vector(coef(tempregglm))  
res$CoeffCFull <- cbind(res$CoeffCFull,c(tempCoeffC,rep(NA,nt-kk)))
tempCoeffConstante <- tempCoeffC[1]
res$CoeffConstante <- cbind(res$CoeffConstante,tempCoeffConstante)
tempCoeffC <- tempCoeffC[-1]
}
else
{
tt<-res$tt
tempregglm <- glm(YwotNA~tt,family=family)
rm(tt)
res$AIC <- cbind(res$AIC,AIC(tempregglm))
res$BIC <- cbind(res$BIC,AIC(tempregglm, k = log(res$nr)))
res$Coeffsmodel_vals <- cbind(res$Coeffsmodel_vals,rbind(summary(tempregglm)$coefficients,matrix(rep(NA,4*(nt-kk)),ncol=4)))
res$ChisqPearson <- c(res$ChisqPearson,crossprod(residuals.glm(tempregglm,type="pearson")))
if (modele %in% c("pls-glm-logistic")) {
res$MissClassed <- cbind(res$MissClassed,sum(unclass(res$RepY)!=ifelse(predict(tempregglm,type="response") < 0.5, 0,1)))
}
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
            Varyy <- function(piVaryy) {
            diag(piVaryy[-length(piVaryy)])-piVaryy[-length(piVaryy)]%*%t(piVaryy[-length(piVaryy)])
            }
            Chisqcomp <- function(yichisq,pichisq) {
            t(yichisq[-length(yichisq)]-pichisq[-length(pichisq)])%*%solve(Varyy(pichisq))%*%(yichisq[-length(yichisq)]-pichisq[-length(pichisq)])
            }
            Chiscompmatrix <- function(rowspi,rowsyi) {
            sum(mapply(Chisqcomp,rowsyi,rowspi))
            }
if (kk==1) {
tempconstpolr <- polr(YwotNA~1,na.action=na.exclude,Hess=TRUE)
res$AIC <- AIC(tempconstpolr)
res$BIC <- AIC(tempconstpolr, k = log(res$nr))
res$MissClassed <- sum(!(unclass(predict(tempconstpolr,type="class"))==unclass(res$RepY)))
res$Coeffsmodel_vals <- rbind(summary(tempconstpolr)$coefficients,matrix(rep(NA,3*nt),ncol=3))
tempmodord <- predict(tempconstpolr,type="class")
tempfff <- ~tempmodord-1
tempm <- model.frame(tempfff, tempmodord)
tempmat <- model.matrix(tempfff, model.frame(tempfff, tempmodord))
res$ChisqPearson <- sum(Chiscompmatrix(as.list(as.data.frame(t(predict(tempconstpolr,type="probs")))),as.list(as.data.frame(t(tempmat)))))
rm(tempconstpolr)
tts<-res$tt
tempregpolr <- polr(YwotNA~tts,na.action=na.exclude,Hess=TRUE)
rm(tts)
res$AIC <- cbind(res$AIC,AIC(tempregpolr))
res$BIC <- cbind(res$BIC,AIC(tempregpolr, k = log(res$nr)))
res$MissClassed <- cbind(res$MissClassed,sum(!(unclass(predict(tempregpolr,type="class"))==unclass(res$RepY))))
res$Coeffsmodel_vals <- cbind(res$Coeffsmodel_vals,rbind(summary(tempregpolr)$coefficients,matrix(rep(NA,3*(nt-kk)),ncol=3)))
tempmodord <- predict(tempregpolr,type="class")
tempfff <- ~tempmodord-1
tempm <- model.frame(tempfff, tempmodord)
tempmat <- model.matrix(tempfff, model.frame(tempfff, tempmodord))
res$ChisqPearson <- c(res$ChisqPearson,sum(Chiscompmatrix(as.list(as.data.frame(t(predict(tempregpolr,type="probs")))),as.list(as.data.frame(t(tempmat))))))
tempCoeffC <- -1*as.vector(tempregpolr$coef)
tempCoeffConstante <- as.vector(tempregpolr$zeta)
res$CoeffCFull <- matrix(c(tempCoeffConstante,tempCoeffC,rep(NA,nt-kk)),ncol=1)
res$CoeffConstante <- tempCoeffConstante
} else {
if (!(na.miss.X | na.miss.Y)) {
tts <- res$tt
tempregpolr <- polr(YwotNA~tts,na.action=na.exclude,Hess=TRUE)
rm(tts)
res$AIC <- cbind(res$AIC,AIC(tempregpolr))
res$BIC <- cbind(res$BIC,AIC(tempregpolr, k = log(res$nr)))
res$MissClassed <- cbind(res$MissClassed,sum(!(unclass(predict(tempregpolr,type="class"))==unclass(res$RepY))))
res$Coeffsmodel_vals <- cbind(res$Coeffsmodel_vals,rbind(summary(tempregpolr)$coefficients,matrix(rep(NA,3*(nt-kk)),ncol=3)))
tempmodord <- predict(tempregpolr,type="class")
tempfff <- ~tempmodord-1
tempm <- model.frame(tempfff, tempmodord)
tempmat <- model.matrix(tempfff, model.frame(tempfff, tempmodord))
res$ChisqPearson <- c(res$ChisqPearson,sum(Chiscompmatrix(as.list(as.data.frame(t(predict(tempregpolr,type="probs")))),as.list(as.data.frame(t(tempmat))))))
tempCoeffC <- -1*as.vector(tempregpolr$coef)  
tempCoeffConstante <- as.vector(tempregpolr$zeta)
res$CoeffCFull <- cbind(res$CoeffCFull,c(tempCoeffConstante,tempCoeffC,rep(NA,nt-kk)))
res$CoeffConstante <- cbind(res$CoeffConstante,tempCoeffConstante)
}
else
{
tts<-res$tt
tempregpolr <- polr(YwotNA~tts,na.action=na.exclude,Hess=TRUE)
rm(tts)
res$AIC <- cbind(res$AIC,AIC(tempregpolr))
res$BIC <- cbind(res$BIC,AIC(tempregpolr, k = log(res$nr)))
res$MissClassed <- cbind(res$MissClassed,sum(!(unclass(predict(tempregpolr,type="class"))==unclass(res$RepY))))
res$Coeffsmodel_vals <- cbind(res$Coeffsmodel_vals,rbind(summary(tempregpolr)$coefficients,matrix(rep(NA,3*(nt-kk)),ncol=3)))
tempmodord <- predict(tempregpolr,type="class")
tempfff <- ~tempmodord-1
tempm <- model.frame(tempfff, tempmodord)
tempmat <- model.matrix(tempfff, model.frame(tempfff, tempmodord))
res$ChisqPearson <- c(res$ChisqPearson,sum(Chiscompmatrix(as.list(as.data.frame(t(predict(tempregpolr,type="probs")))),as.list(as.data.frame(t(tempmat))))))
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
res$residYChapeau <- res$tt%*%tempCoeffC
if (kk==1) {
res$RSSresidY <- crossprod(RepY-mean(RepY))
}
res$RSSresidY <- cbind(res$RSSresidY,crossprod(res$residY-res$residYChapeau))


tempCoeffs <- res$wwetoile%*%res$CoeffC*attr(res$RepY,"scaled:scale")/attr(res$ExpliX,"scaled:scale")
tempConstante <- attr(res$RepY,"scaled:center")-sum(tempCoeffs*attr(res$ExpliX,"scaled:center"))
res$Coeffs <- rbind(tempConstante,tempCoeffs)

res$YChapeau <- attr(res$RepY,"scaled:center")+attr(res$RepY,"scaled:scale")*res$tt%*%res$CoeffC             
res$Yresidus <- dataY-res$YChapeau
if (kk==1) {
res$RSS <- crossprod(dataY-mean(dataY))
}
res$RSS <- cbind(res$RSS,crossprod(res$Yresidus))
}
##############################################


##############################################
######              PLS-GLM             ######
##############################################
if (modele %in% c("pls-glm-gaussian","pls-glm-logistic")) {
res$residYChapeau <- tempregglm$linear.predictors
if (kk==1) {
res$RSSresidY <- crossprod(RepY-mean(RepY))
}
res$RSSresidY <- cbind(res$RSSresidY,crossprod(res$residY-res$residYChapeau))


tempCoeffs <- res$wwetoile%*%res$CoeffC*attr(res$RepY,"scaled:scale")/attr(res$ExpliX,"scaled:scale")
tempConstante <- attr(res$RepY,"scaled:center")-sum(tempCoeffs*attr(res$ExpliX,"scaled:center"))+attr(res$RepY,"scaled:scale")*res$Std.Coeffs[1]
res$Coeffs <- rbind(tempConstante,tempCoeffs)

res$YChapeau <- tempregglm$fitted.values          
res$Yresidus <- dataY-res$YChapeau
if (kk==1) {
res$RSS <- crossprod(dataY-mean(dataY))
}
res$RSS <- cbind(res$RSS,crossprod(res$Yresidus))
}

if (modele %in% c("pls-glm-polr")) {
tempCoeffs <- res$wwetoile%*%res$CoeffC*attr(res$RepY,"scaled:scale")/attr(res$ExpliX,"scaled:scale")
tempConstante <- attr(res$RepY,"scaled:center")-sum(tempCoeffs*attr(res$ExpliX,"scaled:center"))+attr(res$RepY,"scaled:scale")*tempCoeffConstante
res$Coeffs <- rbind(as.matrix(tempConstante),tempCoeffs)
rownames(res$Coeffs) <- rownames(res$Std.Coeffs)
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
res$residYChapeau <- res$tt%*%tempCoeffC
if (kk==1) {
res$RSSresidY <- crossprod(RepY-mean(RepY))
}
res$RSSresidY <- cbind(res$RSSresidY,crossprod(res$residY-res$residYChapeau))


tempCoeffs <- res$wwetoile%*%res$CoeffC*attr(res$RepY,"scaled:scale")/attr(res$ExpliX,"scaled:scale")
tempConstante <- attr(res$RepY,"scaled:center")-sum(tempCoeffs*attr(res$ExpliX,"scaled:center"))
res$Coeffs <- rbind(tempConstante,tempCoeffs)

res$YChapeau <- attr(res$RepY,"scaled:center")+attr(res$RepY,"scaled:scale")*res$tt%*%res$CoeffC            
res$Yresidus <- dataY-res$YChapeau
if (kk==1) {
res$RSS <- crossprod(dataY-mean(dataY))
}
res$RSS <- cbind(res$RSS,crossprod(res$Yresidus))
}
##############################################



##############################################
######              PLS-GLM             ######
##############################################
if (modele %in% c("pls-glm-gaussian","pls-glm-logistic")) {
res$residYChapeau <- tempregglm$linear.predictors
if (kk==1) {
res$RSSresidY <- crossprod(RepY-mean(RepY))
}
res$RSSresidY <- cbind(res$RSSresidY,crossprod(res$residY-res$residYChapeau))

tempCoeffs <- res$wwetoile%*%res$CoeffC*attr(res$RepY,"scaled:scale")/attr(res$ExpliX,"scaled:scale")
tempConstante <- attr(res$RepY,"scaled:center")-sum(tempCoeffs*attr(res$ExpliX,"scaled:center"))+attr(res$RepY,"scaled:scale")*res$Std.Coeffs[1]
res$Coeffs <- rbind(tempConstante,tempCoeffs)

res$YChapeau <- tempregglm$fitted.values                      
res$Yresidus <- dataY-res$YChapeau
if (kk==1) {
res$RSS <- crossprod(dataY-mean(dataY))
}
res$RSS <- cbind(res$RSS,crossprod(res$Yresidus))
}

if (modele %in% c("pls-glm-polr")) {
tempCoeffs <- res$wwetoile%*%res$CoeffC*attr(res$RepY,"scaled:scale")/attr(res$ExpliX,"scaled:scale")
tempConstante <- attr(res$RepY,"scaled:center")-sum(tempCoeffs*attr(res$ExpliX,"scaled:center"))+attr(res$RepY,"scaled:scale")* tempCoeffConstante
res$Coeffs <- rbind(as.matrix(tempConstante),tempCoeffs)
rownames(res$Coeffs) <- rownames(res$Std.Coeffs)
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


if (kk==1) {
res$AIC.std <- AIC(lm(res$RepY~1))
res$AIC.std <- cbind(res$AIC.std,AICpls(res$nr,kk,res$residY))
res$AIC <- AIC(lm(dataY~1))
res$AIC <- cbind(res$AIC,AICpls2(kk,dataY,res$YChapeau,res$Yresidus))
if (MClassed) {
res$MissClassed <- sum(unclass(dataY)!=ifelse(predict(lm(dataY~1)) < 0.5, 0,1))
res$MissClassed <- cbind(res$MissClassed,sum(unclass(dataY)!=ifelse(res$YChapeau < 0.5, 0,1)))
tempprob <- res$Probs <- predict(lm(dataY~1))
tempprob <- ifelse(tempprob<0,0,tempprob)
res$Probs.trc <- ifelse(tempprob>1,1,tempprob)
res$Probs <- cbind(res$Probs,res$YChapeau)
tempprob <- ifelse(res$YChapeau<0,0,res$YChapeau)
tempprob <- ifelse(tempprob>1,1,tempprob)
res$Probs.trc <- cbind(res$Probs.trc,tempprob)
}
} else {
res$AIC.std <- cbind(res$AIC.std,AICpls(res$nr,kk,res$residY))
res$AIC <- cbind(res$AIC,AICpls2(kk,dataY,res$YChapeau,res$Yresidus))
if (MClassed) {
res$MissClassed <- cbind(res$MissClassed,sum(unclass(dataY)!=ifelse(res$YChapeau < 0.5, 0,1)))
res$Probs <- cbind(res$Probs,res$YChapeau)
tempprob <- ifelse(res$YChapeau<0,0,res$YChapeau)
tempprob <- ifelse(tempprob>1,1,tempprob)
res$Probs.trc <- cbind(res$Probs.trc,tempprob)
}
}


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
rm(tempCoeffs)
rm(tempConstante)
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


if (pvals.expli) {
res$Coeffsmodel_vals<-res$Coeffsmodel_vals[1:(dim(res$Coeffsmodel_vals)[1]-(nt-res$computed_nt)),]
}


##############################################
#                                            #
#           Predicting components            #
#                                            #
##############################################

if (!(na.miss.X | na.miss.Y)) {
if(kk==1){
cat("____Predicting X without NA neither in X or Y____\n")
}
res$ttPredictY <- PredictYwotNA%*%res$wwetoile 
colnames(res$ttPredictY) <- paste("tt",1:res$computed_nt,sep="")
}
else {
if (na.miss.X & !na.miss.Y) {
if(kk==1){
cat("____Predicting X with NA in X and not in Y____\n")
}
for (ii in 1:nrow(PredictYwotNA)) {  
      res$ttPredictY <- rbind(res$ttPredictY,t(solve(t(res$pp[PredictYNA[ii,],])%*%res$pp[PredictYNA[ii,],])%*%t(res$pp[PredictYNA[ii,],])%*%(PredictYwotNA[ii,])[PredictYNA[ii,]]))
}
colnames(res$ttPredictY) <- paste("tt",1:res$computed_nt,sep="")
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
res$R2residY <- 1-res$RSSresidY[2:(res$computed_nt+1)]/res$RSSresidY[1]
res$R2 <- 1-res$RSS[2:(res$computed_nt+1)]/res$RSS[1]
if (MClassed==FALSE) {
res$InfCrit <- t(rbind(res$AIC, res$RSS, c(0,res$R2), c(0,res$R2residY), res$RSSresidY, res$AIC.std))
dimnames(res$InfCrit) <- list(paste("Nb_Comp_",0:res$computed_nt), c("AIC", "RSS_Y", "R2_Y", "R2_residY", "RSS_residY", "AIC.std"))
} else {
res$InfCrit <- t(rbind(res$AIC, res$RSS, c(0,res$R2), res$MissClassed, c(0,res$R2residY), res$RSSresidY, res$AIC.std))
dimnames(res$InfCrit) <- list(paste("Nb_Comp_",0:res$computed_nt), c("AIC", "RSS_Y", "R2_Y", "MissClassed", "R2_residY", "RSS_residY", "AIC.std"))
}
}


##############################################
######              PLS-GLM             ######
##############################################
if (modele %in% c("pls-glm-gaussian", "pls-glm-logistic")) {
res$R2residY <- 1-res$RSSresidY[2:(res$computed_nt+1)]/res$RSSresidY[1]
res$R2 <- 1-res$RSS[2:(res$computed_nt+1)]/res$RSS[1]
if (modele %in% c("pls-glm-gaussian")) {
res$InfCrit <- t(rbind(res$AIC, res$BIC, res$ChisqPearson, res$RSS, c(0,res$R2), c(0,res$R2residY), res$RSSresidY))
dimnames(res$InfCrit) <- list(paste("Nb_Comp_",0:res$computed_nt), c("AIC", "BIC", "Chi2_Pearson_Y", "RSS_Y", "R2_Y", "R2_residY", "RSS_residY"))
}
if (modele %in% c("pls-glm-logistic")) {
res$InfCrit <- t(rbind(res$AIC, res$BIC, res$MissClassed, res$ChisqPearson, res$RSS, c(0,res$R2), c(0,res$R2residY), res$RSSresidY))
dimnames(res$InfCrit) <- list(paste("Nb_Comp_",0:res$computed_nt), c("AIC", "BIC", "Missclassed", "Chi2_Pearson_Y", "RSS_Y", "R2_Y", "R2_residY", "RSS_residY"))
}
}


##############################################
######           PLS-GLM-POLR           ######
##############################################
if (modele == "pls-glm-polr") {

res$InfCrit <- t(rbind(res$AIC, res$BIC, res$MissClassed, res$ChisqPearson))
dimnames(res$InfCrit) <- list(paste("Nb_Comp_",0:res$computed_nt), c("AIC", "BIC", "Missclassed", "Chi2_Pearson_Y"))
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
names(res$CoeffC) <- paste("Coeff_Comp_Reg",1:res$computed_nt)
rownames(res$Coeffs) <- c("Intercept",colnames(ExpliX))
}


##############################################
######              PLS-GLM             ######
##############################################
if (modele %in% c("pls-glm-gaussian","pls-glm-logistic")) {
res$YChapeau <- as.matrix(tempregglm$fitted.values)            
rownames(res$YChapeau) <- rownames(ExpliX)

tt <- res$ttPredictY
res$Std.ValsPredictY <- predict(tempregglm,newdata=data.frame(tt))
res$ValsPredictY <- predict(tempregglm,newdata=data.frame(tt),type = "response")

res$Std.XChapeau <- res$tt%*%t(res$pp)
rownames(res$Std.XChapeau) <- rownames(ExpliX)
names(res$CoeffC) <- paste("Coeff_Comp_Reg",1:res$computed_nt)
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

res$ValsPredictY <- predict(tempregpolr, data.frame(tts=I(res$ttPredictY)),type="probs")
res$ValsPredictYCat <- predict(tempregpolr, data.frame(tts=I(res$ttPredictY)),type="class")

res$Std.XChapeau <- res$tt%*%t(res$pp)
rownames(res$Std.XChapeau) <- rownames(ExpliX)
names(res$CoeffC) <- paste("Coeff_Comp_Reg",1:res$computed_nt)
res$FinalModel <- tempregpolr
}


rownames(res$pp) <- colnames(ExpliX)
colnames(res$pp) <- paste("Comp_",1:res$computed_nt)
rownames(res$ww) <- colnames(ExpliX)
colnames(res$ww) <- paste("Comp_",1:res$computed_nt)
rownames(res$wwnorm) <- colnames(ExpliX)
colnames(res$wwnorm) <- paste("Comp_",1:res$computed_nt)
rownames(res$wwetoile) <- colnames(ExpliX)
colnames(res$wwetoile) <- paste("Coord_Comp_",1:res$computed_nt)
rownames(res$tt) <- rownames(ExpliX)
colnames(res$tt) <- paste("Comp_",1:res$computed_nt)
res$XXwotNA <- XXwotNA
cat("****________________________________________________****\n")
cat("\n")
return(res)
}
