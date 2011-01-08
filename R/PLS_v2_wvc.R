PLS_v2_wvc <- function(dataY,dataX,nt=2,dataPredictY=dataX,modele="pls",family=NULL,scaleX=TRUE,scaleY=NULL,keepcoeffs=FALSE,keepstd.coeffs=FALSE,tol_Xi=10^(-12)) {

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

res <- list(nr=nrow(ExpliX),nc=ncol(ExpliX),ww=NULL,wwnorm=NULL,wwetoile=NULL,tt=NULL,pp=NULL,CoeffC=NULL,uscores=NULL,YChapeau=NULL,residYChapeau=NULL,RepY=RepY,na.miss.Y=na.miss.Y,YNA=YNA,residY=RepY,ExpliX=ExpliX,na.miss.X=na.miss.X,XXNA=XXNA,residXX=ExpliX,PredictY=PredictYwotNA,RSS=rep(NA,nt),RSSresidY=rep(NA,nt),R2=rep(NA,nt),R2residY=rep(NA,nt),press.ind=NULL,press.tot=NULL,Q2cum=rep(NA, nt),family=family,ttPredictY = NULL,typeVC="none",listValsPredictY=NULL) 
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

temptest <- sqrt(colSums(res$residXX^2, na.rm=TRUE))
if(any(temptest<tol_Xi)) {
break_nt <- TRUE
if (is.null(names(which(temptest<tol_Xi)))) {
cat(paste("Warning : ",paste(names(which(temptest<tol_Xi)),sep="",collapse=" ")," < 10^{-12}\n",sep=""))
} else {
cat(paste("Warning : ",paste((which(temptest<tol_Xi)),sep="",collapse=" ")," < 10^{-12}\n",sep=""))
}
cat(paste("Warning only ",res$computed_nt," components could thus be extracted\n",sep=""))
break
}

res$computed_nt <- kk

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
XXwotNA[!XXNA] <- NA
for (jj in 1:(res$nc)) {
    tempww[jj] <- coef(glm(YwotNA~cbind(res$tt,XXwotNA[,jj]),family=family))[kk+1]
}
XXwotNA[!XXNA] <- 0
rm(jj)}


##############################################
######           PLS-GLM-POLR           ######
##############################################
if (modele %in% c("pls-glm-polr")) {
YwotNA <- as.factor(YwotNA)
XXwotNA[!XXNA] <- NA
library(MASS)
for (jj in 1:(res$nc)) {
    tempww[jj] <- -1*MASS:::polr(YwotNA~cbind(res$tt,XXwotNA[,jj]),na.action=na.exclude)$coef[kk]
}
XXwotNA[!XXNA] <- 0
rm(jj)}




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
res$Coeffsmodel_vals <- rbind(summary(tempconstglm)$coefficients,matrix(rep(NA,4*nt),ncol=4))
rm(tempconstglm)
tt<-res$tt
tempregglm <- glm(YwotNA~tt,family=family)
rm(tt)
res$Coeffsmodel_vals <- cbind(res$Coeffsmodel_vals,rbind(summary(tempregglm)$coefficients,matrix(rep(NA,4*(nt-kk)),ncol=4)))
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
res$Coeffsmodel_vals <- cbind(res$Coeffsmodel_vals,rbind(summary(tempregglm)$coefficients,matrix(rep(NA,4*(nt-kk)),ncol=4)))
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
res$Coeffsmodel_vals <- cbind(res$Coeffsmodel_vals,rbind(summary(tempregglm)$coefficients,matrix(rep(NA,4*(nt-kk)),ncol=4)))
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
tempconstpolr <- MASS:::polr(YwotNA~1,na.action=na.exclude,Hess=TRUE)
res$Coeffsmodel_vals <- rbind(summary(tempconstpolr)$coefficients,matrix(rep(NA,3*nt),ncol=3))
rm(tempconstpolr)
tempregpolr <- MASS:::polr(YwotNA~res$tt,na.action=na.exclude,Hess=TRUE)
res$Coeffsmodel_vals <- cbind(res$Coeffsmodel_vals,rbind(summary(tempregpolr)$coefficients,matrix(rep(NA,3*(nt-kk)),ncol=3)))
tempCoeffC <- -1*as.vector(tempregpolr$coef)
tempCoeffConstante <- as.vector(tempregpolr$zeta)
res$CoeffCFull <- matrix(c(tempCoeffConstante,tempCoeffC,rep(NA,nt-kk)),ncol=1)
res$CoeffConstante <- tempCoeffConstante
} else {
if (!(na.miss.X | na.miss.Y)) {
tempregpolr <- MASS:::polr(YwotNA~res$tt,na.action=na.exclude,Hess=TRUE)
res$Coeffsmodel_vals <- cbind(res$Coeffsmodel_vals,rbind(summary(tempregpolr)$coefficients,matrix(rep(NA,3*(nt-kk)),ncol=3)))
tempCoeffC <- -1*as.vector(tempregpolr$coef)  
tempCoeffConstante <- as.vector(tempregpolr$zeta)
res$CoeffCFull <- cbind(res$CoeffCFull,c(tempCoeffConstante,tempCoeffC,rep(NA,nt-kk)))
res$CoeffConstante <- cbind(res$CoeffConstante,tempCoeffConstante)
}
else
{
tempregpolr <- MASS:::polr(YwotNA~res$tt,na.action=na.exclude,Hess=TRUE)
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


tempCoeffs <- res$wwetoile%*%res$CoeffC*attr(res$RepY,"scaled:scale")/attr(res$ExpliX,"scaled:scale")
tempConstante <- attr(res$RepY,"scaled:center")-sum(tempCoeffs*attr(res$ExpliX,"scaled:center"))
res$Coeffs <- rbind(tempConstante,tempCoeffs)

res$YChapeau <- attr(res$RepY,"scaled:center")+attr(res$RepY,"scaled:scale")*res$tt%*%res$CoeffC             
res$Yresidus <- dataY-res$YChapeau
}
##############################################


##############################################
######              PLS-GLM             ######
##############################################
if (modele %in% c("pls-glm-gaussian","pls-glm-logistic")) {
res$residYChapeau <- tempregglm$linear.predictors


tempCoeffs <- res$wwetoile%*%res$CoeffC*attr(res$RepY,"scaled:scale")/attr(res$ExpliX,"scaled:scale")
tempConstante <- attr(res$RepY,"scaled:center")-sum(tempCoeffs*attr(res$ExpliX,"scaled:center"))+attr(res$RepY,"scaled:scale")*res$Std.Coeffs[1]
res$Coeffs <- rbind(tempConstante,tempCoeffs)

res$YChapeau <- tempregglm$fitted.values          
res$Yresidus <- dataY-res$YChapeau
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


tempCoeffs <- res$wwetoile%*%res$CoeffC*attr(res$RepY,"scaled:scale")/attr(res$ExpliX,"scaled:scale")
tempConstante <- attr(res$RepY,"scaled:center")-sum(tempCoeffs*attr(res$ExpliX,"scaled:center"))
res$Coeffs <- rbind(tempConstante,tempCoeffs)

res$YChapeau <- attr(res$RepY,"scaled:center")+attr(res$RepY,"scaled:scale")*res$tt%*%res$CoeffC            
res$Yresidus <- dataY-res$YChapeau
}
##############################################



##############################################
######              PLS-GLM             ######
##############################################
if (modele %in% c("pls-glm-gaussian","pls-glm-logistic")) {
res$residYChapeau <- tempregglm$linear.predictors


tempCoeffs <- res$wwetoile%*%res$CoeffC*attr(res$RepY,"scaled:scale")/attr(res$ExpliX,"scaled:scale")
tempConstante <- attr(res$RepY,"scaled:center")-sum(tempCoeffs*attr(res$ExpliX,"scaled:center"))
res$Coeffs <- rbind(tempConstante,tempCoeffs)

res$YChapeau <- tempregglm$fitted.values                      
res$Yresidus <- dataY-res$YChapeau
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



##############################################
#                                            #
#           Predicting components            #
#                                            #
##############################################

if (!(na.miss.X | na.miss.Y)) {
if(kk==1){
cat("____Predicting X without NA neither in X nor in Y____\n")
}
res$ttPredictY <- PredictYwotNA%*%res$wwetoile 
colnames(res$ttPredictY) <- paste("tt",1:kk,sep="")
}
else {
if (na.miss.X & !na.miss.Y) {
if(kk==1){
cat("____Predicting X with NA in X and not in Y____\n")
}
res$ttPredictY <- NULL
for (ii in 1:nrow(PredictYwotNA)) {  
      res$ttPredictY <- rbind(res$ttPredictY,t(solve(t(res$pp[PredictYNA[ii,],])%*%res$pp[PredictYNA[ii,],])%*%t(res$pp[PredictYNA[ii,],])%*%(PredictYwotNA[ii,])[PredictYNA[ii,]]))
}
colnames(res$ttPredictY) <- paste("tt",1:kk,sep="")
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


##############################################
######              PLS-GLM             ######
##############################################


##############################################
######           PLS-GLM-POLR           ######
##############################################





##########################################
#                                        #
#          Predicting responses          #
#                                        #
##########################################


##############################################
######               PLS                ######
##############################################
if (modele == "pls") {
res$listValsPredictY <- cbind(res$listValsPredictY,attr(res$RepY,"scaled:center")+attr(res$RepY,"scaled:scale")*res$ttPredictY%*%res$CoeffC)
}


##############################################
######              PLS-GLM             ######
##############################################
if (modele %in% c("pls-glm-gaussian","pls-glm-logistic")) {
tt <- res$ttPredictY
res$listValsPredictY <- cbind(res$listValsPredictY,predict(object=tempregglm,newdata=data.frame(tt),type = "response")) 
}


##############################################
######           PLS-GLM-POLR           ######
##############################################
if (modele %in% c("pls-glm-polr")) {
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


cat("****________________________________________________****\n")
cat("\n")
if (!keepcoeffs) {
if (!keepstd.coeffs) {return(list(valsPredict=res$listValsPredictY))} else {return(list(valsPredict=res$listValsPredictY, std.coeffs=res$Std.Coeffs))}}
else {
if (!keepstd.coeffs) {return(list(valsPredict=res$listValsPredictY, coeffs=res$Coeffs))} else {return(list(valsPredict=res$listValsPredictY, coeffs=res$Coeffs, std.coeffs=res$Std.Coeffs))}
}
}
