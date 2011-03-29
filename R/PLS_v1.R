PLS_v1 <- function(dataY,dataX,nt=2,limQ2set=.0975,dataPredictY=dataX,modele="pls",
                    family=NULL,typeVC="none",EstimXNA=FALSE,scaleX=TRUE,scaleY=NULL,pvals.expli=FALSE,alpha.pvals.expli=.05,MClassed=FALSE,tol_Xi=10^(-12)) {

##################################################
#                                                #
#    Initialization and formatting the inputs    #
#                                                #
##################################################

cat("____************************************************____\n")
if (!is.data.frame(dataX)) {dataX <- data.frame(dataX)}
if (!(modele %in% c("pls"))) {break}
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

res <- list(nr=nrow(ExpliX),nc=ncol(ExpliX),nt=nt,ww=NULL,wwnorm=NULL,wwetoile=NULL,tt=NULL,pp=NULL,CoeffC=NULL,uscores=NULL,YChapeau=NULL,residYChapeau=NULL,RepY=RepY,na.miss.Y=na.miss.Y,YNA=YNA,residY=RepY,ExpliX=ExpliX,na.miss.X=na.miss.X,XXNA=XXNA,residXX=ExpliX,PredictY=PredictYwotNA,press.ind=NULL,press.tot=NULL,family=family,ttPredictY = NULL,typeVC=typeVC,dataY=dataY)
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
res$pp <- cbind(res$pp,temppp)   
res$tt <- cbind(res$tt,temptt)       




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
if (typeVC == "none") {} else {
if (typeVC %in% c("standard","adaptative")) {
if (kk==1) {
cat(paste("____TypeVC____",typeVC,"____\n"))
}
temppred <- rep(0, res$nr)
for (i in 1:(res$nr)) {    
                tempww.cv <- t(XXwotNA[-i, ])%*%YwotNA[-i]/(t(XXNA[-i, ])%*%YwotNA[-i]^2)
                tempwwnorm.cv <- tempww.cv/sqrt(drop(crossprod(tempww.cv)))
                temptt.cv <- XXwotNA[-i, ]%*%tempwwnorm.cv/(XXNA[-i, ]%*%(tempwwnorm.cv^2))
                tempc.cv <- solve(t(temptt.cv)%*%temptt.cv)%*%t(temptt.cv)%*%(YwotNA[-i])
                tempc.cv <- as.vector(tempc.cv)
                temppred[i] <- tempc.cv * XXwotNA[i, ] %*% tempwwnorm.cv
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
                if(det(t(temppp.cv[XXNA[i,],])%*%temppp.cv[XXNA[i,],])<tol_Xi) {
                break_nt_vc <- TRUE
                res$computed_nt_vc <- kk-1
                cat(paste("Warning : determinant of t(temppp.cv[XXNA[",i,",],])%*%temppp.cv[XXNA[",i,",],]) < 10^{-12}\n",sep=""))
                cat(paste("Please choose a smaller number of component to cross validate\n",sep=""))
                break
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
if (kk==1) {
res$RSSresidY <- crossprod(RepY)
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
                if(det(t(temppp.cv[XXNA[i,],])%*%temppp.cv[XXNA[i,],])<tol_Xi) {
                break_nt_vc <- TRUE
                res$computed_nt_vc <- kk-1
                cat(paste("Warning : determinant of t(temppp.cv[XXNA[",i,",],])%*%temppp.cv[XXNA[",i,",],]) < 10^{-12}\n",sep=""))
                cat(paste("Please choose a smaller number of component to cross validate\n",sep=""))
                break
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
                if(det(t(temppp.cv[XXNA[i,],])%*%temppp.cv[XXNA[i,],])<tol_Xi) {
                break_nt_vc <- TRUE
                res$computed_nt_vc <- kk-1
                cat(paste("Warning : determinant of t(temppp.cv[XXNA[",i,",],])%*%temppp.cv[XXNA[",i,",],]) < 10^{-12}\n",sep=""))
                cat(paste("Please choose a smaller number of component to cross validate\n",sep=""))
                break
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
if (kk==1) {
res$RSSresidY <- crossprod(RepY)
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
res$AIC.std <- cbind(res$AIC.std,AICpls(kk,res$residY))
res$AIC <- AIC(lm(dataY~1))
res$AIC <- cbind(res$AIC,AICpls(kk,res$Yresidus))
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
res$AIC.std <- cbind(res$AIC.std,AICpls(kk,res$residY))
res$AIC <- cbind(res$AIC,AICpls(kk,res$Yresidus))
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

cat("____Component____",kk,"____\n")
if(break_nt_vc==TRUE) {break}
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
if(kk==1){
cat("____Predicting X without NA neither in X nor in Y____\n")
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
if (typeVC %in% c("standard","missingdata","adaptative")) {
res$press.tot <- colSums(res$press.ind)
res$press.tot2 <- colSums(res$press.ind2)
res$Q2 <- 1-res$press.tot/res$RSSresidY[-(res$computed_nt+1)]
res$limQ2 <- rep(limQ2set,res$computed_nt)
res$Q2_2 <- 1-res$press.tot2/res$RSS[-(res$computed_nt+1)]


for (k in 1:res$computed_nt) {res$Q2cum[k] <- prod(res$press.tot[1:k])/prod(res$RSSresidY[1:k])}
res$Q2cum <- 1 - res$Q2cum
for (k in 1:res$computed_nt) {res$Q2cum_2[k] <- prod(res$press.tot2[1:k])/prod(res$RSS[1:k])}
res$Q2cum_2 <- 1 - res$Q2cum_2
if (MClassed==FALSE) {
res$CVinfos <- t(rbind(res$AIC,c(0,res$Q2cum_2), c(NA,res$limQ2), c(0,res$Q2_2[1:res$computed_nt]), c(0,res$press.tot2[1:res$computed_nt]), res$RSS, c(0,res$R2), c(0,res$R2residY), res$RSSresidY, c(0,res$press.tot), c(0,res$Q2), c(NA,res$limQ2), c(0,res$Q2cum), res$AIC.std))
dimnames(res$CVinfos) <- list(paste("Nb_Comp_",0:res$computed_nt), c("AIC", "Q2cum_Y", "LimQ2_Y", "Q2_Y", "PRESS_Y", "RSS_Y", "R2_Y", "R2_residY", "RSS_residY", "PRESS_residY", "Q2_residY", "LimQ2", "Q2cum_residY", "AIC.std"))
} else {
res$CVinfos <- t(rbind(res$AIC,c(0,res$Q2cum_2), c(NA,res$limQ2), c(0,res$Q2_2[1:res$computed_nt]), c(0,res$press.tot2[1:res$computed_nt]), res$RSS, c(0,res$R2), res$MissClassed, c(0,res$R2residY), res$RSSresidY, c(0,res$press.tot), c(0,res$Q2), c(NA,res$limQ2), c(0,res$Q2cum), res$AIC.std))
dimnames(res$CVinfos) <- list(paste("Nb_Comp_",0:res$computed_nt), c("AIC", "Q2cum_Y", "LimQ2_Y", "Q2_Y", "PRESS_Y", "RSS_Y", "R2_Y", "MissClassed", "R2_residY", "RSS_residY", "PRESS_residY", "Q2_residY", "LimQ2", "Q2cum_residY", "AIC.std"))
}


} else {
if (MClassed==FALSE) {
res$InfCrit <- t(rbind(res$AIC, res$RSS, c(0,res$R2), c(0,res$R2residY), res$RSSresidY, res$AIC.std))
dimnames(res$InfCrit) <- list(paste("Nb_Comp_",0:res$computed_nt), c("AIC", "RSS_Y", "R2_Y", "R2_residY", "RSS_residY", "AIC.std"))
} else {
res$InfCrit <- t(rbind(res$AIC, res$RSS, c(0,res$R2), res$MissClassed, c(0,res$R2residY), res$RSSresidY, res$AIC.std))
dimnames(res$InfCrit) <- list(paste("Nb_Comp_",0:res$computed_nt), c("AIC", "RSS_Y", "R2_Y", "MissClassed", "R2_residY", "RSS_residY", "AIC.std"))
}
}
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
