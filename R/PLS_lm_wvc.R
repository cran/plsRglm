PLS_lm_wvc <- function(dataY,dataX,nt=2,dataPredictY=dataX,modele="pls",scaleX=TRUE,scaleY=NULL,keepcoeffs=FALSE,keepstd.coeffs=FALSE,tol_Xi=10^(-12)) {



##################################################
#                                                #
#    Initialization and formatting the inputs    #
#                                                #
##################################################

cat("____************************************************____\n")
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

res <- list(nr=nrow(ExpliX),nc=ncol(ExpliX),ww=NULL,wwnorm=NULL,wwetoile=NULL,tt=NULL,pp=NULL,CoeffC=NULL,uscores=NULL,YChapeau=NULL,residYChapeau=NULL,RepY=RepY,na.miss.Y=na.miss.Y,YNA=YNA,residY=RepY,ExpliX=ExpliX,PredictY=PredictYwotNA,ttPredictY = NULL,residXX=ExpliX,listValsPredictY=NULL)
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

#Search for singularities
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

#Search for singularities in t(pp)%*%pp if there is any missing value in dataX
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
if (!(na.miss.X | na.miss.Y)) {# Cette distinction n'est pas necessaire a priori
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
res$residYChapeau <- res$tt%*%tempCoeffC
#res$RSSresidY[kk] <- crossprod(res$residY-res$residYChapeau)


tempCoeffs <- res$wwetoile%*%res$CoeffC*attr(res$RepY,"scaled:scale")/attr(res$ExpliX,"scaled:scale")
tempConstante <- attr(res$RepY,"scaled:center")-sum(tempCoeffs*attr(res$ExpliX,"scaled:center"))
res$Coeffs <- rbind(tempConstante,tempCoeffs)

res$YChapeau <- attr(res$RepY,"scaled:center")+attr(res$RepY,"scaled:scale")*res$tt%*%res$CoeffC             
res$Yresidus <- dataY-res$YChapeau
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
cat("____Il y a des NA dans X et pas dans Y____\n")
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
}

else {
if (kk==1) {
cat("____Il y a des NA dans X et aussi dans Y____\n")
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
#                                            #
#           Predicting components            #
#                                            #
##############################################

if (!(na.miss.X | na.miss.Y)) {
if(kk==1){
cat("____Predicting X without NA neither in X or Y____\n")
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
if(kk==1){
cat("____Il y a des NA dans X et aussi dans Y____\n")
}
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

