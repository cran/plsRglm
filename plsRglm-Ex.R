pkgname <- "plsRglm"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('plsRglm')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("AICpls")
### * AICpls

flush(stderr()); flush(stdout())

### Name: AICpls
### Title: AIC function for plsR models
### Aliases: AICpls
### Keywords: models regression utilities

### ** Examples

data(pine)
ypine <- pine[,11]
Xpine <- pine[,1:10]
(Pinscaled <- as.data.frame(cbind(scale(log(ypine)),scale(as.matrix(Xpine)))))
colnames(Pinscaled)[1] <- "yy"

lm(yy~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10,data=Pinscaled)

modpls <- plsR(log(ypine),Xpine,10)
modpls$Std.Coeffs
lm(yy~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10,data=Pinscaled)

AIC(lm(yy~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10,data=Pinscaled))
print(logLik(lm(yy~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10,data=Pinscaled)))

sum(dnorm(modpls$RepY, modpls$Std.ValsPredictY, sqrt(mean(modpls$residY^2)), log=TRUE))
sum(dnorm(Pinscaled$yy,fitted(lm(yy~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10,data=Pinscaled)),sqrt(mean(residuals(lm(yy~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10,data=Pinscaled))^2)), log=TRUE))
loglikpls(modpls$residY)
loglikpls(residuals(lm(yy~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10,data=Pinscaled)))
AICpls(10,residuals(lm(yy~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10,data=Pinscaled)))
AICpls(10,modpls$residY)



cleanEx()
nameEx("CorMat")
### * CorMat

flush(stderr()); flush(stdout())

### Name: CorMat
### Title: Correlation matrix for simulating plsR datasets
### Aliases: CorMat
### Keywords: datasets

### ** Examples

data(CorMat)
## maybe str(CorMat) ; plot(CorMat) ...



cleanEx()
nameEx("Cornell")
### * Cornell

flush(stderr()); flush(stdout())

### Name: Cornell
### Title: Cornell dataset
### Aliases: Cornell
### Keywords: datasets

### ** Examples

data(Cornell)
## maybe str(Cornell) ; plot(Cornell) ...



cleanEx()
nameEx("PLS_glm")
### * PLS_glm

flush(stderr()); flush(stdout())

### Name: PLS_glm
### Title: Partial least squares Regression generalized linear models
### Aliases: PLS_glm
### Keywords: models regression

### ** Examples

data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
PLS_glm(yCornell,XCornell,3)$uscores
PLS_glm(yCornell,XCornell,3)$pp
PLS_glm(yCornell,XCornell,3)$Coeffs
PLS_glm(yCornell,XCornell,10)$InfCrit
PLS_glm(yCornell,XCornell,10,modele="pls-glm-gaussian")$InfCrit
data.frame(pls=PLS_glm(yCornell,XCornell,3)$Coeffs,pls_glm_gaussian=PLS_glm(yCornell,XCornell,10,modele="pls-glm-gaussian")$Coeffs,pls_glm_family_gaussian=PLS_glm(yCornell,XCornell,10,modele="pls-glm-family",family=gaussian())$Coeffs)

mod <- PLS_glm(yCornell,XCornell,10,pvals.expli =TRUE)
mod2 <- PLS_glm(yCornell,XCornell,10,sparse=TRUE)
mod3 <- PLS_glm(yCornell,XCornell,10,sparse=TRUE,sparseStop=FALSE)
rm(list=c("XCornell","yCornell"))


## User specified links can be used.
## Example of user-specified link, a logit model for p^days
## See Shaffer, T.  2004. Auk 121(2): 526-540 and ?family.
logexp <- function(days = 1)
{
    linkfun <- function(mu) qlogis(mu^(1/days))
    linkinv <- function(eta) plogis(eta)^days
    mu.eta <- function(eta) days * plogis(eta)^(days-1) *
      .Call("logit_mu_eta", eta, PACKAGE = "stats")
    valideta <- function(eta) TRUE
    link <- paste("logexp(", days, ")", sep="")
    structure(list(linkfun = linkfun, linkinv = linkinv,
                   mu.eta = mu.eta, valideta = valideta, name = link),
              class = "link-glm")
}
binomial(logexp(3))

data(aze_compl)
Xaze_compl<-aze_compl[,2:34]
yaze_compl<-aze_compl$y
modpls <- PLS_glm(yaze_compl,Xaze_compl,nt=10,modele="pls-glm-family",family=binomial(link=logexp(3)),MClassed=TRUE,pvals.expli=TRUE)
modpls$InfCrit
rm(list=c("Xaze_compl","yaze_compl","modpls","logexp","aze_compl"))





cleanEx()
nameEx("PLS_glm_formula")
### * PLS_glm_formula

flush(stderr()); flush(stdout())

### Name: PLS_glm_formula
### Title: Partial least squares Regression generalized linear models
### Aliases: PLS_glm_formula
### Keywords: models regression

### ** Examples

data(Cornell)
PLS_glm_formula(Y~.,data=Cornell,3)$uscores
PLS_glm_formula(Y~.,data=Cornell,3)$pp
PLS_glm_formula(Y~.,data=Cornell,3)$Coeffs
PLS_glm_formula(Y~.,data=Cornell,10)$InfCrit
PLS_glm_formula(Y~.,data=Cornell,10,modele="pls-glm-gaussian")$InfCrit
data.frame(pls=PLS_glm_formula(Y~.,data=Cornell,3)$Coeffs,PLS_glm_formula_gaussian=PLS_glm_formula(Y~.,data=Cornell,10,modele="pls-glm-gaussian")$Coeffs,PLS_glm_formula_family_gaussian=PLS_glm_formula(Y~.,data=Cornell,10,modele="pls-glm-family",family=gaussian())$Coeffs)

mod <- PLS_glm_formula(Y~.,data=Cornell,10,pvals.expli =TRUE)
mod2 <- PLS_glm_formula(Y~.,data=Cornell,10,sparse=TRUE)
mod3 <- PLS_glm_formula(Y~.,data=Cornell,10,sparse=TRUE,sparseStop=FALSE)


## User specified links can be used.
## Example of user-specified link, a logit model for p^days
## See Shaffer, T.  2004. Auk 121(2): 526-540 and ?family.
logexp <- function(days = 1)
{
    linkfun <- function(mu) qlogis(mu^(1/days))
    linkinv <- function(eta) plogis(eta)^days
    mu.eta <- function(eta) days * plogis(eta)^(days-1) *
      .Call("logit_mu_eta", eta, PACKAGE = "stats")
    valideta <- function(eta) TRUE
    link <- paste("logexp(", days, ")", sep="")
    structure(list(linkfun = linkfun, linkinv = linkinv,
                   mu.eta = mu.eta, valideta = valideta, name = link),
              class = "link-glm")
}
binomial(logexp(3))

data(aze_compl)
modpls <- PLS_glm_formula(y~.,data=aze_compl,nt=10,modele="pls-glm-family",family=binomial(link=logexp(3)),MClassed=TRUE,pvals.expli=TRUE)
modpls$InfCrit
rm(list=c("modpls","logexp","aze_compl"))





cleanEx()
nameEx("PLS_glm_kfoldcv")
### * PLS_glm_kfoldcv

flush(stderr()); flush(stdout())

### Name: PLS_glm_kfoldcv
### Title: Partial least squares regression glm models with kfold cross
###   validation
### Aliases: PLS_glm_kfoldcv
### Keywords: models regression

### ** Examples

data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
bbb <- PLS_glm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=10,NK=1,modele="pls")
kfolds2CVinfos_glm(bbb)

PLS_glm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=3,modele="pls-glm-gaussian",K=12)
PLS_glm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=3,modele="pls-glm-gaussian",K=6,NK=2,random=TRUE,keepfolds=TRUE)$results_kfolds

#Different ways of model specifications
PLS_glm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=3,modele="pls-glm-gaussian",K=6,NK=2,random=FALSE,keepfolds=TRUE)$results_kfolds
PLS_glm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=3,modele="pls-glm-family",family=gaussian,K=6,NK=2,random=FALSE,keepfolds=TRUE)$results_kfolds
PLS_glm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=3,modele="pls-glm-family",family=gaussian(),K=6,NK=2,random=FALSE,keepfolds=TRUE)$results_kfolds
PLS_glm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=3,modele="pls-glm-family",family=gaussian(link=log),K=6,NK=2,random=FALSE,keepfolds=TRUE)$results_kfolds




cleanEx()
nameEx("PLS_glm_kfoldcv_formula")
### * PLS_glm_kfoldcv_formula

flush(stderr()); flush(stdout())

### Name: PLS_glm_kfoldcv_formula
### Title: Partial least squares regression glm models with kfold cross
###   validation
### Aliases: PLS_glm_kfoldcv_formula
### Keywords: models regression

### ** Examples

data(Cornell)
bbb <- PLS_glm_kfoldcv_formula(Y~.,data=Cornell,nt=10,NK=1,modele="pls")
kfolds2CVinfos_glm(bbb)

PLS_glm_kfoldcv_formula(Y~.,data=Cornell,nt=3,modele="pls-glm-gaussian",K=12)
PLS_glm_kfoldcv_formula(Y~.,data=Cornell,nt=3,modele="pls-glm-gaussian",K=6,NK=2,random=TRUE,keepfolds=TRUE)$results_kfolds

#Different ways of model specifications
PLS_glm_kfoldcv_formula(Y~.,data=Cornell,nt=3,modele="pls-glm-gaussian",K=6,NK=2,random=FALSE,keepfolds=TRUE)$results_kfolds
PLS_glm_kfoldcv_formula(Y~.,data=Cornell,nt=3,modele="pls-glm-family",family=gaussian,K=6,NK=2,random=FALSE,keepfolds=TRUE)$results_kfolds
PLS_glm_kfoldcv_formula(Y~.,data=Cornell,nt=3,modele="pls-glm-family",family=gaussian(),K=6,NK=2,random=FALSE,keepfolds=TRUE)$results_kfolds
PLS_glm_kfoldcv_formula(Y~.,data=Cornell,nt=3,modele="pls-glm-family",family=gaussian(link=log),K=6,NK=2,random=FALSE,keepfolds=TRUE)$results_kfolds




cleanEx()
nameEx("PLS_glm_wvc")
### * PLS_glm_wvc

flush(stderr()); flush(stdout())

### Name: PLS_glm_wvc
### Title: Light version of PLS\_glm for cross validation purposes
### Aliases: PLS_glm_wvc
### Keywords: models regression

### ** Examples

data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
PLS_glm_wvc(dataY=yCornell,dataX=XCornell,nt=3,modele="pls-glm-gaussian",dataPredictY=XCornell[1,])
PLS_glm_wvc(dataY=yCornell,dataX=XCornell,nt=3,modele="pls-glm-family",family=gaussian(),dataPredictY=XCornell[1,])
PLS_glm_wvc(dataY=yCornell[-1],dataX=XCornell[-1,],nt=3,modele="pls-glm-gaussian",dataPredictY=XCornell[1,])
PLS_glm_wvc(dataY=yCornell[-1],dataX=XCornell[-1,],nt=3,modele="pls-glm-family",family=gaussian(),dataPredictY=XCornell[1,])
rm("XCornell","yCornell")




cleanEx()
nameEx("PLS_lm")
### * PLS_lm

flush(stderr()); flush(stdout())

### Name: PLS_lm
### Title: Partial least squares Regression models with leave one out cross
###   validation
### Aliases: PLS_lm
### Keywords: models regression

### ** Examples

data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
PLS_lm(yCornell,XCornell,10)$InfCrit
PLS_lm(yCornell,XCornell,10,typeVC="standard")$CVinfos
PLS_lm(yCornell,XCornell,6)$AIC 
PLS_lm(yCornell,XCornell,6)$AIC.std   

modpls <- PLS_lm(yCornell,XCornell,6,pvals.expli =TRUE)
modpls2 <- PLS_lm(yCornell,XCornell,6,sparse=TRUE)
modpls3 <- PLS_lm(yCornell,XCornell,6,sparse=TRUE,sparseStop=FALSE)
rm(list=c("XCornell","yCornell"))





cleanEx()
nameEx("PLS_lm_formula")
### * PLS_lm_formula

flush(stderr()); flush(stdout())

### Name: PLS_lm_formula
### Title: Partial least squares Regression models with leave one out cross
###   validation
### Aliases: PLS_lm_formula
### Keywords: models regression

### ** Examples

data(Cornell)
PLS_lm_formula(Y~.,data=Cornell,nt=10)$InfCrit
PLS_lm_formula(Y~.,data=Cornell,nt=10,typeVC="standard")$CVinfos
PLS_lm_formula(Y~.,data=Cornell,nt=6)$AIC 
PLS_lm_formula(Y~.,data=Cornell,nt=6)$AIC.std    

modpls <- PLS_lm_formula(Y~.,data=Cornell,nt=6,pvals.expli =TRUE) 
modpls2 <- PLS_lm_formula(Y~.,data=Cornell,nt=6,sparse=TRUE)    
modpls3 <- PLS_lm_formula(Y~.,data=Cornell,nt=6,sparse=TRUE,sparseStop=FALSE)    





cleanEx()
nameEx("PLS_lm_kfoldcv")
### * PLS_lm_kfoldcv

flush(stderr()); flush(stdout())

### Name: PLS_lm_kfoldcv
### Title: Partial least squares regression models with kfold cross
###   validation
### Aliases: PLS_lm_kfoldcv
### Keywords: models regression

### ** Examples

data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
PLS_lm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=3,K=12,keepfolds=TRUE)
PLS_lm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=3,K=12,keepfolds=FALSE)
PLS_lm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=3,K=6,NK=2,random=FALSE,keepfolds=TRUE)
PLS_lm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=3,K=6,NK=2,random=TRUE,keepfolds=TRUE)
PLS_lm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=3,keepcoeffs=TRUE,keepfolds=TRUE)
PLS_lm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=3,keepcoeffs=TRUE,keepfolds=FALSE)

bbb <- PLS_lm_kfoldcv(dataY=yCornell,dataX=data.frame(scale(as.matrix(XCornell))[,]),nt=6,K=12,NK=1)
bbb2 <- PLS_lm_kfoldcv(dataY=yCornell,dataX=data.frame(scale(as.matrix(XCornell))[,]),nt=6,K=6,NK=1)
kfolds2CVinfos_lm(bbb)
kfolds2CVinfos_lm(bbb2)
PLS_lm(yCornell,XCornell,6,typeVC="standard")$CVinfos
rm(list=c("XCornell","yCornell","bbb","bbb2"))





cleanEx()
nameEx("PLS_lm_kfoldcv_formula")
### * PLS_lm_kfoldcv_formula

flush(stderr()); flush(stdout())

### Name: PLS_lm_kfoldcv_formula
### Title: Partial least squares regression models with kfold cross
###   validation
### Aliases: PLS_lm_kfoldcv_formula
### Keywords: models regression

### ** Examples

data(Cornell)
PLS_lm_kfoldcv_formula(Y~.,Cornell,nt=3,K=12,keepfolds=TRUE)
PLS_lm_kfoldcv_formula(Y~.,Cornell,nt=3,K=12,keepfolds=FALSE)
PLS_lm_kfoldcv_formula(Y~.,Cornell,nt=3,K=6,NK=2,random=FALSE,keepfolds=TRUE)
PLS_lm_kfoldcv_formula(Y~.,Cornell,nt=3,K=6,NK=2,random=TRUE,keepfolds=TRUE)
PLS_lm_kfoldcv_formula(Y~.,Cornell,nt=3,keepcoeffs=TRUE,keepfolds=TRUE)
PLS_lm_kfoldcv_formula(Y~.,Cornell,nt=3,keepcoeffs=TRUE,keepfolds=FALSE)

bbb <- PLS_lm_kfoldcv_formula(Y~scale(as.matrix(Cornell))[,-8],Cornell,nt=6,K=12,NK=1)
bbb2 <- PLS_lm_kfoldcv_formula(Y~scale(as.matrix(Cornell))[,-8],Cornell,nt=6,K=6,NK=1)
kfolds2CVinfos_lm(bbb)
kfolds2CVinfos_lm(bbb2)
PLS_lm_formula(Y~.,Cornell,6,typeVC="standard")$CVinfos
rm(list=c("bbb","bbb2"))





cleanEx()
nameEx("PLS_lm_wvc")
### * PLS_lm_wvc

flush(stderr()); flush(stdout())

### Name: PLS_lm_wvc
### Title: Light version of PLS\_lm for cross validation purposes
### Aliases: PLS_lm_wvc
### Keywords: models regression

### ** Examples

data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
PLS_lm_wvc(dataY=yCornell,dataX=XCornell,nt=3,dataPredictY=XCornell[1,])
PLS_lm_wvc(dataY=yCornell[-c(1,2)],dataX=XCornell[-c(1,2),],nt=3,dataPredictY=XCornell[c(1,2),])
PLS_lm_wvc(dataY=yCornell[-c(1,2)],dataX=XCornell[-c(1,2),],nt=3,dataPredictY=XCornell[c(1,2),],keepcoeffs=TRUE)
rm("XCornell","yCornell")

## With an incomplete dataset (X[1,2] is NA)
data(pine)
ypine <- pine[,11]
data(XpineNAX21)
PLS_lm_wvc(dataY=log(ypine)[-1],dataX=XpineNAX21[-1,],nt=3)
PLS_lm_wvc(dataY=log(ypine)[-1],dataX=XpineNAX21[-1,],nt=3,dataPredictY=XpineNAX21[1,])
PLS_lm_wvc(dataY=log(ypine)[-2],dataX=XpineNAX21[-2,],nt=3,dataPredictY=XpineNAX21[2,])
PLS_lm_wvc(dataY=log(ypine),dataX=XpineNAX21,nt=3)
rm("XpineNAX21","ypine")



cleanEx()
nameEx("XbordeauxNA")
### * XbordeauxNA

flush(stderr()); flush(stdout())

### Name: XbordeauxNA
### Title: Missing data analysis for the quality of wine dataset
### Aliases: XbordeauxNA
### Keywords: datasets

### ** Examples

data(XbordeauxNA)
## maybe str(XbordeauxNA) ; plot(XbordeauxNA) ...



cleanEx()
nameEx("XpineNAX21")
### * XpineNAX21

flush(stderr()); flush(stdout())

### Name: XpineNAX21
### Title: Missing data analysis for the pine dataset
### Aliases: XpineNAX21
### Keywords: datasets

### ** Examples

data(XpineNAX21)
## maybe str(XpineNAX21) ; plot(XpineNAX21) ...



cleanEx()
nameEx("aic.dof")
### * aic.dof

flush(stderr()); flush(stdout())

### Name: aic.dof
### Title: Akaike and Bayesian Information Criteria and Generalized minimum
###   description length
### Aliases: aic.dof bic.dof gmdl.dof
### Keywords: models regression utilities

### ** Examples

data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
modpls <- plsR(yCornell,XCornell,4)
dof.object <- plsR.dof(modpls)
aic.dof(modpls$RSS,modpls$nr,dof.object$DoF,dof.object$sigmahat)
bic.dof(modpls$RSS,modpls$nr,dof.object$DoF,dof.object$sigmahat)
gmdl.dof(dof.object$sigmahat,modpls$nr,dof.object$DoF,dof.object$yhat)
naive.object <- plsR.dof(modpls,naive=TRUE)
aic.dof(modpls$RSS,modpls$nr,naive.object$DoF,naive.object$sigmahat)
bic.dof(modpls$RSS,modpls$nr,naive.object$DoF,naive.object$sigmahat)
gmdl.dof(naive.object$sigmahat,modpls$nr,naive.object$DoF,naive.object$yhat)



cleanEx()
nameEx("aze")
### * aze

flush(stderr()); flush(stdout())

### Name: aze
### Title: Microsat Dataset
### Aliases: aze
### Keywords: datasets

### ** Examples

data(aze)
## maybe str(aze) ; plot(aze) ...



cleanEx()
nameEx("aze_compl")
### * aze_compl

flush(stderr()); flush(stdout())

### Name: aze_compl
### Title: As aze without missing values
### Aliases: aze_compl
### Keywords: datasets

### ** Examples

data(aze_compl)
## maybe str(aze_compl) ; plot(aze_compl) ...



cleanEx()
nameEx("bootpls")
### * bootpls

flush(stderr()); flush(stdout())

### Name: bootpls
### Title: Non-parametric Bootstrap for PLS models
### Aliases: bootpls
### Keywords: models

### ** Examples

data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]

# Lazraq-Cléroux PLS ordinary bootstrap

set.seed(250)
Cornell.boot <- bootpls(plsR(yCornell,XCornell,3), sim="ordinary", stype="i", R=250)
boot::boot.array(Cornell.boot, indices=TRUE)

# Graph similar to the one of Bastien et al. in CSDA 2005
boxplot(as.vector(Cornell.boot$t[,-1])~factor(rep(1:7,rep(250,7))), main="Bootstrap distributions of standardised bj (j = 1, ..., 7).")
points(c(1:7),Cornell.boot$t0[-1],col="red",pch=19)
# Using the boxplots.bootpls function
boxplots.bootpls(Cornell.boot,indices=2:8)
# Confidence intervals plotting
confints.bootpls(Cornell.boot,indices=2:8)
plots.confints.bootpls(confints.bootpls(Cornell.boot,indices=2:8))





cleanEx()
nameEx("bootplsglm")
### * bootplsglm

flush(stderr()); flush(stdout())

### Name: bootplsglm
### Title: Non-parametric Bootstrap for PLS generalized linear models
### Aliases: bootplsglm
### Keywords: models

### ** Examples



options(contrasts = c(unordered = "contr.treatment",ordered = "contr.poly"))
cleanEx()
nameEx("bordeaux")
### * bordeaux

flush(stderr()); flush(stdout())

### Name: bordeaux
### Title: Quality of wine dataset
### Aliases: bordeaux
### Keywords: datasets

### ** Examples

data(bordeaux)
## maybe str(bordeaux) ; plot(bordeaux) ...



cleanEx()
nameEx("boxplots.bootpls")
### * boxplots.bootpls

flush(stderr()); flush(stdout())

### Name: boxplots.bootpls
### Title: Boxplot bootstrap distributions
### Aliases: boxplots.bootpls
### Keywords: regression models

### ** Examples

data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]

# Lazraq-Cléroux PLS ordinary bootstrap

set.seed(250)
Cornell.boot <- bootpls(plsR(yCornell,XCornell,3), sim="ordinary", stype="i", R=250)

# Graph similar to the one of Bastien et al. in CSDA 2005
boxplots.bootpls(Cornell.boot,indices=2:8)




cleanEx()
nameEx("coefs.plsR")
### * coefs.plsR

flush(stderr()); flush(stdout())

### Name: coefs.plsR
### Title: Coefficients for bootstrap computations
### Aliases: coefs.plsR
### Keywords: models

### ** Examples

data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]



cleanEx()
nameEx("coefs.plsRglm")
### * coefs.plsRglm

flush(stderr()); flush(stdout())

### Name: coefs.plsRglm
### Title: Coefficients for bootstrap computations
### Aliases: coefs.plsRglm
### Keywords: models

### ** Examples

data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]



cleanEx()
nameEx("confints.bootpls")
### * confints.bootpls

flush(stderr()); flush(stdout())

### Name: confints.bootpls
### Title: Bootstrap confidence intervals
### Aliases: confints.bootpls
### Keywords: regression models

### ** Examples

data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]

# Lazraq-Cléroux PLS ordinary bootstrap

set.seed(250)
Cornell.boot <- bootpls(plsR(yCornell,XCornell,3), sim="ordinary", stype="i", R=250)

(temp.ci <- confints.bootpls(Cornell.boot,2:8))
plots.confints.bootpls(temp.ci)
(temp.ci <- confints.bootpls(Cornell.boot,2:8,typeBCa=FALSE))
plots.confints.bootpls(temp.ci)
(temp.ci <- confints.bootpls(Cornell.boot,c(2,4,6)))
plots.confints.bootpls(temp.ci)




cleanEx()
nameEx("dicho")
### * dicho

flush(stderr()); flush(stdout())

### Name: dicho
### Title: Dichotomization
### Aliases: dicho
### Keywords: utilities

### ** Examples

dimX <- 6
Astar <- 4
(dataAstar4 <- t(replicate(10,simul_data_YX(dimX,Astar))))

dicho(dataAstar4)

rm(list=c("dimX","Astar"))



cleanEx()
nameEx("fowlkes")
### * fowlkes

flush(stderr()); flush(stdout())

### Name: fowlkes
### Title: Fowlkes dataset
### Aliases: fowlkes
### Keywords: datasets

### ** Examples

data(fowlkes)
## maybe str(fowlkes) ; plot(fowlkes) ...



cleanEx()
nameEx("infcrit.dof")
### * infcrit.dof

flush(stderr()); flush(stdout())

### Name: infcrit.dof
### Title: Information criteria
### Aliases: infcrit.dof
### Keywords: models regression utilities

### ** Examples

data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
modpls <- plsR(yCornell,XCornell,4)
infcrit.dof(modpls)



cleanEx()
nameEx("kfolds2CVinfos_glm")
### * kfolds2CVinfos_glm

flush(stderr()); flush(stdout())

### Name: kfolds2CVinfos_glm
### Title: Extracts and computes information criteria and fits statistics
###   for kfold cross validated partial least squares glm models
### Aliases: kfolds2CVinfos_glm
### Keywords: models regression

### ** Examples



options(contrasts = c(unordered = "contr.treatment",ordered = "contr.poly"))
cleanEx()
nameEx("kfolds2CVinfos_lm")
### * kfolds2CVinfos_lm

flush(stderr()); flush(stdout())

### Name: kfolds2CVinfos_lm
### Title: Extracts and computes information criteria and fits statistics
###   for kfold cross validated partial least squares models
### Aliases: kfolds2CVinfos_lm
### Keywords: internal models regression

### ** Examples



cleanEx()
nameEx("kfolds2Chisq")
### * kfolds2Chisq

flush(stderr()); flush(stdout())

### Name: kfolds2Chisq
### Title: Computes Predicted Chisquare for kfold cross validated partial
###   least squares regression models.
### Aliases: kfolds2Chisq
### Keywords: models regression

### ** Examples



cleanEx()
nameEx("kfolds2Chisqind")
### * kfolds2Chisqind

flush(stderr()); flush(stdout())

### Name: kfolds2Chisqind
### Title: Computes individual Predicted Chisquare for kfold cross
###   validated partial least squares regression models.
### Aliases: kfolds2Chisqind
### Keywords: models regression

### ** Examples



cleanEx()
nameEx("kfolds2Mclassed")
### * kfolds2Mclassed

flush(stderr()); flush(stdout())

### Name: kfolds2Mclassed
### Title: Number of missclassified individuals for kfold cross validated
###   partial least squares regression models.
### Aliases: kfolds2Mclassed
### Keywords: models regression

### ** Examples




cleanEx()
nameEx("kfolds2Mclassedind")
### * kfolds2Mclassedind

flush(stderr()); flush(stdout())

### Name: kfolds2Mclassedind
### Title: Number of missclassified individuals per group for kfold cross
###   validated partial least squares regression models.
### Aliases: kfolds2Mclassedind
### Keywords: models regression

### ** Examples




cleanEx()
nameEx("kfolds2Press")
### * kfolds2Press

flush(stderr()); flush(stdout())

### Name: kfolds2Press
### Title: Computes PRESS for kfold cross validated partial least squares
###   regression models.
### Aliases: kfolds2Press
### Keywords: models regression

### ** Examples

data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
bbb <- PLS_lm_kfoldcv(dataY=yCornell,dataX=data.frame(scale(as.matrix(XCornell))[,]),nt=6,K=12,NK=1)
bbb2 <- PLS_lm_kfoldcv(dataY=yCornell,dataX=data.frame(scale(as.matrix(XCornell))[,]),nt=6,K=6,NK=1)
kfolds2Press(bbb)
kfolds2Press(bbb2)
rm(list=c("XCornell","yCornell","bbb","bbb2"))





cleanEx()
nameEx("kfolds2Pressind")
### * kfolds2Pressind

flush(stderr()); flush(stdout())

### Name: kfolds2Pressind
### Title: Computes individual PRESS for kfold cross validated partial
###   least squares regression models.
### Aliases: kfolds2Pressind
### Keywords: models regression

### ** Examples

data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
bbb <- PLS_lm_kfoldcv(dataY=yCornell,dataX=data.frame(scale(as.matrix(XCornell))[,]),nt=6,K=12,NK=1)
bbb2 <- PLS_lm_kfoldcv(dataY=yCornell,dataX=data.frame(scale(as.matrix(XCornell))[,]),nt=6,K=6,NK=1)
kfolds2Pressind(bbb)
kfolds2Pressind(bbb2)
rm(list=c("XCornell","yCornell","bbb","bbb2"))




cleanEx()
nameEx("kfolds2coeff")
### * kfolds2coeff

flush(stderr()); flush(stdout())

### Name: kfolds2coeff
### Title: Extracts coefficients from kfold cross validated partial least
###   squares regression models
### Aliases: kfolds2coeff
### Keywords: models regression

### ** Examples

data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
data(aze_compl)
Xaze_compl<-aze_compl[,2:34]
yaze_compl<-aze_compl$y
data(pine)
Xpine<-pine[,1:10]
ypine<-pine[,11]
bbb <- PLS_lm_kfoldcv(dataY=yCornell,dataX=XCornell,nt=3,keepcoeffs=TRUE)
bbb2 <- PLS_lm_kfoldcv(dataY=log(ypine),dataX=Xpine,nt=4,keepcoeffs=TRUE)
kfolds2coeff(bbb)
boxplot(kfolds2coeff(bbb)[,1])
kfolds2coeff(bbb2)
boxplot(kfolds2coeff(bbb2)[,2])
rm(list=c("XCornell","yCornell","Xpine","ypine","bbb","bbb2"))



cleanEx()
nameEx("loglikpls")
### * loglikpls

flush(stderr()); flush(stdout())

### Name: loglikpls
### Title: loglikelihood function for plsR models
### Aliases: loglikpls
### Keywords: models regression utilities

### ** Examples

data(pine)
ypine <- pine[,11]
Xpine <- pine[,1:10]
(Pinscaled <- as.data.frame(cbind(scale(log(ypine)),scale(as.matrix(Xpine)))))
colnames(Pinscaled)[1] <- "yy"

lm(yy~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10,data=Pinscaled)

modpls <- plsR(log(ypine),Xpine,10)
modpls$Std.Coeffs
lm(yy~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10,data=Pinscaled)

AIC(lm(yy~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10,data=Pinscaled))
print(logLik(lm(yy~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10,data=Pinscaled)))

sum(dnorm(modpls$RepY, modpls$Std.ValsPredictY, sqrt(mean(modpls$residY^2)), log=TRUE))
sum(dnorm(Pinscaled$yy,fitted(lm(yy~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10,data=Pinscaled)),sqrt(mean(residuals(lm(yy~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10,data=Pinscaled))^2)), log=TRUE))
loglikpls(modpls$residY)
loglikpls(residuals(lm(yy~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10,data=Pinscaled)))
AICpls(10,residuals(lm(yy~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10,data=Pinscaled)))
AICpls(10,modpls$residY)



cleanEx()
nameEx("permcoefs.plsR")
### * permcoefs.plsR

flush(stderr()); flush(stdout())

### Name: permcoefs.plsR
### Title: Coefficients computation for permutation bootstrap
### Aliases: permcoefs.plsR
### Keywords: models

### ** Examples



cleanEx()
nameEx("permcoefs.plsRglm")
### * permcoefs.plsRglm

flush(stderr()); flush(stdout())

### Name: permcoefs.plsRglm
### Title: Coefficients computation for permutation bootstrap
### Aliases: permcoefs.plsRglm
### Keywords: models

### ** Examples



cleanEx()
nameEx("pine")
### * pine

flush(stderr()); flush(stdout())

### Name: pine
### Title: Pine dataset
### Aliases: pine
### Keywords: datasets

### ** Examples

data(pine)
## maybe str(pine) ; plot(pine) ...



cleanEx()
nameEx("pine_full")
### * pine_full

flush(stderr()); flush(stdout())

### Name: pine_full
### Title: Full pine dataset
### Aliases: pine_full
### Keywords: datasets

### ** Examples

data(pine_full)
## maybe str(pine_full) ; plot(pine_full) ...



cleanEx()
nameEx("pine_sup")
### * pine_sup

flush(stderr()); flush(stdout())

### Name: pine_sup
### Title: Supplementary data for pine dataset
### Aliases: pine_sup
### Keywords: datasets

### ** Examples

data(pine_sup)
## maybe str(pine_sup) ; plot(pine_sup) ...



cleanEx()
nameEx("plots.confints.bootpls")
### * plots.confints.bootpls

flush(stderr()); flush(stdout())

### Name: plots.confints.bootpls
### Title: Plot bootstrap confidence intervals
### Aliases: plots.confints.bootpls
### Keywords: regression models

### ** Examples

data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]

# Lazraq-Cléroux PLS ordinary bootstrap

set.seed(250)
Cornell.boot <- bootpls(plsR(yCornell,XCornell,3), sim="ordinary", stype="i", R=250)
temp.ci <- confints.bootpls(Cornell.boot,2:8)

plots.confints.bootpls(temp.ci)
plots.confints.bootpls(temp.ci,prednames=FALSE)
plots.confints.bootpls(temp.ci,prednames=FALSE,articlestyle=FALSE,main="Bootstrap confidence intervals for the bj")
plots.confints.bootpls(temp.ci,indices=1:3,prednames=FALSE)
plots.confints.bootpls(temp.ci,c(2,4,6),"bottomright")
plots.confints.bootpls(temp.ci,c(2,4,6),articlestyle=FALSE,main="Bootstrap confidence intervals for some of the bj")

temp.ci <- confints.bootpls(Cornell.boot,typeBCa=FALSE)
plots.confints.bootpls(temp.ci,prednames=TRUE)
plots.confints.bootpls(temp.ci,prednames=FALSE)





cleanEx()
nameEx("plsR")
### * plsR

flush(stderr()); flush(stdout())

### Name: plsR
### Title: Partial least squares Regression models with leave one out cross
###   validation
### Aliases: plsR plsRmodel.default plsRmodel.formula
### Keywords: models regression

### ** Examples

data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
plsR(yCornell,XCornell,10)$InfCrit
plsR(yCornell,XCornell,10,typeVC="standard")$CVinfos
plsR(yCornell,XCornell,6)$AIC 
plsR(yCornell,XCornell,6)$AIC.std    
rm(list=c("XCornell","yCornell"))





cleanEx()
nameEx("plsR.dof")
### * plsR.dof

flush(stderr()); flush(stdout())

### Name: plsR.dof
### Title: Computation of the Degrees of Freedom
### Aliases: plsR.dof
### Keywords: models regression utilities

### ** Examples

data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
modpls <- plsR(yCornell,XCornell,4)
plsR.dof(modpls) 
plsR.dof(modpls,naive=TRUE) 



cleanEx()
nameEx("plsRglm-package")
### * plsRglm-package

flush(stderr()); flush(stdout())

### Name: plsRglm-package
### Title: Partial least squares Regression for generalized linear models
### Aliases: plsRglm-package
### Keywords: package

### ** Examples

data(pine)
ypine <- pine[,11]
Xpine <- pine[,1:10]
(Pinscaled <- as.data.frame(cbind(scale(log(ypine)),scale(as.matrix(Xpine)))))
colnames(Pinscaled)[1] <- "yy"

modpls <- plsR(log(ypine),Xpine,10)
modpls$Std.Coeffs



cleanEx()
nameEx("plsRglm")
### * plsRglm

flush(stderr()); flush(stdout())

### Name: plsRglm
### Title: Partial least squares Regression generalized linear models
### Aliases: plsRglm plsRglmmodel.default plsRglmmodel.formula
### Keywords: models regression

### ** Examples

data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
plsRglm(yCornell,XCornell,3)$uscores
plsRglm(yCornell,XCornell,3)$pp
plsRglm(yCornell,XCornell,3)$Coeffs
plsRglm(yCornell,XCornell,10)$InfCrit
plsRglm(yCornell,XCornell,10,modele="pls-glm-gaussian")$InfCrit
rm(list=c("XCornell","yCornell"))


data(aze_compl)
Xaze_compl<-aze_compl[,2:34]
yaze_compl<-aze_compl$y
plsRglm(yaze_compl,Xaze_compl,nt=10,modele="pls",MClassed=TRUE)$InfCrit
modpls <- plsRglm(yaze_compl,Xaze_compl,nt=10,modele="pls-glm-logistic",MClassed=TRUE,pvals.expli=TRUE)
modpls$InfCrit
modpls$valpvalstep
modpls$Coeffsmodel_vals

plot(plsRglm(yaze_compl,Xaze_compl,4,modele="pls-glm-logistic")$FinalModel)
plsRglm(yaze_compl[-c(99,72)],Xaze_compl[-c(99,72),],4,modele="pls-glm-logistic",pvals.expli=TRUE)$pvalstep
plot(plsRglm(yaze_compl[-c(99,72)],Xaze_compl[-c(99,72),],4,modele="pls-glm-logistic",pvals.expli=TRUE)$FinalModel)
rm(list=c("Xaze_compl","yaze_compl","modpls"))


data(bordeaux)
Xbordeaux<-bordeaux[,1:4]
ybordeaux<-factor(bordeaux$Quality,ordered=TRUE)
modpls <- plsRglm(ybordeaux,Xbordeaux,10,modele="pls-glm-polr")
modpls$Coeffsmodel_vals
modpls$InfCrit

XbordeauxNA<-Xbordeaux
XbordeauxNA[1,1] <- NA
modplsNA <- plsRglm(ybordeaux,XbordeauxNA,10,modele="pls-glm-polr")
modplsNA$Coeffsmodel_vals
modplsNA$InfCrit
rm(list=c("Xbordeaux","XbordeauxNA","ybordeaux","modplsNA"))





options(contrasts = c(unordered = "contr.treatment",ordered = "contr.poly"))
cleanEx()
nameEx("print.plsRglmmodel")
### * print.plsRglmmodel

flush(stderr()); flush(stdout())

### Name: print.plsRglmmodel
### Title: Print method for plsRglm models
### Aliases: print.plsRglmmodel
### Keywords: methods print

### ** Examples

data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
modplsglm <- plsRglm(yCornell,XCornell,3,modele="pls-glm-gaussian")
class(modplsglm)
print(modplsglm)
rm(list=c("XCornell","yCornell","modplsglm"))



cleanEx()
nameEx("print.plsRmodel")
### * print.plsRmodel

flush(stderr()); flush(stdout())

### Name: print.plsRmodel
### Title: Print method for plsR models
### Aliases: print.plsRmodel
### Keywords: methods print

### ** Examples

data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
modpls <- plsRglm(yCornell,XCornell,3,modele="pls")
class(modpls)
print(modpls)
rm(list=c("XCornell","yCornell","modpls"))



cleanEx()
nameEx("print.summary.plsRglmmodel")
### * print.summary.plsRglmmodel

flush(stderr()); flush(stdout())

### Name: print.summary.plsRglmmodel
### Title: Print method for summaries of plsRglm models
### Aliases: print.summary.plsRglmmodel
### Keywords: methods print

### ** Examples

data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
modplsglm <- plsRglm(yCornell,XCornell,3,modele="pls-glm-gaussian")
class(modplsglm)
print(summary(modplsglm))
rm(list=c("XCornell","yCornell","modplsglm"))



cleanEx()
nameEx("print.summary.plsRmodel")
### * print.summary.plsRmodel

flush(stderr()); flush(stdout())

### Name: print.summary.plsRmodel
### Title: Print method for summaries of plsR models
### Aliases: print.summary.plsRmodel
### Keywords: methods print

### ** Examples

data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
modpls <- plsRglm(yCornell,XCornell,3,modele="pls")
class(modpls)
print(summary(modpls))
rm(list=c("XCornell","yCornell","modpls"))



cleanEx()
nameEx("simul_data_UniYX")
### * simul_data_UniYX

flush(stderr()); flush(stdout())

### Name: simul_data_UniYX
### Title: Data generating function for univariate plsR models
### Aliases: simul_data_UniYX
### Keywords: datagen utilities

### ** Examples

simul_data_UniYX(20,6)                          




cleanEx()
nameEx("simul_data_UniYX_binom")
### * simul_data_UniYX_binom

flush(stderr()); flush(stdout())

### Name: simul_data_UniYX_binom
### Title: Data generating function for univariate binomial plsR models
### Aliases: simul_data_UniYX_binom
### Keywords: datagen utilities

### ** Examples




cleanEx()
nameEx("simul_data_YX")
### * simul_data_YX

flush(stderr()); flush(stdout())

### Name: simul_data_YX
### Title: Data generating function for multivariate plsR models
### Aliases: simul_data_YX
### Keywords: datagen utilities

### ** Examples

simul_data_YX(20,6)                          




cleanEx()
nameEx("simul_data_complete")
### * simul_data_complete

flush(stderr()); flush(stdout())

### Name: simul_data_complete
### Title: Data generating detailed process for multivariate plsR models
### Aliases: simul_data_complete
### Keywords: datagen utilities

### ** Examples

simul_data_complete(20,6)                          

dimX <- 6
Astar <- 2
simul_data_complete(dimX,Astar)


dimX <- 6
Astar <- 3
simul_data_complete(dimX,Astar)


dimX <- 6
Astar <- 4
simul_data_complete(dimX,Astar)

rm(list=c("dimX","Astar"))



cleanEx()
nameEx("summary.plsRglmmodel")
### * summary.plsRglmmodel

flush(stderr()); flush(stdout())

### Name: summary.plsRglmmodel
### Title: Summary method for plsRglm models
### Aliases: summary.plsRglmmodel
### Keywords: methods print

### ** Examples

data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
modplsglm <- plsRglm(yCornell,XCornell,3,modele="pls-glm-gaussian")
class(modplsglm)
summary(modplsglm)
rm(list=c("XCornell","yCornell","modplsglm"))



cleanEx()
nameEx("summary.plsRmodel")
### * summary.plsRmodel

flush(stderr()); flush(stdout())

### Name: summary.plsRmodel
### Title: Summary method for plsR models
### Aliases: summary.plsRmodel
### Keywords: methods print

### ** Examples

data(Cornell)
XCornell<-Cornell[,1:7]
yCornell<-Cornell[,8]
modpls <- plsR(yCornell,XCornell,4)
class(modpls)
summary(modpls)
rm(list=c("XCornell","yCornell","modpls"))



cleanEx()
nameEx("tilt.bootpls")
### * tilt.bootpls

flush(stderr()); flush(stdout())

### Name: tilt.bootpls
### Title: Tilted bootstrap for PLS models
### Aliases: tilt.bootpls
### Keywords: models

### ** Examples




cleanEx()
nameEx("tilt.bootplsglm")
### * tilt.bootplsglm

flush(stderr()); flush(stdout())

### Name: tilt.bootplsglm
### Title: Tilted bootstrap for PLS models
### Aliases: tilt.bootplsglm
### Keywords: models

### ** Examples



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
