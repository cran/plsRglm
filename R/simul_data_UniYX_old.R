simul_data_UniYX_old <- function(totdim,ncomp) {
varsR <- c(10,5,2,1/2)
varepsilon <- .01
varsF <- c(.25,.125,.05,.0125)

ksi1 <- c(1,1,1,1,1,1)/sqrt(6)
ksi2 <- c(1/2,1/2,-1,1/2,1/2,-1)/sqrt(3)
ksi3 <- c(1,1,1,-1,-1,-1)/sqrt(6)
ksi4 <- c(-1,1,0,-1,1,0)/2
ksi <- cbind(ksi1,ksi2,ksi3,ksi4)[1:totdim,1:ncomp]

epsilon <- rnorm(totdim,mean=rep(0,totdim),sd=varepsilon)

r <- rnorm(ncomp,mean=rep(0,ncomp),sd=varsR[1:ncomp])

simX <- r%*%t(ksi)+epsilon


if(ncomp==2) {
HH <- 3
eta21 <- c(1,1,1)/sqrt(3)
eta22 <- c(1,1,1)/sqrt(3)
eta <- cbind(eta21,eta22)
}

if(ncomp==3) {
HH <- 3
eta31 <- c(1,1,1)/sqrt(3)
eta32 <- c(1,1,1)/sqrt(3)
eta33 <- c(1,1,1)/sqrt(3)
eta <- cbind(eta31,eta32,eta33)
}

if(ncomp==4) {
HH <- 4
eta41 <- c(1,1,1,1)/2
eta42 <- c(1,1,1,1)/2
eta43 <- c(1,1,1,1)/2
eta44 <- c(1,1,1,1)/2
eta <- cbind(eta41,eta42,eta43,eta44)
}

f <- rnorm(ncomp,mean=rep(0,ncomp),sd=varsF[1:ncomp])
z <- f+r


sigmaScarre <- .001
lambda <- .6

sigmaPsi <- sigmaScarre*((1-lambda)*diag(rep(1,HH))+lambda*rep(1,HH)%*%t(rep(1,HH)))
Psi <- rmvnorm(1,mean=rep(0,HH),sigma=sigmaPsi)


Y <- z%*%t(eta)+Psi

res <- c(Y[1],simX)
names(res) <- c("Y",paste("X",1:totdim,sep=""))

return(res)
}
