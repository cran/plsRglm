boxplots.bootpls <- function(bootobject,indices=NULL,prednames=TRUE,articlestyle=TRUE,...){
nr <- length(bootobject$t0)
nboot <- dim(bootobject$t)[1]
if(is.null(indices)){indices <- 1:nr}
plotpos <- (1:nr)[1:length(indices)]
if(articlestyle){oldpar <- par();par(mar = c(2, 2, 1, 1) + 0.1, mgp=c(2,1,0))}
boxplot(as.vector(bootobject$t[,indices])~factor(rep(1:length(indices),rep(nboot,length(indices)))),ylim=c(max(-5,min(as.vector(bootobject$t[,indices]))),min(5,max(as.vector(bootobject$t[,indices])))),xaxt="n",...)
#if(prednames){axis(1, at = plotpos+.225, labels = rownames(bootobject$t0)[indices])} else {axis(1, at = plotpos+.225, labels = paste("x",(1:nr)[indices],sep=""))} 
if(prednames){axis(1, at = plotpos, labels = rownames(bootobject$t0)[indices])} else {axis(1, at = plotpos, labels = paste("x",(1:nr)[indices],sep=""))} 
abline(h=0,lty=2,col="blue",lwd=2)
points(plotpos,bootobject$t0[indices],col="red",pch=19)
if(articlestyle){par(oldpar)}
}
