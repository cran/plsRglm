AICpls2 <- function(ncomp,dataYpls,YChapeaupls,residpls) {
return(-2*loglikpls2(dataYpls,YChapeaupls,residpls)+2*(ncomp+1+1))}
