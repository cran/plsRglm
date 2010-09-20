AICpls <- function(nn,ncomp,residpls) {
return(-2*loglikpls(nn,ncomp,residpls)+2*(ncomp+1+1))}
