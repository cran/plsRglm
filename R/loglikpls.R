loglikpls <- function(nn,ncomp,residpls) {
return(-nn/2*(log(2*pi) + 1 - log(nn) + log(crossprod(residpls))))}
