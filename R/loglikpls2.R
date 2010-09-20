loglikpls2 <- function(dataYpls,YChapeaupls,residpls) {
return(sum(dnorm(dataYpls, YChapeaupls, sqrt(mean(residpls^2)), log=T)))}
