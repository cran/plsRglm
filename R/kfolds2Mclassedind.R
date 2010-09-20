kfolds2Mclassedind <- function(pls_kfolds) {
    if (!is.null(pls_kfolds$call$family)) {
        fam_var <- pls_kfolds$call$family$variance
        cat(paste("glm-family : ",as.character(pls_kfolds$call$family$family)))
        cat("\n")
    } else {
        if (pls_kfolds$call$modele=="pls") {
            fam_var <- function(vals) {return(1)}
            pls_kfolds$call$family$family <- "pls"
        }
        if (pls_kfolds$call$modele=="pls-glm-polr") {
            pls_kfolds$call$family$family <- "pls-glm-polr"
            Varyy <- function(piVaryy) {
            diag(piVaryy[-length(piVaryy)])-piVaryy[-length(piVaryy)]%*%t(piVaryy[-length(piVaryy)])
            }
            Chisqcomp <- function(yichisq,pichisq) {
            t(yichisq[-length(yichisq)]-pichisq[-length(pichisq)])%*%solve(Varyy(pichisq))%*%(yichisq[-length(yichisq)]-pichisq[-length(pichisq)])
            }
            Chiscompmatrix <- function(rowspi,rowsyi) {
            sum(mapply(Chisqcomp,rowsyi,rowspi))
            }
        }
    }
    cat(paste("glm-family : ",as.character(pls_kfolds$call$family$family)))
    cat("\n")

if (!(pls_kfolds$call$modele=="pls-glm-polr")) {
    if (length(pls_kfolds$results_kfolds)==1) {Mclassedind_kfolds <- list(vector("list", length(pls_kfolds$results_kfolds[[1]])))}
    else
    {
      if (length(pls_kfolds$results_kfolds)>1)
      {
      Mclassedind_kfolds <-vector("list",length(pls_kfolds$results_kfolds))
        for (jj in 1:length(pls_kfolds$results_kfolds))
        {
          Mclassedind_kfolds[[jj]] <-vector("list",length(pls_kfolds$results_kfolds[[jj]]))
        }
      rm(jj)
      }
    }
    for (nnkk in 1:length(pls_kfolds$results_kfolds))
    {
        for (ii in 1:length(pls_kfolds$results_kfolds[[1]]))
        {
            if (dim(pls_kfolds$results_kfolds[[nnkk]][[ii]])[1]==1)
            {
                Mclassedind_kfolds[[nnkk]][[ii]] <- as.numeric(ifelse(pls_kfolds$results_kfolds[[nnkk]][[ii]] < 0.5, 0, 1) != unclass(pls_kfolds$dataY_kfolds[[nnkk]][[ii]]))
            }
            else
            {
                Mclassedind_kfolds[[nnkk]][[ii]] <- colSums(apply(ifelse(pls_kfolds$results_kfolds[[nnkk]][[ii]] < 0.5, 0, 1),2,'!=',unclass(pls_kfolds$dataY_kfolds[[nnkk]][[ii]])))
            }
        }
    }
}
if (pls_kfolds$call$modele=="pls-glm-polr") {
    if (length(pls_kfolds$results_kfolds)==1) {Mclassedind_kfolds <- list(vector("list", length(pls_kfolds$results_kfolds[[1]])))}
    else
    {
      if (length(pls_kfolds$results_kfolds)>1)
      {
      Mclassedind_kfolds <-vector("list",length(pls_kfolds$results_kfolds))
        for (jj in 1:length(pls_kfolds$results_kfolds))
        {
          Mclassedind_kfolds[[jj]] <-vector("list",length(pls_kfolds$results_kfolds[[jj]]))
        }
      rm(jj)
      }
    }
    for (nnkk in 1:length(pls_kfolds$results_kfolds))
    {
        for (ii in 1:length(pls_kfolds$results_kfolds[[1]]))
        {
            Mclassedind_kfolds[[nnkk]][[ii]] <- unlist(lapply(lapply(lapply(pls_kfolds$results_kfolds[[nnkk]],lapply,apply,1,which.max)[[ii]],"!=",unclass(pls_kfolds$dataY_kfolds[[nnkk]][[ii]])),sum))
        }
    }
rm(ii)
rm(nnkk)
}
return(Mclassedind_kfolds)
}
