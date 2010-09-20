kfolds2Chisqind <- function(pls_kfolds) {
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


    if (length(pls_kfolds$results_kfolds)==1) {preChisqind_kfolds <- list(vector("list", length(pls_kfolds$results_kfolds[[1]])))}
    else
    {
      if (length(pls_kfolds$results_kfolds)>1)
      {
      preChisqind_kfolds <-vector("list",length(pls_kfolds$results_kfolds))
        for (jj in 1:length(pls_kfolds$results_kfolds))
        {
          preChisqind_kfolds[[jj]] <-vector("list",length(pls_kfolds$results_kfolds[[jj]]))
        }
      rm(jj)
      }
    }

    for (nnkk in 1:length(pls_kfolds$results_kfolds))
    {
        for (ii in 1:length(pls_kfolds$results_kfolds[[1]]))
        {
            if (pls_kfolds$call$modele=="pls-glm-polr") {
                    fff <- ~pls_kfolds$dataY_kfolds[[nnkk]][[ii]]-1
                    m <- model.frame(fff, pls_kfolds$dataY_kfolds[[nnkk]][[ii]])
                    mat <- model.matrix(fff, model.frame(fff, pls_kfolds$dataY_kfolds[[nnkk]][[ii]]))
                    preChisqind_kfolds[[nnkk]][[ii]] <- (unlist(lapply(lapply(pls_kfolds$results_kfolds[[nnkk]][[ii]],function(xxx) {as.list(as.data.frame(t(xxx)))}),
                    Chiscompmatrix,as.list(as.data.frame(t(mat))))))
                    rm(fff)
                    rm(m)
                    rm(mat)
            }
            else {
                    if (dim(pls_kfolds$results_kfolds[[nnkk]][[ii]])[1]==1)
                    {
                        preChisqind_kfolds[[nnkk]][[ii]] <- (pls_kfolds$dataY_kfolds[[nnkk]][[ii]]-pls_kfolds$results_kfolds[[nnkk]][[ii]])^2/(fam_var(pls_kfolds$results_kfolds[[nnkk]][[ii]]))
                    }                
                    else {
                        preChisqind_kfolds[[nnkk]][[ii]] <- colSums((apply(pls_kfolds$results_kfolds[[nnkk]][[ii]],2,'-',pls_kfolds$dataY_kfolds[[nnkk]][[ii]]))^2/(fam_var(pls_kfolds$results_kfolds[[nnkk]][[ii]]))) 
                    }
            }
        }
    }
rm(ii)
rm(nnkk)
return(preChisqind_kfolds)
}
