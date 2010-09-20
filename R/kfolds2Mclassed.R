kfolds2Mclassed <- function(pls_kfolds) {
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
        }
    }
    cat(paste("glm-family : ",as.character(pls_kfolds$call$family$family)))
    cat("\n")

if (!(pls_kfolds$call$modele=="pls-glm-polr")) {
    max_nt <- rep(NA,length(pls_kfolds$results_kfolds))
    if (length(pls_kfolds$results_kfolds)==1) {
      max_nt[1] <- min(unlist(lapply(pls_kfolds$results_kfolds[[1]],ncol)))
      Mclassed_kfolds <- list(rep(0, max_nt[1]))
    }
    else
    {
      if (length(pls_kfolds$results_kfolds)>1)
      {
      Mclassed_kfolds <-vector("list",length(pls_kfolds$results_kfolds))
        for (jj in 1:length(pls_kfolds$results_kfolds))
        {
          max_nt[jj] <- min(unlist(lapply(pls_kfolds$results_kfolds[[jj]],ncol)))
          Mclassed_kfolds[[jj]] <- rep(0,max_nt[jj])
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
                Mclassed_kfolds[[nnkk]] <- Mclassed_kfolds[[nnkk]]+as.numeric(ifelse(pls_kfolds$results_kfolds[[nnkk]][[ii]][1:max_nt[nnkk]] < 0.5, 0, 1) != unclass(pls_kfolds$dataY_kfolds[[nnkk]][[ii]]))
            }
            else
            {
                Mclassed_kfolds[[nnkk]] <- Mclassed_kfolds[[nnkk]]+colSums(apply(ifelse(pls_kfolds$results_kfolds[[nnkk]][[ii]][,1:max_nt[nnkk]] < 0.5, 0, 1),2,'!=',unclass(pls_kfolds$dataY_kfolds[[nnkk]][[ii]]))) 
            }
        }
    }
rm(ii)
rm(nnkk)
}

if (pls_kfolds$call$modele=="pls-glm-polr") {
    max_nt <- rep(NA,length(pls_kfolds$results_kfolds))
    if (length(pls_kfolds$results_kfolds)==1) {
        max_nt[1] <- min(unlist(lapply(pls_kfolds$results_kfolds[[1]],length)))
        Mclassed_kfolds <- list(rep(0, max_nt[1]))
    }
    else
    {
      if (length(pls_kfolds$results_kfolds)>1)
      {
      Mclassed_kfolds <-vector("list",length(pls_kfolds$results_kfolds))
        for (jj in 1:length(pls_kfolds$results_kfolds))
        {
          max_nt[jj] <- min(unlist(lapply(pls_kfolds$results_kfolds[[jj]],length)))
          Mclassed_kfolds[[jj]] <- rep(0,max_nt[jj])
        }
      rm(jj)
      }
    }

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
            Mclassedind_kfolds[[nnkk]][[ii]] <- unlist(lapply(lapply(lapply(pls_kfolds$results_kfolds[[nnkk]],lapply,apply,1,which.max)[[ii]],"!=",unclass(pls_kfolds$dataY_kfolds[[nnkk]][[ii]])),sum))[1:max_nt[nnkk]]
        }
    }

    for (nnkk in 1:length(pls_kfolds$results_kfolds))
    {
            Mclassed_kfolds[[nnkk]] <- colSums(matrix(unlist(Mclassedind_kfolds[[nnkk]]),ncol=max_nt[nnkk],byrow=TRUE))
    }

rm(ii)
rm(nnkk)
}
return(Mclassed_kfolds)
}
