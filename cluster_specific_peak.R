
getbetaMLE <- function(data, cluster, offset, beta_ini, maxiter, thres, quiet){
  ## data is n by p
  ## cluster is the cluster membership matrix, length n. Take entries 1 to nCluster
  ## offset corresponds to h_i, here it is a vector of length n
  p <- ncol(data)
  x <- getdesign(cluster)
  beta <- beta_ini
  betaDiffs <- c() 
  percConverged <- c()
  converged_flag <- rep(0, p)
  
  converged <- 0
  iter <- 1
  while (!converged & iter<maxiter){
    mu <- offset*exp(x%*%beta) 
    u <- crossprod(x, data-mu)
    updateBeta <- function(r){
      if (converged_flag[r]==0 & !is.na(beta[1,r]) ){
        Jr <- crossprod(x, mu[,r]*x)
        betaDiff <- try(solve(Jr, u[,r]), silent=T)
        if (class(betaDiff)=="try-error"){
          return( rep(NA, length(beta[,r])) )    
        } else {
          return(beta[,r] + betaDiff)  
        }
      } else {
        return(beta[,r])
      }
    }
    
    betaPre <- beta
    beta <- sapply(1:p, updateBeta)
    betaDiffs <- rbind(betaDiffs, colSums(abs(beta - betaPre)) )
    converged_flag <- (betaDiffs[iter,] <= thres) + 0
    converged_flag[which(is.na(converged_flag))] <- 0
    singular_flag <- is.na(beta[1,]) + 0
    percConverged <- c(percConverged, sum(converged_flag[which(singular_flag==0)])/sum(singular_flag==0) )
    
    if (!quiet){
        cat("\r", round(percConverged[length(percConverged)]*100), "%", "converged")  
    }
    if (percConverged[length(percConverged)]==1){
      converged <- 1
    }
    iter <- iter + 1
  }
  return(beta)
}


getSigmaPrior <- function(beta, q=0.05){
  ## get the empirical prior estimate
  ## q is the quantile to match
  
  ## center each column in beta
  betaC <- t(t(beta) - colMeans(beta))
  ## 
  sigmas <- rep(0, nrow(betaC))
  for (i in 1:nrow(betaC)){
    tmp <- betaC[i,]
    tmp <- tmp[which(!is.na(tmp))]
    sigmas[i] <- quantile(abs(tmp), 1-q)/qnorm(1-q/2)
  }
  return(sigmas)
}



getdesign <- function(cluster){
  nCluster <- length(unique(cluster))
  n <- length(cluster)
  x <- matrix(0, nrow=n, ncol=nCluster)
  for (i in 1:nCluster){
    x[which(cluster==i), i] <- 1
  }
  return(x)
}

getContrast <- function(nCluster){
  ## get all the contrasts that is needed to calculate the p-value
  constrasts <- matrix(0, nrow=(nCluster-1)*nCluster, ncol=nCluster)
  constrastsCluster <- rep(1:nCluster, each=nCluster-1)
  for (i in 1:nCluster){
    constrasttmp <- matrix(0, nrow=nCluster-1, ncol=nCluster)
    constrasttmp[, i] <- 1
    count <- 1
    for (j in 1:(nCluster-1)){
      constrasttmp[count, (c(1:nCluster)[-i])[j]] <- -1
      count <- count + 1
    }
    constrasts[(i*(nCluster-1)-(nCluster-1)+1):(i*(nCluster-1)), ] <- constrasttmp
  }
  return(list(constrasts=constrasts, constrastsCluster=constrastsCluster))
}




getbetaMAP <- function(data, cluster, offset, sigmas, beta_ini, maxiter, thres, quiet){
  ## data is n by p
  ## beta_ini includes the intercept
  ## cluster is the cluster membership matrix, length n. Take entries 1 to nCluster
  ## offset corresponds to h_i, here it is a vector of length n
  
  p <- ncol(data)
  x <- getdesign(cluster)
  ## add the intercept in x
  x <- cbind(1, x)
  ## get lambda
  lambda <- c(0, 1/sigmas^2)
  
  beta <- beta_ini
  betaDiffs <- c()
  percConverged <- c()
  converged_flag <- rep(0, p)
  
  converged <- 0
  iter <- 1
  while (!converged & iter<maxiter){
    mu <- offset*exp(x%*%beta) 
    z <- log(mu/offset) + (data-mu)/mu
    
    updateBetaRidge <- function(r){
      if (converged_flag[r]==0){
        JrRidge <- crossprod(x, mu[,r]*x) + diag(lambda)
        betar <- solve(JrRidge, crossprod(x, matrix(mu[,r]*z[,r], ncol=1)) )
        return( as.vector(betar) )
      } else {
        return(beta[,r])
      }
    }
    
    betaPre <- beta
    beta <- sapply(1:p, updateBetaRidge)
    betaDiffs <- rbind(betaDiffs, colSums(abs(beta - betaPre)) )
    converged_flag <- (betaDiffs[iter,] <= thres) + 0
    percConverged <- c(percConverged, sum(converged_flag)/p )
    
    if (!quiet){
        cat("\r", round(percConverged[length(percConverged)]*100), "%", "converged")  
    }
    
    if (percConverged[length(percConverged)]==1){
      converged <- 1
    }
    iter <- iter + 1
  }
  mu <- offset*exp(x%*%beta) 
  
  ## get all the contrast matrices and cluster label for each contrast
  nCluster <- length(unique(cluster))
  tmp <- getContrast(nCluster)
  constrasts <- tmp$constrasts
  
  ## pad the intercept with 0 in the constrast
  constrasts <- cbind(0, constrasts)
  constrastsCluster <- tmp$constrastsCluster
  
  calRidgePvalueA <- function(r){
    ## calculates the p-value for a small hypothesis
    if ( !is.na(beta[1,r]) ){
      Jr <- crossprod(x, mu[,r]*x)
      JrRidgeInv <- solve(Jr + diag(lambda))
      covRidge <- JrRidgeInv%*%Jr%*%JrRidgeInv
      ##
      betaC <- constrasts%*%matrix(beta[,r], ncol=1)  
      CcovC <- constrasts%*%covRidge%*%t(constrasts)
      SEbetaC <- sqrt(diag(CcovC))
      ##
      pvs <- pnorm(betaC/SEbetaC, lower.tail = FALSE)
      return(pvs)
    } else {
      return(rep(NA, length(constrastsCluster)))
    }
  }
  if (!quiet){
    cat("\nCalculating the p-values\n")
  }
  pvsA <- sapply(1:p, calRidgePvalueA) 
  ## get the p-value for the large null hypothesis by taking maximum
  if (nCluster>2){
    pvs <- c()
    for (i in 1:nCluster){
      pvs <- cbind(pvs, apply(pvsA[which(constrastsCluster==i),], 2, max))
    }  
  } else {
    pvs <- t(pvsA)
  }
  return(list(beta=beta, pvalue=pvs))
}


getClusterSpecificPvalue <- function(data, cluster, offset, landmark=NULL, maxiter=1000, thresMLE=10^-3, thresMAP=10^-5, quiet=FALSE){
  ## the main function for peak selection
  ## data is nPeaks(p) by Cells(n)
  ## cluster is the cluster membership matrix, length n. Take entries 1 to nCluster
  ## offset corresponds to h_i, here it is a vector of length n
  ## landmark is optional. If landmark is provided, we only do hypothesis testing on the union of landmark peaks
  if (!is.null(landmark)){
    ## take union of the landmark peaks
    unionlandmark <- which(rowSums(landmark)!=0)
    p <- nrow(data)
    data <- data[unionlandmark,]
  } else {
  }
  
  data <- t(data) # make data n by p
  ## remove the offset=0 samples
  data <- data[which(offset>0),]
  cluster <- cluster[which(offset>0)]
  offset <- offset[which(offset>0)]
  ## get beta_MLE
  if (!quiet){
    cat("\nEstimating beta MLE\n")  
  }
  betaMLE_ini <- matrix(0, nrow=length(unique(cluster)), ncol=ncol(data))
  betaMLE <- getbetaMLE(data=data, cluster=cluster, offset=offset, beta_ini=betaMLE_ini, maxiter=maxiter, thres=thresMLE, quiet=quiet)
  
  ## get the empirical prior sigma
  sigmas <- getSigmaPrior(betaMLE)
  
  ## get beta_MAP and the p-value
  if (!quiet){
    cat("\nEstimating beta MAP\n")
  }
  betaMAP_ini <- rbind(0, matrix(0, nrow=length(unique(cluster)), ncol=ncol(data)))
  result <- getbetaMAP(data=data, cluster=cluster, offset=offset, sigmas=sigmas, beta_ini=betaMAP_ini, maxiter=maxiter, thres=thresMAP, quiet=quiet)
  
  if (!is.null(landmark)){
    betaMAP <- matrix( nrow=length(unique(cluster))+1, ncol=p )
    pvalue <- matrix( nrow=p, ncol=length(unique(cluster)) )
    betaMAP[, unionlandmark] <- result$beta
    pvalue[unionlandmark,] <- result$pvalue
  } else {
    betaMAP <- result$beta
    pvalue <- result$pvalue
  }
  colnames(pvalue) <- paste("Cluster", 1:length(unique(cluster)) )
  return(list(betaMAP=betaMAP, sigmas=sigmas, pvalue=pvalue))
}




data_type <- 'blood'
data_name <- 'donor_BM0828'
multi_donor <- F
bulk_name <- ''
peak_selection_k <- 3
K2 <- 5

if (multi_donor == T) {
  task_name <- sprintf('Our_%s_peak%03d_dim%d_multidonor%s',data_name,100*peak_selection_k,K2,bulk_name);
}else{
  task_name <- sprintf('Our_%s_peak%03d_dim%d%s',data_name,100*peak_selection_k,K2,bulk_name);
}
data <- readMat(sprintf('~/data/dataForComparison_0314/%s_peak%03d.mat',data_name,100*peak_selection_k))

label_file_name <- sprintf('/home/zhangwenyu/dataForComparison_clustering/clusters_our/%s_clusters.tsv',task_name)
cluster_assign <- read.csv(label_file_name, sep='\t')
data$cell.labels <- cluster_assign$louvain
data$cell.labels <- unlist(data$cell.labels) + 1

resultPeakSelection <- getClusterSpecificPvalue(data=data$count.mat, cluster=data$cell.labels, offset=colSums(data$count.mat))

save(resultPeakSelection, file = sprintf("%s_peak%03d.RData",data_name,100*peak_selection_k))

