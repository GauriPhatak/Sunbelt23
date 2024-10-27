library(tidyverse)
library(dplyr)
library(igraph)
library(network)
library(intergraph)
library(ergm)
library(RColorBrewer)
source("Code/HelperFuncs.R")

prob2logit <- function(x) {
  return(log(x / (1 - x)))
}

MSE <- function(Wtmat,Z,beta,N){
  #cbind(1, Wtmat)
  pred <- cbind(1, Wtmat) %*% t(beta)
  mse <- colSums((pred - Z) ^ 2) / N
  return(mse)
}

NWSimBin <- function(nc, k,  pC, N, pClust, B, o, dist,dir,covTypes = c("binary", "continuous"), 
                     CovNamesLin = c() ,CovNamesLP = c(), missing = NULL) {
  C <- 0
  while (length(table(C)) < nc) {
    C <- sample(
      x = letters[1:nc],
      size = N,
      replace = TRUE,
      prob = pClust
    )
  }
  ## Create an empty network i.e. no cluster assignments or category assignments or edges
  net <- network(N, directed = dir, density = 0)
  ## Assign clusters to the nodes
  net %v% 'Cluster' <- C
  if (c("binary") %in% covTypes) {
    ## Based on the probability matrix pC assign indicates the binary assignment of covariates to a community
    binVal <- matrix(ncol = k, nrow = N, NA)
    for (i in 1:nc) {
      for (j in 1:k) {
        binVal[which(as.numeric(net %v% 'Cluster') == i), j] <- rbinom(length(which(as.numeric(
          net %v% 'Cluster'
        ) == i)), 1, pC[i, j])
      }
    }
    for (j in 1:k) {
      binVal[sample(1:N, round(missing * N / 100)), j] <- NA
      net %v% CovNamesLP[j] <- binVal[, j]
    }
  }
  if (c("continuous") %in% covTypes) {
    ## Based on the mean and variance indicated by the dist matrix we assign continuous values to the network covariates.
    contVal <- matrix(ncol = o, nrow = N, NA)
    
    for (i in 1:nc) {
      for (j in 1:o) {
        contVal[which(as.numeric(net %v% 'Cluster') == i), j] <- rnorm(length(which(as.numeric(
          net %v% 'Cluster'
        ) == i)), dist[i, j][[1]][[1]], dist[i, j][[1]][[2]])
      }
    }
    for (j in 1:o) {
      contVal[sample(1:N, round(missing * N / 100)), j] <- NA
      net %v% CovNamesLin[j] <- contVal[, j]
    }
  }
  ## Use the probability of connection values defined earlier for both cluster and category for fitting the data.
  g.sim <- simulate(
    net ~ nodemix("Cluster", levels = TRUE, levels2 = TRUE),
    nsim = 1,
    coef = prob2logit(c(B)),
    control = ergm::control.simulate(MCMC.burnin = 10000, MCMC.interval =
                                       1000)
  )
  return(g.sim)
}

flip <- function(p) {
  return(sample(
    x = c(T, F),
    size = 1,
    prob = c(p, 1 - p)
  ))
}

SetBinarycov <- function(G, k, N,nc, k_start, Cnew,connType, pC,CovNamesLP, missing) {
  ## Based on the probability matrix pC assign indicates the binary assignment of covariates to a community
  ## setting covariates for outgoing communities. hence Cnew will be -1
  binVal <- matrix(ncol = k, nrow = N, 0)
  for (i in 1:nc) {
    pc_i <- k_start+i
    for (j in 1:k) {
      pc_j <- k_start+j
      idx <- which(Cnew[, i] == connType)
      binVal[idx, j] <- rbinom(length(idx), 1, pC[pc_i, pc_j])
    }
  }
  if (missing > 0) {
    
    for (j in 1:k) {
      binVal[sample(1:N, round(missing * N / 100)), j] <- NA
      G <- G %>%
        set_vertex_attr(name = CovNamesLP[j], value = c(binVal[, j]))
    }
  } else{
    for (j in 1:k) {
      G <- G %>%
        set_vertex_attr(name = CovNamesLP[j], value = c(binVal[, j]))
    }
  }
  
  return(G)
}

SetContinuousCov <- function(G,o, N,nc, o_start, Cnew, connType, dist,CovNamesLin,missing){
  ## Based on the mean and variance indicated by the dist matrix we assign continuous values to the network covariates.
  contVal <- matrix(ncol = o, nrow = N, rnorm(N,1,1)) #rnorm(N,0,1)
  
  for (i in 1:nc) {
    pc_i <- o_start+i
    for (j in 1:o) {
      pc_j <- o_start+j
      idx <- which(Cnew[, i] == connType)
      contVal[idx, j] <- rnorm(length(idx),
                               dist[[pc_j]][[i]][1],
                               dist[[pc_j]][[i]][2])
    }
  }
  
  if(missing >0 ){
    for (j in 1:o) {
      contVal[sample(1:N, round(missing * N / 100)), j] <- NA
      G <- G %>% set_vertex_attr(name = CovNamesLin[j],
                                 value = c(contVal[, j]))
    }
  } else{
    for (j in 1:o) {
      G <- G %>%set_vertex_attr(name = CovNamesLin[j], 
                                value = c(contVal[, j]))
    }
  }
  
  return(G)
  
}

genBipartite <- function(N,nc,pClust,k_in,k_out,o_in,o_out,pC,dist,covTypes,
                         CovNamesLPin,CovNamesLPout,CovNamesLinin,CovNamesLinout,
                         pConn, dir,dirPct, epsilon,missing) {
  C <- 0
  while (length(table(C)) < nc) {
    C <- sample(
      x = letters[1:nc],
      size = N,
      replace = TRUE,
      prob = pClust[1:nc]
    )
  }
  op <- model.matrix( ~ C - 1)
  colnames(op) <- letters[1:nc]
  
  ## overlapping various combinations of communities.
  Cnew <- data.frame(op)
  fin <- c()
  for (i in 2:(nc - 1)) {
    c <- t(combn(letters[1:nc], i))
    for (j in 1:dim(c)[1]) {
      com <- which(C %in% c[j, ])
      fin <- c(fin, com)
      (idx <- sample(com , size = length(com) * pClust[nc + i + j - 2], replace =
                       FALSE))
      Cnew[idx, c[j, ]] <- 1
    }
  }
  
  ## percentage overlapping all the communities
  idx <- sample(unique(fin),
                size = length(unique(fin)) * pClust[length(pClust)],
                replace = FALSE)
  Cnew[idx, ] <- 1
  
  ## Determine the direction of the community membership. i.e. if it is is incoming or outgoing community affiliation.
  ## -1 for outgoing membership and +1 for incoming membership
  Cold <- Cnew
  if (dir == TRUE) {
    for (i in 1:nc) {
      ## sampled outgoing membership by percentage
      sampOut <- sample(which(Cnew[, i] == 1),
                        size = length(which(Cnew[, i] == 1)) * dirPct[i],
                        replace = FALSE)
      Cnew[sampOut, i] <- -1
    }
  }
  
  A <- matrix(0, nrow = N, ncol = N)
  Apuv <- A
  for (i in 1:N) {
    for (j in 1:N) {
      p_uv <- 1 - prod(1 - ((as.numeric(Cnew[i, ] - Cnew[j, ] > 1)) * pConn)) + epsilon
      Apuv[i, j] <- p_uv
      A[i, j] <- purrr::rbernoulli(1, p_uv)
    }
  }
  
  diag(A) <- 0
  g.sim <- graph_from_adjacency_matrix(A, mode = "directed")
  
  for (i in 1:nc) {
    g.sim <- g.sim %>%
      set_vertex_attr(name = letters[i], value = c(Cnew[, i]))
  }
  
  ## Creating binary covariates in the network.
  if (c("binary") %in% covTypes) {
    if (length(CovNamesLPin) > 0) {
      g.sim <- SetBinarycov(g.sim, k_in, N, nc, k_start = 0, Cold, 
                            connType = 1,pC, CovNamesLPin, missing)
    }
    if (length(CovNamesLPout) > 0) {
      g.sim <- SetBinarycov(g.sim, k = k_out, N, nc, k_start = k_in, Cold, 
                            connType = 1, pC, CovNamesLPout, missing)
    }
  }
  
  if (c("continuous") %in% covTypes) {
    if (length(CovNamesLinin) > 0) {
      g.sim <- SetContinuousCov(g.sim, o_in, N, nc, o_start = 0, Cnew, 
                                connType = 1, dist,CovNamesLinin,missing)
    }
    if (length(CovNamesLinout) > 0) {
      g.sim <- SetContinuousCov(g.sim, o_out, N,nc, o_start = o_in, Cnew, 
                                connType = -1, dist, CovNamesLinout, missing)
    }
  }
  
  g.sim <- g.sim  %>%
    set_vertex_attr(name = "Cluster", value = C)
  return(list(g.sim, Apuv))
}

sigmoid <- function(W, Ftotn) {
  Q <- 1 / (1 + exp(-1 * (W %*% t(Ftotn))))
  return(Q)
}

LRParamUpdt <- function(W, X, Ftot,Htot, alpha, lambda, missing) {
  Ftotn <- cbind(1, Ftot, Htot)
  Q <- sigmoid(W, Ftotn)
  
  W_grad <- t(X - t(Q)) %*% Ftotn
  W_new <- W + alpha * (W_grad - lambda * (sign(W)))
  return(W_new)
}

PredictCovLR <- function(X, W, Ftot, missing) {
  for (i in 1:dim(X)[2]) {
    X[which(missing[, i]), i] <- as.numeric(sigmoid(W[i, ], cbind(1, Ftot[which(missing[, i]), ])) > 0.5)
  }
  return(X)
}

SigmaSqCalc <- function(Z, beta, Ftot, missVals) {
  sigmaSq <- rep(0, dim(beta)[1])
  Fm <- cbind(1, Ftot)
  for (i in 1:dim(beta)[1]) {
    sigmaSq[i] <- sum((Z[!missVals[, i], i] - (t(beta[i, ]) %*% t(Fm[!missVals[, i], ]))) ^ 2, na.rm = TRUE) / sum(!missVals[, i])
  }
  return(sigmaSq)
}

PredictCovLin <- function(Ftot, Z, beta, missVals) {
  Fm <-  Ftot#cbind(1,Ftot)
  for (i in 1:dim(Z)[2]) {
    if(sum(missVals) >1 ) {
      Z[which(missVals[, i]), i] <- Fm[which(missVals[, i]),] %*% as.matrix(beta[i, ])
      
    } else{
      Z[, i] <- Fm %*% as.matrix(beta[i,])
    }
  }
  return(Z)
}

LinRParamUpdt <- function(beta, Z, Mat, alpha, lambda, N, missVals) {
  if(printFlg == TRUE){
    print("In Linparamupdate Func")
  }
  beta_new <- matrix(0, nrow = dim(beta)[1], ncol = dim(beta)[2])
  Fm_n <- cbind(1, Mat)
  #for (i in 1:dim(beta)[1]) {
  #  beta_new <- t(MASS::ginv(t(Fm_n) %*% Fm_n) %*% (t(Fm_n) %*% Z))
  #print(beta_new)
  #}
  
  ## Calculating update to the beta using gradient descent
  y <-  PredictCovLin(Fm_n, Z, beta, missVals)
  gradient <- (t(y - Z) %*% Fm_n)/N
  #print(paste0(gradient))
  #alpha <- 0.01
  beta_new <-  beta - (alpha*gradient)
  return(beta_new)
}

CalcQuk <- function(Fmat, Hmat, W){
  Q_uk <- 1 / (1 + exp(-1 * (cbind(1, Fmat, Hmat) %*% t(W))))
  return(Q_uk)
}

Lx_cal <- function(W, N, Fmat, Hmat, X){
  S <- 0
  
  if (length(W) > 0) {
    Q <- CalcQuk(Fmat, Hmat, W)
    #Win %*% t(cbind(1, Fmat))
    S <- sum(log(Q) * X + log(1-Q) * (1- X))
    # for (i in 1:N) {
    #   Q_uk <- CalcQuk(Fmat[i,], W)#1 / (1 + exp(-1 * (c(1, Fmat[i, ]) %*% t(W))))#1/(1+exp(-1*c(1, Ftot[i,]) %*% t(W)))
    #   #print(Q_uk)
    #   su <- sum((log(Q_uk) * X[i, ]) + (log(1 - Q_uk) * (1 - X[i, ])))
    #   #print(su)
    #   if(is.nan(su)){su <- 0}
    #   S <- S + su
    # }
  }
  
  return(S)
}

Lz_cal <- function(beta,N,Z,Fmat,sigmaSq ){
  if(printFlg == TRUE){
    print("In Lz_Calc Func")
  }
  S <- 0
  if (length(beta) > 0) {
    for (j in 1:dim(beta)[1]) {
      S_tmp <- sum(( Z[,j]  - cbind(1,Fmat) %*% beta[j, ]) ^ 2)
      #for (i in 1:N) {
      #c(1, Fmat[i, ])
      #  S_tmp <- sum(c(S_tmp, (Z[i, j] - sum(c(1, Fmat[i, ]) * beta[j, ], na.rm = TRUE)) ^ 2),na.rm = TRUE)
      #}
      S <- S - ((N * log(sigmaSq[j]) / 2) - (S_tmp / (2 * sigmaSq[j])))
    }
  }
  return(S)
}

Lg_cal <- function(G, Ftot, Htot,Eneg){
  scale <- 1
  E <- as_edgelist(G)
  Fv <- Ftot[E[, 1], ]
  Hv <- Htot[E[, 2], ]
  N <-  gorder(G)
  ## Variable to save the loglik contributed by nodes with edges between them
  ## Implementing rowsums Hopefully this is a little faster
  S11 <- sum(log(1-exp(-1 * rowSums(Fv * Hv))))
  
  ## Creating a graph that is negative of the given graph to get edges not present in G
  E2 <-  Eneg 
  Fv2 <- Ftot[E2[, 1], ]
  Hv2 <- Htot[E2[, 2], ]
  
  ## Variable to save the loglik contributed by nodes without edges between them
  
  S12 <- sum(rowSums(Fv2 *Hv2))
  S1 <- S11 - S12
  
  return(S1)
}

#G,Ftot,Htot,Win, Wout,X_in, X_out,betain, betaout,Z_in, Z_out, sigmaSq,alphaLL
findLLDir <- function(G,  Ftot, Htot, Win = NA, Wout = NA, X_in = NA,X_out = NA,
                      betain = NULL, betaout = NULL,Z_in = NA, Z_out = NA,
                      sigmaSqin = NA,sigmaSqout = NA, alphaLL, Eneg) {
  if(printFlg == TRUE){
    print("In Find LLDir Func")
  }
  # scale <- 1
  # E <- as_edgelist(G)
  # Fv <- Ftot[E[, 1], ]
  # Hv <- Htot[E[, 2], ]
  # N <-  gorder(G)
  # ## Calculating the L_g part of the log likelihood for structure of the network
  # 
  # ## Variable to save the loglik contributed by nodes with edges between them
  # ## Implementing rowsums Hopefully this is a little faster
  # S11 <- sum(log(1-exp(-1 * rowSums(Fv * Hv))))
  # ## Older implementation
  # #S11 <- 0
  # #for (i in 1:dim(E)[1]) {
  # #  S11 <- S11 + log(1 - exp(-1 * sum(Fv[i, ] * Hv[i, ], na.rm=TRUE)))
  # #}
  # 
  # ## Creating a graph that is negative of the given graph to get edges not present in G
  # #A <- 1 - (1 * as.matrix(as_adjacency_matrix(G)))
  # E2 <-  Eneg #as_edgelist(graph_from_adjacency_matrix(A))
  # Fv2 <- Ftot[E2[, 1], ]
  # Hv2 <- Htot[E2[, 2], ]
  # 
  # ## Variable to save the loglik contributed by nodes without edges between them
  # #S12 <- 0
  # #for (i in 1:dim(E2)[1]) {
  # #  S12 <- S12 + sum(Fv2[i, ] * Hv2[i, ], na.rm = TRUE)
  # #}
  # S12 <- sum(rowSums(Fv2 *Hv2))
  # S1 <- S11 - S12
  ## Calculating the L_g part of the log likelihood for structure of the network
  S1 <- Lg_cal(G, Ftot, Htot, Eneg)
  
  ## Calculating the L_xin part of the log likelihood for binary covariates
  S2_in <- 0 #Lx_cal(Win, N, Ftot, X_in)
  
  ## Calculating the L_xout part of the log likelihood for binary covariates
  S2_out <- Lx_cal(W=Wout, N = N, Fmat= Ftot, Hmat =Htot, X=X_out)
  
  ## adding the loglikelihood from the covariates
  S2 <- S2_in + S2_out
  
  ## Calculating the L_zin part of the log likelihood for continuous covariates
  S3_in <- Lz_cal(beta = betain,N = N,Z = Z_in,Fmat = Ftot,sigmaSq = sigmaSqin )
  
  ## Calculating the L_zout part of the log likelihood for continuous covariates
  S3_out <- Lz_cal(betaout,N,Z_out,Htot,sigmaSqout )
  S3 <- S3_in + S3_out
  #print(paste0(S1, " ", S2, " ", S3))
  
  ## Calculating the final log likelihood
  #if (is.null(alphaLL)) {
  ll <- S1 + (S2 + S3)
  #} else{
  #  ll <- alphaLL * S1 + (1 - alphaLL) * (S2 + S3)
  #}
  return(ll)
}

mode <- function(x) {
  which.max(tabulate(x))
}

initCov <- function(covtmp, CovNames) {
  X_in <- as.matrix(covtmp %>%
                      dplyr::select(all_of(CovNames)))
  missVals <- is.na(X_in)
  numMissVal <- sum(missVals)
  
  return(list(X_in, numMissVal))
}

initWtmat <- function(G, mode, N, nc){
  #Wtmat <- igraph::degree(G, mode = mode) / sum(igraph::degree(G, mode = mode))
  Wtmat <- matrix(nrow = N, ncol = nc, runif(nc * N))
  #Wtmat <- matrix(nrow = N, ncol = nc, 0.001)
  return(Wtmat)
}

initW <- function(k, nc,missVals, X, Fmat, lambda, alpha){
  W <- NULL
  if (k > 0) {
    W <- matrix(nrow = k, ncol = (nc*2) + 1, 0)
    ###*** For future updates ***######
    if (sum(missVals) > 0) {
      X_in <- PredictCovLR(X, W, Fmat, missVals)
    }
    #Win <- LRParamUpdt(W, X, Fmat, alpha, lambda, missVals)
  }
  return(W)
}

initbeta <- function(o, nc,missVals, Z, Fmat, alpha,N, lambda){
  beta <- NULL
  sigmaSq <- NULL
  if (o > 0) {
    beta <- matrix(nrow = o, ncol = (nc) + 1, 0)
    
    beta <- LinRParamUpdt(beta, Z, Fmat, alpha, lambda, N, missVals)
    sigmaSq <- SigmaSqCalc(Z, beta, Fmat, missVals)
    if (sum(missVals) > 0) {
      Z <- PredictCovLin(Fmat, Z, beta, missVals)
    }
  }
  return(list(beta, sigmaSq))
}

GraphComntyWtUpdt <- function(f_u,h_u,h_neigh, h_sum){
  ## First part of log likelihood of G based on the structure of the network
  #a <- exp(-1 * (f_u %*% t(Fvmat)))
  a <- exp(-1 * (f_u %*% t(h_neigh)))
  b <- a / (1 - a)
  llG_1 <- t(h_neigh) %*% t(b)
  ## 2nd part of log likelihood of G
  llG_2 <-  as.matrix(h_sum - t(h_u) - as.matrix(colSums(h_neigh)))#as.matrix(colSums(Fvnot))
  ## Total llG: This math has been verified
  llG <- llG_1 - llG_2
  ## experimental scaled version of the above llog lik
  scale <- 1
  llG <- scale * llG
  
  return(llG)
}

BincovCmntyWtUpdt <- function(nc, f_u,h_u, W, X_u){
  llX <- rep(0, nc*2 + 1)
  if (length(W) > 0) {
    ## Updating the logistic regression weights
    ## Gradient function to update log likelihood for covariates
    ## adding intercept term
    Q_u <- CalcQuk(f_u,h_u,W)
    ## This math has been verified
    llX <- t((X_u - Q_u) %*% W)
  }
  return(llX)
}

ContcovCmntyWtUpdt <- function(nc, Z_u, beta, f_u, sigmaSq){
  
  llZ <- rep(0, nc*2 + 1)
  f_u <- matrix(c(1, f_u), ncol = nc + 1)
  if (length(beta) > 0) {
    ## Adding the gradient based on continuous covariates
    llZ <-  ((Z_u - beta %*% t(f_u)) / sigmaSq) %*% beta
  }
  
  return(llZ)
}


#f_u,  Fvmat, Fvnot,
CmntyWtUpdt <- function(f_u ,h_u ,h_neigh ,h_sum ,X_u = NULL ,W = NULL ,Z_u = NULL,
                        beta = NULL ,sigmaSq ,alpha ,alphaLL ,nc ,start ,end ,mode ) {
  if(printFlg == TRUE){
    print("In CmntyWtupdate Func")
  }
  tv <- 0 #0.000001
  
  ## Gradient function to update log likelihood of G
  #First part of Likelihood based on graph only.
  llG <- GraphComntyWtUpdt(f_u,h_u,h_neigh,h_sum)#Fvmat, Fvnot)
  
  #Second part of Likelihood based on binary covariates
  
  if(mode == "out"){
    llX <- BincovCmntyWtUpdt(nc, f_u, h_u, W, X_u)
  }else if(mode == "in"){
    W <- cbind(W[,1], W[,start:end], W[,(start-nc):(end-nc)])
    llX <- BincovCmntyWtUpdt(nc, f_u, h_u, W, X_u)
  }
  
  #Third part of likelihood based on binary covariates
  llZ <- ContcovCmntyWtUpdt(nc, Z_u, beta, f_u, sigmaSq)
  
  ## Function to update the community weights vector using non negative matrix factorization
  #if (is.null(alphaLL)) {
  #[2:(nc + 1)]

  f_u_new <- f_u + (alpha * (t(llG) + llX[2:(nc+1)]))
  # } else{
  #   f_u_new <- f_u + t((alpha * (
  #     alphaLL * llG + (1 - alphaLL) * (llX[2:(nc + 1)] + llZ[2:(nc + 1)])
  #   )))
  # }
  f_u_new[f_u_new < 0] <- tv
  return(f_u_new)
}

updateWtmat <- function(G,Wtm1,Wtm2,mode,s,nc,X,Z,k,o,beta,W,alphaLL,missVals,alpha, start, end){
  if(printFlg == TRUE){
    print("In updateWtmat Func")
  }
  Wtmat <- Wtm1
  
  ## Getting sum of all the weights of Fv. Have to modify this for CoDa as opposed to CESNA 
  Wtm2_sum <- as.matrix(colSums(Wtm2))
  
  for (i in s) {
    Wtm1_u <- matrix(Wtmat[i, ], ncol = nc)
    Wtm2_u <- matrix(Wtm2[i, ], ncol = nc)
    
    ##The neighbours of this node
    neigh <- neighbors(G, i, mode = mode)
    Wtm2_neigh <- matrix(Wtm2[neigh, ], ncol = nc)
    
    #neighAll <- neighbors(G, i, mode = "all")
    #Wtm2_vnot<- matrix(Wtm2[-neighAll, ], ncol = nc)
    #Wtm2_u <- matrix(Wtm2[i, ], ncol = nc)
    X_u <- NULL
    if (k > 0) {
      X_u <- matrix(X[i, ], ncol = k)
    }
    
    Z_u <- NULL
    if (o > 0) {
      Z_u <- matrix(Z[i, ], ncol = o)
      sigmaSq <- SigmaSqCalc(Z, beta, Wtm1, missVals)
    }
    
    # Wtmat[i, ] <- CmntyWtUpdt(Wtm1_u,Wtm2_v,Wtm2_vnot,
    #                           X_u,W,Z_u, beta,
    #                           sigmaSq,alpha,alphaLL,nc)
    Wtmat[i, ] <- CmntyWtUpdt(Wtm1_u,Wtm2_u,Wtm2_neigh, Wtm2_sum,
                              X_u,W,Z_u, beta,sigmaSq,alpha,alphaLL,
                              nc, start, end, mode)
    
  }
  
  return(Wtmat)
}

updateLogisticParam <- function(W,X,Wtm,Wtm2, missVals, lambda,alpha){
  Wret <- W
  Xret <- X
  if (length(W) > 0) {
    if (sum(missVals) > 0) {
      Wret <-  LRParamUpdt(W, X, Wtm,Wtm2 ,alpha, lambda, missVals)
      Xret <-  PredictCovLR(X, W, Wtm, missVals)
    } else{
      Wret <-  LRParamUpdt(W, X, Wtm,Wtm2, alpha, lambda, missVals)
    }
  }
  
  return(list(Wret,Xret))
}

updateLinearRegParam <- function(beta,missVals,Z,Wtm, alpha,lambda,N){
  if(printFlg == TRUE){
    print("In updateLinearRegParam Func")
  }
  betaret <- beta
  Zret <-  Z
  sigmaSq <- NULL
  if (length(beta) > 0) {
    betaret <- LinRParamUpdt(beta, Z, Wtm, alpha, lambda, N, missVals)
    sigmaSq <-  SigmaSqCalc(Z, betaret, Wtm, missVals)
    
    if (sum(missVals) > 0) {
      Zret <- PredictCovLin(Wtm, Z, beta, missVals)
    }
  }
  return(list(betaret,Zret,sigmaSq))
}

accuCalc <- function(k,W,Wtm,X){
  if(printFlg == TRUE){
    print("In accuCalc Func")
  }
  accuracy <- matrix(nrow = 0, ncol = k)
  if (k > 0) {
    predLR <- apply(sigmoid(W, cbind(1, Wtm)) > 0.5, 1, as.numeric)
    acc <- rep(0, k)
    for (i in 1:k) {
      cm <-  confusionMatrix(factor(predLR[, i], levels = c(0, 1)), factor(X[, i]))
      acc[i] <- round(cm$overall[1], 2)
    }
    accuracy <- rbind(accuracy, acc)
  }
  
  return(accuracy)
}

getNWcov <- function(G,CovNamesLinin, CovNamesLinout,CovNamesLPin,CovNamesLPout,
                     k_in, k_out, o_in, o_out ){
  covtmp <-  as.data.frame(vertex_attr(G)) %>%
    dplyr::select(all_of(
      c(CovNamesLinin, CovNamesLinout, CovNamesLPin, CovNamesLPout)
    ))
  X_in <- 0
  X_out <- 0
  Z_in <- 0
  Z_out <- 0
  
  missValsin <- 0
  numMissValin <- 0
  missValsout <- 0
  numMissValout <- 0
  
  if (k_in > 0) {
    X_in <- as.matrix(covtmp %>%
                        dplyr::select(all_of(CovNamesLPin)))
    missValsin <- is.na(X_in)
    numMissValin <- sum(missValsin)
  }
  if (k_out > 0) {
    X_out <- as.matrix(covtmp %>%
                         dplyr::select(all_of(CovNamesLPout)))
    missValsout <- is.na(X_out)
    numMissValout <- sum(missValsout)
  }
  
  if (o_in > 0) {
    Z_in <- as.matrix(covtmp %>%
                        dplyr::select(all_of(CovNamesLinin)))
    ## Find the missing values
    missValsin <- is.na(Z_in)
    numMissValin <- sum(missValsin)
  }
  if (o_out > 0) {
    Z_out <- as.matrix(covtmp %>%
                         dplyr::select(all_of(CovNamesLinout)))
    ## Find the missing values
    missValsout <- is.na(Z_out)
    numMissValout <- sum(missValsout)
  }
  
  return(list(X_in, X_out, Z_in, Z_out, 
              missValsin, missValsout, 
              numMissValin, numMissValout))
}

getDelta <- function(N){
  delta <- sqrt(-1*log(1-(1/N)))
  #delta <- sqrt(-1*log(1-epsilon))
  return(delta)
}

memOverlapCalc <- function(Fm, Hm, delta, N, nc){
  ol <- as.data.frame(matrix(0, ncol = nc, nrow = N))
  for (i in 1:N) {
    ol[i, ] <- (as.numeric(Fm[i, ] > delta) + as.numeric(Hm[i, ] > delta)) > 0
    #memoverlap[i, ] <- as.numeric(memoverlap[i, ] > 0)
  }
  return(ol)
}

OmegaIdx <- function(G, Fm, Hm, N, delta, nc) {
  if(printFlg == TRUE){
    print("In OmegaIdx Func")
  }
  memoverlap <- memOverlapCalc(Fm, Hm, delta, N, nc)
  OrigVal <-  as.data.frame(vertex_attr(G)) %>%
    dplyr::select(all_of(c(letters[1:nc])))
  OrigVal[OrigVal == -1] <- 1
  
  ol <- list()
  for(i in 1:nc){
    ol[[i]] <- which(OrigVal[,i] == 1)
  }
  
  posCom <- rbind(expand.grid(ol[[1]],ol[[1]]),
                  expand.grid(ol[[2]],ol[[2]]),
                  expand.grid(ol[[3]],ol[[3]]))
  posCom <- paste0(posCom[,1],"-",posCom[,2])
  valuesOrig <- as.data.frame(table(posCom))
  
  ol <- list()
  for(i in 1:nc){
    ol[[i]] <- which(memoverlap[,i] == 1)
  }
  
  posCom <- rbind(expand.grid(ol[[1]],ol[[1]]),
                  expand.grid(ol[[2]],ol[[2]]),
                  expand.grid(ol[[3]],ol[[3]]))
  posCom <- paste0(posCom[,1],"-",posCom[,2])
  valuesMem <- as.data.frame(table(posCom))
  
  values <- full_join(valuesOrig, valuesMem, by = "posCom")
  tot <- expand.grid(1:N, 1:N)
  tot <- as.data.frame( paste0(tot[,1],"-",tot[,2]))
  colnames(tot) <- c("posCom")
  gbg <- left_join(tot,values, by = "posCom")
  gbg[is.na(gbg)] <- 0
  gbg$matched <- gbg$Freq.x == gbg$Freq.y
  oi <- sum(gbg$matched)
  
  oi <- oi / (N ^ 2)
  return(oi)
}

GTLogLik <- function(G,nc,pConn,alphaLL,CovNamesLinin,CovNamesLinout, 
                     CovNamesLPin,CovNamesLPout,k_in,k_out,o_in,o_out,epsilon){
  alpha_c <- -1*log((1-pConn))
  comAff <-  as.data.frame(vertex_attr(G)) %>%
    dplyr::select(all_of(letters[1:nc]))
  
  b <- getNWcov(G, CovNamesLinin, CovNamesLinout, 
                CovNamesLPin, CovNamesLPout, 
                k_in,k_out,o_in,o_out)
  
  X_in <- b[[1]]
  X_out <- b[[2]]
  
  Z_in <- b[[3]]
  Z_out <- b[[4]]
  
  missValsin <- b[[5]]
  missValsout <- b[[6]]
  numMissValin <- b[[7]]
  numMissValout <- b[[8]]
  
  t <- 0 #0.000001
  Fm <- matrix(t, nrow = dim(comAff)[1],ncol = dim(comAff)[2])
  Fm[comAff == 1] <- 1
  Fm <- sweep(Fm, MARGIN=2, alpha_c, '*')
  Fm <- Fm + (-(1/2)*log((1-epsilon)))
  
  Hm <- matrix(t, nrow = dim(comAff)[1],ncol = dim(comAff)[2])
  Hm[comAff == -1] <- 1
  Hm <- sweep(Hm, MARGIN=2, alpha_c, '*')
  Hm <- Hm + (-(1/2)*log((1-epsilon)))
  
  
  Win <- NULL
  # if(k_in >0 ){
  #   Win <- matrix(0, nrow = k_in, ncol = nc*2+1 )
  #   for(i in 1:dim(X_in)[2]){
  #     Win[i,] <- coef(glm(X_in[,i] ~ 1+Fm[,1]+Fm[,2]+Fm[,3], 
  #                         family = "binomial"))
  #   }
  # }
  
  Wout <- NULL
  binLL <- 0
  if(k_out >0){
    Wout <- matrix(0, nrow = k_out, ncol = nc*2+1 )
    for(i in 1:dim(X_out)[2]){
      mod <- glm(X_out[,i] ~ Fm[,1]+Fm[,2]+Fm[,3]+Hm[,1]+Hm[,2]+Hm[,3],family = "binomial")
      binLL <- binLL + as.numeric(logLik(mod))
      Wout[i,] <- coef(mod)
    }
  }
  betain <- NULL
  if(o_in >0 ){
    betain <- matrix(0, nrow = o_in, ncol = nc+1)
    for(i in 1:dim(Z_in)[2]){
      betain[i,] <- coef(lm(Z_in[,i] ~ 1+Fm[,1]+Fm[,2]+Fm[,3]))
    }
    sigmaSqin <- SigmaSqCalc(Z_in,betain,Fm,missValsin)
  }
  
  betaout <- NULL
  if(o_out >0 ){
    betaout <- matrix(0, nrow = o_out, ncol = nc+1 )
    for(i in 1:dim(Z_out)[2]){
      betaout[i,] <- coef(lm(Z_out[,i] ~ 1+Hm[,1]+Hm[,2]+Hm[,3]))
    }
    sigmaSqout <- SigmaSqCalc(Z_out,betaout,Hm,missValsout)
  }
  
  A <- 1 - (1 * as.matrix(as_adjacency_matrix(G)))
  Eneg  <- as_edgelist(graph_from_adjacency_matrix(A))
  
  ll <- binLL + Lg_cal(G = G, Ftot = Fm, Htot = Hm,Eneg = Eneg)
  
  # ll <- findLLDir(G, Ftot= Fm, Htot= Hm,Win,Wout,X_in,X_out,
  #                 betain,betaout,Z_in,Z_out,
  #                 sigmaSqin,sigmaSqout,alphaLL,Eneg)
  return(ll)
}

CoDA <- function(G,nc, k = c(0, 0) ,o = c(0, 0) , N,  alpha, lambda, thresh, 
                 nitermax, orig,randomize = TRUE, CovNamesLinin = c(),
                 CovNamesLinout = c(),CovNamesLPin = c(),CovNamesLPout = c(), 
                 alphaLL = NULL,test = TRUE,missing = 0) {
  if(printFlg == TRUE){
    print("In CoDa Func")
  }
  ## Community weights outgoing connections
  
  Ftot <- initWtmat(G,"all",N,nc)
  Finit <- Ftot
  
  ## Community weights outgoing connections
  Htot <- initWtmat(G,"all",N,nc)
  Hinit <- Htot
  
  iter <- 0
  k_in <- k[1]
  k_out <- k[2]
  o_in <- o[1]
  o_out <- o[2]
  
  ########## Get the covariate matrix for binary data and continuous data ##########
  b <- getNWcov(G, CovNamesLinin, CovNamesLinout,CovNamesLPin, CovNamesLPout,k_in,k_out,o_in,o_out)
  
  X_in <- b[[1]]
  X_out <- b[[2]]
  
  Z_in_orig <- b[[3]]
  Z_out_orig <- b[[4]]
  
  missValsin <- b[[5]]
  missValsout <- b[[6]]
  numMissValin <- b[[7]]
  numMissValout <- b[[8]]
  
  ## Degree condition for incoming and outgoing edges
  in_E <-  igraph::degree(G, mode = "in")
  out_E <- igraph::degree(G, mode = "out")
  
  ############ Initialize the coef matrix for covariates ##############
  ############## Setting initial values of binary covariates coefs weights ###############
  Win <- NULL#initW(k_in, nc,missValsin, X_in, Ftot, lambda, alpha)
  Win_orig <- Win
   
  Wout <- initW(k_out, nc,missValsout, X_out, Htot, lambda, alpha)
  Wout_orig <- Wout
  
  #################### Setting initial values of linear regression weights
  b <- initbeta(o_in, nc,missValsin, Z_in_orig, Ftot, alpha, N, lambda)
  betain <- b[[1]]
  sigmaSqin <- b[[2]]

  b <- initbeta(o_out, nc,missValsout, Z_out_orig, Htot, alpha,N, lambda)
  betaout <- b[[1]]
  sigmaSqout <- b[[2]]
  
  ## Getting the negative adjacency matrix
  Aneg <- 1 - (1 * as.matrix(as_adjacency_matrix(G)))
  Eneg  <- as_edgelist(graph_from_adjacency_matrix(Aneg))
  
  ##Initial log likelihood value
  LLold <- -10000
  # LLold <- findLLDir(G,Ftot,Htot, 
  #                    Win, Wout,X_in,X_out,
  #                    betain,betaout,Z_in,Z_out,
  #                    sigmaSqin, sigmaSqout,alphaLL, Eneg)
  lllst <- c()
  #lllst <- c(lllst, LLold)
  
  arilst <- c()
  OmegaVal <- c()
  mse <- matrix(nrow = 0, ncol = 2*(o_in+o_out))
  accuracy <- matrix(nrow = 0, ncol = k_in + k_out)
  
  delta <- getDelta(N)
  continue <- TRUE
  s <- 1:N
  
  while (continue) {
    
    ## Update the log likelihood.
    LLnew <- findLLDir(G,Ftot,Htot, 
                       Win,Wout,X_in,X_out, 
                       betain,betaout,Z_in_orig,Z_out_orig,
                       sigmaSqin, sigmaSqout,alphaLL, Eneg)
    #print(LLnew)
    #if(any(c(((LLnew - LLold) < thresh), (iter > nitermax)))){
    pctInc <- ((LLnew - LLold) / abs(LLold)) * 100
    #print(pctInc)
    #if (any((pctInc < thresh), (iter > nitermax))) {
    if(pctInc < thresh){
      print(paste0("The final percent change ", pctInc, " ,total iterations ", iter))
      break#continue <- FALSE
    }
    
    LLold <- LLnew
    lllst <- c(lllst, LLnew)
    iter <- iter + 1
    Ffin <- Ftot
    Hfin <- Htot
    Winfin <- Win
    Woutfin <- Wout
    betainfin <- betain
    betaoutfin <- betaout
    
    
    ## Updating the community weight parameter
    if (randomize == TRUE) {
      s <- sample(1:N, size = N, replace = FALSE)
    }
    ## Updating F i.e. all outgoing connections
    ## We look for neighbors our node is connected to with the edges directed outward from our node.
    Ftot <- updateWtmat(G,Ftot,Htot,"out",s,nc,
                        X_out,Z_in_orig,k_out,o_in,betain,Wout,
                        alphaLL,missValsin,alpha, start = 2, end = nc+1)
    
    ## Updating H i.e. all incoming connections
    ## We look for neighbors our node is connected to with the edges directed outward from our node.
    Htot <- updateWtmat(G,Htot,Ftot,"in",s,nc,
                        X_out,Z_out_orig,k_out,o_out,betaout,Wout,
                        alphaLL,missValsout,alpha, start = nc+2, end = nc*2+1)
    
    ## Updating logistic paramters
    ## Incoming correlated covariates
    # b1 <- updateLogisticParam(Win,X_in,Ftot, missValsin,lambda,alpha)
    # Win <- b1[[1]]
    # X_in <- b1[[2]]
    ## outgoing correlated covariates
    b2 <- updateLogisticParam(Wout,X_out,Ftot,Htot, missValsout, lambda,alpha)
    Wout <- b2[[1]]
    X_out <- b2[[2]]
    
    ### Updating the beta matrix for continuous covariates
    #print("Gradient for out")
    c1 <- updateLinearRegParam(betaout,missValsout,Z_out_orig,Htot,alpha,lambda,N)
    betaout <- c1[[1]]
    Z_out <- c1[[2]]
    sigmaSqout <- c1[[3]]
    
    #print("Gradient for in")
    c2 <- updateLinearRegParam(betain, missValsin, Z_in_orig, Ftot,alpha,lambda,N)
    betain <- c2[[1]]
    Z_in <- c2[[2]]
    sigmaSqin <- c2[[3]]
    
    ## Look for Communities based on the F matrix
    if(test == TRUE){
      mem <- rep(NA, N)
      for (i in 1:N) {
        m <- data.frame(com = rep(letters[1:nc], times = 2),
                        val = c(Ftot[i,], Htot[i,]))
        mem[i] <-  m$com[which.max(m$val)]
      }
      arilst <- c(arilst, mclust::adjustedRandIndex(orig, mem))
      OmegaVal <- c(OmegaVal, OmegaIdx(G, Ftot, Htot, N, delta,nc))
      
      if((o_in + o_out) >0){
        msein  <- 0
        mseout <- 0
        mseintot  <- 0
        mseouttot <- 0
        if (o_in > 0) {
          sbst <- (Z_in == 0)
          msein <- MSE(Ftot[!sbst,],Z_in[!sbst],betain,N)
          mseintot <- MSE(Ftot,Z_in,betain,N)
        }
        if (o_out > 0) {
          sbst <- (Z_out == 0)
          mseout <- MSE(Htot[!sbst,],Z_out[!sbst],betaout,N)
          mseouttot <- MSE(Htot,Z_out,betaout,N)
        }
        
        #print(paste0("msein ", msein," mseout ", mseout))
        mse <- rbind(mse, c(mseintot, mseouttot, msein, mseout))
      }
      #### Calculating the accuracy for the binary covariates
      #accuracy_in <- accuCalc(k = k_in,W = Win,Wtm = Ftot,X = X_in)
      accuracy_out <- accuCalc(k_out,Wout,cbind(Ftot,Htot),X_out)
      accuracy <- rbind(accuracy, cbind(accuracy_out))
    }
  }
  #print("Done with CoDA")
  
  return(
    list(
      Ffin = Ffin,
      Forig = Finit,
      Hfin = Hfin,
      Horig =  Hinit,
      Win_orig = Win_orig,
      Win_fin = Winfin,
      Wout_orig = Wout_orig,
      Wout_fin = Woutfin,
      Loglik = lllst,
      ARI = arilst,
      betain = betainfin,
      Zin_cov = Z_in,
      betaout = betaoutfin,
      Zout_cov = Z_out,
      MSE = mse,
      Acc = accuracy,
      OmegaIndex = OmegaVal
    )
  )
}


