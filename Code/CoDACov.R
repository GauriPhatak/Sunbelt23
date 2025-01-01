library(tidyverse)
library(dplyr)
library(igraph)
library(network)
library(intergraph)
library(ergm)
library(RColorBrewer)
#source("Code/HelperFuncs.R")
source(paste0(getwd(),"/Code/HelperFuncs.R"))

prob2logit <- function(x) {
  return(log(x / (1 - x)))
}

MSE <- function(Wtmat1,Wtmat2,Z,beta,N){
  #cbind(1, Wtmat)
  pred <- cbind(1, Wtmat1, Wtmat2) %*% t(beta)
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
  binVal_orig <- binVal
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
  
  return(list(G, binVal_orig))
}

SetContinuousCov <- function(G,o, N,nc, o_start, C, connType, dist,CovNamesLin,missing){
  ## Based on the mean and variance indicated by the dist matrix we assign continuous values to the network covariates.
  contVal <- matrix(ncol = o, nrow = N, 0) #rnorm(N,0,1)
  connType <- c(-1,1)
  s <- 1#c(1:(length(connType)*nc))
  for(k in 1:o){
    for(j in 1:length(connType)){
      for (i in 1:nc) {
        idx <- which(C[, i] == connType[j])
        contVal[idx, k] <- rnorm(length(idx), dist[[s]][1],dist[[s]][2])
        s <- s + 1
      }
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

genNested <- function(N,nc,pClust,k_in,k_out,o_in,o_out,pC,dist,covTypes,
                      CovNamesLPin,CovNamesLPout,CovNamesLinin,CovNamesLinout,
                      pConn, dir,dirPct, epsilon,missing, alpha, beta){
 
  ## Assigning the first network
  
   
}


genBipartite <- function(N,nc,pClust,k_in,k_out,o_in,o_out,pC,dist,covTypes,
                         CovNamesLPin,CovNamesLPout,CovNamesLinin,CovNamesLinout,
                         pConn, dir,dirPct, epsilon,missing, alpha, beta) {
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
  #Cold <- Cnew
  fin <- c()
  for (i in 2:(nc - 1)) {
    c <- t(combn(letters[1:nc], i))
    for (j in 1:dim(c)[1]) {
      com <- which(C %in% c[j, ])
      fin <- c(fin, com)
      print(c[j, ])
      print(pClust[nc + i + j - 2])
      idx <- sample(com , size = length(com) * pClust[nc + i + j - 2], replace = FALSE)
      print(idx)
      Cnew[idx, c[j, ]] <- 1
      
      ## sampled outgoing membership by percentage
      #sampOut <- sample(idx,
      #                  size = length(idx) * dirPctInt[j],
      #                  replace = FALSE)
      #Cnew[sampOut, c[j,1]] <- -1
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
  if (dir == "directed") {
    for (i in 1:nc) {
      ## sampled outgoing membership by percentage
      sampOut <- sample(which(Cnew[, i] == 1),
                        size = length(which(Cnew[, i] == 1)) * dirPct[i],
                        replace = FALSE)
      Cnew[sampOut, i] <- -1
    }
  }
  
  ## calculate the weight matrices
  F_u <- matrix(0,nrow = N, ncol = nc)
  H_u <- matrix(0,nrow = N, ncol = nc)
  if( dir == "directed" ){
    F_u[Cnew > 0] <- rbeta(sum(Cnew > 0),alpha,beta)
    H_u[Cnew < 0] <- rbeta(sum(Cnew < 0),alpha,beta)
  }else {
    F_u[Cnew > 0] <- rbeta(sum(Cnew > 0),alpha,beta)
    H_u <- F_u
  }
  
  A <- matrix(0, nrow = N, ncol = N)
  Apuv <- A
  for (i in 1:N) {
    for (j in 1:N) {
      ## probability of edge from node u to v
      #p_uv <- 1 - prod(1 - ((as.numeric(Cnew[i, ] - Cnew[j, ] > 1)) * pConn)) + epsilon
      p_uv <- 1 - exp(-1 * sum(F_u[i, ] * H_u[j, ])) + epsilon
      Apuv[i, j] <- p_uv
      A[i, j] <- rbinom(1,1,p_uv)#purrr::rbernoulli(1, p_uv)
    }
  }
  
  diag(A) <- 0
  g.sim <- graph_from_adjacency_matrix(A, mode = "directed")
  
  for (i in 1:nc) {
    g.sim <- g.sim %>%
      set_vertex_attr(name = letters[i], 
                      value = c(Cnew[, i]))
  }
  
  ## Creating binary covariates in the network.
  covVal_orig <- 0
  if (c("binary") %in% covTypes) {
    if (length(CovNamesLPin) > 0) {
      op <- SetBinarycov(g.sim, k_in, N, nc, k_start = 0, Cold, 
                            connType = 1,pC, CovNamesLPin, missing)
    }
    if (length(CovNamesLPout) > 0) {
      op <- SetBinarycov(g.sim, k = k_out, N, nc, k_start = k_in, Cold, 
                            connType = 1, pC, CovNamesLPout, missing)
    }
    g.sim <- op[[1]]
    covVal_orig <- op[[2]]
  }
  
  if (c("continuous") %in% covTypes) {
    if (length(CovNamesLinin) > 0) {
      op <- SetContinuousCov(g.sim, o_in, N, nc, o_start = 0, Cnew, 
                                connType = 1, dist,CovNamesLinin,missing)
    }
    if (length(CovNamesLinout) > 0) {
      op <- SetContinuousCov(g.sim, o_out, N,nc, o_start = o_in, Cnew, 
                                connType = 1, dist, CovNamesLinout, missing)
    }
    g.sim <- op[[1]]
    covVal_orig <- op[[2]]
  }
  
  g.sim <- g.sim  %>%
    set_vertex_attr(name = "Cluster", value = C)
  return(list(g.sim, covVal_orig, F_u, H_u))
}

sigmoid <- function(W, Ftotn) {
  Q <- 1 / (1 + exp(-1 * (W %*% t(Ftotn))))
  return(Q)
}

LRParamUpdt <- function(W, X, Ftot,Htot, alpha, lambda, missing, dir) {
  if(dir =="directed"){
    Ftotn <- cbind(1, Ftot, Htot)
  } else{
    Ftotn <- cbind(1, Ftot)
  }
  Q <- sigmoid(W, Ftotn)
  
  W_grad <- t(X - t(Q)) %*% Ftotn
  W_new <- W + alpha * (W_grad - lambda * (sign(W)))
  return(W_new)
}

PredictCovLR <- function(X, W, Fmat, Hmat, missing,dir) {
  if(dir =="directed"){
    Ftotn <- cbind(1, Fmat, Hmat)
  } else{
    Ftotn <- cbind(1, Fmat)
  }
  
  for (i in 1:dim(X)[2]) {
    X[which(missing[, i]), i] <- as.numeric(sigmoid(W[i, ], Ftotn[which(missing[, i]), ]) > 0.5)
  }
  return(X)
}

updateLogisticParam <- function(W,X,Wtm,Wtm2, missVals, lambda,alpha, dir){
  Wret <- W
  Xret <- X
  if (length(W) > 0) {
    if (sum(missVals) > 0) {
      Wret <-  LRParamUpdt(W, X, Wtm,Wtm2 ,alpha, lambda, missVals, dir)
      Xret <-  PredictCovLR(X, W, Wtm,Wtm2, missVals, dir)
    } else{
      Wret <-  LRParamUpdt(W, X, Wtm,Wtm2, alpha, lambda, missVals, dir)
    }
  }
  
  return(list(Wret,Xret))
}

SigmaSqCalc <- function(Z, beta, Ftot, Htot, missVals) {
  sigmaSq <- rep(0, dim(beta)[1])
  Fm <- cbind(1,Ftot, Htot)
  for (i in 1:dim(beta)[1]) {
    sigmaSq[i] <- sum((Z[!missVals[, i], i] - (t(beta[i, ]) %*% t(Fm[!missVals[, i], ]))) ^ 2, na.rm = TRUE) / sum(!missVals[, i])
  }
  return(sigmaSq)
}

PredictCovLin <- function(Ftot,Htot, Z, beta, missVals) {
  Fm <-  cbind(1, Ftot, Htot)#cbind(1,Ftot)
  for (i in 1:dim(Z)[2]) {
    if(sum(missVals) >1 ) {
      Z[which(missVals[, i]), i] <- Fm[which(missVals[, i]),] %*% as.matrix(beta[i, ])
      
    } else{
      Z[, i] <- Fm %*% as.matrix(beta[i,])
    }
  }
  return(Z)
}

LinRParamUpdt <- function(beta, Z, Fmat,Hmat, alpha, lambda, N, missVals) {
  if(printFlg == TRUE){
    print("In Linparamupdate Func")
  }
  beta_new <- matrix(0, nrow = dim(beta)[1], ncol = dim(beta)[2])
  Fm_n <- cbind(1, Fmat, Hmat)
  #for (i in 1:dim(beta)[1]) {
  #  beta_new <- t(MASS::ginv(t(Fm_n) %*% Fm_n) %*% (t(Fm_n) %*% Z))
  #print(beta_new)
  #}
  
  ## Calculating update to the beta using gradient descent
  y <-  PredictCovLin(Fmat,Hmat, Z, beta, missVals)
  gradient <- (t(y - Z) %*% Fm_n)/N
  #print(paste0(gradient))
  #alpha <- 0.01
  beta_new <-  beta - (alpha*(gradient- (lambda * sign(beta))))
  return(beta_new)
}

CalcQuk <- function(Fmat, Hmat, W, dir){
  
  if(dir == "directed"){
    Q_uk <- 1 / (1 + exp(-1 * (cbind(1, Fmat, Hmat) %*% t(W))))
    
  } else{
    Q_uk <- 1 / (1 + exp(-1 * (cbind(1, Fmat) %*% t(W))))
  }
  
  return(Q_uk)
}

Lx_cal <- function(W, N, Fmat, Hmat, X, dir){
  S <- 0
  
  ## Filter the NA values
  #Vals <- is.na(X)
  #Fmat <- Fmat[!naVals]
  
  if (length(W) > 0) {
    cases <- complete.cases(X)
    Q <- CalcQuk(Fmat, Hmat, W, dir)
    S <- sum((log(Q[cases, ]) * X[cases,]) + (log(1-Q[cases,]) * (1- X[cases,])), na.rm = TRUE)
  }
  
  return(S)
}

Lz_cal <- function(beta,N,Z,Fmat, Hmat,sigmaSq ){
  if(printFlg == TRUE){
    print("In Lz_Calc Func")
  }
  S <- 0
  if (length(beta) > 0) {
    for (j in 1:dim(beta)[1]) {
      S_tmp <- sum(( Z[,j]  - cbind(1,Fmat, Hmat) %*% beta[j, ]) ^ 2)
      #for (i in 1:N) {
      #c(1, Fmat[i, ])
      #  S_tmp <- sum(c(S_tmp, (Z[i, j] - sum(c(1, Fmat[i, ]) * beta[j, ], na.rm = TRUE)) ^ 2),na.rm = TRUE)
      #}
      S <- S - ((N * log(sigmaSq[j]) / 2) - (S_tmp / (2 * sigmaSq[j])))
    }
  }
  return(S)
}

Lg_cal <- function(G, Ftot, Htot,Eneg,epsilon, dir){
  scale <- 1
  E <- as_edgelist(G)
  Fv <- Ftot[E[, 1], ]
  Hv <- Htot[E[, 2], ]
  N <-  gorder(G)
  ## Variable to save the loglik contributed by nodes with edges between them
  ## Implementing rowsums Hopefully this is a little faster
  mulop <- log(1 - exp(-1 * rowSums(Fv * Hv)))
  S11 <- sum(is.finite(mulop) * mulop,
             na.rm = TRUE)
  
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
                      sigmaSqin = NA,sigmaSqout = NA, alphaLL, Eneg, dir, epsilon) {
  if(printFlg == TRUE){
    print("In Find LLDir Func")
  }
  
  N <-  gorder(G)
  
  ## Calculating the L_g part of the log likelihood for structure of the network
  S1 <- Lg_cal(G, Ftot, Htot, Eneg,epsilon, dir)
  
  ## Calculating the L_xin part of the log likelihood for binary covariates
  S2_in <- 0 #Lx_cal(Win, N, Ftot, X_in)
  
  ## Calculating the L_xout part of the log likelihood for binary covariates
  S2_out <- Lx_cal(W=Wout, N = N, Fmat= Ftot, Hmat =Htot, X=X_out, dir)
  
  ## adding the loglikelihood from the covariates
  S2 <- S2_in + S2_out
  
  ## Calculating the L_zin part of the log likelihood for continuous covariates
  S3_in <- 0# Lz_cal(beta = betain,N = N,Z = Z_in,Fmat = Ftot,sigmaSq = sigmaSqin )
  
  ## Calculating the L_zout part of the log likelihood for continuous covariates
  S3_out <- Lz_cal(beta = betaout,N = N,Z = Z_out,
                   Fmat = Ftot, Hmat = Htot,sigmaSq = sigmaSqout)
  S3 <- S3_in + S3_out
  #print(paste0(S1, " ", S2, " ", S3))
  
  ## Calculating the final log likelihood
  #if (is.null(alphaLL)) {
  ll <- S1 + (S2 + S3)
  #} else{
  #  ll <- alphaLL * S1 + (1 - alphaLL) * (S2 + S3)
  #}
  return(c(ll, S1, S2, S3))
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
  Wtmat <- matrix(nrow = N, ncol = nc, runif(nc * N, 0, 1))
  #Wtmat <- matrix(nrow = N, ncol = nc, 0.001)
  return(Wtmat)
}

initW <- function(k, nc,missVals, X, Fmat, lambda, alpha, ncoef){
  W <- NULL
  if (k > 0) {
    W <- matrix(nrow = k, ncol = (ncoef) + 1, 0)
    ###*** For future updates ***######
    #if (sum(missVals) > 0) {
    #  X_in <- PredictCovLR(X, W, Fmat, missVals)
    #}
    #Win <- LRParamUpdt(W, X, Fmat, alpha, lambda, missVals)
  }
  return(W)
}

initbeta <- function(o, ncoef,missVals, Z, Fmat, Hmat, alpha,N, lambda){
  beta <- NULL
  sigmaSq <- NULL
  if (o > 0) {
    beta <- matrix(nrow = o, ncol = (ncoef) + 1, 0)
    
    beta <- LinRParamUpdt(beta, Z, Fmat,Hmat, alpha, lambda, N, missVals)
    sigmaSq <- SigmaSqCalc(Z, beta, Fmat,Hmat, missVals)
    if (sum(missVals) > 0) {
      Z <- PredictCovLin(Fmat, Hmat, Z, beta, missVals)
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

BincovCmntyWtUpdt <- function(ncoef, f_u,h_u, W, X_u, dir){
  llX <- rep(0, ncoef + 1)
  if (length(W) > 0) {
    ## Updating the logistic regression weights
    ## Gradient function to update log likelihood for covariates
    ## adding intercept term
    Q_u <- CalcQuk(f_u,h_u,W, dir)
    ## This math has been verified
    llX <- t((X_u - Q_u) %*% W)
  }
  return(llX)
}

ContcovCmntyWtUpdt <- function(nc, Z_u, beta, f_u,h_u, sigmaSq){
  
  llZ <- rep(0, nc*2 + 1)
  f <- matrix(c(1, f_u,h_u), ncol = (nc*2) + 1)
  
  #print(Z_u - t(beta %*% t(f)))
  #print(t(f))
  if (length(beta) > 0) {
    ## Adding the gradient based on continuous covariates
    llZ <-  ((Z_u - t(beta %*% t(f))) / sigmaSq) %*% beta
  }
  
  return(llZ)
}

CmntyWtUpdt <- function(f_u ,h_u ,h_neigh ,h_sum ,X_u = NULL ,W = NULL ,Z_u = NULL,
                        beta = NULL ,sigmaSq ,alpha ,alphaLL ,nc ,start ,end ,mode,dir ) {
  if(printFlg == TRUE){
    print("In CmntyWtupdate Func")
  }
  tv <- 0 #0.000001
  
  ## Gradient function to update log likelihood of G
  #First part of Likelihood based on graph only.
  llG <- GraphComntyWtUpdt(f_u,h_u,h_neigh,h_sum)#Fvmat, Fvnot)
  
  #Second part of Likelihood based on binary covariates
  if(dir == "directed"){
    if(mode == "out"){
      llX <- BincovCmntyWtUpdt(nc*2, f_u, h_u, W, X_u, dir)
    }else if(mode == "in"){
      W <- cbind(W[,1], W[,start:end], W[,(start-nc):(end-nc)])
      llX <- BincovCmntyWtUpdt(nc*2, f_u, h_u, W, X_u, dir)
    }
  } else{
    llX <- BincovCmntyWtUpdt(nc, f_u, h_u, W, X_u, dir)
  }
  
  #Third part of likelihood based on continuous covariates
  if(dir == "directed"){
    if(mode == "out"){
      llZ <- ContcovCmntyWtUpdt(nc, Z_u, beta, f_u, h_u, sigmaSq)
    }else if(mode == "in"){
      #beta <- as.matrix(beta, ncol = nc*2+1)
      beta <- cbind(beta[,1], beta[,start:end], beta[,(start-nc):(end-nc)])
      llZ <- ContcovCmntyWtUpdt(nc, Z_u, beta, f_u,h_u, sigmaSq)
    }
  } else{
    llZ <- ContcovCmntyWtUpdt(nc, Z_u, beta, f_u, h_u, sigmaSq)
  }
  
  
  ## Function to update the community weights vector using non negative matrix factorization
  #if (is.null(alphaLL)) {
  #[2:(nc + 1)]

  f_u_new <- f_u + (alpha * (t(llG) + llX[2:(nc+1)] + llZ[2:(nc+1)]))
  # } else{
  #   f_u_new <- f_u + t((alpha * (
  #     alphaLL * llG + (1 - alphaLL) * (llX[2:(nc + 1)] + llZ[2:(nc + 1)])
  #   )))
  # }
  f_u_new[f_u_new < 0] <- tv
  return(f_u_new)
}

updateWtmat <- function(G,Wtm1,Wtm2,mode,s,nc,X,Z,k,o,beta,W,alphaLL,missVals,alpha, start, end, dir){
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
      if(mode == "in"){
        #beta <- as.matrix(beta, ncol = nc*2+1)
        #beta <- cbind(beta[,1], beta[,start:end], beta[,(start-nc):(end-nc)])
        sigmaSq <- SigmaSqCalc(Z, beta, Wtm2,Wtm1, missVals)
        
      }else{
        sigmaSq <- SigmaSqCalc(Z, beta, Wtm1,Wtm2, missVals)
        
      }
    }
    
    # Wtmat[i, ] <- CmntyWtUpdt(Wtm1_u,Wtm2_v,Wtm2_vnot,
    #                           X_u,W,Z_u, beta,
    #                           sigmaSq,alpha,alphaLL,nc)
    Wtmat[i, ] <- CmntyWtUpdt(Wtm1_u,Wtm2_u,Wtm2_neigh, Wtm2_sum,
                              X_u,W,Z_u, beta,sigmaSq,alpha,alphaLL,
                              nc, start, end, mode, dir)
    
  }
  
  return(Wtmat)
}

updateLinearRegParam <- function(beta,missVals,Z,Wtm1,Wtm2, alpha,lambda,N){
  if(printFlg == TRUE){
    print("In updateLinearRegParam Func")
  }
  betaret <- beta
  Zret <-  Z
  sigmaSq <- NULL
  if (length(beta) > 0) {
    betaret <- LinRParamUpdt(beta, Z, Wtm1,Wtm2, alpha, lambda, N, missVals)
    sigmaSq <-  SigmaSqCalc(Z, betaret, Wtm1,Wtm2, missVals)
    
    if (sum(missVals) > 0) {
      Zret <- PredictCovLin(Wtm1, Wtm2, Z, beta, missVals)
    }
  }
  return(list(betaret,Zret,sigmaSq))
}

accuCalc <- function(k,W,Wtm1,Wtm2,X, dir){
  
  if(printFlg == TRUE){
    print("In accuCalc Func")
  }
  if(dir == "directed"){
    Wtm <- cbind(Wtm1, Wtm2)
  } else{
    Wtm <- Wtm1
  }
  accuracy <- matrix(nrow = 0, ncol = k)
  if (k > 0) {
    predLR <- apply(sigmoid(W, cbind(1, Wtm)) > 0.5, 1, as.numeric)
    acc <- rep(0, k)
    for (i in 1:k) {
      cm <-  confusionMatrix(factor(predLR[, i], levels = c(0, 1)), factor(X[, i]))
      acc[i] <- round(cm$overall[1], 2)
    }
    accuracy <- rbind(accuracy,acc)
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

convertToUndir <- function(G, nc){
  G <- igraph::as_undirected(G, mode = "collapse")
  comAff <- as.data.frame(vertex_attr(G)) %>%
    dplyr::select(all_of(letters[1:nc]))
  for(i in 1:nc){
    comAff[ !(comAff[,i] == 0) ,i] <- 1
    G <- set_vertex_attr(G, letters[i], index= V(G), comAff[,i])
  }
  return(G)
}

GTLogLik <- function(G,nc,pConn,alphaLL,CovNamesLinin,CovNamesLinout, 
                     CovNamesLPin,CovNamesLPout,
                     k_in,k_out,o_in,o_out,epsilon, dir,
                     orig_Cov, Fm, Hm){
  
  ## Trying to use beta distribution instead of hardcoding the probability of connection
  #alpha_c <- -1*log((1-pConn))
  comAff <-  as.data.frame(vertex_attr(G)) %>%
    dplyr::select(all_of(letters[1:nc]))
  
  b <- getNWcov(G, CovNamesLinin, CovNamesLinout, 
                CovNamesLPin, CovNamesLPout, 
                k_in,k_out,o_in,o_out)
  
  X_in <- b[[1]]
  X_out <- orig_Cov#b[[2]]
  
  Z_in <- b[[3]]
  Z_out <- orig_Cov#b[[4]]
  
  missValsin <- b[[5]]
  missValsout <- b[[6]]
  numMissValin <- b[[7]]
  numMissValout <- b[[8]]
  
  Win <- NULL
  
  if(dir =="undirected"){
    ncoef <- nc
  } else{ 
    ncoef = nc*2
  }
  
  
  Wout <- NULL
  binLL <- 0
  if(k_out >0){
    
    Wout <- matrix(0, nrow = k_out, ncol = ncoef+1 )
    for(i in 1:dim(X_out)[2]){
      if(dir == "undirected"){
        mod <- glm(X_out[,i] ~ Fm[,1]+Fm[,2]+Fm[,3],family = "binomial")
      } else{
        mod <- glm(X_out[,i] ~ Fm[,1]+Fm[,2]+Fm[,3]+Hm[,1]+Hm[,2]+Hm[,3],family = "binomial")
      }
      #print(summary(mod))
      binLL <- binLL + as.numeric(logLik(mod))
      Wout[i,] <- coef(mod)
    }
  }
  
  
  betaout <- NULL
  contLL <- 0
  if(o_out >0 ){
    betaout <- matrix(0, nrow = o_out, ncol = (ncoef)+1 )
    for(i in 1:dim(Z_out)[2]){
      if(dir =="undirected"){
        mod <- glm(Z_out[,i] ~ Fm[,1]+Fm[,2]+Fm[,3],family = "binomial")
      } else{
        mod <- glm(Z_out[,i] ~ Fm[,1]+Fm[,2]+Fm[,3]+Hm[,1]+Hm[,2]+Hm[,3],family = "binomial")
      }      
      betaout[i,] <- coef(mod)
      contLL <- contLL + as.numeric(logLik(mod))
      
    }
    #sigmaSqout <- sigma(mod)^2#SigmaSqCalc(Z_out,betaout,Fm,Hm,missValsout)
  }
  
  A <- 1 - (1 * as.matrix(as_adjacency_matrix(G)))
  Eneg  <- as_edgelist(igraph::graph_from_adjacency_matrix(A, mode = dir))
  
  GrphLL <- Lg_cal(G = G, Ftot = Fm, Htot = Hm,Eneg = Eneg, epsilon = epsilon, dir = dir)
  ll <- binLL + contLL + GrphLL
  
  findLLDir(G, Ftot = Fm, Htot= Hm,Win,Wout,X_in,X_out,
            betain,betaout,Z_in,Z_out,
            sigmaSqin,sigmaSqout,alphaLL,Eneg, 
            dir, epsilon)
  #Lx_cal(W = Wout, N=100, Fmat = Fm, Hmat = Hm, X = X_out, dir= "directed")
  return(list(ll, c(binLL,contLL,GrphLL)))
}

CoDA <- function(G,nc, k = c(0, 0) ,o = c(0, 0) , N,  alpha, lambda, thresh, 
                 nitermax, orig,randomize = TRUE, CovNamesLinin = c(),
                 CovNamesLinout = c(),CovNamesLPin = c(),CovNamesLPout = c(), dir, 
                 alphaLL = NULL,test = TRUE,missing = 0, covOrig, epsilon) {
  if(printFlg == TRUE){
    print("In CoDA Func")
  }
  
  if(dir == "directed"){
    ncoef = nc*2
  }else{ncoef = nc}
  
  ## Community weights outgoing connections
  
  Ftot <- initWtmat(G,"all",N,nc)
  Finit <- Ftot
  
  ## Community weights outgoing connections
  if(dir == "directed"){
    Htot <- initWtmat(G,"all",N,nc)
    Hinit <- Htot
  } else{
    Htot <- Ftot
    Hinit <- Htot
  }
  
  
  iter <- 0
  k_in <- k[1]
  k_out <- k[2]
  o_in <- o[1]
  o_out <- o[2]
  
  ########## Get the covariate matrix for binary data and continuous data ##########
  b <- getNWcov(G, CovNamesLinin,CovNamesLinout,
                CovNamesLPin,CovNamesLPout,
                k_in,k_out,o_in,o_out)
  
  X_in <- b[[1]]
  X_out <- b[[2]]
  
  Z_in_orig <- b[[3]]
  Z_out_orig <- b[[4]]
  
  missValsin <- b[[5]]
  missValsout <- b[[6]]
  numMissValin <- b[[7]]
  numMissValout <- b[[8]]
  
  ## For covariates where the value is missing, fill it in initially with random values
  if(k_out >0){
    for(i in 1:dim(X_out)[2]){
      X_out[missValsout[,i],i] <- rbinom(n = sum(missValsout[,i]), 
                                         size = 1, prob = 0.5)
    }
  }
  
  ## Degree condition for incoming and outgoing edges
  in_E <-  igraph::degree(G, mode = "in")
  out_E <- igraph::degree(G, mode = "out")
  
  ############ Initialize the coef matrix for covariates ##############
  ############## Setting initial values of binary covariates coefs weights ###############
  Win <- NULL#initW(k_in, nc,missValsin, X_in, Ftot, lambda, alpha)
  Win_orig <- Win
   
  Wout <- initW(k_out, nc,missValsout, X_out, Htot, lambda, alpha, ncoef)
  Wout_orig <- Wout
  
  #################### Setting initial values of linear regression weights
  b <- NULL#initbeta(o_in, nc,missValsin, Z_in_orig, Ftot, alpha, N, lambda)
  betain <- b#b[[1]]
  sigmaSqin <- b#b[[2]]

  b <- initbeta(o_out, ncoef,missValsout, Z_out_orig, Ftot, Htot, alpha,N, lambda)
  betaout <- b[[1]]
  sigmaSqout <- b[[2]]
  
  ## Getting the negative adjacency matrix
  Aneg <- 1 - (1 * as.matrix(as_adjacency_matrix(G)))
  Eneg  <- as_edgelist(graph_from_adjacency_matrix(Aneg, mode= dir))
  
  ##Initial log likelihood value
  LLold <- -10000
  # LLold <- findLLDir(G,Ftot,Htot, 
  #                    Win, Wout,X_in,X_out,
  #                    betain,betaout,Z_in,Z_out,
  #                    sigmaSqin, sigmaSqout,alphaLL, Eneg)
  lllst <- matrix(nrow = 0, ncol = 4)
  #lllst <- c(lllst, LLold)
  
  arilst <- c()
  OmegaVal <- c()
  mse <- matrix(nrow = 0, ncol = (o_in+o_out))
  accuracy <- matrix(nrow = 0, ncol = k_in + k_out)
  accuracyMD <- matrix(nrow = 0, ncol = k_in + k_out)
  
  delta <- getDelta(N)
  continue <- TRUE
  s <- 1:N
  
  while (continue) {
    
    ## Update the log likelihood.
    LLvec <- findLLDir(G,Ftot,Htot, 
                       Win,Wout,X_in,X_out, 
                       betain,betaout,Z_in_orig,Z_out_orig,
                       sigmaSqin, sigmaSqout,alphaLL, 
                       Eneg, dir, epsilon)
    LLnew <- LLvec[1]
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
    lllst <- rbind(lllst, LLvec)
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
    
    if(dir =="directed"){
      ## Updating F i.e. all outgoing connections
      ## We look for neighbors our node is connected to with the edges directed outward from our node.
      Ftot <- updateWtmat(G,Ftot,Htot,"out",s,nc,
                          X_out,Z_in_orig,k_out,o_in,betain,Wout,
                          alphaLL,missValsin,alpha, start = 2, end = nc + 1, dir)
      
      
      ## Updating H i.e. all incoming connections
      ## We look for neighbors our node is connected to with the edges directed outward from our node.
      Htot <- updateWtmat(G,Htot,Ftot,"in",s,nc,
                          X_out,Z_out_orig,k_out,o_out,betaout,Wout,
                          alphaLL,missValsout,alpha, start = nc+2, end = nc*2+1, dir)
    } else{
      Ftot <- updateWtmat(G,Ftot,Htot,"all",s,nc,
                          X_out,Z_in_orig,k_out,o_in,betain,Wout,
                          alphaLL,missValsin,alpha, start = 2, end = nc+1, dir)
      Htot <- Ftot
    }
    
    
    ## Updating logistic paramters
    ## Incoming correlated covariates
    # b1 <- updateLogisticParam(Win,X_in,Ftot, missValsin,lambda,alpha)
    # Win <- b1[[1]]
    # X_in <- b1[[2]]
    ## outgoing correlated covariates
    b2 <- updateLogisticParam(W = Wout,X = X_out , Wtm = Ftot ,Wtm2 = Htot,missVals = missValsout,
                              lambda = lambda,alpha = alpha, dir =dir)
    Wout <- b2[[1]]
    X_out <- b2[[2]]
    
    ### Updating the beta matrix for continuous covariates
    #print("Gradient for out")
    c1 <- updateLinearRegParam(betaout,missValsout,Z_out_orig,Ftot,Htot,alpha,lambda,N)
    betaout <- c1[[1]]
    Z_out <- c1[[2]]
    sigmaSqout <- c1[[3]]
    
    #print("Gradient for in")
    # c2 <- updateLinearRegParam(betain, missValsin, Z_in_orig, Ftot,alpha,lambda,N)
    # betain <- c2[[1]]
     Z_in <- 0#c2[[2]]
    # sigmaSqin <- c2[[3]]
    
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
      
      if((o_in + o_out) > 0){
        #msein  <- 0
        mseout <- 0
        #mseintot  <- 0
        mseouttot <- 0
        # if (o_in > 0) {
        #   sbst <- (Z_in == 0)
        #   msein <- MSE(Ftot[!sbst,],Htot[!sbst,],Z_in[!sbst],betain,N)
        #   mseintot <- MSE(Ftot,Htot,Z_in,betain,N)
        # }
        if (o_out > 0) {
          #sbst <- (Z_out == 0)
          #mseout <- MSE(Ftot[!sbst], Htot[!sbst,],Z_out[!sbst],betaout,N)
          mseouttot <- MSE(Ftot,Htot,Z_out,betaout,N)
        }
        
        #print(paste0("msein ", msein," mseout ", mseout))
        mse <- rbind(mse, c( mseouttot))
      }
      
      if((k_in + k_out) > 0){
        #### Calculating the accuracy for the binary covariates
        
        ## Overall Accuracy
        accuracy_out <- accuCalc(k_out,Wout,Ftot,Htot,covOrig,dir)
        ## Missing data prediction accuracy
        #accuracy_outMD <- c(rep(0,  k_in + k_out))
        if(missing > 0){
        accuracy_outMD <- accuCalc(k_out,Wout,
                                   Ftot[as.logical(rowSums(missValsout)),],
                                   Htot[as.logical(rowSums(missValsout)),],
                                   covOrig[as.logical(rowSums(missValsout)),],dir)
        accuracyMD <- rbind(accuracyMD, cbind(accuracy_outMD))
        }
        accuracy <- rbind(accuracy, cbind(accuracy_out))
      }
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
      OmegaIndex = OmegaVal,
      AccMD = accuracyMD
    )
  )
}


