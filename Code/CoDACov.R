library(tidyverse)
library(dplyr)
library(igraph)
library(glmnet)
library(RColorBrewer)

source(paste0(getwd(),"/Code/NetworkMetrics.R"))
source(paste0(getwd(),"/Code/HelperFuncs.R"))


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
  missValsout_bin <- 0
  numMissValout_bin <- 0
  missValsout_cont <- 0
  numMissValout_cont <- 0
  
  if (k_in > 0) {
    X_in <- as.matrix(covtmp %>%
                        dplyr::select(all_of(CovNamesLPin)))
    missValsin <- is.na(X_in)
    numMissValin <- sum(missValsin)
  }
  if (k_out > 0) {
    X_out <- as.matrix(covtmp %>%
                         dplyr::select(all_of(CovNamesLPout)))
    missValsout_bin <- is.na(X_out)
    numMissValout_bin <- sum(missValsout_bin)
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
    missValsout_cont <- is.na(Z_out)
    numMissValout_cont <- sum(missValsout_cont)
  }
  
  return(list(X_in, X_out, Z_in, Z_out, 
              missValsin, missValsout_bin,missValsout_cont, 
              numMissValin, numMissValout_bin, numMissValout_cont))
}

convertToUndir <- function(G, nc){
  G <- igraph::as_undirected(G, mode = "collapse")
  ColNames <- make_letter_names(nc)
  comAff <- as.data.frame(vertex_attr(G)) %>%
    dplyr::select(any_of(ColNames))
  for(i in 1:dim(comAff)[2]){
    comAff[ !(comAff[,i] == 0) ,i] <- 1
    G <- set_vertex_attr(G, ColNames[i], index = V(G), comAff[,i])
  }
  return(G)
}

GTLogLik <- function(G,nc,pConn,alphaLL,CovNamesLinin,CovNamesLinout, 
                     CovNamesLPin,CovNamesLPout,
                     k_in,k_out,o_in,o_out,epsilon, dir,
                     orig_Cov, Fm, Hm){
  
  ## Trying to use beta distribution instead of hardcoding the probability of connection
  ColNames <- make_letter_names(nc)
  comAff <-  as.data.frame(vertex_attr(G)) %>%
    dplyr::select(any_of(ColNames))
  
  b <- getNWcov(G, CovNamesLinin, CovNamesLinout, 
                CovNamesLPin, CovNamesLPout, 
                k_in,k_out,o_in,o_out)
  #X_in, X_out, Z_in, Z_out, missValsin, missValsout_bin, missValsout_cont, numMissValin, numMissValout_bin, numMissValout_cont
  X_in <- b[[1]]
  X_out <- orig_Cov
  
  Z_in <- b[[3]]
  Z_out <- orig_Cov
  
  missValsin <- b[[5]]
  missValsout_bin <- b[[6]]
  missValsout_cont <- b[[7]]
  numMissValin <- b[[8]]
  numMissValout_bin  <- b[[9]]
  numMissValout_cont <- b[[10]]
  
  Win <- NULL
  VIFval <- 0
  
  if(dir =="undirected"){
    ncoef <- nc
  } else{ 
    ncoef = nc*2
  }
  
  ## Calculating the VIF value
  Wout <- NULL
  binLL <- 0
  if(k_out >0){
    
    Wout <- matrix(0, nrow = k_out, ncol = ncoef+1 )
    for(i in 1:dim(X_out)[2]){
      if(dir == "undirected"){
        mod <- glm(X_out[,i] ~ Fm,family = "binomial")
      } else{
        mod <- glm(X_out[,i] ~ Fm + Hm,family = "binomial")
      }
      binLL <- binLL + as.numeric(logLik(mod))
      Wout[i,] <- coef(mod)
    }
    
  }
  
  betaout <- NULL
  contLL <- 0
  sigmaSqout <- rep(0, o_out)
  if(o_out >0 ){
    betaout <- matrix(0, nrow = o_out, ncol = (ncoef)+1 )
    for(i in 1:dim(Z_out)[2]){
      if(dir =="undirected"){
        mod <- lm(Z_out[,i] ~ Fm)
      }else{
        mod <- lm(Z_out[,i] ~ Fm + Hm )
      }
      
      betaout[i,] <- coef(mod)
      contLL <- contLL + as.numeric(logLik(mod))
      sigmaSqout[i] <- sigma(mod)^2
    }
  }
  
  A <- 1 - (1 * as.matrix(as_adjacency_matrix(G)))
  Eneg  <- as_edgelist(igraph::graph_from_adjacency_matrix(A, mode = dir))
  E <- igraph::as_edgelist(G)
  GrphLL <- Lg_cal(G = G, Ftot = Fm, Htot = Hm,Eneg = Eneg, epsilon = epsilon, dir = dir, E)#(G, Ftot, Htot,Eneg,epsilon, dir, E)
  ll <- binLL + contLL + GrphLL
  
  return(list(ll, c(binLL,contLL,GrphLL), VIFval))
}

prob2logit <- function(x) {
  return(log(x / (1 - x)))
}

flip <- function(p) {
  return(sample(
    x = c(T, F),
    size = 1,
    prob = c(p, 1 - p)
  ))
}

mode <- function(x) {
  which.max(tabulate(x))
}

## Network generation code
SetBinarycov <- function(G, k, N,nc, k_start, Cnew,connType, pC,
                         CovNamesLP, missing,missType) {
  ## Based on the probability matrix pC assign indicates the binary assignment of covariates to a community
  ## setting covariates for outgoing communities. hence Cnew will be -1
  binVal <- matrix(ncol = k, nrow = N, 0)
  
  ## generate community to covariate correlation assignments.
  if(nc < k){
    assgnments <- c( sample(1:nc), sample(1:nc, k - nc, replace = TRUE))
  }else{
    assgnments <- sample(1:nc, k,replace  = FALSE)
  }
  
  for (i in 1:k) {
  #  if(i <= dim(pC)[1]){
    j <- assgnments[i]
        idx <- which(!(Cnew[, j] == 0))
        binVal[idx, i] <- rbinom(length(idx), 1, pC[i, 1])
        idx_n <- which(Cnew[, j] == 0)
        binVal[idx_n, i] <- rbinom(length(idx_n), 1, pC[i, 2])
        
  #  }
  }
 
  binVal_orig <- binVal
  if(!is.null(missing)) {
    for (j in 1:k) {
      if(missType[j] == "Random"){ # 1) Random (MCAR)
        binVal[sample(1:N, round(missing[j] * N / 100)), j] <- NA
      }else if(missType[j] == "CovDep_MAR"){ # 4) if the value is 1 then we have a certain percent missing. No missing in 0
        idx <- which(binVal[, j] == 1)
        idx <- sample(idx, round(missing[j] * length(idx) / 100) )
        binVal[idx, j] <- NA
      }else if(missType[j] == "CommDep_MAR"){
        binVal[sample(1:N, round(missing[j] * N / 100)), j] <- NA
      }
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

SetContinuousCov <- function(G,o, N,nc, o_start, C, connType, 
                             dist,CovNamesLin,missing,missType, MARparam=NULL, k){
  ## Based on the mean and variance indicated by the dist matrix we assign continuous values to the network covariates.
  #print("Stop here")
  contVal <- matrix(ncol = o, nrow = N, 0) 
  s <- 1 
  
  ## generate community to covariate correlation assignments.
  if(nc < o){
    assgnments <- c( sample(1:nc), sample(1:nc, o - nc, replace = TRUE))
  }else{
    assgnments <- sample(1:nc, o,replace  = FALSE)
  }
  
  ## for each community set the covariate values. 
  for (i in 1:o) {
    j <- assgnments[i]
    if(s < length(dist)){
      ## set incoming community covariate
      idx_in <- which(C[,j] == -1)
      contVal[idx_in, i] <- rnorm(length(idx_in), dist[[s]][1], dist[[s]][2])
      #Set outgoing community covariate
      idx_out <- which(C[,j] == 1)
      contVal[idx_out, i] <- rnorm(length(idx_out), dist[[s+1]][1], dist[[s+1]][2])
      ## setting covariates for the non community nodes
      idx_n <- which(!(1:N %in% c(idx_in, idx_out)))
      contVal[idx_n, i] <- rnorm(length(idx_n), dist[[s+2]][1], dist[[s+2]][2])
      s <- s+3
    }
  }
  
  contVal_orig <- contVal
  if(!is.null(missing)) {
    for (j in 1:o) {
      ## Create missingness based on type of missingness. 
      if(missType[j] == "Random"){ # 1) Random (MCAR)
        contVal[sample(1:N, round(missing[j] * N / 100)), j] <- NA
      }else if(missType[j] == "GT_MAR"){ # 2) Greater Than a particular value (MAR)
        sampTot <- which(floor(contVal[,j]) > MARparam[1]) ## Indices to be sampled from based on restriction
        contVal[sample(sampTot,round(missing[j]*MARparam[2]/100)),j] <- NA
        contVal[sample(which(!((1:N) %in% sampTot)), round(missing[j]*(100-MARparam[2])/100)), j] <- NA
      }else if(missType[j] == "LT_MAR"){ # 3) Less than a particular value (MAR)
        sampTot <- which(contVal[,j] < MARparam[1]) ## Indices to be sampled from based on restriction
        contVal[sample(sampTot,round(missing[j]*MARparam[2]/100)),j] <- NA
        contVal[sample(which(!((1:N) %in% sampTot)), round(missing[j]*(100-MARparam[2])/100)), j] <- NA
      }else if(missType[j] == "CovDep_MAR"){ # 4) Dependent on another covariate that is not missing (MAR)
        
      }
      
      G <- G %>% set_vertex_attr(name = CovNamesLin[j], value = c(contVal[, j]))
    }
  } else{
    for (j in 1:o) {
      G <- G %>% set_vertex_attr(name = CovNamesLin[j], value = c(contVal[, j]))
    }
  }
  
  return(list(G, contVal_orig))
  
}

genBipartite <- function(N,nc,pClust,k_in,k_out,o_in,o_out,pC,dist,covTypes,
                         CovNamesLPin,CovNamesLPout,CovNamesLinin,CovNamesLinout,
                         pConn,dir,dirPct,epsilon,missing,a,b,Type,pClustOL,
                         missType= rep("Random", sum(k_in,k_out,o_in,o_out)), 
                         MARparam=NULL, seed){
  
  set.seed(seed)
  C <- 0
  ## make column names
  ColNames <- make_letter_names(nc)
  while (length(table(C)) < nc) {
    C <- sample(
      x = ColNames, #letters[1:nc],
      size = N,
      replace = TRUE,
      prob = pClust[1:nc]
    )
  }
  op <- model.matrix( ~ C - 1)
  colnames(op) <- ColNames #letters[1:nc]
  
  ## overlapping various combinations of communities.
  Cnew <- as.data.frame(op, check.names = FALSE)
  fin <- c()
  #3 Right now just set this to do 2 group overlap . need to add 3 group and higher later.
  for (i in 2) {
    c <- t(combn(ColNames, i))
    if(Type == "Nested"){ ## For nested network currently only 3 communities are supported.
      l <- 1
      for(j in 1:dim(c)[1]){
        com <- which(C %in% c[j, ])
        fin <- c(fin, com)
        idx <- sample(com , size = length(com) * pClustOL[l], replace = FALSE)
        Cnew[idx, c[j,1 ]] <- 1
        l <- l + 1
      }
    }else if(Type == "2-Mode"){ ## For 2-Mode network currently only 3 communities are supported.
      ## For this type there should be two cohesive communities and one community that is only outgoing
      l <- 1
      for (j in seq(from =1, to = (dim(c)[1]*2), by = 2)){
        com1 <- which(C %in% c[l,1]) ## Finding the original assigned nodes to community 1
        com2 <- which(C %in% c[l,2]) ## Finding the original assigned nodes to community 2
        
        ## Number of nodes to be assigned to community 2 from community 1
        idx1 <- sample(com1 , size = length(com1) * pClustOL[j], replace = FALSE)
        
        ## Number of nodes to be assigned to community 1 from community 2
        idx2 <- sample(com2 , size = length(com2) * pClustOL[j+1], replace = FALSE)
        
        Cnew[idx1, c[l,2 ]] <- 1 ## Assiging values in community 2 based on overlap with community 1
        Cnew[idx2, c[l,1 ]] <- 1 ## Assiging values in community 1 based on overlap with community 0
        l <- l + 1
      }
    }else if(Type == "Cohesive"){
      l <- 1
      id <- which(pClustOL != 0)
      for(j in id){#1:dim(c)[1]){
        com <- which(C %in% c[j, ])
        fin <- c(fin, com)
        idx <- sample(com , size = length(com) * pClustOL[j], replace = FALSE)
        Cnew[idx, c[j, ]] <- 1
      }
    }
  }
  ## percentage overlapping all the communities
  idx <- sample(unique(fin),
                size = length(unique(fin)) * pClustOL[length(pClustOL)],
                replace = FALSE)
  Cnew[idx, ] <- 1
  
  ## Determine the direction of the community membership. i.e. if it is is incoming or outgoing community affiliation.
  ## -1 for outgoing membership and +1 for incoming membership
  
  Cold <- Cnew
  if (dir == "directed") {
    if(Type == "2-Mode"){
      c <- t(combn(ColNames, 2))
      idx <- which(Cnew[,c[1,1]] == 1 & Cnew[,c[1,2]] == 1)
      Cnew[idx, c[1,1]] <- -1
      sampOut <- sample( which(Cnew[,c[1,2]] == 1 ) ,
                         size = length(which(Cnew[, c[1,2]] == 1)) * dirPct[1], 
                         replace = FALSE)
      Cnew[sampOut, c[1,2]] <- -1
      
      sampOut <- sample( which(Cnew[,c[2,2]] == 1 ) ,
                         size = length(which(Cnew[, c[2,2]] == 1)) * dirPct[2], 
                         replace = FALSE)
      Cnew[sampOut, c[2,2]] <- -1
      
    } else{
      for (i in 1:nc) {
        ## sampled outgoing membership by percentage
        sampOut <- sample(which(Cnew[, i] == 1),
                          size = length(which(Cnew[, i] == 1)) * dirPct[i],
                          replace = FALSE)
        Cnew[sampOut, i] <- -1
      }
    }
  }
  
  ## calculate the weight matrices
  F_u <- matrix(0,nrow = N, ncol = nc)
  H_u <- matrix(0,nrow = N, ncol = nc)
  if( dir == "directed" ){
    F_u[Cnew > 0] <- rbeta(sum(Cnew > 0),a,b)
    H_u[Cnew < 0] <- rbeta(sum(Cnew < 0),a,b)
  }else {
    F_u[Cnew > 0] <- 0.4
    H_u <- F_u
  }
  
  A <- matrix(0, nrow = N, ncol = N)
  Apuv <- A
  for (i in 1:N) {
    for (j in 1:N) {
      ## probability of edge from node u to v
      p_uv <- 1 - exp(-1 * sum(F_u[i, ] * H_u[j, ])) + epsilon
      Apuv[i, j] <- p_uv
      A[i, j] <- rbinom(1,1,p_uv)
    }
  }
  
  diag(A) <- 0
  g.sim <- graph_from_adjacency_matrix(A, mode = "directed")
  for (i in 1:nc) {
    g.sim <- g.sim %>%
      set_vertex_attr(name = ColNames[i],#letters[i], 
                      value = c(Cnew[, i]))
  }
  
  ## Setting covariates in the network
  covVal_orig <- 0
  covVal_orig_bin<-0
  covVal_orig_cont <- 0
  ## Creating binary covariates in the network.
  if (c("binary") %in% covTypes) {
    if (length(CovNamesLPin) > 0) {
      op <- SetBinarycov(g.sim,k_in,N,nc,k_start=0,Cnew,connType=1,pC,
                         CovNamesLPin,missing[1:k_in],missType[1:k_out])
    }
    if (length(CovNamesLPout) > 0) {
      op <- SetBinarycov(g.sim,k=k_out,N,nc,k_start=k_in,Cnew,connType=1,pC,
                         CovNamesLPout,missing[1:k_out],missType[1:k_out])
      
    }
    g.sim <- op[[1]]
    covVal_orig_bin <- op[[2]]
  }
  
  if (c("continuous") %in% covTypes) {
    if (length(CovNamesLinin) > 0) {
      op <- SetContinuousCov(g.sim, o_in, N, nc, o_start = 0, Cnew, 
                             connType = 1, dist,CovNamesLinin,missing[(k_out+1):(o_out+k_out)],
                             missType[(k_out+1):(o_out+k_out)], MARparam, k_out)
    }
    if (length(CovNamesLinout) > 0) {
      
      op <- SetContinuousCov(g.sim, o_out, N,nc, o_start = o_in, Cnew, 
                             connType = 1, dist, CovNamesLinout, missing[(k_out+1):(o_out+k_out)],
                             missType[(k_out+1):(o_out+k_out)], MARparam,k_out)
    }
    g.sim <- op[[1]]
    covVal_orig_cont <- op[[2]]
  }
  g.sim <- g.sim  %>%
    set_vertex_attr(name = "Cluster", value = C)
  return(list(g.sim, list(covVal_orig_bin, covVal_orig_cont), F_u, H_u, Apuv))
}

updateWtmat <- function(G, Ftot, Htot, mode, s, nc, X, Z, k, o, beta, W, alphaLL, 
                        missVals_bin,missVals_cont , alpha, start, end, dir, 
                        inNeigh, outNeigh, lambda_bin, lambda_lin, penalty, lambda_grph){
  if(printFlg == TRUE){
    print("in updateWtmat")
  }
  ## Getting sum of all the weights of Fv. Have to modify this for CoDa as opposed to CESNA 
  
  # calculate the sigma sq value
  if (o > 0) {
    sigmaSq <- SigmaSqCalc(Z, beta, Ftot ,Htot, missVals_cont,dir)
  }
  
  ## Calculating the sum at the beginning of the loop
  Htot_sum <- as.matrix(colSums(Htot))
  Ftot_sum <- as.matrix(colSums(Ftot))
  for (i in s) {
    ## If there are covariates set the variables here
   
    X_u <- NULL
    if (k > 0) {
      X_u <- matrix(X[i, ], ncol = k)
    }
    
    Z_u <- NULL
    if(o > 0 ){
      Z_u <- matrix(Z[i, ], ncol = o)
    }
    
    if(dir == "directed"){
      mode = "in"
    }else{
      mode = "all"
    }
    neigh_Ftot <- inNeigh[[i]] 
    
    ## Htot neighbours
    Htot_neigh <- matrix(Htot[neigh_Ftot, ], ncol = nc)
    
    Ftot[i, ] <- CmntyWtUpdt(matrix(Ftot[i,], ncol = nc), matrix(Htot[i,], ncol = nc),
                             Htot_neigh, Htot_sum,X_u,W,Z_u, beta,sigmaSq,alpha,
                             alphaLL,nc, 2, nc+1,mode, dir,lambda_bin, lambda_lin, penalty,lambda_grph)
    if(dir == "directed"){
      ## Ftot neighbours
      neigh_Htot <- outNeigh[[i]] 
      
      Ftot_neigh <- matrix(Ftot[neigh_Htot, ], ncol = nc)
      
      Htot[i, ] <- CmntyWtUpdt(matrix(Htot[i,],ncol = nc), matrix(Ftot[i,],ncol = nc),
                               Ftot_neigh,Ftot_sum,X_u,W,Z_u,beta,sigmaSq,alpha,
                               alphaLL,nc, nc + 2, (nc*2)+1,"out",dir,lambda_bin, lambda_lin, penalty,lambda_grph)
    }else{
      Htot[i, ] <- Ftot[i, ]
    }
    
  }
 
  return(list(Ftot,Htot))
}



CmntyWtUpdt <- function(W1_u ,W2_u ,W2_neigh ,W2_sum ,X_u = NULL ,W = NULL ,
                        Z_u = NULL, beta = NULL ,sigmaSq ,alpha ,alphaLL ,
                        nc ,start ,end ,mode,dir , lambda_bin, lambda_lin, penalty, lambda_grph ) {
  
  llG <- rep(0,nc)
  
  ## Gradient function to update log likelihood of G
  #First part of Likelihood based on graph only.
  llG <- GraphComntyWtUpdt(W1_u,W2_u,W2_neigh,W2_sum)
  
  #Second part of Likelihood based on binary covariates
  if(dir == "directed"){
    if(mode == "in"){
      llX <- BincovCmntyWtUpdt(nc, W1_u, W2_u, W, X_u, dir)
    }else if(mode == "out"){
      W <- cbind(W[,1], W[,start:end], W[,(start-nc):(end-nc)])
      llX <- BincovCmntyWtUpdt(nc, W1_u, W2_u, W, X_u, dir)
    }
  }else{
    llX <- BincovCmntyWtUpdt(nc, W1_u, W2_u, W, X_u, dir)
  }
  #llX <- BincovCmntyWtUpdt(nc, W1_u, W2_u, W, X_u, dir)
  
  #Third part of likelihood based on continuous covariates
  if(dir == "directed"){
    if(mode == "in"){
      llZ <- ContcovCmntyWtUpdt(nc, Z_u, beta, W1_u, W2_u, sigmaSq, dir)
    }else if(mode == "out"){
      beta <- cbind(beta[,1], beta[,start:end], beta[,(start-nc):(end-nc)])
      llZ <- ContcovCmntyWtUpdt(nc, Z_u, beta, W1_u,W2_u, sigmaSq, dir)
    }
  }else{
    llZ <- ContcovCmntyWtUpdt(nc, Z_u, beta, W1_u, W2_u, sigmaSq, dir)
  }
  
  ## Function to update the community weights vector using non negative matrix factorization
  #Updated value
  UpVal <- (c(llG) +  alphaLL *(llX[1:nc] + llZ[1:nc]))
  ## experimenting with adding L1 enalty to the weight matrix
   #L1_Wtpenalty <- 0.00005 lambda_grph
  #  for(i in 1:length(UpVal)){
  #   UpVal[i] <- sign(UpVal[i]) * max(abs(UpVal[i]) - lambda_grph, 0 )
  # }
  W1_u_new <- W1_u + (alpha * UpVal)
  
  W1_u_new[(W1_u_new < 0) | (W1_u_new == 0)] <- 0
  W1_u_new[(W1_u_new > 1)] <- 1
  
  return(W1_u_new)
  
}

## Graph Community Weight updates
GraphComntyWtUpdt <- function(W1_u,W2_u,W2_neigh, W2_sum){
  ## First part of log likelihood of G based on the structure of the network
  ## Check if there are neighbours for this node and proceed only if there are neighbours
  llG_1 <- rep(0, dim(W1_u)[2])
  a <- exp(-1 * (W1_u %*% t(W2_neigh)))
  b <- a / (1 - a)
  b[is.infinite(b)] <- 0
  llG_1 <-  b %*% W2_neigh 
  ## 2nd part of log likelihood of G
  llG_2 <-  as.matrix(W2_sum - t(W2_u) - as.matrix(colSums(W2_neigh)))
  ## Total llG: This math has been verified
  llG <- c(llG_1) - c(llG_2)
  ## experimental scaled version of the above llog lik
  scale <- 1
  ## experimenting with adding L1 enalty to the weight matrix
  # L1_Wtpenalty <- 0.001
  # for(i in 1:length(llG)){
  #   llG[i] <- sign(llG[i]) * max(abs(llG[i]) - L1_Wtpenalty, 0 )
  # }
  return(llG)
}

## Binary covariates community weight
BincovCmntyWtUpdt <- function(nc, f_u,h_u, W, X_u, dir){
  llX <- rep(0, nc )
  if(dir =="directed"){
    WtMat <- matrix(c(f_u,h_u), nrow=1)
    nc_n <- nc*2

  }else{
    WtMat <- matrix(f_u, nrow=1)
    nc_n <- nc
  }
  if (length(W) > 0) {
    ## Updating the logistic regression weights
    ## Gradient function to update log likelihood for covariates
    ## adding intercept term
    Q_u <- CalcQuk(WtMat,W)
    ## This math has been verified
    llX <- t((X_u - Q_u) %*% W[,2:(nc+1)])
  }
  return(llX)
}

## Continuous covariates community weight
ContcovCmntyWtUpdt <- function(nc, Z_u, beta, W1_u,W2_u, sigmaSq, dir){
  
  if(dir =="directed"){
    llZ <- rep(0, nc )
    WW <- matrix(c(1,W1_u,W2_u), ncol = 1+(nc*2))
    
  }else{
    llZ <- rep(0, nc)
    WW <- matrix(c(1,W1_u), ncol = (1+nc))
    
  }
  if (length(beta) > 0) {
    ## Adding the gradient based on continuous covariates
    llZ <- 1 * ((Z_u - t(beta %*% t(WW))) %*% beta[,2:(nc+1)])
  }
  return(llZ)
}

CalcQuk <- function(wtmat, W){
  Q_uk <- rep(0, dim(W)[1])
  Q_uk <- 1 / (1 + exp(-1 * (cbind(1, wtmat) %*% t(W))))
  
  return(Q_uk)
}

updateLogisticParam_glmnet <- function(W,BC,Wtm1,Wtm2,missVals,lambda,alphaLR,dir,impType,seed){
  if(printFlg == TRUE){
    print("In update linear param")
  }
  #set.seed(seed)
  
  if(dir == "directed"){
    cm <- cbind( Wtm1 , Wtm2)
  }else{
    cm <- Wtm1
  }
  
  Wret <- W
  BCret <-  BC
  nc <- ncol(Wtm1)
  sigmaSq <- NULL
  mod <- list()
  logLik <- 0
  p <- list()
  total_loss <- rep(0,dim(BC)[2])
  if (length(W) > 0) {
    
    for(i in 1:dim(BC)[2]){
      ##"Ridge","LASSO","ElasticNet" 
      ## filter out the missing data so we update only based on the avaiable data
      if(sum(missVals[,i]) > 0){
        mIdx <- missVals[,i]
        X <-  cm[!mIdx,]
        y <- BC[!mIdx,i]
      }else{
        mIdx <- rep(FALSE, dim(BC)[1])
        X <- cm
        y <- BC[,i]
      }
      suppressWarnings({

        mod[[i]] <- glmnet(X, y, family = "binomial", lambda = lambda, alpha = 1,maxit = 50000)
        
      })
      #W_old <- W[i,] ## Temp
      W[i,] <- as.matrix(coef(mod[[i]]))[,1]
      #W[i,] <- alphaLR * W[i,] + (1 - alphaLR) * W_old ##temp
      
      p[[i]] <- predict(mod[[i]], newx = cm, s = lambda, type = "response")
      
      ## predicting to calculate logistic loss
      eta <- predict(mod[[i]], newx = X, s = lambda, type = "link")
      log_loss <- mean(log(1 + exp(eta)) - y * eta)
      beta <- as.vector(coef(mod[[i]], s = lambda))[-1]  # remove intercept
      lasso_penalty <- lambda * sum(abs(beta))
      total_loss[i] <- log_loss + lasso_penalty
      
      
      # MANUALLY CALCULATE LOG-LIKELIHOOD
      # For binomial regression (Logistic Model)
      eps <- 1e-15 ## to make sure the predicted probability does not go to 1 or 0 and cause numerical error
      pt <- pmin(pmax(p[[i]], eps), 1 - eps)
      logLik <- logLik + sum(y * log(pt[!mIdx]) + (1 - y) * log(1 - pt[!mIdx]), na.rm = TRUE)
      if(is.infinite(logLik)){
        print("Infinite value here")
      }
      
    }
    
    Wret <- W
    if (sum(missVals) > 0) {
      for(i in 1:dim(BC)[2]){
        idx <- which(missVals[, i])
        BCret[idx, i] <- ifelse(p[[i]][idx] > 0.5, 1, 0)
      }
    }
  }
  
  return(list(Wret, BCret, logLik, total_loss))
}

##logistic regression proximal gradient
updateLogisticParam_ProximalGradient <- function(W,BC,Wtm1,Wtm2,missVals,lambda,alpha,dir,impType,seed){
  if(printFlg == TRUE){
    print("In update linear param")
  }
  #set.seed(seed)
  
  if(dir == "directed"){
    cm <- cbind( Wtm1 , Wtm2)
  }else{
    cm <- Wtm1
  }
  
  Wret <- W
  BCret <-  BC
  nc <- ncol(Wtm1)
  sigmaSq <- NULL
  mod <- list()
  logLik <- 0
  p <- list()
  total_loss <- rep(0,dim(BC)[2])
  prox_steps <- 10
  if (length(W) > 0) {
    
    for(i in 1:dim(BC)[2]){
      ##"Ridge","LASSO","ElasticNet" 
      ## filter out the missing data so we update only based on the avaiable data
      if(sum(missVals[,i]) > 0){
        mIdx <- missVals[,i]
        X <-  cbind(1,cm[!mIdx,])
        y <- BC[!mIdx,i]
      }else{
        mIdx <- rep(FALSE, dim(BC)[1])
        X <- cbind(1,cm)
        y <- BC[,i]
      }
      for(k in 1:prox_steps){
        eta <- X %*% W[i,]
        p_tmp <- 1 / (1 + exp(-eta))
        grad <- t(X) %*% (y - p_tmp) / nrow(X)
        # gradient ascent
        W[i, ] <- W[i,] + alpha * grad       
        # proximal step
        W[i,-1] <- sign(W[i,-1]) * pmax(abs(W[i,-1]) - alpha * lambda, 0)  
      }
      
      #Prediction full data
      eta_full <- cbind(1,cm) %*% W[i,]
      p[[i]] <- 1 / (1 + exp(-eta_full))
      # Loss + logLik
      eta_train <- X %*% W[i,]
      
      log_loss <- mean(log1p(exp(eta_train)) - y * eta_train)
      l1_penalty <- lambda * sum(abs(W[i,-1]))
      total_loss[i] <- log_loss + l1_penalty
      
      logLik <- logLik + sum(y * eta_train - log1p(exp(eta_train)))
      
    }
    
    Wret <- W
    if (sum(missVals) > 0) {
      for(i in 1:dim(BC)[2]){
        idx <- which(missVals[, i])
        BCret[idx, i] <- ifelse(p[[i]][idx] > 0.5, 1, 0)
      }
    }
  }
  
  return(list(Wret, BCret, logLik, total_loss))
}

## logsitic regression coordinate descent
updateLogisticParam <- function(W,BC,Wtm1,Wtm2,missVals,lambda,alpha,dir,impType,seed) {
  
  if (printFlg == TRUE) {
    print("In update logistic param (Coordinate Descent)")
  }
  
  # Build covariates
  if (dir == "directed") {
    cm <- cbind(Wtm1, Wtm2)
  } else {
    cm <- Wtm1
  }
  
  Wret <- W
  BCret <- BC
  logLik <- 0
  if(!is.null(ncol(BC))){
    total_loss <- rep(0, ncol(BC))
  }else {total_loss  <- 0}
  p_list <- list()
  
  max_outer <- 25   # IRLS iterations
  max_inner <- 50   # coordinate descent sweeps
  tol <- 1e-6
  
  if (length(W) > 0) {
    
    for (i in 1:ncol(BC)) {
      
      # ---- Handle missing ----
      if (sum(missVals[, i]) > 0) {
        mIdx <- missVals[, i]
        X <- cbind(1, cm[!mIdx, , drop = FALSE])
        y <- BC[!mIdx, i]
      } else {
        mIdx <- rep(FALSE, nrow(BC))
        X <- cbind(1, cm)
        y <- BC[, i]
      }
      
      n <- nrow(X)
      p_dim <- ncol(X)
      
      # Initialize
      b <- W[i, ]
      
      for (outer in 1:max_outer) {
        
        eta <- X %*% b
        p <-  ifelse(eta >= 0,1 / (1 + exp(-eta)), exp(eta) / (1 + exp(eta)))#1 / (1 + exp(-eta))
        
        # Avoid numerical issues
        w <- p * (1 - p)
        w <- pmax(w, 1e-5)
        
        z <- eta + (y - p) / w
        z <- pmin(pmax(z, -1e5), 1e5)
        
        
        # Precompute weighted norms
        col_norms <- t(X^2) %*% w /n#colSums(w * X^2) / n
        col_norms[col_norms < 1e-8] <- 1e-8
        
        
        # Residual for weighted LS
        r <- z - X %*% b
        
        # ---- Coordinate Descent ----
        for (inner in 1:max_inner) {
          b_old <- b
          
          for (j in 1:p_dim) {
            
            # Add back contribution
            r <- r + X[, j] * b[j]
            
            # Weighted rho
            rho <- sum(w * X[, j] * r) / n
            
            if (j == 1) {
              # Intercept (no penalty)
              b[j] <- rho / col_norms[j]
            } else {
              b[j] <- sign(rho) * max(abs(rho) - lambda, 0) / col_norms[j]
            }
            
            # Remove updated contribution
            r <- r - X[, j] * b[j]
          }
          
          if (max(abs(b - b_old)) < tol) break
        }
        
        # Check outer convergence
        if (max(abs(X %*% b - eta)) < tol) break
      }
      
      # Save coefficients
      W[i, ] <- b
      
      # ---- Predictions (FULL data) ----
      eta_full <- cbind(1, cm) %*% b
      p_full <- 1 / (1 + exp(-eta_full))
      p_list[[i]] <- p_full
      
      # ---- Loss + logLik ----
      eta_train <- X %*% b
      
      log_loss <- mean(log1pexp(eta_train) - y * eta_train)
      l1_penalty <- lambda * sum(abs(b[-1]))
      total_loss[i] <- log_loss + l1_penalty
      
      logLik <- logLik + sum(y * eta_train - log1pexp(eta_train))
      if(is.infinite(logLik)){
        print("stop")
      }
    }
    
    Wret <- W
    
    # ---- Impute missing ----
    if (sum(missVals) > 0) {
      for (i in 1:ncol(BC)) {
        idx <- which(missVals[, i])
        BCret[idx, i] <- ifelse(p_list[[i]][idx] > 0.5, 1, 0)
      }
    }
  }
  
  return(list(Wret, BCret, logLik, total_loss))
}

log1pexp <- function(x) {
  ifelse(x > 0,
         x + log1p(exp(-x)),
         log1p(exp(x)))
}
## Continuous covariates parameter updates
SigmaSqCalc <- function(Z, beta, Ftot, Htot, missVals,dir) {
  sigmaSq <- rep(0, dim(beta)[1])
  if(dir == "directed"){
    Fm <- cbind(1,Ftot, Htot)
  }else{  
    Fm <- cbind(1,Ftot)
  }
  
  for (i in 1:dim(beta)[1]) {
    z_sum <- Z[!missVals[, i], i] - (t(beta[i, ]) %*% t(Fm[!missVals[, i], ]))
    sigmaSq[i] <- sum(z_sum ^ 2, na.rm = TRUE) / sum(!missVals[, i])
    if(is.infinite(sigmaSq[i])){
      print("Null Sigma")
    }
  }
  return(sigmaSq)
}


updateLinearRegParam_glmnet <- function(beta,missVals,Z,Wtm1,Wtm2, alpha,lambda,N,dir, 
                                 impType, alphaLin, penalty,seed){
  if(printFlg == TRUE){
    print("In update linear param")
  }
  #set.seed(seed)
  
  if(dir == "directed"){
    cm <- cbind( Wtm1 , Wtm2)
  }else{
    cm <- Wtm1
  }
  
  betaret <- beta
  Zret <-  Z
  nc <- ncol(Wtm1)
  sigmaSq <- NULL
  mod <- list()
  logLik <- 0
  predictions <- list()
  if (length(beta) > 0) {
    
    for(i in 1:dim(Z)[2]){
      ##"Ridge","LASSO","ElasticNet" 
      ## filter out the missing data so we update only based on the avaiable data
      if(sum(missVals[,i]) > 0){
        mIdx <- missVals[,i]
        X <-  cm[!mIdx,]
        y <- Z[!mIdx,i]
      }else{
        mIdx <- rep(FALSE, dim(Z)[1])
        X <- cm
        y <- Z[,i]
      }
      
      if(penalty == "LASSO"){
        suppressWarnings({
          mod[[i]] <- glmnet(X, y, family = "gaussian" , lambda = lambda, alpha = 1, maxit =1000)
        })

        #beta_old <- beta[i,] ## Temp
        beta[i,] <- as.matrix(coef(mod[[i]]))[,1]
        #beta[i,] <- alpha * beta[i,] + (1 - alpha) * beta_old ##temp
        
      }else if(penalty =="ElasticNet"){
        
        mod[[i]] <- glmnet(X, y, family = "gaussian",lambda = lambda, alpha = 0.5 , maxit = 1000)
        beta[i,] <- as.matrix(coef(mod[[i]]))[,1]
        
      }else if(penalty == "Ridge"){
        mod[[i]] <- glmnet(X, y, family = "gaussian", lambda = lambda, alpha = 0, maxit = 2000)
        beta[i,] <- as.matrix(coef(mod[[i]]))[,1]
      }
      
      predictions[[i]] <- predict(mod[[i]], newx = cm, s = lambda, type = "response")
      
      # MANUALLY CALCULATE LOG-LIKELIHOOD
      # For Gaussian regression (Linear Model)
      logLik <- logLik + sum(dnorm(y, 
                                   mean = predictions[[i]][!mIdx], 
                                   sd = sqrt(mean((y - predictions[[i]][!mIdx])^2)), log = TRUE))
    }
    
    betaret <- beta
    sigmaSq <-  SigmaSqCalc(Z, betaret, Wtm1,Wtm2, missVals,dir)
    
    if (sum(missVals) > 0) {
      for(i in 1:dim(Z)[2]){
        idx <- which(missVals[, i])
        Zret[idx, i] <- predictions[[i]][idx]
      }
    }
  }
  return( list( betaret, Zret, sigmaSq, logLik) )
}

updateLinearRegParam_ProximalGradient <- function(beta,missVals,Z,Wtm1,Wtm2, alpha,lambda,N,dir, 
                                        impType, alphaLin, penalty,seed){
  if(printFlg == TRUE){
    print("In update linear param")
  }
  #set.seed(seed)
  
  if(dir == "directed"){
    cm <- cbind( Wtm1 , Wtm2)
  }else{
    cm <- Wtm1
  }
  alpha <- alpha /4 ## learning rate for linear regression needs to be smaller compared to learning rate for binary covariates.
  betaret <- beta
  Zret <-  Z
  nc <- ncol(Wtm1)
  sigmaSq <- NULL
  mod <- list()
  logLik <- 0
  predictions <- list()
  prox_steps <- 10
  if (length(beta) > 0) {
    
    for(i in 1:dim(Z)[2]){
      # Filter missing
      if(sum(missVals[,i]) > 0){
        mIdx <- missVals[,i]
        X <- cm[!mIdx,]
        y <- Z[!mIdx,i]
      } else {
        mIdx <- rep(FALSE, nrow(Z))
        X <- cm
        y <- Z[,i]
      }
      
      # Add intercept
      X_design <- cbind(1, X)
      
      # Proximal gradient steps
      for(k in 1:prox_steps){
        # Residual
        r <- X_design %*% beta[i,] - y
        # Gradient of squared loss
        grad <- t(X_design) %*% r / nrow(X_design)
        # Gradient step
        beta_temp <- beta[i,] - alpha * grad
        beta[i,] <- beta_temp
        beta[i,-1] <- sign(beta_temp[-1]) * pmax(abs(beta_temp[-1]) - alpha * lambda, 0)
        }
      
      # Predictions (FULL data)
      X_full <- cbind(1, cm)
      predictions[[i]] <- X_full %*% beta[i,]
      # Log-likelihood (TRAINING ONLY)
      resid <- y - (X_design %*% beta[i,])
      sigma2 <- mean(resid^2)
      
      logLik <- logLik + sum( -0.5 * log(2 * pi * sigma2) - (resid^2) / (2 * sigma2))
    }
    
    betaret <- beta
    sigmaSq <-  SigmaSqCalc(Z, betaret, Wtm1,Wtm2, missVals,dir)
    
    if (sum(missVals) > 0) {
      for(i in 1:dim(Z)[2]){
        idx <- which(missVals[, i])
        Zret[idx, i] <- predictions[[i]][idx]
      }
    }
  }
  return( list( betaret, Zret, sigmaSq, logLik) )
}

## Using coordinate descent so we don't have to worry about alpha.
updateLinearRegParam <- function(beta,missVals,Z,Wtm1,Wtm2, alpha,lambda,N,dir, 
                                    impType, alphaLin, penalty,seed) {
  
  if (printFlg == TRUE) {
    print("In update linear param (Coordinate Descent)")
  }
  
  # Build covariate matrix
  if (dir == "directed") {
    cm <- cbind(Wtm1, Wtm2)
  } else {
    cm <- Wtm1
  }
  
  betaret <- beta
  Zret <- Z
  logLik <- 0
  sigmaSq <- NULL
  predictions <- list()
  
  max_iter <- 50   # number of coordinate descent sweeps
  tol <- 1e-6      # convergence tolerance
  
  if (length(beta) > 0) {
    
    for (i in 1:ncol(Z)) {
      
      # ---- Handle missing ----
      if (sum(missVals[, i]) > 0) {
        mIdx <- missVals[, i]
        X <- cm[!mIdx, , drop = FALSE]
        y <- Z[!mIdx, i]
      } else {
        mIdx <- rep(FALSE, nrow(Z))
        X <- cm
        y <- Z[, i]
      }
      
      # ---- Design matrix ----
      X_design <- cbind(1, X)
      n <- nrow(X_design)
      p <- ncol(X_design)
      
      # ---- Initialize beta ----
      b <- beta[i, ]
      
      # ---- Precompute column norms ----
      col_norms <- colSums(X_design^2) / n
      
      # ---- Initialize residual ----
      r <- y - X_design %*% b
      
      # ---- Coordinate Descent ----
      for (iter in 1:max_iter) {
        b_old <- b
        
        for (j in 1:p) {
          
          # Add back current feature contribution
          r <- r + X_design[, j] * b[j]
          
          # Compute rho
          rho <- sum(X_design[, j] * r) / n
          
          # Update coefficient
          if (j == 1) {
            # Intercept (no penalty)
            b[j] <- rho
          } else {
            b[j] <- sign(rho) * max(abs(rho) - lambda, 0) / col_norms[j]
          }
          
          # Remove updated contribution
          r <- r - X_design[, j] * b[j]
        }
        
        # ---- Convergence check ----
        if (max(abs(b - b_old)) < tol) break
      }
      
      # Save updated beta
      beta[i, ] <- b
      
      # ---- Predictions (FULL data) ----
      X_full <- cbind(1, cm)
      predictions[[i]] <- X_full %*% b
      
      # ---- Log-likelihood (training only) ----
      resid <- y - (X_design %*% b)
      sigma2 <- mean(resid^2)
      
      logLik <- logLik + sum(
        -0.5 * log(2 * pi * sigma2) - (resid^2) / (2 * sigma2)
      )
    }
    
    betaret <- beta
    
    sigmaSq <- SigmaSqCalc(Z, betaret, Wtm1, Wtm2, missVals, dir)
    
    # ---- Impute missing ----
    if (sum(missVals) > 0) {
      for (i in 1:ncol(Z)) {
        idx <- which(missVals[, i])
        Zret[idx, i] <- predictions[[i]][idx]
      }
    }
  }
  
  return(list(betaret, Zret, sigmaSq, logLik))
}

sdErr <- function(Z_i, beta, Fm, N, p, k){
  sdErr <- sqrt(sum((Z_i -  Fm %*% as.matrix(beta))^2)/ (N-p-k))
  return(sdErr)
}

## Soft thresholding function
soft_thresh<- function(z_j,lambda){
  if(z_j < lambda){
    return(z_j - lambda)
  }else if(z_j < (-1*lambda)){
    return(z_j + lambda)
  }else{
    return(0)
  }
}

## Log likleihood calculations
Lx_cal <- function(W, N, Fmat, Hmat, X, dir){
  S <- 0
  nc <- dim(Fmat)[2]
  if(dir =="directed"){
    wtmat <- cbind(Fmat,Hmat)
  }else{
    wtmat <- Fmat
  }
  if (length(W) > 0) {
    cases <- complete.cases(X)
    Q <- CalcQuk(wtmat, W)
    S <- sum((log(Q[cases, ]) * X[cases,]) + (log(1-Q[cases,]) * (1- X[cases,])), na.rm = TRUE)
  }
  
  return(S)
}

Lz_cal <- function(beta,N,Z,Fmat, Hmat,sigmaSq,dir,lambda,missVals, penalty, alphaLin ){
  if(printFlg == TRUE){
    print("In Lz_Calc Func")
  }
  if(dir == "directed"){
    Fin <- cbind(1,Fmat, Hmat)
  }else {
    Fin <- cbind(1,Fmat)
  }
  S <- 0
  
  if (!is.null(beta)) {
    sigmaSq <- rep(1, dim(beta)[1])
    
    for (j in 1:dim(beta)[1]) {
      df <- sum(beta[j,] != 0)  
      
      RSS <- sum(( Z[,j]  - (Fin %*% beta[j, ])) ^ 2)
      S <- S - RSS/2
      
    }
  }
  return(S)
}

Lg_cal <- function(G, Ftot, Htot,Eneg,epsilon, dir, E){
  scale <- 1
  Fv <- Ftot[E[, 1], ]
  Hv <- Htot[E[, 2], ]
  N <-  gorder(G)
  ## Variable to save the loglik contributed by nodes with edges between them
  ## Implementing rowsums Hopefully this is a little faster
  mulop <- log(1 - exp(-1 * rowSums(Fv * Hv)))
  S11 <- sum(is.finite(mulop) * mulop, na.rm = TRUE)
  
  ## Creating a graph that is negative of the given graph to get edges not present in G
  E2 <-  Eneg 
  Fv2 <- Ftot[E2[, 1], ]
  Hv2 <- Htot[E2[, 2], ]
  
  ## Variable to save the loglik contributed by nodes without edges between them
  S12 <- sum(rowSums(Fv2 *Hv2))
  S1 <- S11 - S12
  return(S1)
}

findLLDir <- function(G,  Ftot, Htot, Win = NA, Wout = NA, X_in = NA,X_out = NA,
                      betain = NULL, betaout = NULL,Z_in = NA, Z_out = NA,
                      sigmaSqin = NA,sigmaSqout = NA,  alphaLL, Eneg, dir, 
                      epsilon, lambda_bin, lambda_lin,missVals_bin,missVals_cont , penalty, alphaLin, E, LinLL =0, BinLL =0) {
  if(printFlg == TRUE){
    print("In Find LLDir Func")
  }
  
  N <-  gorder(G)
  
  ## Calculating the L_g part of the log likelihood for structure of the network
  S1 <- Lg_cal(G, Ftot, Htot, Eneg,epsilon, dir, E)
  
  ## Calculating the L_xin part of the log likelihood for binary covariates
  S2_in <- 0
  
  ## Calculating the L_xout part of the log likelihood for binary covariates
  S2_out <- 0#Lx_cal(W=Wout, N = N, Fmat= Ftot, Hmat =Htot, X=X_out, dir)
  
  ## adding the loglikelihood from the covariates
  S2 <- BinLL#S2_in + S2_out
  
  ## Calculating the L_zin part of the log likelihood for continuous covariates
  S3_in <- 0
  
  ## Calculating the L_zout part of the log likelihood for continuous covariates
  S3_out <- 0#Lz_cal(beta = betaout,N = N,Z = Z_out,
            #       Fmat = Ftot, Hmat = Htot,
            #       sigmaSq = sigmaSqout,dir,lambda, missVals, penalty, alphaLin)
  S3 <- LinLL
  
  ## Calculating the final log likelihood
  ll <- (S1 +  (S2 + S3))
  if(is.nan(ll)){
    print(ll)
    print("stop here in find ll")
  }
  return(c(ll, S1, S2, S3))
}

## Initialization functions
initCov <- function(covtmp, CovNames) {
  X_in <- as.matrix(covtmp %>%
                      dplyr::select(all_of(CovNames)))
  missVals <- is.na(X_in)
  numMissVal <- sum(missVals)
  
  return(list(X_in, numMissVal))
}

initWtmat <- function(G, mode, N, nc, seed,specOP,epsilon){
  #set.seed(seed)
  Wtmat <- matrix(nrow = N, ncol = nc,runif(nc * N, 0.0001,0.0005))
  rownames(Wtmat) <- names(igraph::degree(G, mode = mode))
  return(Wtmat)
}

initW <- function(k, nc,missVals, X, Fmat, lambda, alpha, ncoef){
  W <- NULL
  if (k > 0) {
    W <- matrix(nrow = k, ncol = (ncoef) + 1, 0)
    for(i in 1:k){
      y_bar <- mean(X[,i], na.rm =TRUE)
      W[i,1] <- log((1-y_bar)/y_bar)
    }
  }
  return(W)
}

initbeta <- function(o, ncoef,missVals, Z, Fmat, Hmat, alpha,N, lambda,dir,
                     impType, alphaLin, penalty){
  beta <- NULL
  sigmaSq <- NULL
  if (o > 0) {
    beta <- matrix(nrow = o, ncol = (ncoef) + 1, 0)
    
    for(i in 1:o){
      beta[i,1] <- mean(Z[,i], na.rm =TRUE)
    }
    sigmaSq <- 1
  }
  return(list(beta, sigmaSq))
}

randomizeIdx <- function(N, randomize){
  ## Updating the community weight parameter
  if (randomize == TRUE) {
    s <- sample(1:N, size = N, replace = FALSE)
  }else{s <- 1:N}
  return(s)
}

CoDA <- function(G,nc, k = c(0, 0) ,o = c(0, 0) , N,  alpha, lambda_lin, lambda_bin, thresh, 
                 nitermax, orig,randomize = TRUE, CovNamesLinin = c(),
                 CovNamesLinout = c(),CovNamesLPin = c(),CovNamesLPout = c(), 
                 dir, alphaLL = NULL,test = TRUE,missing = NULL, covOrig, 
                 epsilon,  impType, alphaLin, penalty, seed, covInit, specOP,nc_sim,lambda_grph) {
  if(printFlg == TRUE){
    print("In CoDA Func")
  }
  
  if(dir == "directed"){
    ncoef = nc*2
  }else{ncoef = nc}
  
  ## Setting the original covariates
  if(!is.null(covOrig[[1]])){
    covOrig_bin <- covOrig[[1]]
  }
  if(!is.null(covOrig[[2]])){
    covOrig_cont <- covOrig[[2]]
  }
  
  ## Community weights outgoing connections
  Ftot <- initWtmat(G,"out",N,nc, seed, specOP,epsilon)
  Finit <- Ftot
  
  ## Community weights outgoing connections
  if(dir == "directed"){
    Htot <- initWtmat(G,"in",N,nc, (seed+10), specOP,epsilon)
    Hinit <- Htot
  } else{
    Htot <- Ftot
    Hinit <- Htot
  }
  
  #set.seed(seed)
  
  iter <- 0
  k_in <- k[1]
  k_out <- k[2]
  o_in <- o[1]
  o_out <- o[2]
  
  ########## Get the covariate matrix for binary data and continuous data ##########
  b <- getNWcov(G, CovNamesLinin,CovNamesLinout,CovNamesLPin,CovNamesLPout,k_in,k_out,o_in,o_out)
  
  #X_in, X_out, Z_in, Z_out, missValsin, missValsout_bin, missValsout_cont, numMissValin, numMissValout_bin, numMissValout_cont
  X_in <- b[[1]]
  X_out <- b[[2]]
  
  Z_in <- b[[3]]
  Z_out <- b[[4]]
  
  missValsin <- b[[5]]
  missValsout_bin <- b[[6]]
  missValsout_cont <- b[[7]]
  numMissValin <- b[[8]]
  numMissValout_bin  <- b[[9]]
  numMissValout_cont <- b[[10]]
  
  ## For covariates where the value is missing, fill it in initially with random values
  if(k_out >0){
    for(i in 1:dim(X_out)[2]){
      X_out[missValsout_bin[,i],i] <- rbinom(n = sum(missValsout_bin[,i]), size = 1, prob = (N- sum(missValsout_bin[,i]))/N )
    }
  }
  
  backTransParams <- list()
  if(o_out > 0){
    
    for(i in 1:dim(Z_out)[2]){
      
      if(sum(missValsout_cont[,i]) > 0){
        
        if(covInit == "mean"){
          ## Mean imputation of initial values
          Z_out[missValsout_cont[,i],i] <- mean(Z_out[,i], na.rm =TRUE) 
        }else if(covInit == "Nmode"){
          ## Imputing values based on the neighbours. Using mode of the neighbour values
          mv <- which(missValsout_cont[,i])
          for(v in mv){
            neigh <- igraph::neighbors(G, v, mode = "all")
            if(length(neigh) > 0){
              if(all(is.na(Z_out[c(neigh), i]))){
                Z_out[v,i] <- mode(Z_out[,i])
              }else{
                Z_out[v, i] <- mode(Z_out[c(neigh), i])
              }
            } else{
              ## Mean imputation of initial values if no neighbours
              Z_out[v,i] <- mode(Z_out[,i])
            }
          }
        }else if(covInit == "Nmean"){
          ## Imputing values based on the neighbours. Using mode of the neighbour values
          mv <- which(missValsout_cont[,i])
          for(v in mv){
            neigh <- igraph::neighbors(G, v, mode = "all")
            if(length(neigh) > 0){
              if(all(is.na(Z_out[c(neigh), i]))){
                Z_out[v,i] <- mean(Z_out[,i], na.rm =TRUE)
              }else{
                Z_out[v, i] <- mean(Z_out[c(neigh), i], na.rm =TRUE)
              }
            } else{
              ## Mean imputation of initial values if no neighbours
              Z_out[v,i] <- mean(Z_out[,i], na.rm =TRUE)
            }
          }
        }else if(covInit == "Nmedian"){
          ## Imputing values based on the neighbours. Using mode of the neighbour values
          mv <- which(missValsout_cont[,i])
          for(v in mv){
            neigh <- igraph::neighbors(G, v, mode = "all")
            if(length(neigh) > 0){
              if(all(is.na(Z_out[c(neigh), i]))){
                Z_out[v,i] <- median(Z_out[,i], na.rm = TRUE)
              }else{
                Z_out[v, i] <- median(Z_out[c(neigh), i], na.rm = TRUE)
              }
            } else{
              ## Mean imputation of initial values if no neighbours
              Z_out[v,i] <- median(Z_out[,i], na.rm = TRUE)
            }
          }
        }else if(covInit == "NNmedian"){
          ## Imputing values based on the neighbours. Using mode of the neighbour values
          mv <- which(missValsout_cont[,i])
          for(v in mv){
            n1 <- igraph::neighbors(G, v, mode = "all")
            n_of_n <- unique(unlist(lapply(n1, function(x) neighbors(G, v = x))))
            # Remove the original vertex and its direct neighbors if desired
            neigh <- setdiff(n_of_n, c(1, n1))
            if(length(neigh) > 0){
              if(all(is.na(Z_out[c(neigh), i]))){
                Z_out[v,i] <- median(Z_out[,i], na.rm = TRUE)
              }else{
                Z_out[v, i] <- median(Z_out[c(neigh), i], na.rm = TRUE)
              }
            } else{
              ## Mean imputation of initial values if no neighbours
              Z_out[v,i] <- median(Z_out[,i], na.rm = TRUE)
            }
          }
        }
      }
    }
  }
  
  
  ## scaling the continuous covariats
  if(o_out > 0 ){
    Z_out <- scale(Z_out)
    covOrig_cont <- scale(covOrig_cont)  
  }
  
  
  ############ Initialize the coef matrix for covariates ##############
  
  ############## Setting initial values of binary covariates coefs weights ###############
  Win <- NULL
  Win_orig <- Win
  
  Wout <- initW(k_out, nc,missValsout_bin, X_out, Htot, lambda_bin, alpha, ncoef)
  Wout_orig <- Wout
  
  #################### Setting initial values of linear regression weights
  b <- NULL
  betain <- b
  sigmaSqin <- b
  
  b <- initbeta(o_out, ncoef, missValsout_cont, Z_out, Ftot, Htot, alpha, N, lambda_lin, dir, impType, alphaLin, penalty)
  betaout <- b[[1]]
  sigmaSqout <- b[[2]]
  
  ## List of neighbour for the nodes
  if(dir == "directed"){
    inNeigh <- lapply(V(G), function(x) neighbors(G, x, "in"))
    outNeigh <- lapply(V(G), function(x) neighbors(G, x, "out"))
  }else{
    inNeigh <- lapply(V(G), function(x) neighbors(G, x, "all"))
    outNeigh <- inNeigh
  }
  
  ## Getting the negative adjacency matrix
  Aneg <- 1 - (1 * as.matrix(as_adjacency_matrix(G)))
  Eneg  <- as_edgelist(simplify(graph_from_adjacency_matrix(Aneg, mode = dir)))
  
  ##Initial log likelihood value
  LLold <- -10000
  lllst <- matrix(nrow = 0, ncol = 4)
  
  arilst <- c()
  OmegaVal <- c()
  mse <- matrix(nrow = 0, ncol = (o_in+o_out))
  mseMD <- matrix(nrow = 0, ncol = (o_in + o_out))
  
  accuracy <- matrix(nrow = 0, ncol = k_in + k_out)
  accuracyMD <- matrix(nrow = 0, ncol = k_in + k_out)
  
  ##logistic regression loss
  bin_loss <- matrix(nrow=0,ncol=k_in+k_out)
  
  corMat <- matrix(0, nrow =0 , ncol = 15)
  
  delta <- getDelta(N,epsilon)
  continue <- TRUE
  s <- 1:N
  beta_old <- betaout
  
  ## get the edges in the graph
  E <- igraph::as_edgelist(G)
  LinLL = 0
  BinLL = 0 
  repeat{  
    ## Update the log likelihood.
    LLvec <- findLLDir(G, Ftot, Htot, Win, Wout, X_in, X_out, betain, betaout,
                       Z_in, Z_out, sigmaSqin, sigmaSqout, alphaLL, Eneg, 
                       dir, epsilon, lambda_bin, lambda_lin, missValsout_bin,missValsout_cont,
                       penalty, alphaLin, E, LinLL, BinLL)
    LLnew <- LLvec[1]
    pctInc <- abs((LLnew - LLold)/(LLold)) * 100
    if(is.nan(pctInc)){
      C <- as.matrix(memOverlapCalc(Ftot,Htot,delta, N, nc),ncol = nc )
      if(test == FALSE){
        
        arilst <- c(arilst, ARIop(Ftot,Htot,orig,nc,N))
        OmegaVal <- c(OmegaVal, OmegaIdx(G, C, N, nc, nc_sim) )
        
        MSEtmp <- MSEop(Ftot, Htot, covOrig_cont, betaout,N, dir, o_in,o_out, missValsout_cont)
        mseMD <- rbind(mseMD, MSEtmp[[1]])
        mse <- rbind(mse, MSEtmp[[2]])
        
        accutmp <- AccuracyCalc(k_out,Wout,Ftot, Htot,covOrig_bin,dir,missValsout_bin)
        accuracyMD <- rbind(accuracyMD, accutmp[[1]])
        accuracy <- rbind(accuracy, accutmp[[2]])
        bin_loss <- total_loss
        
      }
      FAILURE <- TRUE
      omgIdx <- round(OmegaIdx(G, C, N, nc, nc_sim) ,3)
      print(paste0("ERROR in Pct increase ,total iterations ", 
                   iter, " with seed value ", seed, " Final OI ",omgIdx))
      break
    }
    
    if((pctInc < thresh) | (dim(lllst)[1] == nitermax)){
      C <- as.matrix(memOverlapCalc(Ftot,Htot,delta, N, nc),ncol = nc )
      if(test == FALSE){
        lllst <- rbind(lllst, LLvec)
        
        arilst <- c(arilst, ARIop(Ftot,Htot,orig,nc,N))
        OmegaVal <- c(OmegaVal, OmegaIdx(G, C, N, nc, nc_sim))
        if(o_out > 0 ){
          MSEtmp <- MSEop(Ftot, Htot, covOrig_cont, betaout,N, dir, o_in,o_out, missValsout_cont)
          mseMD <- rbind(mseMD, MSEtmp[[1]])
          mse <- rbind(mse, MSEtmp[[2]])  
        }
        
        if(k_out > 0 ){
          accutmp <- AccuracyCalc(k_out,Wout,Ftot, Htot,covOrig_bin,dir,missValsout_bin)
          accuracyMD <- rbind(accuracyMD, accutmp[[2]])
          accuracy <- rbind(accuracy, accutmp[[1]])
          bin_loss <- total_loss  
        }
      }
      FAILURE <- FALSE
      print(paste0("The final percent change ", round(pctInc,6), 
                   " ,total iterations ", iter, " with seed value ", seed, 
                   " Final LogLik ", round(LLnew,3),
                   " Final OI ",round(OmegaIdx(G, C, N, nc, nc_sim),3), 
                   " MSE ", paste0(round(tail(mse,1),3), collapse = ", "),
                   " MSE MD ", paste0(round(tail(mseMD,1),3), collapse = ", "),
                   " Acc ", paste0(round(tail(accuracy,1),3), collapse = ", "),
                   " Acc MD ", paste0(round(tail(accuracyMD,1),3), collapse = ", ")))
      break
    }
    
    LLold <- LLnew
    lllst <- rbind(lllst, LLvec)
    iter <- iter + 1
    
    ## Randomize the community weights update
    s <- randomizeIdx(N, randomize)
    
    ## Updating F i.e. all outgoing connections
    ## We look for neighbors our node is connected to with the edges directed outward from our node.
    opMat <- updateWtmat(G,Ftot,Htot,"all",s,nc,X_out,Z_out,k_out,o_out,
                         betaout,Wout, alphaLL,missValsout_bin, missValsout_cont,alpha, start = 2, 
                         end = nc + 1, dir, inNeigh, outNeigh, lambda_bin, lambda_lin, penalty,lambda_grph)
    Ftot <- opMat[[1]]
    Htot <- opMat[[2]]
    
    ## Updating logistic paramters
    b2 <- updateLogisticParam(W = Wout, BC = X_out , Wtm1 = Ftot, Wtm2 = Htot, 
                              missVals = missValsout_bin, lambda = lambda_bin, 
                              alpha = alpha, dir = dir, impType, seed)
    Wout <- b2[[1]]
    X_out <- b2[[2]]
    BinLL <- b2[[3]]
    total_loss <- b2[[4]]
    
    ### Updating the beta matrix for continuous covariates
    c1 <- updateLinearRegParam(betaout,missValsout_cont,Z_out,Ftot,Htot,alpha,
                               lambda_lin, N, dir, impType, alphaLin, penalty,seed)
    betaout <- c1[[1]]
    Z_out <- c1[[2]]
    sigmaSqout <- c1[[3]]
    LinLL <- c1[[4]]
    Z_in <- 0
   
    
    if(test == TRUE){
      
      C <- as.matrix(memOverlapCalc(Ftot,Htot,delta, N, nc),ncol = nc )
      OIn <- OmegaIdx(G, C, N, nc, nc_sim)
      
      OmegaVal <- c(OmegaVal, OIn)
      MSEtmp <- MSEop(Ftot, Htot, covOrig_cont, betaout,N, dir, o_in,o_out, missValsout_cont)
      mseMD <- rbind(mseMD, MSEtmp[[1]])
      mse <- rbind(mse, MSEtmp[[2]])
      accutmp <- AccuracyCalc(k_out,Wout,Ftot, Htot,covOrig_bin,dir,missValsout_bin)
      accuracyMD <- rbind(accuracyMD, accutmp[[1]])
      accuracy <- rbind(accuracy, accutmp[[2]])
      
      bin_loss <- rbind(bin_loss, total_loss)
      
    }
  }
  
  if(dir == "directed"){
    WtMat <- cbind(Ftot, Htot)
  }else{
    WtMat <- Ftot
  }
  BICv <- c(BIC(LLnew, nc, N, ecount(G)), 
            ICL(LLnew, WtMat, nc, N,ecount(G)), 
            LLnew, nc, N, ecount(G))
  
  
  
  Ffin <- Ftot
  Hfin <- Htot
  Winfin <- Win
  Woutfin <- Wout
  betainfin <- betain
  betaoutfin <- betaout
  backTransParams <- backTransParams
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
      Xin_cov = X_in,
      betaout = betaoutfin,
      Zout_cov = Z_out,
      Xout_cov = X_out,
      MSE = mse,
      Acc = accuracy,
      OmegaIndex = OmegaVal,
      AccMD = accuracyMD,
      MSEMD = mseMD,
      FAILURE,
      backTransParams,
      corMat,
      bic = BICv,
      bin_loss = bin_loss
    )
  )
}
