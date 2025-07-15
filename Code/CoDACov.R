library(tidyverse)
library(dplyr)
library(igraph)
library(network)
library(intergraph)
library(ergm)
library(glmnet)
library(gglasso)
library(RColorBrewer)
library(clustAnalytics)

source(paste0(getwd(),"/Code/NetworkMetrics.R"))
source(paste0(getwd(),"/Code/HelperFuncs.R"))

## Helper, error and accuracy calculations
## Set the values to be between range 0 and 1 
range01 <- function(x){
  if(!all(x == x[1])){
    (x-min(x))/(max(x)-min(x))
  }else{
    x
  }
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
  comAff <-  as.data.frame(vertex_attr(G)) %>%
    dplyr::select(all_of(letters[1:nc]))
  
  b <- getNWcov(G, CovNamesLinin, CovNamesLinout, 
                CovNamesLPin, CovNamesLPout, 
                k_in,k_out,o_in,o_out)
  
  X_in <- b[[1]]
  X_out <- orig_Cov
  
  Z_in <- b[[3]]
  Z_out <- orig_Cov
  
  missValsin <- b[[5]]
  missValsout <- b[[6]]
  numMissValin <- b[[7]]
  numMissValout <- b[[8]]
  
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
        mod <- glm(X_out[,i] ~ Fm[,1]+Fm[,2]+Fm[,3],family = "binomial")
      } else{
        mod <- glm(X_out[,i] ~ Fm[,1]+Fm[,2]+Fm[,3]+Hm[,1]+Hm[,2]+Hm[,3],family = "binomial")
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
        mod <- lm(Z_out[,i] ~ Fm[,1]+Fm[,2]+Fm[,3])
      }else{
        mod <- lm(Z_out[,i] ~ Fm[,1]+Fm[,2]+Fm[,3]+Hm[,1]+Hm[,2]+Hm[,3])
      }
      
      betaout[i,] <- coef(mod)
      contLL <- contLL + as.numeric(logLik(mod))
      sigmaSqout[i] <- sigma(mod)^2
    }
  }
  
  A <- 1 - (1 * as.matrix(as_adjacency_matrix(G)))
  Eneg  <- as_edgelist(igraph::graph_from_adjacency_matrix(A, mode = dir))
  
  GrphLL <- Lg_cal(G = G, Ftot = Fm, Htot = Hm,Eneg = Eneg, epsilon = epsilon, dir = dir)
  ll <- binLL + contLL + GrphLL
  
  return(list(ll, c(binLL,contLL,GrphLL), VIFval))
}

prob2logit <- function(x) {
  return(log(x / (1 - x)))
}

NWSimBin <- function(nc, k,  pC, N, pClust, B, o, dist,dir,
                     covTypes = c("binary", "continuous"), 
                     CovNamesLin = c(),CovNamesLP = c(), missing = NULL) {
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
        binVal[which(as.numeric(net %v% 'Cluster') == i), j] <- rbinom(length(which(as.numeric(net %v% 'Cluster') == i)), 
                                                                       1, 
                                                                       pC[i, j])
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
        contVal[which(as.numeric(net %v% 'Cluster') == i), j] <- rnorm(length(which(as.numeric(net %v% 'Cluster') == i)), 
                                                                       dist[i, j][[1]][[1]], 
                                                                       dist[i, j][[1]][[2]])
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
    control = ergm::control.simulate(MCMC.burnin = 10000, 
                                     MCMC.interval = 1000)
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

mode <- function(x) {
  which.max(tabulate(x))
}

## Network generation code
SetBinarycov <- function(G, k, N,nc, k_start, Cnew,connType, pC,
                         CovNamesLP, missing) {
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
  if(!is.null(missing)) {
    for (j in 1:k) {
      binVal[sample(1:N, round(missing[j] * N / 100)), j] <- NA
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

OldSetContinuousCov <- function(G,o, N,nc, o_start, C, connType, 
                                dist,CovNamesLin,missing,missType, MARparam=NULL){
  ## Based on the mean and variance indicated by the dist matrix we assign continuous values to the network covariates.
  contVal <- matrix(ncol = o, nrow = N, 0) 
  s <- 1 
  
  for (i in 1:nc) {
    idx <- which(C[,i] == connType)
    contVal[idx, i] <- rnorm(length(idx), dist[[s]][1], dist[[s]][2])
    idx_n <- which(!(1:N %in% idx))
    contVal[idx_n, i] <- rnorm(length(idx_n), dist[[s+1]][1], dist[[s+1]][2])
    s <- s+2
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

SetContinuousCov <- function(G,o, N,nc, o_start, C, connType, 
                             dist,CovNamesLin,missing,missType, MARparam=NULL){
  ## Based on the mean and variance indicated by the dist matrix we assign continuous values to the network covariates.
  contVal <- matrix(ncol = o, nrow = N, 0) 
  s <- 1 
  
  for (i in 1:nc) {
    ## set incoming community covariate
    idx_in <- which(C[,i] == -1)
    contVal[idx_in, i] <- rnorm(length(idx_in), dist[[s]][1], dist[[s]][2])
    #Set outgoing community covariate
    idx_out <- which(C[,i] == 1)
    contVal[idx_out, i] <- rnorm(length(idx_out), dist[[s+1]][1], dist[[s+1]][2])
    ## setting covariates for the non community nodes
    idx_n <- which(!(1:N %in% c(idx_in, idx_out)))
    contVal[idx_n, i] <- rnorm(length(idx_n), dist[[s+2]][1], dist[[s+2]][2])
    s <- s+3
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
      for(j in 1:dim(c)[1]){
        com <- which(C %in% c[j, ])
        fin <- c(fin, com)
        idx <- sample(com , size = length(com) * pClustOL[l], replace = FALSE)
        Cnew[idx, c[j, ]] <- 1
        l <- l + 1
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
      c <- t(combn(letters[1:nc], 2))
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
    F_u[Cnew > 0] <- 0.5#rbeta(sum(Cnew > 0),a,b)
    H_u[Cnew < 0] <- 0.5#rbeta(sum(Cnew < 0),a,b)
  }else {
    F_u[Cnew > 0] <- 0.5#rbeta(sum(Cnew > 0),a,b)
    H_u <- F_u
  }
  
  A <- matrix(0, nrow = N, ncol = N)
  Apuv <- A
  for (i in 1:N) {
    for (j in 1:N) {
      ## probability of edge from node u to v
      p_uv <- 1 - exp(-1 * sum(F_u[i, ] * H_u[j, ])) + epsilon
      Apuv[i, j] <- p_uv
      A[i, j] <- rbinom(1,1,p_uv)#median(rbinom(25,1,p_uv))#purrr::rbernoulli(1, p_uv)
    }
  }
  
  diag(A) <- 0
  g.sim <- graph_from_adjacency_matrix(A, mode = "directed")
  for (i in 1:nc) {
    g.sim <- g.sim %>%
      set_vertex_attr(name = letters[i], 
                      value = c(Cnew[, i]))
  }
  
  ## Setting covariates in the network
  covVal_orig <- 0
  
  ## Creating binary covariates in the network.
  if (c("binary") %in% covTypes) {
    if (length(CovNamesLPin) > 0) {
      op <- SetBinarycov(g.sim,k_in,N,nc,k_start=0,Cold,connType=1,pC,
                         CovNamesLPin,missing)
    }
    if (length(CovNamesLPout) > 0) {
      op <- SetBinarycov(g.sim,k=k_out,N,nc,k_start=k_in,Cold,connType=1,pC,
                         CovNamesLPout,missing)
      
    }
    g.sim <- op[[1]]
    covVal_orig <- op[[2]]
  }
  
  if (c("continuous") %in% covTypes) {
    if (length(CovNamesLinin) > 0) {
      op <- SetContinuousCov(g.sim, o_in, N, nc, o_start = 0, Cnew, 
                             connType = 1, dist,CovNamesLinin,missing,
                             missType, MARparam)
    }
    if (length(CovNamesLinout) > 0) {
      
      op <- SetContinuousCov(g.sim, o_out, N,nc, o_start = o_in, Cnew, 
                             connType = 1, dist, CovNamesLinout, missing,
                             missType, MARparam)
    }
    g.sim <- op[[1]]
    covVal_orig <- op[[2]]
  }
  g.sim <- g.sim  %>%
    set_vertex_attr(name = "Cluster", value = C)
  return(list(g.sim, covVal_orig, F_u, H_u, Apuv))
}

updateWtmat <- function(G, Ftot, Htot, mode, s, nc, X, Z, k, o, beta, W, alphaLL, 
                        missVals, alpha, start, end, dir, inNeigh, outNeigh, lambda, penalty){
  
  ## Getting sum of all the weights of Fv. Have to modify this for CoDa as opposed to CESNA 
  
  # calculate the sigma sq value
  if (o > 0) {
    sigmaSq <- SigmaSqCalc(Z, beta, Ftot ,Htot, missVals,dir)
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
                             alphaLL,nc, 2, nc+1,mode, dir,lambda, penalty)
    if(dir == "directed"){
      ## Ftot neighbours
      neigh_Htot <- outNeigh[[i]] 
      
      Ftot_neigh <- matrix(Ftot[neigh_Htot, ], ncol = nc)
      
      Htot[i, ] <- CmntyWtUpdt(matrix(Htot[i,],ncol = nc), matrix(Ftot[i,],ncol = nc),
                               Ftot_neigh,Ftot_sum,X_u,W,Z_u,beta,sigmaSq,alpha,
                               alphaLL,nc, nc + 2, (nc*2)+1,"out",dir,lambda, penalty)
    }else{
      Htot[i, ] <- Ftot[i, ]
    }
    
  }
  
  return(list(Ftot,Htot))
}

bb_step_size <- function(x_prev, x_current, grad_prev, grad_current, variant = 1) {
  s <- x_current - x_prev
  y <- grad_current - grad_prev
  
  if (variant == 1) {
    # BB1 variant
    as.numeric((t(s) %*% s) / (t(s) %*% y))
  } else {
    # BB2 variant
    as.numeric((t(s) %*% y) / (t(y) %*% y))
  }
}


CmntyWtUpdt <- function(W1_u ,W2_u ,W2_neigh ,W2_sum ,X_u = NULL ,W = NULL ,
                        Z_u = NULL, beta = NULL ,sigmaSq ,alpha ,alphaLL ,
                        nc ,start ,end ,mode,dir , lambda, penalty ) {
  
  llG <- rep(0,nc)
  
  ## Gradient function to update log likelihood of G
  #First part of Likelihood based on graph only.
  llG <- GraphComntyWtUpdt(W1_u,W2_u,W2_neigh,W2_sum)
  
  #Second part of Likelihood based on binary covariates
  llX <- BincovCmntyWtUpdt(nc*2, W1_u, W2_u, W, X_u, dir)
  
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
  W1_u_new <- W1_u + (alpha * (c(llG) +  alphaLL *(llX[2:(nc+1)] + llZ[1:nc])))
  W1_u_new[(W1_u_new < 0) | (W1_u_new == 0)] <- 0
  
  return(W1_u_new)
  
}


WithBacktracking_CmntyWtUpdt <- function(W1_u ,W2_u ,W2_neigh ,W2_sum ,X_u = NULL ,W = NULL ,
                                         Z_u = NULL, beta = NULL ,sigmaSq ,alpha ,alphaLL ,
                                         nc ,start ,end ,mode,dir, lambda, penalty ) {
  rho <- 0.8
  cv <- 0.1
  llG <- rep(0,nc)
  
  ## Gradient function to update log likelihood of G
  #First part of Likelihood based on graph only.
  llG <- GraphComntyWtUpdt(W1_u,W2_u,W2_neigh,W2_sum)
  
  #Second part of Likelihood based on binary covariates
  llX <- BincovCmntyWtUpdt(nc*2, W1_u, W2_u, W, X_u, dir)
  
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
  
  f_curr <- objective_f(W1_u,W2_u,W2_neigh,W2_sum, beta, Z_u, lambda, penalty)
  grad <- c(llG) +  llX[2:(nc+1)] + llZ[1:nc]
  W_u_prop <- W1_u + (alpha * grad)
  W_u_prop[W_u_prop < 0] <- 0
  f_prop <- objective_f(W_u_prop,W2_u,W2_neigh,W2_sum, beta, Z_u, lambda, penalty)
  
  while(f_prop < (f_curr + (cv * alpha * sum(grad^2)) )){
    alpha <- rho *alpha
    W_u_prop <- W1_u + (alpha * grad)
    W_u_prop[W_u_prop < 0] <- 0
    f_prop <- objective_f(W_u_prop,W2_u,W2_neigh,W2_sum, beta, Z_u, lambda, penalty)
  }
  W1_u_new <- W_u_prop
  return(W1_u_new)
  
}
## future implementation
objective_f <- function(W1_u,W2_u,W2_neigh, W2_sum, beta, Z, lambda, penalty){
  ## Graph section of loglikelihood
  ##part one
  a <- exp(-1 * (W1_u %*% t(W2_neigh)))
  p1 <- sum(log(1 - a))
  
  ##part 2
  p2 <- W1_u %*% as.matrix(W2_sum - t(W2_u) - as.matrix(colSums(W2_neigh)))
  
  ## regression part of loglikelihood
  S <- S_tmp <- 0
  for (j in 1:dim(beta)[1]) {
    S_tmp <- ( Z[,j]  - (c(1,W1_u,W2_u) %*% beta[j, ])) ^ 2
    S <- S - (S_tmp) - lambda * penLL(penalty , beta[j,])
  }
  return(p1 + p2 + S)
}

OldCmntyWtUpdt <- function(W1_u ,W2_u ,W2_neigh ,W2_sum ,X_u = NULL ,W = NULL ,
                           Z_u = NULL, beta = NULL ,sigmaSq ,alpha ,alphaLL ,
                           nc ,start ,end ,mode,dir ) {
  
  llG <- rep(0,nc)
  
  ## Gradient function to update log likelihood of G
  #First part of Likelihood based on graph only.
  llG <- GraphComntyWtUpdt(W1_u,W2_u,W2_neigh,W2_sum)
  
  #Second part of Likelihood based on binary covariates
  llX <- BincovCmntyWtUpdt(nc, W1_u, W2_u, W, X_u, dir)
  #Third part of likelihood based on continuous covariates
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
  W1_u_new <- W1_u + (alpha * (c(llG) +  alphaLL *(llX[1:nc] + llZ[1:nc])))
  W1_u_new[(W1_u_new < 0) | (W1_u_new == 0)] <- 0
  
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
  llG <- scale * llG
  
  return(llG)
}

## Binary covariates community weight
BincovCmntyWtUpdt <- function(nc, f_u,h_u, W, X_u, dir){
  # if(dir == "directed"){
  #   ncoef <- (nc*2)+1
  # }else{
  #   ncoef <- nc+1
  # }
  llX <- rep(0, nc)
  if (length(W) > 0) {
    ## Updating the logistic regression weights
    ## Gradient function to update log likelihood for covariates
    ## adding intercept term
    Q_u <- CalcQuk(f_u,h_u,W,nc, dir)
    ## This math has been verified
    llX <- t((X_u - Q_u) %*% W[,2:(nc+1)])
  }
  return(llX)
}

## Continuous covariates community weight
OldContcovCmntyWtUpdt <- function(nc, Z_u, beta, W1_u,W2_u, sigmaSq, dir){
  
  if(dir =="directed"){
    llZ <- rep(0, nc*2 )
    WW <- matrix(c(1,W1_u,W2_u), ncol = 1+(nc*2))
    
  }else{
    llZ <- rep(0, nc )
    WW <- matrix(c(1,W1_u), ncol = (1+nc))
    
  }
  if (length(beta) > 0) {
    ## Adding the gradient based on continuous covariates
    llZ <- 1 * ((Z_u - t(beta %*% t(WW))) %*% beta[,2:(nc+1)])/sigmaSq
  }
  return(llZ)
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
    llZ <- 1 * ((Z_u - t(beta %*% t(WW))) %*% beta[,2:(nc+1)])#/sigmaSq
  }
  return(llZ)
}

CalcQuk <- function(Fmat, Hmat, W, nc, dir){
  Q_uk <- rep(0, dim(W)[1])
  if(dir == "directed"){
    Q_uk <- 1 / (1 + exp(-1 * (cbind(1, Fmat, Hmat) %*% t(W))))
    
  } else{
    Q_uk <- 1 / (1 + exp(-1 * (cbind(1, Fmat) %*% t(W))))
  }
  
  return(Q_uk)
}
## Logistic covariates parameter updates
sigmoid <- function(W, Ftotn) {
  Q <- 1 / (1 + exp(-1 * (W %*% t(Ftotn))))
  return(Q)
}

LRParamUpdt <- function(W, X, Ftot,Htot, alphaLR, lambda, missing, dir) {
  
  ## Do not add penalty to the intercept hence removing the 1 from the Ftotn matrix
  if(dir =="directed"){
    Ftotn <- cbind(1,Ftot, Htot)
  } else{
    Ftotn <- cbind(1,Ftot)
  }
  N <- dim(X)[1]
  for(i in 1:dim(X)[2]){
    p_i <- sigmoid(W[i,], Ftotn)
    W_grad <- ((p_i - X[,i]) %*% Ftotn[,-c(1)])/N
    W_0 <- sum(p_i - X[,i])/N
    W[i,1] <- W[i,1] - alphaLR*W_0
    W[i, -c(1)] <- sign(W[i,-c(1)] - alphaLR * W_grad) * max(abs(W[i,-c(1)] - alphaLR * W_grad) - alphaLR * lambda, 0)
  }
  
  return(W)
}

PredictCovLR <- function(X, W, Fmat, Hmat, missing,dir, impType) {
  if(dir =="directed"){
    Ftotn <- cbind(1, Fmat, Hmat)
  } else{
    Ftotn <- cbind(1, Fmat)
  }
  size <- 100
  for (i in 1:dim(X)[2]) {
    if(impType == "StochasticReg"){
      X[which(missing[, i]), i] <- as.numeric((rbinom(1,size,sigmoid(W[i, ], Ftotn[which(missing[, i]), ]))/size) > 0.5)
    }else{
      X[which(missing[, i]), i] <- as.numeric(sigmoid(W[i, ], Ftotn[which(missing[, i]), ]) > 0.5)## This uses sigmoid function 
    }
  }
  return(X)
}

updateLogisticParam <- function(W,X,Wtm,Wtm2,missVals,lambda,alphaLR,dir,impType){
  Wret <- W
  Xret <- X
  if (length(W) > 0) {
    if (sum(missVals) > 0) {
      Wret <-  LRParamUpdt(W, X, Wtm,Wtm2 ,alphaLR, lambda, missVals, dir)
      Xret <-  PredictCovLR(X, W, Wtm,Wtm2, missVals, dir, impType)
    } else{
      Wret <-  LRParamUpdt(W, X, Wtm,Wtm2, alphaLR, lambda, missVals, dir)
    }
  }
  
  return(list(Wret,Xret))
}

## Continuous covariates parameter updates
SigmaSqCalc <- function(Z, beta, Ftot, Htot, missVals,dir) {
  sigmaSq <- rep(0, dim(beta)[1])
  if(dir == "directed"){
    Fm <- cbind(1,Ftot, Htot)
  }else{  Fm <- cbind(1,Ftot)}
  
  for (i in 1:dim(beta)[1]) {
    z_sum <- Z[!missVals[, i], i] - (t(beta[i, ]) %*% t(Fm[!missVals[, i], ]))
    sigmaSq[i] <- sum(z_sum ^ 2, na.rm = TRUE) / sum(!missVals[, i])
    if(is.infinite(sigmaSq[i])){
      print("Null Sigma")
    }
  }
  return(sigmaSq)
}

sdErr <- function(Z_i, beta, Fm, N, p, k){
  sdErr <- sqrt(sum((Z_i -  Fm %*% as.matrix(beta))^2)/ (N-p-k))
  return(sdErr)
}

PredictCovLin <- function(Ftot,Htot, Z, beta, missVals, dir, impType) {
  Z_out <- Z
  if(dir == "directed"){
    Fm <-  cbind(1, Ftot, Htot)
  }else{  
    Fm <-  cbind(1, Ftot)
  }
  for (i in 1:dim(Z_out)[2]) {
    if(sum(missVals) > 0 ) {
      idx <- which(missVals[, i])
      if(impType == "StochasticReg"){
        s <- sdErr(Z_out[-idx, i], beta[i, ], Fm[-idx,], dim(Z_out)[1], ncol(Z_out), dim(Z_out[idx,])[1])
        err <- mean(rnorm(10,mean = 0,sd = s))
        Z_out[idx, i] <- (Fm[idx,] %*% as.matrix(beta[i, ])) + err
      }else if(impType == "Reg"){
        Z_out[idx, i] <- Fm[idx,] %*% as.matrix(beta[i, ])
      }
    }
  }
  return(Z_out)
}

safe_sign <- function(x) ifelse(x == 0, 0, sign(x))

grad_norm <- function(grad) {
  sqrt(sum(grad^2))  # L2 norm
}

# Soft-thresholding function
soft_threshold <- function(z, lambda) {
  sign(z) * pmax(abs(z) - lambda, 0)
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

#AllUpdatemethods
LinRParamUpdt <- function(beta, Z, Fmat,Hmat, alpha, lambda, N, missVals,dir, 
                          impType, alphaLin, penalty ) {
  if(printFlg == TRUE){
    print("In Linparamupdate Func")
  }
  
  beta_new <- beta
  if(dir =="directed"){
    Fm_n <- cbind(1,Fmat, Hmat)
    X <- cbind(Fmat,Hmat)
  }else{  
    Fm_n <- cbind(1, Fmat)
  }
  ## setting group size 
  p_g <- 2
  
  sigmaSq <- SigmaSqCalc(Z, beta, Fmat, Hmat, missVals,dir)
  for (idx in 1:ncol(Z)) {
    if(penalty == "LASSO"){
      
      beta_new[idx, 1] <- mean(Z[, idx] - (Fm_n[, -1] %*% beta_new[idx, -1]))
      
      for (j in 2:ncol(Fm_n)) {
        ## residual and soft thresholding
        r_j <- Z[,idx] - (Fm_n[,-j] %*% beta_new[idx ,-j])
        z_j <- sum(Fm_n[, j] * r_j)
        d_j <- sum(Fm_n[,j]^2)
        s_j <- soft_thresh(z_j,lambda)
        beta_new[idx,j] <- s_j/d_j
      }
      
    }else if(penalty =="ElasticNet"){
      mixP <- 0.5
      beta_new[idx, 1] <- mean(Z[, idx] - (Fm_n[, -1] %*% beta_new[idx, -1]))
      
      for (j in 2:ncol(Fm_n)) {
        # partial residual w/o feature j
        r_j <- Z[,idx] - (Fm_n[,-j] %*% beta_new[idx,-j]) 
        z_j <- sum(Fm_n[, j] * r_j)
        
        # Update beta_j with soft-thresholding
        denom <- sum(Fm_n[, j]^2) + lambda * (1 - mixP) * N
        beta_new[idx,j] <- sign(z_j) * pmax(abs(z_j) - lambda * mixP * N, 0)/denom
      }
      
    }else if(penalty  == "GroupLASSO"){
      beta_new[idx, 1] <- mean(Z[, idx] - (Fm_n[, -1] %*% beta_new[idx, -1]))
      r <- Z[,idx] - Fm_n[,-1] %*% beta_new[idx,-1] - beta_new[idx,1]
      for (g in 2:(ncol(Fmat)+1)) {
        i <- c(g, g+3)
        X_g <- Fm_n[, i]
        beta_g <- beta_new[idx, i]
        
        # Residual excluding group g (using current intercept)
        r <- r + X_g %*% beta_g
        
        # Least squares sol for g
        z_g <- crossprod(X_g,r)/N
        norm_z <- sqrt(sum(z_g^2))
        if(norm_z >lambda){
          shirnkage <- 1-lambda/norm_z
          beta_g_new <- shirnkage * z_g
        }else{
          beta_g_new <- rep(0, length(i))
        }
        beta_new[idx,i] <-beta_g_new
        r <- r - X_g %*% beta_g_new
      }
      
      
    }else if(penalty == "Ridge"){
      beta_new[idx, 1] <- mean(Z[, idx] - Fm_n[, -1] %*% beta_new[idx, -1])
      
      for (j in 2:ncol(Fm_n)) {
        # Partial residual w/o  j
        r_j <- Z[,idx] - Fm_n[,-j] %*% beta_new[idx,-j]
        #L2 penalty
        xj <- Fm_n[, j]
        beta_new[idx, j] <- sum(xj * r_j) / (sum(xj^2) + N * lambda)
      }
      
    }else if(penalty == "GroupLASSOProxGrad"){
      #Current parameters
      beta_current <- beta[idx,]
      intercept <- beta_current[1]
      weights <- beta_current[-1]
      
      # Compute gradient (only for the smooth part)
      residual <- Z[,idx] - intercept - X %*% weights
      grad_intercept <- -mean(residual)
      grad_weights <- (-1/N) * crossprod(X, residual)
      
      # Gradient step
      beta_temp <- beta_current - alpha * c(grad_intercept, grad_weights)
      
      # Proximal operator (group soft thresholding)
      beta_new <- beta_temp
      # Don't penalize intercept
      for (g in 2:(ncol(Fmat)+1)) {
        ind <- c(g, g+3)
        #for (g in seq_along(group_indices)) {
        # ind <- group_indices[[g]] + 1  # +1 because of intercept
        group_norm <- sqrt(sum(beta_temp[ind]^2))
        if (group_norm <= lambda * alpha) {
          beta_new[ind] <- 0
        } else {
          beta_new[ind] <- beta_temp[ind] * (1 - (lambda * alpha)/group_norm)
        }
      }
      ## update the beta
      beta[idx, ] <- beta_new
      
      
      # Transform coefficients back to original scale
      #intercept <- beta_new[1] + mean(Z[,idx]) - sum(beta_new[-1] * X_means / X_sds)
      #scaled_coef <- beta_new[-1] / X_sds
      #beta[idx,] <- c(intercept, scaled_coef)
      
    }else{
      print("Please provide the penalty method")
      break}
  }
  
  if(sum(is.infinite(beta_new)) > 0){
    print("Beta_new is inifinite")
  }
  return(beta_new)
}


updateLinearRegParam <- function(beta,missVals,Z,Wtm1,Wtm2, alpha,lambda,N,dir, 
                                 impType, alphaLin, penalty){
  
  if(dir == "directed"){
    cm <- cbind( Wtm1 , Wtm2)
  }else{
    cm <- Wtm1
  }
  ##FOR GRUOP LASSO
  ## GROUPLASSO CODE IS DISABLED
  #cm <- matrix(0, nrow = N, ncol =0)
  #for(i in 1:ncol(Wtm1)){ cm <- cbind(cm, Wtm1[,i], Wtm2[,i])}
  
  betaret <- beta
  Zret <-  Z
  nc <- ncol(Wtm1)
  sigmaSq <- NULL
  mod <- list()
  if (length(beta) > 0) {
    for(i in 1:dim(Z)[2]){
      ## "GroupLASSO","Ridge","LASSO","ElasticNet")
      ## filter out the missing data so we update only based on the avaiable data
      if(sum(missVals) > 0){
        mIdx <- missVals[,i]
        X <-  cm[-mIdx,]
        y <- Z[-mIdx,i]
      }else{
        X <- cm
        y <- Z[,i]
      }
      
      if(penalty =="GroupLASSO"){
        ## GROUPLASSO CODE IS DISABLED
        ##mod[[i]] <- gglasso(cm[-mIdx,], Z[-mIdx,i], group = rep(1:3, each= 2), lambda = lambda)
        ##beta[i,] <- unname(as.matrix(coef(mod[[i]]))[,1])
        
        {##Find the R-sq value
          # Get predictions
          #preds <- predict(mod[[i]], newx = cm)
          # For regression (continuous y)
          #mse <- mean((Z[,i] - preds)^2)
          #rsq <- 1 - sum((Z[,i]-preds)^2)/sum((Z[,i]-mean(Z[,i]))^2)
          # print(paste("r-squared for ", i, " :",rsq, " mse ", mse))
        }
      }else if(penalty =="LASSO"){
        #mod[[i]] <- cv.glmnet(X , y, alpha = 1, maxit = 1000) 
        #beta[i,] <- coef(mod[[i]], s = "lambda.min")
        
        mod[[i]] <- glmnet(X ,y,lambda = lambda, alpha = 1, maxit = 1000) 
        beta[i,] <- as.matrix(coef(mod[[i]]))[,1]
        
      }else if(penalty =="ElasticNet"){
        #mod[[i]] <- cv.glmnet(X , y, alpha = 0.5, maxit = 1000) 
        #beta[i,] <- coef(mod[[i]], s = "lambda.min")
        
        mod[[i]] <- glmnet(X ,y ,lambda = lambda, alpha = 0.5, maxit = 1000) 
        beta[i,] <- as.matrix(coef(mod[[i]]))[,1]
        
      }else if(penalty == "Ridge"){
        #mod[[i]] <- cv.glmnet(X , y, alpha = 0, maxit = 1000) 
        #beta[i,] <- as.matrix(coef(mod[[i]], s = "lambda.min"))
        
        mod[[i]] <- glmnet(X ,y , family = "gaussian",lambda = lambda, alpha = 0, maxit = 1000) 
        beta[i,] <- as.matrix(coef(mod[[i]]))
        
      }
    }
    
    ## GROUP LASSO CODE IS DISABLED 
    if(penalty == "GroupLASSO"){
      betaret <- cbind(beta[,1], beta[, seq(2,((nc*2)+1),by=2)], beta[,seq(3,((nc*2)+1), by=2)])
    } else{ betaret <- beta}
    
    sigmaSq <-  SigmaSqCalc(Z, betaret, Wtm1,Wtm2, missVals,dir)
    
    if (sum(missVals) > 0) {
      Zret <- PredictCovLin(Wtm1, Wtm2, Z, betaret, missVals,dir, impType)
    }
  }
  return(list(betaret,Zret,sigmaSq))
}



OldupdateLinearRegParam <- function(beta,missVals,Z,Wtm1,Wtm2, alpha,lambda,N,dir, 
                                    impType, alphaLin, penalty){
  if(printFlg == TRUE){
    print("In updateLinearRegParam Func")
  }
  betaret <- beta
  Zret <-  Z
  sigmaSq <- NULL
  if (length(beta) > 0) {
    betaret <- LinRParamUpdt(beta, Z, Wtm1,Wtm2, alpha, lambda, N, missVals,dir, impType, alphaLin, penalty)
    sigmaSq <-  SigmaSqCalc(Z, betaret, Wtm1,Wtm2, missVals,dir)
    
    if (sum(missVals) > 0) {
      Zret <- PredictCovLin(Wtm1, Wtm2, Z, beta, missVals,dir, impType)
    }
  }
  return(list(betaret,Zret,sigmaSq))
}

## Log likleihood calculations
Lx_cal <- function(W, N, Fmat, Hmat, X, dir){
  S <- 0
  
  if (length(W) > 0) {
    cases <- complete.cases(X)
    Q <- CalcQuk(Fmat, Hmat, W, dir)
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
    sigmaSq <- c(1,1,1)
    
    for (j in 1:dim(beta)[1]) {
      S_tmp <- sum(( Z[,j]  - (Fin %*% beta[j, ])) ^ 2)
      S <- S - (N*log(2*pi*sigmaSq[j])/2) - (S_tmp / (2 * sigmaSq[j])) - lambda * penLL(penalty , beta[j,])
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
                      epsilon, lambda,missVals, penalty, alphaLin) {
  if(printFlg == TRUE){
    print("In Find LLDir Func")
  }
  
  N <-  gorder(G)
  
  ## Calculating the L_g part of the log likelihood for structure of the network
  S1 <- Lg_cal(G, Ftot, Htot, Eneg,epsilon, dir)
  
  ## Calculating the L_xin part of the log likelihood for binary covariates
  S2_in <- 0
  
  ## Calculating the L_xout part of the log likelihood for binary covariates
  S2_out <- Lx_cal(W=Wout, N = N, Fmat= Ftot, Hmat =Htot, X=X_out, dir)
  
  ## adding the loglikelihood from the covariates
  S2 <- S2_in + S2_out
  
  ## Calculating the L_zin part of the log likelihood for continuous covariates
  S3_in <- 0
  
  ## Calculating the L_zout part of the log likelihood for continuous covariates
  S3_out <- Lz_cal(beta = betaout,N = N,Z = Z_out,
                   Fmat = Ftot, Hmat = Htot,
                   sigmaSq = sigmaSqout,dir,lambda, missVals, penalty, alphaLin)
  S3 <- S3_in + S3_out
  
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

initWtmat <- function(G, mode, N, nc, seed,specOP){
  set.seed(seed)
  colwts <- igraph::degree(G, mode = mode) / sum(igraph::degree(G, mode = mode))
  Wtmat <- replicate(nc, colwts) + matrix(nrow = N, ncol = nc, runif(nc * N, 0.0001,0.001))
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

CoDA <- function(G,nc, k = c(0, 0) ,o = c(0, 0) , N,  alpha, lambda, thresh, 
                 nitermax, orig,randomize = TRUE, CovNamesLinin = c(),
                 CovNamesLinout = c(),CovNamesLPin = c(),CovNamesLPout = c(), 
                 dir, alphaLL = NULL,test = TRUE,missing = NULL, covOrig, 
                 epsilon,  impType, alphaLin, penalty, seed, covInit, specOP) {
  if(printFlg == TRUE){
    print("In CoDA Func")
  }
  
  if(dir == "directed"){
    ncoef = nc*2
  }else{ncoef = nc}
  
  ## Community weights outgoing connections
  Ftot <- initWtmat(G,"out",N,nc, seed, specOP)
  Finit <- Ftot
  
  ## Community weights outgoing connections
  if(dir == "directed"){
    Htot <- initWtmat(G,"in",N,nc, (seed+10), specOP)
    Hinit <- Htot
  } else{
    Htot <- Ftot
    Hinit <- Htot
  }
  
  set.seed(seed)
  
  iter <- 0
  k_in <- k[1]
  k_out <- k[2]
  o_in <- o[1]
  o_out <- o[2]
  
  ########## Get the covariate matrix for binary data and continuous data ##########
  b <- getNWcov(G, CovNamesLinin,CovNamesLinout,CovNamesLPin,CovNamesLPout,k_in,k_out,o_in,o_out)
  
  X_in <- b[[1]]
  X_out <- b[[2]]
  
  Z_in <- b[[3]]
  Z_out <- b[[4]]
  
  missValsin <- b[[5]]
  missValsout <- b[[6]]
  numMissValin <- b[[7]]
  numMissValout <- b[[8]]
  
  ## For covariates where the value is missing, fill it in initially with random values
  if(k_out >0){
    for(i in 1:dim(X_out)[2]){
      X_out[missValsout[,i],i] <- rbinom(n = sum(missValsout[,i]), size = 1, prob = (N- sum(missValsout[,i]))/N )
    }
  }
  
  backTransParams <- list()
  if(o_out > 0){
    
    for(i in 1:dim(Z_out)[2]){
      
      if(sum(missValsout[,i]) > 0){
        
        if(covInit == "mean"){
          ## Mean imputation of initial values
          Z_out[missValsout[,i],i] <- mean(Z_out[,i], na.rm =TRUE) 
        }else if(covInit == "Nmode"){
          ## Imputing values based on the neighbours. Using mode of the neighbour values
          mv <- which(missValsout[,i])
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
          mv <- which(missValsout[,i])
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
          mv <- which(missValsout[,i])
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
          mv <- which(missValsout[,i])
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
  
  ############ Initialize the coef matrix for covariates ##############
  
  ############## Setting initial values of binary covariates coefs weights ###############
  Win <- NULL
  Win_orig <- Win
  
  Wout <- initW(k_out, nc,missValsout, X_out, Htot, lambda, alpha, ncoef)
  Wout_orig <- Wout
  
  #################### Setting initial values of linear regression weights
  b <- NULL
  betain <- b
  sigmaSqin <- b
  
  b <- initbeta(o_out, ncoef, missValsout, Z_out, Ftot, Htot, alpha, N, lambda, dir, impType, alphaLin, penalty)
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
  Eneg  <- as_edgelist(graph_from_adjacency_matrix(Aneg, mode = dir))
  
  ##Initial log likelihood value
  LLold <- -10000
  lllst <- matrix(nrow = 0, ncol = 4)
  
  arilst <- c()
  OmegaVal <- c()
  mse <- matrix(nrow = 0, ncol = (o_in+o_out))
  mseMD <- matrix(nrow = 0, ncol = (o_in + o_out))
  
  accuracy <- matrix(nrow = 0, ncol = k_in + k_out)
  accuracyMD <- matrix(nrow = 0, ncol = k_in + k_out)
  
  corMat <- matrix(0, nrow =0 , ncol = 15)
  
  delta <- getDelta(N)
  continue <- TRUE
  s <- 1:N
  beta_old <- betaout
  
  repeat{  
    ## Update the log likelihood.
    LLvec <- findLLDir(G,Ftot,Htot, Win,Wout,X_in,X_out, betain,betaout,
                       Z_in,Z_out, sigmaSqin, sigmaSqout,alphaLL, Eneg, 
                       dir, epsilon,lambda, missValsout, penalty, alphaLin)
    LLnew <- LLvec[1]
    pctInc <- abs((LLnew - LLold) /(LLold)) * 100
    if(is.nan(pctInc)){
      if(test == FALSE){
        
        arilst <- c(arilst, ARIop(Ftot,Htot,orig,nc,N))
        OmegaVal <- c(OmegaVal, OmegaIdx(G, Ftot, Htot, N, delta,nc))
        
        MSEtmp <- MSEop(Ftot, Htot, covOrig, betaout,N, dir, o_in,o_out, missValsout)
        mseMD <- rbind(mseMD, MSEtmp[[1]])
        mse <- rbind(mse, MSEtmp[[2]])
        
        accutmp <- accuOP(k_in, k_out,Wout,Ftot, Htot,covOrig,dir,missValsout)
        accuracyMD <- rbind(accuracyMD, accutmp[[1]])
        accuracy <- rbind(accuracy, accutmp[[2]])
        
      }
      FAILURE <- TRUE
      print(paste0("ERROR in Pct increase ,total iterations ", 
                   iter, " with seed value ", seed, " Final OI ",OmegaIdx(G, Ftot, Htot, N, delta, nc)))
      break
    }
    
    if((pctInc < thresh) | (dim(lllst)[1] == nitermax)){
      if(test == FALSE){
        lllst <- rbind(lllst, LLvec)
        
        arilst <- c(arilst, ARIop(Ftot,Htot,orig,nc,N))
        OmegaVal <- c(OmegaVal, OmegaIdx(G, Ftot, Htot, N, delta,nc))
        
        MSEtmp <- MSEop(Ftot, Htot, covOrig, betaout,N, dir, o_in,o_out, missValsout)
        mseMD <- rbind(mseMD, MSEtmp[[1]])
        mse <- rbind(mse, MSEtmp[[2]])
        
        accutmp <- accuOP(k_in, k_out,Wout,Ftot, Htot,covOrig,dir,missValsout)
        accuracyMD <- rbind(accuracyMD, accutmp[[1]])
        accuracy <- rbind(accuracy, accutmp[[2]])
        
      }
      FAILURE <- FALSE
      print(paste0("The final percent change ", round(pctInc,6), " ,total iterations ", 
                   iter, " with seed value ", seed, " Final OI ",OmegaIdx(G, Ftot, Htot, N, delta, nc), 
                   " MSE ", paste0(round(tail(mse,1),3), collapse = ", "),
                   " MSE MD ", paste0(round(tail(mseMD,1),3), collapse = ", ")))
      break
    }
    
    LLold <- LLnew
    lllst <- rbind(lllst, LLvec)
    iter <- iter + 1
    if(dir == "directed"){
      corMat <- rbind(corMat, findCor(Ftot,Htot,"directed"))
    }
    
    ## Randomize the community weights update
    s <- randomizeIdx(N, randomize)
    
    ## Updating F i.e. all outgoing connections
    ## We look for neighbors our node is connected to with the edges directed outward from our node.
    opMat <- updateWtmat(G,Ftot,Htot,"all",s,nc,X_out,Z_out,k_out,o_out,
                         betaout,Wout, alphaLL,missValsout,alpha, start = 2, 
                         end = nc + 1, dir, inNeigh, outNeigh, lambda, penalty)
    Ftot <- opMat[[1]]
    Htot <- opMat[[2]]
    
    ## Updating logistic paramters
    
    b2 <- updateLogisticParam(W = Wout, X = X_out , Wtm = Ftot, Wtm2 = Htot, 
                              missVals = missValsout, lambda = lambda, 
                              alphaLR = alpha, dir = dir, impType)
    Wout <- b2[[1]]
    X_out <- b2[[2]]
    
    ### Updating the beta matrix for continuous covariates
    c1 <- updateLinearRegParam(betaout,missValsout,Z_out,Ftot,Htot,alpha,
                               lambda,N,dir, impType, alphaLin, penalty)
    betaout <- c1[[1]]
    Z_out <- c1[[2]]
    sigmaSqout <- c1[[3]]
    
    Z_in <- 0
    
    if(test == TRUE){
      
      arilst <- c(arilst, ARIop(Ftot,Htot,orig,nc,N))
      OIn <- OmegaIdx(G, Ftot, Htot, N, delta,nc)
      OmegaVal <- c(OmegaVal, OIn)
      MSEtmp <- MSEop(Ftot, Htot, covOrig, betaout,N, dir, o_in,o_out, missValsout)
      mseMD <- rbind(mseMD, MSEtmp[[1]])
      mse <- rbind(mse, MSEtmp[[2]])
      
      accutmp <- accuOP(k_in, k_out,Wout,Ftot, Htot,covOrig,dir,missValsout)
      accuracyMD <- rbind(accuracyMD, accutmp[[1]])
      accuracy <- rbind(accuracy, accutmp[[2]])
      
    }
  }
  
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
      betaout = betaoutfin,
      Zout_cov = Z_out,
      MSE = mse,
      Acc = accuracy,
      OmegaIndex = OmegaVal,
      AccMD = accuracyMD,
      MSEMD = mseMD,
      FAILURE,
      backTransParams,
      corMat
    )
  )
}
