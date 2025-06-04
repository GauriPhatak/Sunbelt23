library(tidyverse)
library(dplyr)
library(igraph)
library(network)
library(intergraph)
library(ergm)
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
  #VIFval <- matrix(0, nrow = 0, ncol = ncoef+1)
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
      #VIFval <- rbind(VIFval, c(i, c(vif(mod))))
      binLL <- binLL + as.numeric(logLik(mod))
      Wout[i,] <- coef(mod)
    }
    # VIFval <- vif(mod)
    
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
    #VIFval <- vif(mod)
  }
  
  A <- 1 - (1 * as.matrix(as_adjacency_matrix(G)))
  Eneg  <- as_edgelist(igraph::graph_from_adjacency_matrix(A, mode = dir))
  
  GrphLL <- Lg_cal(G = G, Ftot = Fm, Htot = Hm,Eneg = Eneg, epsilon = epsilon, dir = dir)
  ll <- binLL + contLL + GrphLL
  
  # findLLDir(G, Ftot = Fm, Htot= Hm,Win,Wout,X_in,X_out,
  #           betain,betaout,Z_in,Z_out,
  #           sigmaSqin,sigmaSqout,alphaLL,Eneg, 
  #           dir, epsilon,lambda,missVals)
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

SetContinuousCov <- function(G,o, N,nc, o_start, C, connType, 
                             dist,CovNamesLin,missing,missType, MARparam=NULL){
  ## Based on the mean and variance indicated by the dist matrix we assign continuous values to the network covariates.
  contVal <- matrix(ncol = o, nrow = N, 0) 
  s <- 1 
  # for(k in 1:o){
  #   for (i in 1:nc) {
  #     idx <- which(C[, i] == connType)
  #     contVal[idx, k] <- rnorm(length(idx), dist[[s]][1],dist[[s]][2])
  #     s <- s + 1
  #   }
  # }
  
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


genBipartite <- function(N,nc,pClust,k_in,k_out,o_in,o_out,pC,dist,covTypes,
                         CovNamesLPin,CovNamesLPout,CovNamesLinin,CovNamesLinout,
                         pConn,dir,dirPct,epsilon,missing,a,b,Type,pClustOL,
                         missType= rep("Random", sum(k_in,k_out,o_in,o_out)), MARparam=NULL, seed){
  
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
  
  ## Creating continuous covariates in the network.
  if (c("continuous") %in% covTypes) {
    if (length(CovNamesLinin) > 0) {
      op <- SetContinuousCov(g.sim, o_in, N, nc, o_start = 0, Cold, 
                             connType = 1, dist,CovNamesLinin,missing,
                             missType, MARparam)
    }
    if (length(CovNamesLinout) > 0) {
      
      op <- SetContinuousCov(g.sim, o_out, N,nc, o_start = o_in, Cold, 
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

updateWtmatBothTypes <- function(G, Ftot, Htot, mode, s, nc, X, Z, k, o, beta, W, alphaLL, missVals, alpha, start, end, dir){
  
  ## Getting sum of all the weights of Fv. Have to modify this for CoDa as opposed to CESNA 
  
  # calculate the sigma sq value
  if (o > 0) {
    if(dir == "directed"){
      #beta_t <- cbind(beta[,1], beta[,start:end], beta[,(start-nc):(end-nc)])
      sigmaSq <- SigmaSqCalc(Z, beta, Ftot ,Htot, missVals,dir)
    }else{
      sigmaSq <- SigmaSqCalc(Z, beta, Ftot, Htot , missVals,dir)
    }
  }
  
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
    neigh_Ftot <- igraph::neighbors(G, i,mode)
    
    ## Htot neighbours
    Htot_neigh <- matrix(Htot[neigh_Ftot, ], ncol = nc)
    Htot_sum <- as.matrix(colSums(Htot))
    
    Ftot[i, ] <- CmntyWtUpdt(matrix(Ftot[i,], ncol = nc), matrix(Htot[i,], ncol = nc),
                             Htot_neigh, Htot_sum,X_u,W,Z_u, beta,sigmaSq,alpha,
                             alphaLL,nc, 2, nc+1,mode, dir)
    if(dir == "directed"){
      ## Ftot neighbours
      neigh_Htot <- igraph::neighbors(G, i,"out")
      
      Ftot_neigh <- matrix(Ftot[neigh_Htot, ], ncol = nc)
      Ftot_sum <- as.matrix(colSums(Ftot))
      
      Htot[i, ] <- CmntyWtUpdt(matrix(Htot[i,],ncol = nc),matrix(Ftot[i,],ncol = nc),
                               Ftot_neigh,Ftot_sum,X_u,W,Z_u,beta,sigmaSq,alpha,
                               alphaLL,nc, nc+2, (nc*2)+1,"out",dir)
    }else{
      Htot[i, ] <- Ftot[i, ]
    }
    
  }
  
  return(list(Ftot,Htot))
}

## Community weight updates
updateWtmat <- function(G, Wtm1, Wtm2, mode, s, nc, X, Z, k, o, beta, W, alphaLL, missVals, alpha, start, end, dir){
  if(printFlg == TRUE){
    print("In updateWtmat Func")
  }
  Wtmat <- Wtm1
  
  ## Getting sum of all the weights of Fv. Have to modify this for CoDa as opposed to CESNA 
  Wtm2_sum <- as.matrix(colSums(Wtm2))
  # calculate the sigma sq value
  
  if (o > 0) {
    if(dir == "directed"){
      if(mode == "in"){
        #beta_t <- cbind(beta[,1], beta[,start:end], beta[,(start-nc):(end-nc)])
        sigmaSq <- SigmaSqCalc(Z, beta, Wtm2,Wtm1, missVals,dir)
      }else{
        sigmaSq <- SigmaSqCalc(Z, beta, Wtm1,Wtm2, missVals,dir)
      }
    }else{
      sigmaSq <- SigmaSqCalc(Z, beta, Wtm1,Wtm2, missVals,dir)
    }
    #sigmaSq <- 1
  }
  
  for (i in s) {
    #print(i)
    Wtm1_u <- matrix(Wtm1[i, ], ncol = nc)
    Wtm2_u <- matrix(Wtm2[i, ], ncol = nc)
    
    ##The neighbours of this node
    neigh <- igraph::neighbors(G, i, mode)
    Wtm2_neigh <- matrix(Wtm2[neigh, ], ncol = nc)
    
    X_u <- NULL
    if (k > 0) {
      X_u <- matrix(X[i, ], ncol = k)
    }
    
    Z_u <- NULL
    if(o > 0 ){
      Z_u <- matrix(Z[i, ], ncol = o)
    }
    
    Wtmat[i, ] <- CmntyWtUpdt(Wtm1_u,Wtm2_u,Wtm2_neigh, Wtm2_sum,
                              X_u,W,Z_u, beta,sigmaSq,alpha,alphaLL,
                              nc, start, end, mode, dir)
    
    if(sum(is.infinite(Wtm2_u)) > 0){
      print("infinite h_u")
    }
  }
  #Wtmat <- apply(Wtmat,2,range01) 
  #reg_t_23 <<- rbind(reg_t_23, c(Wtmat[23,]))
  #nc_t_23 <<- rbind(nc_t_23, c(Wtmat[23,]))
  
  #reg_t_92 <<- rbind(reg_t_92, c(Wtmat[92,]))
  #reg_t_97 <<- rbind(reg_t_97, c(Wtmat[97,]))
  
  if(sum(is.nan(Wtmat)) > 0){
    print("some weight value is nan")
  }
  return(Wtmat)
}

CmntyWtUpdt <- function(W1_u ,W2_u ,W2_neigh ,W2_sum ,X_u = NULL ,W = NULL ,
                        Z_u = NULL, beta = NULL ,sigmaSq ,alpha ,alphaLL ,
                        nc ,start ,end ,mode,dir ) {
  
  lb <- 0.00001
  ub <- 1
  
  llG <- rep(0,nc)
  
  ## Gradient function to update log likelihood of G
  #First part of Likelihood based on graph only.
  llG <- GraphComntyWtUpdt(W1_u,W2_u,W2_neigh,W2_sum)
  
  #Second part of Likelihood based on binary covariates
  #if(dir == "directed"){
  #if(mode == "out"){
  llX <- BincovCmntyWtUpdt(nc*2, W1_u, W2_u, W, X_u, dir)
  #}else if(mode == "in"){
  #  W <- cbind(W[,1], W[,start:end], W[,(start-nc):(end-nc)])
  #  llX <- BincovCmntyWtUpdt(nc*2, W1_u, W2_u, W, X_u, dir)
  #}
  #} else{
  #  llX <- BincovCmntyWtUpdt(nc, W1_u, W2_u, W, X_u, dir)
  #}
  
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

## Continuous covariates community weight
ContcovCmntyWtUpdt <- function(nc, Z_u, beta, W1_u,W2_u, sigmaSq, dir){
  if (length(beta) > 0) {
    if(dir =="directed"){
      llZ <- rep(0, nc*2 )
      WW <- matrix(c(1,W1_u,W2_u), ncol = 1+(nc*2))
      
    }else{
      llZ <- rep(0, nc)
      WW <- matrix(c(1,W1_u), ncol = (1+nc))
      ## Adding the gradient based on continuous covariates
      #llZ <- 1 * ((Z_u - t(beta %*% t(WW))) %*% beta[,2:])/sigmaSq
    }
    ## Adding the gradient based on continuous covariates
    llZ <- 1 * ((Z_u - t(beta %*% t(WW))) %*% beta[,2:(nc+1)])/sigmaSq
  }
  return(llZ)
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
  Q <- sigmoid(W, Ftotn)
  
  W_grad <- t(X - t(Q)) %*% Ftotn
  W_new <- W
  W_new[,-c(1)] <- W[,-c(1)] - alphaLR * (W_grad[,-c(1)] + lambda * (sign(W[,-c(1)])))
  return(W_new)
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
  if(dir =="directed"){
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

LinRParamUpdt1.0 <- function(beta, Z, Fmat,Hmat, alpha, lambda, N, missVals,dir, 
                             impType, alphaLin, penalty ) {
  if(printFlg == TRUE){
    print("In Linparamupdate Func")
  }
  beta_new <- beta#matrix(0, nrow = dim(beta)[1], ncol = dim(beta)[2])
  if(dir =="directed"){
    Fm_n <- cbind(1,Fmat, Hmat)
  }else{  
    Fm_n <- cbind(1, Fmat)}
  
  weights <- c(0, 1/apply(Fm_n[, -1], 2 , sd))
  
  # Try rescaling the values just for linear regression
  #Fm_n <- cbind(1, apply(Fm_n[,-1],2,range01) )
  
  ## Calculating update to the beta using gradient descent
  #y <-  Fm_n %*% t(as.matrix(beta)) 
  #gn <- c()
  
  #print("new iter")
  for(idx in 1:dim(beta)[1]){
    lambda <- round(abs(t(Fm_n[,-1]) %*% Z[,idx])/N, 2)
    #print(c(lambda))
    #print(paste("The max lambda value for cov :",idx, " is ", paste(lambda,collapse = ",")))
    lambda_max <- max(lambda)
    # maxLambda <<- c(maxLambda, lambda_max)
    #-1 * t( t(Fm_n) %*% (Z[,idx] - matrix(beta[idx,1],ncol=1, nrow= N ) - (Fm_n[,-1] %*% (beta[idx,-1 ])))) /N
    ## LASSO regression 
    # for(graditer in 1){
    gradient <- -1 * t( t(Fm_n) %*% (Z[,idx] - (Fm_n %*% beta_new[idx, ]) ) ) /N
    #gn <- c( gn, grad_norm(gradient))
    
    if(penalty == "LASSO") {
      #EN_fit <- cv.glmnet(Fm_n[,-1], Z[,idx], alpha = 1, lambda = seq(0.001,1,length = 100))
      #lambda.min <-  EN_fit$lambda.min
      #beta_new[idx,] <- as.numeric(coef(glmnet(Fm_n[,-1], Z[,idx], alpha = 1,lambda = 0.1)))
      beta_new[idx,] <-  beta_new[idx,] - ( alphaLin * (gradient + (lambda_max * weights * safe_sign(beta_new[idx,]) ) ) )
      #beta_new[idx,] - ( alphaLin * (gradient + (lambda * weights * safe_sign(beta_new[idx,]) ) ) )
      #beta_new[idx,] - ( alphaLin * (gradient + (lambda * safe_sign(beta_new[idx,]) ) ) )
      #beta_new[idx,] - ( alphaLin * (gradient + (lambda * weights * safe_sign(beta_new[idx,]) ) ) )
      #beta - ( alphaLin * (gradient + (lambda * sign(beta) ) ) )
      if(is.nan(beta_new[1,1])){
        print("nan beta")}
    }else if(penalty == "Ridge"){
      beta_new <-  beta - (alphaLin * ( gradient + (2*lambda * beta ) ) )
    } else if(penalty == "Proximal"){
      beta_temp <-  beta - ( alphaLin * gradient ) 
      beta_new <- sign(beta_temp) * max( abs(beta_temp) - alphaLin * lambda , 0)
    }
    # }
  }
  #print(gn)
  #print(beta_new)
  #GN <- rbind(GN, gn)
  
  if(sum(is.infinite(beta_new)) > 0){
    print("Beta_new is inifinite")
  }
  return(beta_new)
}

LinRParamUpdt2.0 <- function(beta, Z, Fmat,Hmat, alpha, lambda, N, missVals,dir, 
                             impType, alphaLin, penalty ) {
  if(printFlg == TRUE){
    print("In Linparamupdate Func")
  }
  beta_new <- beta#matrix(0, nrow = dim(beta)[1], ncol = dim(beta)[2])
  if(dir =="directed"){
    Fm_n <- cbind(1,Fmat, Hmat)
  }else{  
    Fm_n <- cbind(1, Fmat)}
  
  weights <- c(0, 1/apply(Fm_n[, -1], 2 , sd))#c(0,rep(1, 6))#
  
  # Try rescaling the values just for linear regression
  #Fm_n <- cbind(1, apply(Fm_n[,-1],2,range01) )
  ## Calculating update to the beta using gradient descent
  #y <-  Fm_n %*% t(as.matrix(beta)) 
  #gn <- c()
  
  
  for(idx in 1:dim(beta)[1]){
    #lambda <- round(abs(t(Fm_n[,-1]) %*% Z[,idx])/N, 2)
    #print(paste("The max lambda value for cov :",idx, " is ", paste(lambda,collapse = ",")))
    #lambda_max <- max(lambda)/10
    # maxLambda <<- c(maxLambda, lambda_max)
    #-1 * t( t(Fm_n) %*% (Z[,idx] - matrix(beta[idx,1],ncol=1, nrow= N ) - (Fm_n[,-1] %*% (beta[idx,-1 ])))) /N
    ## LASSO regression 
    # for(graditer in 1){
    gradient <- -1 * t( t(Fm_n) %*% (Z[,idx] - (Fm_n %*% beta_new[idx, ]) ) ) /N
    #gn <- c( gn, grad_norm(gradient))
    
    if(penalty == "LASSO") {
      #EN_fit <- cv.glmnet(Fm_n[,-1], Z[,idx], alpha = 1, lambda = seq(0.001,1,length = 100))
      #lambda.min <-  EN_fit$lambda.min
      #beta_new[idx,] <- as.numeric(coef(glmnet(Fm_n[,-1], Z[,idx], alpha = 1,lambda = 0.1)))
      beta_new[idx,] <-  beta_new[idx,] - ( alphaLin * (gradient + (lambda * weights * safe_sign(beta_new[idx,]) ) ) )
      #beta_new[idx,] - ( alphaLin * (gradient + (lambda * weights * safe_sign(beta_new[idx,]) ) ) )
      #beta_new[idx,] - ( alphaLin * (gradient + (lambda * safe_sign(beta_new[idx,]) ) ) )
      #beta_new[idx,] - ( alphaLin * (gradient + (lambda * weights * safe_sign(beta_new[idx,]) ) ) )
      #beta - ( alphaLin * (gradient + (lambda * sign(beta) ) ) )
      if(is.nan(beta_new[1,1])){
        print("nan beta")}
    }else if(penalty == "Ridge"){
      beta_new <-  beta - (alphaLin * ( gradient + (2*lambda * beta ) ) )
    } else if(penalty == "Proximal"){
      beta_temp <-  beta - ( alphaLin * gradient ) 
      beta_new <- sign(beta_temp) * max( abs(beta_temp) - alphaLin * lambda , 0)
    }
    # }
  }
  #print(gn)
  #print(beta_new)
  #GN <- rbind(GN, gn)
  
  if(sum(is.infinite(beta_new)) > 0){
    print("Beta_new is inifinite")
  }
  return(beta_new)
}

# Soft-thresholding function
soft_threshold <- function(z, lambda) {
  sign(z) * pmax(abs(z) - lambda, 0)
}

LinRParamUpdt <- function(beta, Z, Fmat,Hmat, alpha, lambda, N, missVals,dir, 
                          impType, alphaLin, penalty ) {
  if(printFlg == TRUE){
    print("In Linparamupdate Func")
  }
  beta_new <- beta
  if(dir =="directed"){
    Fm_n <- cbind(1,Fmat, Hmat)
  }else{  
    Fm_n <- cbind(1, Fmat)
  }
  #Fm_n[, -1] <- scale(Fm_n[, -1], center = TRUE, scale = TRUE)
  
  #weights <- c(0, 1/apply(Fm_n[, -1], 2 , sd))
  sigmaSq <- SigmaSqCalc(Z, beta, Fmat, Hmat, missVals,dir)
  for (idx in 1:ncol(Z)) {  # Loop over responses
    # Update intercept (unpenalized)
    beta_new[idx, 1] <- mean(Z[, idx] - (Fm_n[, -1] %*% beta_new[idx, -1]))
    
    # Update coefficients (j > 1)
    for (j in 2:ncol(Fm_n)) {
      # Partial residual: r = y - intercept - X_{-j} β_{-j}
      #r <- Z[, idx] - beta_new[idx, 1] - (Fm_n[, -c(1,j)] %*% beta_new[idx, -c(1,j)])
      
      # Unpenalized LS estimate (z_j = X_j^T r)
      #z_j <- sum(Fm_n[, j] * r)
      
      # Soft-thresholding
      #beta_new[idx, j] <- sign(z_j) * max(abs(z_j) - N * lambda, 0) / sum(Fm_n[, j]^2)
      # print(cat("beta_j:", z_j, "→", beta_new[idx, j], "\n"))
      
      r_j <- Z[,idx] - (beta_new[idx,1] + Fm_n[,-c(1,j)] %*% beta_new[idx ,-c(1,j)])# + Fm_n[, j] * beta_new[idx, j]
      z_j <- sum(Fm_n[, j] * r_j)# / N
      beta_new[idx, j] <- sign(z_j) * max(abs(z_j) - (N * lambda), 0) / sum(Fm_n[, j]^2)
      #soft_threshold(z_j, lambda*sigmaSq[idx])
    }
  }
  #  print(beta_new)
  
  if(sum(is.infinite(beta_new)) > 0){
    print("Beta_new is inifinite")
  }
  return(beta_new)
}


updateLinearRegParam <- function(beta,missVals,Z,Wtm1,Wtm2, alpha,lambda,N,dir, 
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
  
  if (length(W) > 0) {
    cases <- complete.cases(X)
    Q <- CalcQuk(Fmat, Hmat, W, dir)
    S <- sum((log(Q[cases, ]) * X[cases,]) + (log(1-Q[cases,]) * (1- X[cases,])), na.rm = TRUE)
  }
  
  return(S)
}

Lz_cal <- function(beta,N,Z,Fmat, Hmat,sigmaSq,dir,lambda,missVals ){
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
    sigmaSq <- SigmaSqCalc(Z, beta, Fmat, Hmat, missVals,dir)
    
    for (j in 1:dim(beta)[1]) {
      #S <- S + sum(( Z[,j]  - (Fin %*% beta[j, ])) ^ 2)/ (2*N)
      
      ## Include the penalty term in the cost calculation
      #S <- S - sum(( Z[,j]  - (Fin %*% beta[j, ])) ^ 2)/ (2) #+ lambda * sum(abs(beta[j,]))
      
      S_tmp <- sum(( Z[,j]  - (Fin %*% beta[j, ])) ^ 2)
      S <- S - (N*log(2*pi*sigmaSq[j])/2) - (S_tmp / (2 * sigmaSq[j])) - lambda * sum(abs(beta[j,]))
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

findLLDir <- function(G,  Ftot, Htot, Win = NA, Wout = NA, X_in = NA,
                      X_out = NA,betain = NULL, betaout = NULL,
                      Z_in = NA, Z_out = NA,sigmaSqin = NA,sigmaSqout = NA, 
                      alphaLL, Eneg, dir, epsilon, lambda,missVals) {
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
                   Fmat = Ftot, Hmat = Htot,
                   sigmaSq = sigmaSqout,dir,lambda, missVals)
  S3 <- S3_in + S3_out
  
  ## Calculating the final log likelihood
  ll <- (S1 +  (S2 + S3))#(alphaLL* S1) +  ((1-alphaLL) * (S2 + S3))
  if(is.nan(ll)){
    print(ll)
    print("stop here in find ll")
  }
  #print(c(ll, S1, S2, S3))
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

initWtmat <- function(G, mode, N, nc, seed){
  #set.seed(seed)
  colwts <- igraph::degree(G, mode = mode) / sum(igraph::degree(G, mode = mode))
  Wtmat <- matrix(nrow = N, ncol = nc, runif(nc * N, 0,0.1))
  
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

initbeta <- function(o, ncoef,missVals, Z, Fmat, Hmat, alpha,N, lambda,dir,
                     impType, alphaLin, penalty){
  beta <- NULL
  sigmaSq <- NULL
  if (o > 0) {
    beta <- matrix(nrow = o, ncol = (ncoef) + 1, 0)
    
    for(i in 1:o){
      beta[i,1] <- mean(Z[,i], na.rm =TRUE)
    }
    #LinRParamUpdt(beta, Z, Fmat,Hmat, alpha, lambda, N, missVals,dir, impType, alphaLin, penalty)
    sigmaSq <- 1#SigmaSqCalc(Z, beta, Fmat,Hmat, missVals,dir)
    #if (sum(missVals) > 0) {
    #  Z <- PredictCovLin(Fmat, Hmat, Z, beta, missVals,dir, impType)
    #}
  }
  return(list(beta, sigmaSq))
}

randomizeIdx <- function(N, randomize){
  ## Updating the community weight parameter
  if (randomize == TRUE) {
    s <- sample(1:N, size = N, replace = FALSE)
  }
  return(s)
}

CoDA <- function(G,nc, k = c(0, 0) ,o = c(0, 0) , N,  alpha, lambda, thresh, 
                 nitermax, orig,randomize = TRUE, CovNamesLinin = c(),
                 CovNamesLinout = c(),CovNamesLPin = c(),CovNamesLPout = c(), 
                 dir, alphaLL = NULL,test = TRUE,missing = NULL, covOrig, 
                 epsilon,  impType, alphaLin, penalty, seed, covInit) {
  if(printFlg == TRUE){
    print("In CoDA Func")
  }
  
  if(dir == "directed"){
    ncoef = nc*2
  }else{ncoef = nc}
  
  ## Community weights outgoing connections
  Ftot <- initWtmat(G,"out",N,nc, seed)
  Finit <- Ftot
  
  ## Community weights outgoing connections
  if(dir == "directed"){
    Htot <- initWtmat(G,"in",N,nc, (seed+10))
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
      #rbinom(n = sum(missValsout[,i]), size = 1, prob = 0.5)
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
                Z_out[v,i] <- mean(Z_out[,i], na.rm =TRUE)
              }else{
                Z_out[v, i] <- mode(Z_out[c(neigh), i])
              }
            } else{
              ## Mean imputation of initial values if no neighbours
              Z_out[v,i] <- mean(Z_out[,i], na.rm =TRUE)
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
        }
      }
      backTransParams[[i]] <- summary(Z_out[,i])
      covOrig[,i] <- scale(covOrig[,i])#covOrig[,i] - mean(covOrig[,i], na.rm = TRUE)
      Z_out[,i] <- scale(Z_out[,i])#Z_out[,i] - mean(Z_out[,i], na.rm =TRUE)#range01(Z_out[,i])
      
      #print(summary(Z_out[,i]))
      #print(summary(covOrig[,i]))
    }
  }
  
  ## Degree condition for incoming and outgoing edges
  #in_E <-  igraph::degree(G, mode = "in")
  #out_E <- igraph::degree(G, mode = "out")
  
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
  #print(betaout)
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
  
  delta <- getDelta(N)
  continue <- TRUE
  s <- 1:N
  
  while (continue) {
    
    ## Update the log likelihood.
    LLvec <- findLLDir(G,Ftot,Htot, Win,Wout,X_in,X_out, betain,betaout,Z_in,Z_out,
                       sigmaSqin, sigmaSqout,alphaLL, Eneg, dir, epsilon,lambda, missValsout)
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
                   " MSE ", paste0(round(tail(mse,1),3), collapse = ", ")))
      break
    }
    
    LLold <- LLnew
    lllst <- rbind(lllst, LLvec)
    iter <- iter + 1
    
    
    ## Randomize the community weights update
    s <- randomizeIdx(N, randomize)
    
    #if(dir =="directed"){
    ## Updating F i.e. all outgoing connections
    ## We look for neighbors our node is connected to with the edges directed outward from our node.
    opMat <- updateWtmatBothTypes(G,Ftot,Htot,"all",s,nc,X_out,Z_out,k_out,o_out,
                                  betaout,Wout, alphaLL,missValsout,alpha, start = 2, 
                                  end = nc + 1, dir)
    Ftot <- opMat[[1]]
    Htot <- opMat[[2]]
    
    ## Updating H i.e. all incoming connections
    ## We look for neighbors our node is connected to with the edges directed outward from our node.
    # Htot <- updateWtmat(G,Htot,Ftot,"in",s,nc, X_out,Z_out,k_out,o_out,
    #                     betaout,Wout,alphaLL,missValsout,alpha, start = nc+2, 
    #                     end = nc*2+1, dir)
    #} else{
    #  opMat <- updateWtmatBothTypes(G,Ftot,Htot,"all",s,nc,X_out,Z_out,k_out,o_out,
    #                      betaout,Wout,alphaLL,missValsout,alpha, start = 2, 
    #                      end = nc+1, dir)
    #  Htot <- opMat[[1]]
    #  Ftot <- Htot
    #}
    
    ## Updating logistic paramters
    ## Incoming correlated covariates
    # b1 <- updateLogisticParam(Win,X_in,Ftot, missValsin,lambda,alpha)
    # Win <- b1[[1]]
    # X_in <- b1[[2]]
    ## outgoing correlated covariates
    b2 <- updateLogisticParam(W = Wout, X = X_out , Wtm = Ftot, Wtm2 = Htot, 
                              missVals = missValsout, lambda = lambda, 
                              alphaLR = alphaLin, dir = dir, impType)
    Wout <- b2[[1]]
    X_out <- b2[[2]]
    
    ### Updating the beta matrix for continuous covariates
    #print("Gradient for out")
    c1 <- updateLinearRegParam(betaout,missValsout,Z_out,Ftot,Htot,alpha,
                               lambda,N,dir, impType, alphaLin, penalty)
    betaout <- c1[[1]]
    Z_out <- c1[[2]]
    sigmaSqout <- c1[[3]]
    
    #print("Gradient for in")
    # c2 <- updateLinearRegParam(betain, missValsin, Z_in_orig, Ftot,alpha,lambda,N)
    # betain <- c2[[1]]
    Z_in <- 0#c2[[2]]
    # sigmaSqin <- c2[[3]]
    
    ## Look for Communities based on the F matrix
    #OmegaVal <- c(OmegaVal, OmegaIdx(G, Ftot, Htot, N, delta,nc))
    
    if(test == TRUE){
      
      arilst <- c(arilst, ARIop(Ftot,Htot,orig,nc,N))
      OIn <- OmegaIdx(G, Ftot, Htot, N, delta,nc)
      OmegaVal <- c(OmegaVal, OIn)
      #print(OIn)
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
      backTransParams
    )
  )
}
