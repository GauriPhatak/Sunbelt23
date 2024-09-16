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

## nc : numer of communities
## k : number of covariates
## pClust: Probability of a ndoe belonging to a community
## N : Number of nodes
## pC : matrix of probability for binary assignment for each cluster
NWSimBin <- function(nc,
                     k,
                     pC,
                     N,
                     pClust,
                     B,
                     o,
                     dist,
                     dir,
                     covTypes = c("binary", "continuous"),
                     CovNamesLin = c() ,
                     CovNamesLP = c(),
                     missing = NULL) {
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

NWSimBin2 <- function(nc, N, probedge, pClust) {
  #### create example input data
  nodesSet1 <- letters[1:nc]
  nodesSet2 <- 1:N
  
  ##create edges based on the probability of edge formation
  C <- rep(0, N)
  #for(i in 1:N){
  #  C[i] <- sample(nodesSet1, size = 1, prob = pClust)
  #}
  C <- sample(
    x = letters[1:nc],
    size = N,
    replace = TRUE,
    prob = pClust
  )
  #table(C)
  combo <- data.frame(cbind(expand.grid(grp1 = letters[1:nc], grp2 = letters[1:nc]), B))
  edges <- matrix(nrow = 0, ncol = 2)
  for (c in 1:dim(combo)[1]) {
    print(combo[c, ])
    posEdges <- expand.grid(which(C == combo[c, 1]), which(C == combo[c, 2]))
    print(dim(posEdges))
    edgeStatus <- rep(0, dim(posEdges)[1])
    for (i in 1:dim(posEdges)[1]) {
      edgeStatus[i] <- flip(combo[c, 3])
    }
    edges <- rbind(edges, posEdges[as.logical(edgeStatus), ])
  }
  
  #edgeLst <- as.vector(edges)
  # first we give prefixes to the nodes to discern the two partition
  g <- make_empty_graph()  #graph.empty()
  #g <- add_vertices(g,nv=length(nodesSet1),attr=list(name=letters[1:nc],
  #                                                   type=rep(TRUE,length(nodesSet1))))
  g <- add_vertices(g,
                    nv = length(nodesSet2),
                    attr = list(name = 1:N, type = rep(FALSE, length(nodesSet2))))
  #
  # # we need to turn edgeList into a vector (and using names instead of indexes)
  edgeListVec <- as.vector(t(as.matrix(data.frame(
    S1 = edges[, 1], S2 = edges[, 2]
  ))))
  g <- add_edges(g, edgeListVec)
  
  # check if is recognized as bipartite
  is_bipartite(g)
  
  return(g)
}

genBipartite <- function(N,
                         nc,
                         pClust,
                         k_in,
                         k_out,
                         pC,
                         covTypes,
                         CovNamesLPin,
                         CovNamesLPout,
                         pConn,
                         dir,
                         dirPct,
                         epsilon,
                         missing) {
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
      p_uv <- 1 - prod(1 - ((as.numeric(
        Cnew[i, ] - Cnew[j, ] > 1
      )) * pConn)) + epsilon
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
      g.sim <- SetBinarycov(
        g.sim,
        k_in,
        N,
        nc,
        k_start = 0,
        Cnew,
        connType = 1,
        pC,
        CovNamesLPin
      )
    }
    if (length(CovNamesLPout) > 0) {
      g.sim <- SetBinarycov(
        g.sim,
        k_out,
        N,
        nc,
        k_start = k_in,
        Cnew,
        connType = -1,
        pC,
        CovNamesLPout
      )
    }
  }
  
  g.sim <- g.sim  %>%
    set_vertex_attr(name = "Cluster", value = C)
  return(list(g.sim, Apuv))
}

SetBinarycov <- function(G,
                         k,
                         N,
                         nc,
                         k_start,
                         Cnew,
                         connType,
                         pC,
                         CovNamesLP) {
  ## Based on the probability matrix pC assign indicates the binary assignment of covariates to a community
  ## setting covariates for outgoing communities. hence Cnew will be -1
  binVal <- matrix(ncol = k, nrow = N, 0)
  for (i in 1:nc) {
    pc_i <- k_start+i
    for (j in 1:k) {
      pc_j <- k_start+j
      binVal[which(Cnew[, i] == connType), j] <- rbinom(length(which(Cnew[, i] == connType)), 1, pC[pc_i, pc_j])
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

#f_u : Community weight vector for node u
#Fvmat : Community weight matrix for node neighbours
#Fvnot : Community weight matrix for not neighbours
#X_u: Covariate matrix
#Q_u: logistic probability matrix
#W: Logistic weight matrix
CmntyWtUpdt <- function(f_u,
                        Fvmat,
                        Fvnot,
                        X_u = NULL ,
                        W = NULL ,
                        Z_u = NULL,
                        beta = NULL,
                        sigmaSq,
                        alpha,
                        alphaLL) {
  nc <- ncVal
  ## Gradient function to update log likelihood of G
  ## First part of log likelihood of G based on the structure of the network
  a <- exp(-1 * (f_u %*% t(Fvmat)))
  b <- a / (1 - a)
  llG_1 <- t(Fvmat) %*% t(b)
  ## 2nd part of log likelihood of G
  llG_2 <- as.matrix(colSums(Fvnot))
  ## Total llG: This math has been verified
  llG <- llG_1 - llG_2
  ## experimental scaled version of the above llog lik
  scale <- 1
  llG <- scale * llG
  llX <- rep(0, nc + 1)
  if (length(W) > 0) {
    ## Updating the logistic regression weights
    ## Gradient function to update log likelihood for covariates
    ## adding intercept term
    Q_u <- t(1 / (1 + exp(-1 * (
      W %*% as.matrix(c(1, f_u))
    ))))
    ## This math has been verified
    llX <- t((X_u - Q_u) %*% W)
  }
  llZ <- rep(0, nc + 1)
  if (length(beta) > 0) {
    ## Adding the gradient based on continuous covariates
    llZ <-  ((Z_u - t(beta %*% as.matrix(c(
      1, f_u
    )))) / sigmaSq) %*% beta
  }
  ## Function to update the community weights vector using non negative matrix factorization
  if (is.null(alphaLL)) {
    f_u_new <- f_u + t((alpha * (llG + (llX[2:(nc + 1)] + llZ[2:(nc + 1)]))))
  } else{
    f_u_new <- f_u + t((alpha * (
      alphaLL * llG + (1 - alphaLL) * (llX[2:(nc + 1)] + llZ[2:(nc + 1)])
    )))
  }
  f_u_new[f_u_new < 0] <- 0.000001
  return(f_u_new)
}

sigmoid <- function(W, Ftot) {
  Q <- 1 / (1 + exp(-1 * (W %*% t(Ftot))))
  return(Q)
}

## W: Logistic weight matrix
## X: Covariate values
## Ftot: Community weight params
## alpha: tuning param
## lambda:
LRParamUpdt <- function(W, X, Ftot, alpha, lambda, missing) {
  Ftotn <- cbind(1, Ftot)
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

LinRParamUpdt <- function(beta, Z, Ftot, alpha, lambda, missVals) {
  beta_new <- matrix(0, nrow = dim(beta)[1], ncol = dim(beta)[2])
  Ftot_n <- cbind(1, Ftot)
  for (i in 1:dim(beta)[1]) {
    beta_new[i, ] <- MASS::ginv(t(Ftot_n) %*% Ftot_n) %*% t(Ftot_n) %*% Z[, i]
  }
  return(beta_new)
}

findLL <- function(G,
                   Ftot,
                   W = NA,
                   X = NA,
                   beta = NA ,
                   Z = NA,
                   sigmaSq = NA,
                   alphaLL) {
  scale <- 1
  E <- as_edgelist(G)
  F1 <- Ftot[E[, 1], ]
  F2 <- Ftot[E[, 2], ]
  N <-  gorder(G)
  ## Calculating the L_g part of the log likelihood for structure of the network
  
  ## Variable to save the loglik contributed by nodes with edges between them
  S11 <- 0
  for (i in 1:dim(E)[1]) {
    S11 <- S11 + log(1 - exp(-1 * sum(F1[i, ] * F2[i, ])))
  }
  
  ## Creating a graph that is negative of the given graph to get edges not present in G
  A <- 1 - (1 * as.matrix(as_adjacency_matrix(G)))
  E2 <-  as_edgelist(graph_from_adjacency_matrix(A))
  F1 <- Ftot[E2[, 1], ]
  F2 <- Ftot[E2[, 2], ]
  
  ## Variable to save the loglik contributed by nodes without edges between them
  S12 <- 0
  for (i in 1:dim(E2)[1]) {
    S12 <- S12 + sum(F1[i, ] * F2[i, ])
  }
  S1 <- S11 - S12
  ## Calculating the L_x part of the log likelihood for binary covariates
  S2 <- 0
  if (length(W) > 0) {
    for (i in 1:N) {
      Q_uk <- 1 / (1 + exp(-1 * rowSums(c(1, Ftot[i, ]) * W)))#1/(1+exp(-1*c(1, Ftot[i,]) %*% t(W)))
      S2 <- S2 + sum((log(Q_uk) * X[i, ]) + (log(1 - Q_uk) * (1 - X[i, ])))
    }
  }
  
  
  ## Calculating the L_z part of the log likelihood for continuous covariates
  S3 <- 0
  if (length(beta) > 0) {
    for (j in 1:o) {
      S3_tmp <- 0
      for (i in 1:N) {
        S3_tmp <- sum(c(S3_tmp, (Z[i, j] - sum(
          c(1, Ftot[i, ]) * beta[j, ], na.rm = TRUE
        )) ^ 2), na.rm = TRUE)
      }
      S3 <- S3 - (N * log(sigmaSq[j]) / 2 + S3_tmp / (2 * sigmaSq[j]))
    }
  }
  
  #Loglik[iterationNum, ] <<- c(S1,S2,S3)
  #iterationNum  <<- iterationNum +1
  ## Calculating the final log likelihood
  if (is.null(alphaLL)) {
    ll <- S1 + (S2 + S3)
  } else{
    ll <- alphaLL * S1 + (1 - alphaLL) * (S2 + S3)
  }
  return(ll)
}

#G,Ftot,Htot,Win, Wout,X_in, X_out,betain, betaout,Z_in, Z_out, sigmaSq,alphaLL
findLLDir <- function(G,
                      Ftot,
                      Htot,
                      Win = NA,
                      Wout = NA,
                      X_in = NA,
                      X_out = NA,
                      betain = NA,
                      betaout = NA,
                      Z_in = NA,
                      Z_out = NA,
                      sigmaSq = NA,
                      alphaLL) {
  scale <- 1
  E <- as_edgelist(G)
  Fv <- Ftot[E[, 1], ]
  Hv <- Htot[E[, 2], ]
  N <-  gorder(G)
  ## Calculating the L_g part of the log likelihood for structure of the network
  
  ## Variable to save the loglik contributed by nodes with edges between them
  S11 <- 0
  for (i in 1:dim(E)[1]) {
    S11 <- S11 + log(1 - exp(-1 * sum(Fv[i, ] * Hv[i, ])))
  }
  
  ## Creating a graph that is negative of the given graph to get edges not present in G
  A <- 1 - (1 * as.matrix(as_adjacency_matrix(G)))
  E2 <-  as_edgelist(graph_from_adjacency_matrix(A))
  Fv2 <- Ftot[E2[, 1], ]
  Hv2 <- Htot[E2[, 2], ]
  
  ## Variable to save the loglik contributed by nodes without edges between them
  S12 <- 0
  for (i in 1:dim(E2)[1]) {
    S12 <- S12 + sum(Fv2[i, ] * Hv2[i, ])
  }
  S1 <- S11 - S12
  
  ## Calculating the L_xin part of the log likelihood for binary covariates
  S2_in <- 0
  if (length(Win) > 0) {
    for (i in 1:N) {
      Qin_uk <- 1 / (1 + exp(-1 * rowSums(c(1, Ftot[i, ]) * Win)))#1/(1+exp(-1*c(1, Ftot[i,]) %*% t(W)))
      S2_in <- S2_in + sum((log(Qin_uk) * X_in[i, ]) + (log(1 - Qin_uk) * (1 - X_in[i, ])))
    }
  }
  
  ## Calculating the L_xout part of the log likelihood for binary covariates
  S2_out <- 0
  if (length(Wout) > 0) {
    for (i in 1:N) {
      Qout_uk <- 1 / (1 + exp(-1 * rowSums(c(1, Htot[i, ]) * Wout)))#1/(1+exp(-1*c(1, Ftot[i,]) %*% t(W)))
      S2_out <- S2_out + sum((log(Qout_uk) * X_out[i, ]) + (log(1 - Qout_uk) * (1 - X_out[i, ])))
    }
  }
  
  ## adding the loglikelihood from the covariates
  S2 <- S2_in + S2_out
  
  ## Calculating the L_zin part of the log likelihood for continuous covariates
  S3_in <- 0
  if (length(betain) > 0) {
    for (j in 1:length(betain)) {
      S3in_tmp <- 0
      for (i in 1:N) {
        S3in_tmp <- sum(c(S3in_tmp, (
          Z_in[i, j] - sum(c(1, Ftot[i, ]) * betain[j, ], na.rm = TRUE)
        ) ^ 2), na.rm = TRUE)
      }
      S3_in <- S3_in - (N * log(sigmaSq[j]) / 2 + S3in_tmp / (2 * sigmaSq[j]))
    }
  }
  
  ## Calculating the L_zout part of the log likelihood for continuous covariates
  S3_out <- 0
  if (length(betaout) > 0) {
    for (j in 1:length(betaout)) {
      S3out_tmp <- 0
      for (i in 1:N) {
        S3out_tmp <- sum(c(S3out_tmp, (
          Z_out[i, j] - sum(c(1, Htot[i, ]) * betaout[j, ], na.rm = TRUE)
        ) ^ 2), na.rm = TRUE)
      }
      S3_out <- S3_out - (N * log(sigmaSq[j]) / 2 + S3out_tmp / (2 * sigmaSq[j]))
    }
  }
  
  S3 <- S3_in + S3_out
  
  ## Calculating the final log likelihood
  if (is.null(alphaLL)) {
    ll <- S1 + (S2 + S3)
  } else{
    ll <- alphaLL * S1 + (1 - alphaLL) * (S2 + S3)
  }
  return(ll)
}

SigmaSqCalc <- function(Z, beta, Ftot, missVals) {
  sigmaSq <- rep(0, dim(beta)[1])
  for (i in 1:dim(beta)[1]) {
    sigmaSq[i] <- sum((Z[!missVals[, i], i] - (t(beta[i, ]) %*% t(cbind(
      1, Ftot
    )[!missVals[, i], ]))) ^ 2, na.rm = TRUE) / sum(!missVals[, i])
  }
  return(sigmaSq)
}

PredictCovLin <- function(Ftot, Z, beta, missVals) {
  for (i in 1:dim(Z)[2]) {
    Z[which(missVals[, i]), i] <- cbind(1, Ftot[which(missVals[, i]), ]) %*% as.matrix(beta[i, ])
  }
  return(Z)
}

mode <- function(x) {
  which.max(tabulate(x))
}

CESNA <- function(G,
                  nc,
                  k = 0 ,
                  o = 0 ,
                  N,
                  alpha,
                  lambda,
                  thresh,
                  nitermax,
                  randomize = TRUE,
                  orig,
                  CovNamesLin = c() ,
                  CovNamesLP = c(),
                  alphaLL = NULL) {
  ## Community weights
  Finit <- igraph::degree(G) / sum(igraph::degree(G)) #cbind(Finit, Finit, Finit)
  Ftot <- matrix(nrow = N, ncol = nc, c(rep(Finit, times =  nc) + runif(nc *
                                                                          N)))##sample(seq(0.1,0.5,length.out =10), N*nc, replace = TRUE))
  Ftot_orig <- Ftot
  Ftotold <- Ftot
  iter <- 0
  ## Get the covariate matrix for binary data and continuous data
  covtmp <-  as.data.frame(vertex_attr(G)) %>%
    dplyr::select(all_of(c(CovNamesLin, CovNamesLP)))
  X <- 0
  Z <- 0
  if (k > 0) {
    X <- as.matrix(covtmp %>%
                     dplyr::select(all_of(CovNamesLP)))
    missVals <- is.na(X)
    numMissVal <- sum(missVals)
  }
  if (o > 0) {
    Z <- as.matrix(covtmp %>%
                     dplyr::select(all_of(CovNamesLin)))
    ## Use the beta to predict the missing values
    ## Find the missing values
    missVals <- is.na(Z)
    numMissVal <- sum(missVals)
    ## For each covariate find the missing values and predict the value
  }
  W <- NULL
  W_orig <- NULL
  if (k > 0) {
    W <- matrix(nrow = k, ncol = nc + 1, 0)
    if (numMissVal > 0) {
      X <- PredictCovLR(X, W, Ftot, missVals)
    }
    W <- LRParamUpdt(W, X, Ftot, alpha, lambda, missVals)
    W_orig <- W
  }
  ## Setting initial values of linear regression weights
  beta <- NULL
  if (o > 0) {
    beta <- matrix(nrow = o, ncol = nc + 1, 0)
    if (numMissVal > 0) {
      Z <- PredictCovLin(Ftot, Z, beta, missVals)
    }
    beta <- LinRParamUpdt(beta, Z, Ftot, alpha, lambda, missVals)#matrix(nrow = o, ncol = nc+1, 0)
    beta_orig <- beta
    sigmaSq <- SigmaSqCalc(Z, beta, Ftot, missVals)
  }
  LLold <- findLL(G, Ftot, W, X, beta, Z, sigmaSq, alphaLL)
  lllst <- c()
  lllst <- c(lllst, LLold)
  tempop1 <- 0
  LLG <- 0
  LLX <- 0
  nodeid <- 0
  arilst <- c()
  continue <- TRUE
  mse <- matrix(nrow = 0, ncol = o)
  accuracy <- matrix(nrow = 0, ncol = k)#c()
  while (continue) {
    Ftotold <- Ftot
    ## Updating the commuynity weight parameter
    s <- 1:N
    if (randomize == TRUE) {
      s <- sample(1:N, size = N, replace = FALSE)
    }
    for (i in s) {
      f_u <- matrix(Ftot[i, ], ncol = nc)
      ##The neighbours of this node
      neigh <- neighbors(G, i, mode = "all")
      Fvmat <- matrix(Ftot[neigh, ], ncol = nc)
      Fvnot <- matrix(Ftot[-neigh, ], ncol = nc)
      X_u <- NULL
      if (k > 0) {
        X_u <- matrix(X[i, ], ncol = k)
      }
      Z_u <- NULL
      if (o > 0) {
        Z_u <- matrix(Z[i, ], ncol = o)
        sigmaSq <- SigmaSqCalc(Z, beta, Ftot, missVals)
      }
      Ftot[i, ] <- CmntyWtUpdt(f_u, Fvmat, Fvnot, X_u, W, Z_u, beta, sigmaSq, alpha, alphaLL)
      tempop1 <- append(tempop1, list(Ftot[i, ]))
      nodeid <- append(nodeid, i)
    }
    ## Updating logistic paramters
    if (length(W) > 0) {
      if (numMissVal > 0) {
        W <-  LRParamUpdt(W, X, Ftot, alpha, lambda, missVals)
        X <-  PredictCovLR(X, W, Ftot, missVals)
      } else{
        W <-  LRParamUpdt(W, X, Ftot, alpha, lambda, missVals)
      }
    }
    if (length(beta) > 0) {
      if (numMissVal > 0) {
        beta <- LinRParamUpdt(beta, Z, Ftot, alpha, lambda, missVals)
        sigmaSq <-  SigmaSqCalc(Z, beta, Ftot, missVals)#SigmaSqCalc(as.data.frame(Z[missIdx,]),beta, as.data.frame(Ftot[missIdx,]))
        Z <- PredictCovLin(Ftot, Z, beta, missVals)
      } else{
        beta <- LinRParamUpdt(beta, Z, Ftot, alpha, lambda, missVals)
        sigmaSq <-  SigmaSqCalc(Z, beta, Ftot, missVals)
      }
    }
    ## Look for Communities based on the F matrix
    mem <- rep(NA, N)
    for (i in 1:N) {
      mem[i] <- which.max(Ftot[i, ])
    }
    arilst <- c(arilst, mclust::adjustedRandIndex(orig, mem))
    ## Update the log likelihood.
    LLnew <- findLL(G, Ftot, W, X, beta, Z, sigmaSq , alphaLL)
    #if(any(c(((LLnew - LLold) < thresh), (iter > nitermax)))){
    pctInc <- ((LLnew - LLold) / abs(LLold)) * 100
    if (any(pctInc < thresh, (iter > nitermax))) {
      continue <- FALSE
    }
    LLold <- LLnew
    lllst <- c(lllst, LLnew)
    iter <- iter + 1
    if (o > 0) {
      pred <- cbind(1, Ftot) %*% t(beta)
      mse <- rbind(mse, colSums((pred - Z) ^ 2) / dim(Z)[1])
    }
    if (k > 0) {
      predLR <- apply(sigmoid(W, cbind(1, Ftot)) > 0.5, 1, as.numeric)
      acc <- rep(0, k)
      for (i in 1:k) {
        cm <-  confusionMatrix(factor(predLR[, i], levels = c(0, 1)), factor(X[, i]))
        acc[i] <- round(cm$overall[1], 2)
      }
      accuracy <- rbind(accuracy, acc)
    }
  }
  ## Updating the logistic weight parameters
  return(
    list(
      Ffin = Ftot,
      Forig = Ftot_orig,
      Worig = W_orig,
      ARI = arilst,
      Wfin = W,
      Loglik = lllst,
      beta = beta,
      Zcov = Z,
      MSE = mse,
      Acc = accuracy
    )
  )
}

initCov <- function(covtmp, CovNames) {
  X_in <- as.matrix(covtmp %>%
                      dplyr::select(all_of(CovNames)))
  missVals <- is.na(X_in)
  numMissVal <- sum(missVals)
  
  return(list(X_in, numMissVal))
}

CoDA <- function(G,
                 nc,
                 k = c(0, 0) ,
                 o = c(0, 0) ,
                 N,
                 alpha,
                 lambda,
                 thresh,
                 nitermax,
                 randomize = TRUE,
                 orig,
                 CovNamesLinin = c(),
                 CovNamesLinout = c(),
                 CovNamesLPin = c(),
                 CovNamesLPout = c(),
                 alphaLL = NULL) {
  ## Community weights outgoing connections
  Finit <- igraph::degree(G, mode = "all") / sum(igraph::degree(G, mode = "all")) #cbind(Finit, Finit, Finit)
  Ftot <- matrix(nrow = N, ncol = nc, runif(nc * N))#matrix(nrow = N, ncol = nc, c(rep(Finit, times =  nc) + runif(nc*N)))##sample(seq(0.1,0.5,length.out =10), N*nc, replace = TRUE))
  Finit <- Ftot
  #Ftotold <- Ftot
  
  ## Community weights outgoing connections
  Hinit <- igraph::degree(G, mode = "all") / sum(igraph::degree(G, mode = "all")) #cbind(Finit, Finit, Finit)
  Htot <- matrix(nrow = N, ncol = nc, runif(nc * N))#matrix(nrow = N, ncol = nc, c(rep(Hinit, times =  nc) + runif(nc*N)))##sample(seq(0.1,0.5,length.out =10), N*nc, replace = TRUE))
  Hinit <- Htot
  #Htotold <- Htot
  
  iter <- 0
  
  ########## Get the covariate matrix for binary data and continuous data ##########
  covtmp <-  as.data.frame(vertex_attr(G)) %>%
    dplyr::select(all_of(
      c(CovNamesLinin, CovNamesLinout, CovNamesLPin, CovNamesLPout)
    ))
  X_in <- 0
  X_out <- 0
  
  Z_in <- 0
  Z_out <- 0
  
  ## Setting the number of covariates correlated to incoming or outgoing directions
  k_in <- k[1]
  k_out <- k[2]
  o_in <- o[1]
  o_out <- o[2]
  
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
    ## Use the beta to predict the missing values
    ## Find the missing values
    missValsin <- is.na(Z_in)
    numMissValin <- sum(missValsin)
    ## For each covariate find the missing values and predict the value
  }
  if (o_out > 0) {
    Z_out <- as.matrix(covtmp %>%
                         dplyr::select(all_of(CovNamesLinout)))
    ## Use the beta to predict the missing values
    ## Find the missing values
    missValsout <- is.na(Z_out)
    numMissValout <- sum(missValsout)
    ## For each covariate find the missing values and predict the value
  }
  
  ############ Initialize the coef matrix for covariates ##############
  ############## Setting initial values of binary covariates coefs weights ###############
  Win <- NULL
  Win_orig <- NULL
  if (k_in > 0) {
    Win <- matrix(nrow = k_in, ncol = (nc) + 1, 0)
    if (numMissValin > 0) {
      X_in <- PredictCovLR(X_in, Win, Ftot, missValsin)
    }
    Win <- LRParamUpdt(Win, X_in, Ftot, alpha, lambda, missValsin)
    Win_orig <- Win
  }
  
  Wout <- NULL
  Wout_orig <- NULL
  if (k_out > 0) {
    Wout <- matrix(nrow = k_out, ncol = (nc) + 1, 0)
    if (numMissValout > 0) {
      X_out <- PredictCovLR(X_out, Wout, Htot, missValsout)
    }
    Wout <- LRParamUpdt(Wout, X_out, Htot, alpha, lambda, missValsout)
    Wout_orig <- Wout
  }
  
  
  #################### Setting initial values of linear regression weights
  betain <- NULL
  if (o_in > 0) {
    betain <- matrix(nrow = o_in, ncol = (nc) + 1, 0)
    if (numMissValin > 0) {
      Z_in <- PredictCovLin(Ftot, Z_in, betain, missValsin)
    }
    betain <- LinRParamUpdt(betain, Z_in, Ftot, alpha, lambda, missValsin)#matrix(nrow = o, ncol = nc+1, 0)
    betain_orig <- betain
    sigmaSq <- SigmaSqCalc(Z_in, betain, Ftot, missValsin)
  }
  betaout <- NULL
  if (o_out > 0) {
    betaout <- matrix(nrow = o_out, ncol = (nc) + 1, 0)
    if (numMissValout > 0) {
      Z_out <- PredictCovLin(Htot, Z_out, betaout, missValsout)
    }
    betaout <- LinRParamUpdt(betaout, Z_out, Htot, alpha, lambda, missValsout)#matrix(nrow = o, ncol = nc+1, 0)
    betaout_orig <- betaout
    sigmaSq <- SigmaSqCalc(Z_out, betaout, Htot, missValsout)
  }
  
  
  ##Initial log likelihood value
  LLold <- findLLDir(G,
                     Ftot,
                     Htot,
                     Win,
                     Wout,
                     X_in,
                     X_out,
                     betain,
                     betaout,
                     Z_in,
                     Z_out,
                     sigmaSq,
                     alphaLL)
  lllst <- c()
  lllst <- c(lllst, LLold)
  tempop1 <- 0
  LLG <- 0
  LLX <- 0
  nodeid <- 0
  arilst <- c()
  continue <- TRUE
  mse <- matrix(nrow = 0, ncol = o)
  accuracy <- matrix(nrow = 0, ncol = k)#c()
  while (continue) {
    ## Updating the commuynity weight parameter
    s <- 1:N
    if (randomize == TRUE) {
      s <- sample(1:N, size = N, replace = FALSE)
    }
    
    ## Updating F i.e. all outgoing connections
    ## We look for neighbors our node is connected to with the edges directed outward from our node.
    for (i in s) {
      f_u <- matrix(Ftot[i, ], ncol = nc)
      ##The neighbours of this node
      neigh <- neighbors(G, i, mode = "out")
      Hvmat <- matrix(Htot[neigh, ], ncol = nc)
      neighAll <- neighbors(G, i, mode = "all")
      Hvnot <- matrix(Htot[-neighAll, ], ncol = nc)
      
      Xin_u <- NULL
      #Xout_u <- NULL
      if (k_in > 0) {
        Xin_u <- matrix(X_in[i, ], ncol = k_in)
      }
      # if(k_out > 0){
      #   Xout_u <- matrix(X_out[i,], ncol = k_out)
      # }
      
      Zin_u <- NULL
      if (o_in > 0) {
        Zin_u <- matrix(Z_in[i, ], ncol = o_in)
        sigmaSq <- SigmaSqCalc(Z_in, betain, Ftot, missValsin)
      }
      # Zout_u <- NULL
      # if(o_out > 0){
      #   Zout_u <- matrix(Z_out[i,], ncol = o_out)
      #   sigmaSq <- SigmaSqCalc(Z_out,betaout,Htot, missVals)
      # }
      Ftot[i, ] <- CmntyWtUpdt(f_u,
                               Hvmat,
                               Hvnot,
                               Xin_u,
                               Win,
                               Zin_u,
                               betain,
                               sigmaSq,
                               alpha,
                               alphaLL)
      
    }
    
    ## Updating H i.e. all incoming connections
    ## We look for neighbors our node is connected to with the edges directed outward from our node.
    for (i in s) {
      h_u <- matrix(Htot[i, ], ncol = nc)
      ##The neighbours of this node
      neigh <- neighbors(G, i, mode = "in")
      Fvmat <- matrix(Ftot[neigh, ], ncol = nc)
      neighAll <- neighbors(G, i, mode = "all")
      Fvnot <- matrix(Ftot[-neighAll, ], ncol = nc)
      #Xin_u <- NULL
      Xout_u <- NULL
      # if(k_in > 0){
      #   Xin_u <- matrix(X_in[i,], ncol = k_in)
      # }
      if (k_out > 0) {
        Xout_u <- matrix(X_out[i, ], ncol = k_out)
      }
      
      #Zin_u <- NULL
      Zout_u <- NULL
      # if(o_in > 0){
      #   Zin_u <- matrix(Z_in[i,], ncol = o_in)
      #   sigmaSq <- SigmaSqCalc(Z_in,betain,Htot, missVals)
      # }
      if (o_out > 0) {
        Zout_u <- matrix(Z_out[i, ], ncol = o_out)
        sigmaSq <- SigmaSqCalc(Z_out, betaout, Htot, missValsout)
      }
      
      Htot[i, ] <- CmntyWtUpdt(h_u,
                               Fvmat,
                               Fvnot,
                               Xout_u,
                               Wout,
                               Zout_u,
                               betaout,
                               sigmaSq,
                               alpha,
                               alphaLL)
      
    }
    
    ## Updating logistic paramters
    ## Incoming correlated covariates
    if (length(Win) > 0) {
      if (numMissValin > 0) {
        Win <-  LRParamUpdt(Win, X_in, Ftot, alpha, lambda, missValsin)
        X_in <-  PredictCovLR(X_in, Win, Ftot, missValsin)
      } else{
        Win <-  LRParamUpdt(Win, X_in, Ftot, alpha, lambda, missValsin)
      }
    }
    ## outgoing correlated covariates
    if (length(Wout) > 0) {
      if (numMissValout > 0) {
        Wout <-  LRParamUpdt(Wout, X_out, Htot, alpha, lambda, missValsout)
        X_out <-  PredictCovLR(X_out, Wout, Htot, missValsout)
      } else{
        Wout <-  LRParamUpdt(Wout, X_out, Htot, alpha, lambda, missValsout)
      }
    }
    
    ### Updating the beta matrix for continuous covariates
    if (length(betain) > 0) {
      if (numMissValin > 0) {
        betain <- LinRParamUpdt(betain, Z_in, Ftot, alpha, lambda, missValsin)
        sigmaSq <-  SigmaSqCalc(Z_in, betain, Ftot, missValsin)#SigmaSqCalc(as.data.frame(Z[missIdx,]),beta, as.data.frame(Ftot[missIdx,]))
        Z_in <- PredictCovLin(Ftot, Z_in, betain, missValsin)
      } else{
        betain <- LinRParamUpdt(betain, Z_in, Ftot, alpha, lambda, missValsin)
        sigmaSq <-  SigmaSqCalc(Z_in, betain, Ftot, missValsin)
      }
    }
    
    if (length(betaout) > 0) {
      if (numMissValout > 0) {
        betaout <- LinRParamUpdt(betaout, Z_out, Htot, alpha, lambda, missValsout)
        sigmaSq <-  SigmaSqCalc(Z_out, betaout, Htot, missValsout)#SigmaSqCalc(as.data.frame(Z[missIdx,]),beta, as.data.frame(Ftot[missIdx,]))
        Z_out <- PredictCovLin(Htot, Z_out, betaout, missValsout)
      } else{
        betaout <- LinRParamUpdt(betaout, Z_out, Htot, alpha, lambda, missValsout)
        sigmaSq <-  SigmaSqCalc(Z_out, betaout, Htot, missValsout)
      }
    }
    
    ## Look for Communities based on the F matrix
    mem <- rep(NA, N)
    for (i in 1:N) {
      m <- data.frame(com = rep(letters[1:nc], times = 2),
                      val = c(Ftot[i, ], Htot[i, ]))
      mem[i] <-  m$com[which.max(m$val)]#letters[which.max(Ftot[i,])]
    }
    arilst <- c(arilst, mclust::adjustedRandIndex(orig, mem))
    
    
    ## Update the log likelihood.
    LLnew <- findLLDir(G,
                       Ftot,
                       Htot,
                       Win,
                       Wout,
                       X_in,
                       X_out,
                       betain,
                       betaout,
                       Z_in,
                       Z_out,
                       sigmaSq,
                       alphaLL)
    #if(any(c(((LLnew - LLold) < thresh), (iter > nitermax)))){
    pctInc <- ((LLnew - LLold) / abs(LLold)) * 100
    if (any((pctInc < thresh), (iter > nitermax))) {
      print(paste0("The final percent increase ", pctInc))
      continue <- FALSE
    }
    
    LLold <- LLnew
    lllst <- c(lllst, LLnew)
    iter <- iter + 1
    
    ### Calculating the mean squared error  ************* NEED TO ADD FOR THE OUT COVARIATES ****************
    if (o_in > 0) {
      pred <- cbind(1, Ftot) %*% t(beta)
      mse <- rbind(mse, colSums((pred - Z_in) ^ 2) / dim(Z_in)[1])
    }
    
    
    #### Calculating the accuracy for the binary covariates
    if (k_in > 0) {
      predLR <- apply(sigmoid(Win, cbind(1, Ftot)) > 0.5, 1, as.numeric)
      acc <- rep(0, k_in)
      for (i in 1:k_in) {
        cm <-  confusionMatrix(factor(predLR[, i], levels = c(0, 1)), factor(X_in[, i]))
        acc[i] <- round(cm$overall[1], 2)
      }
      accuracy <- rbind(accuracy, acc)
    }
  }
  ## Updating the logistic weight parameters
  
  return(
    list(
      Ffin = Ftot,
      Forig = Finit,
      Hfin = Htot,
      Horig =  Hinit,
      Win_orig = Win_orig,
      Win_fin = Win,
      Wout_orig = Wout_orig,
      Wout_fin = Wout,
      Loglik = lllst,
      ARI = arilst,
      betain = betain,
      Zin_cov = Z_in,
      betaout = betaout,
      Zout_cov = Z_out,
      MSE = mse,
      Acc = accuracy
    )
  )
}

OmegaIdx <- function(G, op, N, delta, nc) {
  memoverlap <- as.data.frame(matrix(0, ncol = nc, nrow = N))
  
  for (i in 1:N) {
    memoverlap[i, ] <- as.numeric(op[[1]][i, ] > delta)
  }
  
  OrigVal <-  as.data.frame(vertex_attr(G)) %>%
    dplyr::select(all_of(c(letters[1:nc])))
  ## Calculating omega index
  oi <- 0
  for (i in 1:nc) {
    oi <- oi + as.numeric(OrigVal[, i] == memoverlap[, i])
  }
  oi <- oi / (N ^ 2)
}

sim <- TRUE
dir <- TRUE
if (sim == TRUE) {
  ## Creating simulation of the data
  #Types of covariate values
  ## Number of binary covk
  k_out <- 3
  k_in <- 3
  k <- k_out + k_in
  ## Number of continuous cov
  o_out <- 0
  o_in <- 0
  o <- o_in + o_out
  covTypes <- c("binary")#c()#c("continuous")
  CovNamesLinin <- c()#paste0("bv", 1:o)#c()#c("cv1")# c("cv1", "cv2","cv3")
  CovNamesLPin <- paste0("bvin", 1:k_in)#c()#c("bv1", "bv2","bv3")
  CovNamesLinout <- c()#paste0("bv", 1:o)#c()#c("cv1")# c("cv1", "cv2","cv3")
  CovNamesLPout <- paste0("bvout", 1:k_out)#c()#c("bv1", "bv2","bv3")
  ## setting alpha value
  nc <- c(3)#4
  alpha <- c(0.0005)#c(0.0002,0.0005,0.0008,0.001)#c(0.001)#c(0.0002,0.0005,0.0008,0.001)#0.0002
  alphaLL <- c(1)#c(0.7)#c(0, 0.3, 0.5, 0.7, 0.9, 1)#c(0.2,0.5,0.7)#0.5
  lambda <- 0.001
  
  ## Switching to percent increase. Using 0.01% increase
  thresh <- 0.001
  randomize = TRUE
  missing = 0
  randomize = TRUE
  
  ## Number of nodes
  N <- 100
  # Probability of cluster assignment
  pClust <- c(0.3, 0.3, 0.3, 0.0, 0.0, 0.0, 0)
  #c(0.7,0,0,0,0.7,0,0,0,0.7)
  # pC <- matrix(nrow=nc, ncol =k,data = c(1,0,0,
  #                                        0,1,0,
  #                                        0,0,1))
  pC <- matrix(
    nrow = k,
    ncol = k,
    byrow = TRUE,
    data = c(
      1,
      0,
      0,
      0,
      0,
      0,
      0,
      1,
      0,
      0,
      0,
      0,
      0,
      0,
      1,
      0,
      0,
      0,
      0,
      0,
      0,
      1,
      0,
      0,
      0,
      0,
      0,
      0,
      1,
      0,
      0,
      0,
      0,
      0,
      0,
      1
    )
  )
  dist <- data.frame(matrix(
    nrow = nc,
    ncol = o,
    data = list(list(5, 1), list(10, 1), list(15, 1))
  ))
  # dist <- data.frame(matrix(nrow = nc, ncol = o, data=list(list(5,1), list(0,1),list(0,1),
  #                                                          list(0,1), list(10,1),list(0,1),
  #                                                          list(0,1), list(0,1),list(15,1))))
  #dist <- cbind(c(5,10,15), c(1,1,1))
  #B <- c(0.3,0.05,0.3,0.05,0.05,0.3)
  #B <- c(0.1,0.01,0.1,0.01,0.01,0.1)
  ## directed network requires all the combinations of clusters.
  ##      aa  ba    ca    ab    bb   cb   ac   bc    cc
  #B <- c(0.15,0.01, 0.00, 0.00 ,0.15,0.01,0.01,0.00,0.15)
  
  ## this is Pc i.e. the probability of connection between two nodes of the same community.
  pConn <- c(0.5, 0.5, 0.5)
  epsilon <- 0.001
  dirPct <-  c(0.5, 0.5, 0.5)
  #overlapVec <- c(0.1,0.1,0.1)
  #pctInOut <- rep(0.5, k)
  # G <- intergraph::asIgraph(NWSimBin(nc, k, pC, N, pClust, B,
  #                                    o,dist,dir,covTypes,CovNamesLin,
  #                                    CovNamesLP, missing))
  
  ## Creating overlapping clusters
  NWlst <- genBipartite(
    N,
    nc,
    pClust,
    k_in,
    k_out,
    pC,
    covTypes,
    CovNamesLPin,
    pConn,
    dir,
    dirPct,
    epsilon,
    missing,
    k_out
  )
  G <- NWlst[[1]]
  Apuv <- NWlst[[2]]
  plot(
    G,
    vertex.label = NA,
    vertex.size = 5,
    vertex.color = "blue",
    edge.arrow.size = 0.1,
    edge.color = "grey28"
  )
  plot(
    G,
    vertex.label = NA,
    vertex.size = 5,
    vertex.color = as.factor(vertex_attr(G, name = "Cluster")),
    edge.arrow.size = 0.1,
    edge.color = "grey28"
  )
  #plot.igraph(G, vertex.color = V(G)$Cluster, vertex.label=NA)
  orig <- V(G)$Cluster
  Z <-  as.data.frame(vertex_attr(G)) %>%
    dplyr::select(all_of(c(CovNamesLin)))
  X <-  as.data.frame(vertex_attr(G)) %>%
    dplyr::select(all_of(c(CovNamesLP)))
  if (length(Z) > 0) {
    cov <- Z
    if (missing > 0) {
      mns <- colMeans(cov, na.rm  = TRUE)
      for (i in 1:o) {
        cov[is.na(cov[, i]), i] <- mns[i]
      }
    }
    origCov <-  CovAssistedSpecClust(G, cov, nc, alpha = 0.5)
  }
  if (length(X) > 0) {
    if (missing > 0) {
      cov <- X
      for (i in 1:k) {
        cov[is.na(cov[, i]), i] <- mode(cov[, i])
      }
    }
    origCov <-  CovAssistedSpecClust(G, cov, nc, alpha = 0.5)
  }
  nitermax <- 700
  Nsim <- 10
  opf <- list()
  ARIV <- rep(0, Nsim)
  arilst <- list()
  EG <- expand.grid(1:Nsim, alpha, alphaLL, nc)
  delta <- sqrt(-1 * (log(1 - (1 / N))))
  #Loglik <- matrix(0,nrow = nitermax+10, ncol = 3)
  #iterationNum <- 1
  ## saving the final cluster list
  cLst <- c(V(G)$Cluster)
  mse <- matrix(nrow = Nsim, ncol = o + k, 0)
  lvl <- 1
  
  for (a3 in nc) {
    ncVal <- a3
    for (a1 in alphaLL) {
      for (a2 in alpha) {
        for (j in 1:Nsim) {
          print(paste0(
            "iteration ",
            lvl,
            " alpha ",
            a2,
            " alphaLL ",
            a1,
            " num cmnty ",
            ncVal
          ))
          
          if (dir == TRUE) {
            op <- CoDA(
              G,
              ncVal,
              k,
              o,
              N,
              a2,
              lambda,
              thresh,
              nitermax,
              randomize,
              orig,
              CovNamesLin,
              CovNamesLP,
              NULL
            )
          } else{
            op <- CESNA(
              G,
              ncVal,
              k,
              o,
              N,
              a2,
              lambda,
              thresh,
              nitermax,
              randomize,
              orig,
              CovNamesLin,
              CovNamesLP,
              NULL
            )
          }
          #oi <- OmegaIdx(G, op, N, delta, nc)
          #print(paste0("The omega index for simulation Nsim: ", lvl," is ", oi))
          opf[[lvl]] <- op
          lvl <- lvl + 1
        }
      }
    }
  }
}

### Use the igraph for week 50-2020
if (sim == FALSE) {
  iGraph_op <- readRDS("Code/iGraph_NW.rds")
  G <- iGraph_op[["50-2020"]]
  G <- as.undirected(G, mode = "collapse")
  covTypes <- c()#c("continuous")#c("binary")#c("continuous")
  CovNamesLin <- c()#c("covVal")
  CovNamesLP <- c()
  k <- length(CovNamesLP)
  o <- length(CovNamesLin)
  nc <- c(3, 4)#4
  alpha <- c(0.0005, 0.0008, 0.001)#c(0.0002,0.0005,0.0008,0.001)#0.0002
  alphaLL <- c(1)#c(0.7)#c(0, 0.3, 0.5, 0.7, 0.9, 1)#c(0.2,0.5,0.7)#0.5
  lambda <- 0.001
  thresh <- 0.0001
  randomize = TRUE
  N <- gorder(G)
  #RegSpectralClust(G,nc)
  cov <- V(G)$covVal
  cov[is.na(V(G)$covVal)] <- mean(V(G)$covVal, na.rm  = TRUE)
  cov <- as.matrix(cov, ncol = 1)
  # Setting iterations if log likelihood is taking too long
  nitermax <- 800
  Nsim <- 20
  opf <- list()
  ARIV <- rep(0, Nsim)
  arilst <- list()
  ## saving the final cluster list
  cLst <- c(V(G)$Cluster)
  mse <- matrix(nrow = Nsim, ncol = o + k, 0)
  lvl <- 1
  for (a3 in nc) {
    ncVal <- a3
    for (a1 in alphaLL) {
      for (a2 in alpha) {
        for (j in 1:Nsim) {
          print(paste0(
            "iteration ",
            lvl,
            " alpha ",
            a2,
            " alphaLL ",
            a1,
            " num cmnty ",
            ncVal
          ))
          orig <- RegSpectralClust(G, ncVal)
          origCov <-  CovAssistedSpecClust(G, cov, ncVal, alpha = 0.5)
          op <- CESNA(
            G,
            ncVal,
            k,
            o,
            N,
            a2,
            lambda,
            thresh,
            nitermax,
            randomize,
            orig,
            CovNamesLin,
            CovNamesLP,
            a1
          )
          opf[[lvl]] <- op
          lvl <- lvl + 1
        }
      }
    }
  }
  EG <- expand.grid(1:Nsim, alpha, alphaLL, nc)
}

saveRDS(list(opf, op, EG, G), "OutputFile.rds")

## Junk code
group_ids_ol <- lapply(memoverlap %>% split(.$name), function(grp) {
  grp$id
})
group_color_ol <- brewer.pal(length(group_ids_ol), 'Set1')
group_color_fill_ol <- paste0(group_color, '20')
plot.igraph(
  G,
  vertex.label = NA,
  vertex.size = 5,
  vertex.color = as.factor(mem$cmnty),
  edge.arrow.size = 0.1,
  edge.color = "grey28",
  layout = lo,
  mark.groups = group_ids_ol,
  mark.col = group_color_fill_ol,
  mark.border = group_color_ol
)
plot(
  G,
  vertex.label = NA,
  vertex.size = 5,
  vertex.color = "blue",
  edge.arrow.size = 0.1,
  edge.color = "grey28"
)

mem <- rep(NA, N)
for (i in 1:N) {
  mem[i] <- letters[which.max(op[[1]][i, ])]
}
plot(
  G,
  vertex.label = NA,
  vertex.size = 5,
  edge.arrow.size = 0.1,
  edge.color = "grey28",
  vertex.color = as.factor(mem)
)

mem2 <- rep(NA, N)
for (i in 1:N) {
  mem2[i] <- letters[which.max(op[[3]][i, ])]
}
plot(
  G,
  vertex.label = NA,
  vertex.size = 5,
  edge.arrow.size = 0.1,
  edge.color = "grey28",
  vertex.color = as.factor(mem2)
)
nsim <- 1000
coe <- matrix(0, nrow = nsim, ncol = 3)
for (i in 1:nsim) {
  g <- intergraph::asNetwork(genBipartite(N, nc, pClust, B, k, pC, covTypes, CovNamesLP, pConn, epsilon))
  coe[i, ] <- coef(ergm(g ~ nodefactor("a") + nodefactor("b") + nodefactor("c")))
}

as.data.frame(coe) %>%
  mutate(col = 1:dim(coe)[1]) %>%
  pivot_longer(cols = !col) %>%
  ggplot() +
  geom_histogram(aes(value, color = name), fill = NA) +
  theme(legend.position = "none")
