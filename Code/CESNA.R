library(tidyverse)
library(dplyr)
library(igraph)
library(network)
library(intergraph)
library(ergm)


prob2logit <- function(x){
  return(log(x / (1 - x)))
}

## nc : numer of communities
## k : number of covariates
## pClust: Probability of a ndoe belonging to a community
## N : Number of nodes
## pC : matrix of probability for binary assignment for each cluster 
NWSimBin <- function(nc, k, pC,N, pClust,B){
  C <- 0
  while(length(table(C)) < nc){
    C = sample(x = 1:nc,size = N,replace = TRUE, prob = pClust)
  }
  ## Create an empty network i.e. no cluster assignments or category assignments or edges
  net <- network(N, directed = FALSE, density= 0 )
  ## Assign clusters to the nodes 
  net %v% 'Cluster' <- C
  
  ## Based on the probability matrix pC assign indicates the binary assignment of covariates to a community
  binVal <- matrix(ncol =k, nrow =N, NA)
  # for(j in 1:k){
  #   print(pC[,j])
  #   for (i in 1:nc) {
  #     binVal[which(as.numeric(net %v% 'Cluster') == i), j]  <- rbinom(length(which(as.numeric(net %v% 'Cluster') == i)), 1, pC[i,j])
  #   }
  #   net %v% paste0("X",j) <- binVal[,j]
  # }
  
  for(i in 1:nc){
    for(j in 1:k){
      binVal[which(as.numeric(net %v% 'Cluster') == i), j] <- rbinom(length(which(as.numeric(net %v% 'Cluster') == i)), 1, pC[i,j])
    }
  }
  for(j in 1:k){
    net %v% paste0("X",j) <- binVal[,j]
  }
  
  ## Use the probability of connection values defined earlier for both cluster and category for fitting the data. 
  g.sim <- simulate(net ~ nodemix("Cluster", levels = TRUE, levels2 = TRUE),
                    nsim = 1,
                    coef = prob2logit(c(B)),
                    control=ergm::control.simulate( MCMC.burnin=10000, MCMC.interval=1000))
  return(g.sim)
  
}

#f_u : Community weight vector for node u
#Fvmat : Community weight matrix for nodes
#X_u: Covariate matrix
#Q_u: logistic probability matrix
#W: Logistic weight matrix
CmntyWtUpdt <- function(f_u, Fvmat,Fvnot, X_u, Q_u, W, alpha){
  
  ## Gradient function to update log likelihood of G
  ## First part of log likelihood of G
  a <- exp(-1*(f_u %*% t(Fvmat)))
  b <- a/(1-a)
  llG_1 <- t(Fvmat) %*% t(b)
  ## 2nd part of log likelihood of G
  llG_2 <- as.matrix(colSums(Fvnot))
  ## Total llG
  llG <- llG_1 -llG_2
  ## Updating the logistic regression weights
  ## Gradient function to update log likelihood for covariates
  ## adding intercept term
  Q_u <- t(1/(1+exp(-1*(W %*% as.matrix(c(1,f_u))))))
  llX <- t((X_u - Q_u) %*% W)
  
  ## Function to update the community weights vector using non negative matrix factorization
  f_u_new <- f_u + t((alpha *(llG+llX[1:nc])))
  f_u_new[f_u_new < 0] <- 0.0001
  
  return(f_u_new)
}


## W: Logistic weight matrix
## X: Covariate values
## Ftot: Community weight params
## alpha: tuning param
## lambda: 
LRParamUpdt <- function(W, X, Ftot, alpha, lambda){
  Ftotn <- cbind(1,Ftot)
  W_new <- matrix(nrow = dim(W)[1], ncol = dim(W)[2],0)
  for(i in 1:k){
    Q <- t(1/(1+exp(-1* W[i,] %*% t(Ftotn))))
    op <- data.frame(matrix(nrow = length(Q), ncol = dim(W)[2], 0 ))
    for(j in 1:(dim(W)[2])){
      op[,j] <- (X[,i] - Q) * Ftotn[,j]
    }
    
    W_new[i,] <- W[i, ] + alpha*(colSums(op) - lambda*sign(W[i,]))
  }
  
  #print(W_new)
  return(W_new)
}

findLL <- function(G, Ftot, W,X){
  E <- as_edgelist(G)
  F1 <- Ftot[E[,1],]
  F2 <- Ftot[E[,2],]
  S11 <- 0
  for(i in 1:dim(E)[1]){
    S11 <- S11 + log(1- exp(-1*sum(F1[i,] *F2[i,])))
  }
  
  A <- 1- (1 * as.matrix(as_adjacency_matrix(G)))
  E2<-  as_edgelist(graph_from_adjacency_matrix(A))
  F1 <- Ftot[E2[,1],]
  F2 <- Ftot[E2[,2],]
  S12 <- 0
  for(i in 1:dim(E2)[1]){
    S12 <- S12 + sum(F1[i,] * F2[i,])
  }
  S1 <- S11-S12
  S2 <- 0
  for(i in 1:gorder(G)){
    Q_uk <- 1/(1+exp(-1*rowSums(c(1, Ftot[i,]) * W)))#1/(1+exp(-1*c(1, Ftot[i,]) %*% t(W)))
    S2 <- S2 + sum((log(Q_uk) * X[i,]) + (log(1-Q_uk) * (1-X[i,])))
  }
  ll <- S1 + S2
  #print(paste0("Final log lik, ",ll))
  return(ll)
}

CESNA <- function(G, nc, k, N, alpha, lambda, thresh, nitermax, randomize =TRUE, orig){
  
  ## Setting initial values for lgistic weights
  W <- matrix(nrow = k, ncol = nc+1, sample(seq(0.0001,0.005,length.out =10000), k*(nc+1), replace = TRUE))
  
  ## Community weights
  Ftot <- matrix(nrow = N, ncol = nc,sample(seq(0.0001,0.005,length.out =10000), N*nc, replace = TRUE))
  Ftotold <- Ftot
  iter <- 0
  
  ## Get the covariate matrix
  X <- as.matrix(as.data.frame(vertex_attr(G))%>% 
                   dplyr::select(tail(names(.),k)))
  
  ## Updating all the node community weights
  ## Setting the new loglikelihood to 0 
  LLold <- -2000000
  LLnew <- -1000000
  tempop1 <- 0#list()
  LLG <- 0 #list()
  LLX <-0
  nodeid <-0
  
  #LLnew <- findLL(G,Ftot,W,X)
  arilst <- c()
  while (all(c(((LLnew - LLold) > thresh), (iter < nitermax))))  {
    #print(iter)
    Ftotold <- Ftot
    ## Updating the commuynity weight parameter
    s <- 1:N
    if(randomize == TRUE){
      s <- sample(1:N, size = N, replace = FALSE)
    }
    for (i in s) {
      f_u <- t(as.matrix(Ftot[i,]))
      ##The neighbours of this node
      neigh <- neighbors(G, i, mode = "all")
      Fvmat <- as.matrix(Ftot[neigh,])
      Fvnot <- as.matrix(Ftot[-neigh,])
      X_u <- X[i,]
      Ftot[i,] <- CmntyWtUpdt(f_u, Fvmat,Fvnot, X_u, Q_u, W, alpha)      
      tempop1 <- append(tempop1, list(Ftot[i,]))
      nodeid <- append(nodeid, i)
    }
    
    ## Updating logistic paramters
    W <-  LRParamUpdt(W, X, Ftot, alpha, lambda)
    
    ## Look for Communities based on the F matrix
    mem <- rep(NA, N)
    for(i in 1:N){
      mem[i] <- which.max(Ftot[i,])
    }
    arilst <- c(arilst, mclust::adjustedRandIndex(orig, mem))
    
    
    ## Update the log likelihood.
    LLold <- LLnew
    LLnew <- findLL(G,Ftot,W,X)
    iter <- iter+1
  }
  
  ## Updating the logistic weight parameters
  return(list(Ftotold, arilst))
}


## Creating simulation of the data
nc <- 3
k <- 3
N <- 100
pClust <- c(0.3,0.3,0.3)
#c(0.7,0,0,0,0.7,0,0,0,0.7)
pC <- matrix(nrow=nc, ncol =k,data = c(1,0,0,
                                       0,1,0,
                                       0,0,1))
B <- c(0.5,0.01,0.5,0.01,0.01,0.5)

## setting alpha value
alpha <- 0.001

## setting lambda value
lambda <- 0.001

## Setting threshold for log-likelihood difference
thresh <- 0.001
# Setting iterations if log likelihood is taking too long
nitermax <- 500
Nsim <- 20

G <- intergraph::asIgraph(NWSimBin(nc, k, pC, N, pClust, B))
ARIV <- rep(0,Nsim)
arilst <- list()
for(j in 1:Nsim){
  print(j)
  
  op <- CESNA(G, nc, k, N, alpha, lambda, thresh, nitermax, TRUE,c(V(G)$Cluster))
  
  Ffin <- op[[1]]
  arilst <- append(arilst, list(op[[2]]))
  ## Look for Communities based on the F matrix
  mem <- rep(NA, N)
  
  for(i in 1:N){
    mem[i] <- which.max(Ffin[i,])
  }
  Ffin <- as.data.frame(Ffin)
  Ffin$orig <- c(V(G)$Cluster)
  Ffin$cd <- mem
  ARIV[j] <- mclust::adjustedRandIndex(Ffin$orig, Ffin$cd)
  plot(G, vertex.color = mem,vertex.label = NA,vertex.size = 5)
}
ARIV_NewG <- ARIV

### Look for Communities based on the F matrix
plot(G, vertex.color = V(G)$Cluster,vertex.label = NA,vertex.size = 5)

Ffin <- as.data.frame(Ffin)
Ffin$orig <- c(V(G)$Cluster)
Ffin$cd <- mem
Ffin$x1 <- c(V(G)$X1)
Ffin$x2 <- c(V(G)$X2)
Ffin$x3 <- c(V(G)$X3)


mclust::adjustedRandIndex(Ffin$orig, Ffin$cd)

