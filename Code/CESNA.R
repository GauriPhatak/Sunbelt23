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
NWSimBin <- function(nc, k, pC,N, pClust,B,o,dist,covTypes = c("binary", "continuous"),CovNamesLin = c() , CovNamesLP = c(), missing= NULL){
  C <- 0
  while(length(table(C)) < nc){
    C <- sample(x = 1:nc,size = N,replace = TRUE, prob = pClust)
  }
  ## Create an empty network i.e. no cluster assignments or category assignments or edges
  net <- network(N, directed = FALSE, density= 0)
  ## Assign clusters to the nodes 
  net %v% 'Cluster' <- C
  
  if(c("binary") %in% covTypes){
    ## Based on the probability matrix pC assign indicates the binary assignment of covariates to a community
    binVal <- matrix(ncol =k, nrow =N, NA)
    
    for(i in 1:nc){
      for(j in 1:k){
        binVal[which(as.numeric(net %v% 'Cluster') == i), j] <- rbinom(length(which(as.numeric(net %v% 'Cluster') == i)), 1, pC[i,j])
      }
    }
    for(j in 1:k){
      net %v% CovNamesLP[j] <- binVal[,j]
    }
  }
  
  if(c("continuous") %in% covTypes){
    ## Based on the mean and variance indicated by the dist matrix we assign continuous values to the network covariates.
    contVal <- matrix(ncol = o, nrow=N,NA)
    
    for(i in 1:nc){
      for(j in 1:o){
        contVal[which(as.numeric(net %v% 'Cluster') == i), j] <- rnorm(length(which(as.numeric(net %v% 'Cluster') == i)), dist[i,j][[1]][[1]], dist[i,j][[1]][[2]])
      }
    }
    for(j in 1:o){
      contVal[sample(1:N, round(missing*N/100)),j] <- NA
      net %v% CovNamesLin[j] <- contVal[,j]
    }
  }
  
  
  
  ## Use the probability of connection values defined earlier for both cluster and category for fitting the data. 
  g.sim <- simulate(net ~ nodemix("Cluster", levels = TRUE, levels2 = TRUE),
                    nsim = 1,
                    coef = prob2logit(c(B)),
                    control=ergm::control.simulate( MCMC.burnin=10000, MCMC.interval=1000))
  return(g.sim)
  
}

#f_u : Community weight vector for node u
#Fvmat : Community weight matrix for node neighbours
#Fvnot : Community weight matrix for not neighbours
#X_u: Covariate matrix
#Q_u: logistic probability matrix
#W: Logistic weight matrix
CmntyWtUpdt <- function(f_u, Fvmat,Fvnot, X_u = NULL , W = NULL ,Z_u= NULL, beta= NULL, sigmaSq, alpha){
  
  ## Gradient function to update log likelihood of G
  ## First part of log likelihood of G based on the structure of the network
  a <- exp(-1*(f_u %*% t(Fvmat)))
  b <- a/(1-a)
  llG_1 <- t(Fvmat) %*% t(b)
  ## 2nd part of log likelihood of G 
  llG_2 <- as.matrix(colSums(Fvnot))
  ## Total llG: This math has been verified
  llG <- llG_1 -llG_2
  ## experimental scaled version of the above llog lik
  scale <- 1
  llG <- scale *llG
  
  llX <- rep(0,nc+1)
  if(length(W) > 0 ){
    ## Updating the logistic regression weights
    ## Gradient function to update log likelihood for covariates
    ## adding intercept term
    Q_u <- t(1/(1+exp(-1*(W %*% as.matrix(c(1,f_u))))))
    ## This math has been verified
    llX <- t((X_u - Q_u) %*% W)
  }
  
  llZ <- rep(0,nc+1)
  if(length(beta) > 0 ){
    ## Adding the gradient based on continuous covariates
    llZ <-  ((Z_u - t(beta %*% as.matrix(c(1,f_u))))/sigmaSq) %*% beta
  }
  
  ## Function to update the community weights vector using non negative matrix factorization
  f_u_new <- f_u + t((alpha*(llG+llX[2:(nc+1)]+llZ[2:(nc+1)])))
  f_u_new[f_u_new < 0] <- 0.0001
  
  return(f_u_new)
}

sigmoid <- function(W,Ftot){
  Q <- 1/(1+exp(-1*(W %*% t(Ftot))))
  return(Q)
}

## W: Logistic weight matrix
## X: Covariate values
## Ftot: Community weight params
## alpha: tuning param
## lambda: 
LRParamUpdt <- function(W, X, Ftot, alpha, lambda){
  Ftotn <- cbind(1,Ftot)
  Q <- sigmoid(W,Ftotn)
  
  W_grad <- t(X - t(Q)) %*% Ftotn
  alpha <- 0.05
  W_new <- W + alpha*(W_grad - lambda*(sign(W)))
  # W_new <- matrix(nrow = dim(W)[1], ncol = dim(W)[2],0)
  # for(i in 1:k){
  #   Q <- t(1/(1+exp(-1* W[i,] %*% t(Ftotn))))
  #   op <- data.frame(matrix(nrow = length(Q), ncol = dim(W)[2], 0 ))
  #   for(j in 1:(dim(W)[2])){
  #     op[,j] <- (X[,i] - Q) * Ftotn[,j]
  #   }
  #   
  #   W_new[i,] <- W[i, ] + alpha*(colSums(op) - lambda*sign(W[i,]))
  # }
  
  #print(W_new)
  return(W_new)
}

LinRParamUpdt <- function(beta, Z, Ftot, alpha, lambda,missVals){
  #Z <- as.matrix(Z, ncol = dim(beta)[1])
  beta_new <- matrix(0, nrow = dim(beta)[1], ncol = dim(beta)[2])
  Ftot_n <- cbind(1,Ftot)
  for(i in 1:dim(beta)[1]){
    beta_new[i,] <- MASS::ginv(t(Ftot_n) %*% Ftot_n) %*% t(Ftot_n) %*% Z[,i] # (t(X)*X)^-1t(X)y : X : Ftot & y = Z
  }
  return(beta_new)
  
}

findLL <- function(G, Ftot, W= NA, X = NA, beta= NA ,Z= NA, sigmaSq = NA){
  scale <- 1
  E <- as_edgelist(G)
  F1 <- Ftot[E[,1],]
  F2 <- Ftot[E[,2],]
  N <-  gorder(G)
  ## Calculating the L_g part of the log likelihood
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
  
  ## Calculating the L_x part of the log likelihood
  S2 <- 0
  if(length(W) > 0){
    for(i in 1:N){
      Q_uk <- 1/(1+exp(-1*rowSums(c(1, Ftot[i,]) * W)))#1/(1+exp(-1*c(1, Ftot[i,]) %*% t(W)))
      S2 <- S2 + sum((log(Q_uk) * X[i,]) + (log(1-Q_uk) * (1-X[i,])))
    }
  }
  
  
  ## Calculating the L_z part of the log likelihood
  S3 <- 0
  if(length(beta) > 0){
    for(j in 1:o){
      S3_tmp <- 0
      for(i in 1:N){
        #print((Z[i,j] - sum(c(1, Ftot[i,]) * beta[j,], na.rm = TRUE) )^2)
        
        S3_tmp <- sum(c(S3_tmp, (Z[i,j] - sum(c(1, Ftot[i,]) * beta[j,], na.rm = TRUE) )^2),na.rm = TRUE)
      }
      S3 <- S3-(N*log(sigmaSq[j])/2 + S3_tmp/(2*sigmaSq[j]))
    }
  }
  
  
  ## Calculating the final log likelihood 
  #ll <- scale*S1 + (1-scale)* (S2 + S3)
  ll <- S1 +(S2 + S3)
  return(ll)
}

SigmaSqCalc <- function(Z, beta, Ftot, missVals){
  sigmaSq <- rep(0, dim(beta)[1])
  for(i in 1:dim(beta)[1]){
    sigmaSq[i] <- sum((Z[!missVals[,i],i] - (t(beta[i,]) %*% t(cbind(1,Ftot)[!missVals[,i],])) )^2, na.rm = TRUE)/sum(!missVals[,i])
      #colSums((Z - t(beta %*% t(cbind(1,Ftot)))) ^2, na.rm = TRUE)/dim(Z)[1]
  }
  
  return(sigmaSq)
}

PredictCovLin <- function(Ftot, Z, beta, missVals){
  for(i in 1:dim(Z)[2]){
    Z[which(missVals[,i]),i] <- cbind(1,Ftot[which(missVals[,i]),]) %*% as.matrix(beta[i,])
  }
  return(Z)
}

CESNA <- function(G, nc, k =0 , o=0 , N, alpha, lambda, thresh, nitermax, randomize =TRUE, orig, CovNamesLin = c() , CovNamesLP = c()){
  ## Setting initial values for logistic weights
  
  
  ## Community weights
  Finit <- igraph::degree(G)/sum(igraph::degree(G)) #cbind(Finit, Finit, Finit)
  Ftot <- matrix(nrow = N, ncol = nc, c(rep(Finit, times =  nc) + runif(nc*N)))##sample(seq(0.1,0.5,length.out =10), N*nc, replace = TRUE))
  Ftot_orig <- Ftot
  Ftotold <- Ftot
  
  iter <- 0
  
  ## Get the covariate matrix for binary data and continuous data
  covtmp <-  as.data.frame(vertex_attr(G))%>% 
    dplyr::select(all_of(c(CovNamesLin, CovNamesLP)))
  
  X <- 0
  Z <- 0
  
  if(k > 0){
    X <- as.matrix(covtmp %>%
                     dplyr::select(all_of(CovNamesLP)))
    missVals <- is.na(X)
    numMissVal <- sum(missVals)
  }
  if(o > 0){
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
  if(k >0){
    W <- matrix(nrow = k, ncol = nc+1, 0)
    W <- LRParamUpdt(W,X,Ftot,alpha,lambda)#matrix(nrow = k, ncol = nc+1, 0)#sample(seq(0.1,0.5,length.out =10), k*(nc+1), replace = TRUE))
    W_orig <- W
  }
  
  
  ## Setting initial values of linear regression weights
  beta <- NULL
  if(o > 0 ){
    beta <- matrix(nrow = o, ncol = nc+1, 0)
    if(numMissVal > 0){
      Z <- PredictCovLin(Ftot, Z, beta, missVals)
    }
    
    beta <- LinRParamUpdt(beta,Z,Ftot,alpha,lambda,missVals)#matrix(nrow = o, ncol = nc+1, 0)
    beta_orig <- beta
    sigmaSq <- SigmaSqCalc(Z,beta,Ftot, missVals)
    
  }
  
  ## setting the initial sigmaSq value
  #missIdx <-  !apply(missVals,  1, function(x) any(x))
  #if(numMissVal > 0){
  #  sigmaSq <- SigmaSqCalc(as.data.frame(Z[missIdx,]),beta, as.data.frame(Ftot[missIdx,]))
  #}else{
  #sigmaSq <- SigmaSqCalc(Z,beta,Ftot, missVals)
  #}
  
  ## Updating all the node community weights
  ## Setting the new loglikelihood to 0 
  #LLold <- -2000000
  LLold <- findLL(G,Ftot,W,X,beta,Z, sigmaSq)
  lllst <- c()
  lllst <- c(lllst, LLold)
  #LLnew <- -1000000
  tempop1 <- 0#list()
  LLG <- 0 #list()
  LLX <-0
  nodeid <-0
  arilst <- c()
  continue <- TRUE
  mse <- matrix(nrow= 0, ncol =o)
  accuracy <- matrix(nrow= 0, ncol =k)#c()
  #while (all(c(((LLnew - LLold) > thresh), (iter < nitermax))))  {
  while(continue){
    #print(iter)
    #print(iter)
    Ftotold <- Ftot
    ## Updating the commuynity weight parameter
    s <- 1:N
    if(randomize == TRUE){
      #set.seed(42)
      s <- sample(1:N, size = N, replace = FALSE)
      #set.seed(NULL)
    }
    for (i in s) {
      #print(i)
      f_u <- matrix(Ftot[i,], ncol = nc)
      ##The neighbours of this node
      neigh <- neighbors(G, i, mode = "all")
      Fvmat <- matrix(Ftot[neigh,], ncol = nc)
      Fvnot <- matrix(Ftot[-neigh,], ncol = nc )
      X_u <- NULL
      if(k > 0){
        X_u <- X[i,]
      }
      Z_u <- NULL
      if(o > 0){
        Z_u <- matrix(Z[i,], ncol = o)
        #if(numMissVal > 0){
        #  sigmaSq <- SigmaSqCalc(as.data.frame(Z[missIdx,]),beta, as.data.frame(Ftot[missIdx,]))
        #}else{
          sigmaSq <- SigmaSqCalc(Z,beta,Ftot, missVals)
          
        #}        
        }
      
      Ftot[i,] <- CmntyWtUpdt(f_u, Fvmat, Fvnot, X_u, W, Z_u, beta, sigmaSq,alpha)      
      tempop1 <- append(tempop1, list(Ftot[i,]))
      nodeid <- append(nodeid, i)
    }
    
    
    ## Updating logistic paramters
    if(length(W) > 0){
      W <-  LRParamUpdt(W, X, Ftot, alpha, lambda)
    }
    if(length(beta) > 0 ){
      if(numMissVal > 0 ){
        beta <- LinRParamUpdt(beta, Z, Ftot, alpha, lambda, missVals)
        sigmaSq <-  SigmaSqCalc(Z,beta,Ftot, missVals)#SigmaSqCalc(as.data.frame(Z[missIdx,]),beta, as.data.frame(Ftot[missIdx,]))
        Z <- PredictCovLin(Ftot, Z, beta, missVals)
      }else{
        beta <- LinRParamUpdt(beta, Z, Ftot, alpha, lambda, missVals )
        sigmaSq <-  SigmaSqCalc(Z,beta, Ftot,missVals)
      }
      
    }
    
    ## Look for Communities based on the F matrix
    mem <- rep(NA, N)
    for(i in 1:N){
      mem[i] <- which.max(Ftot[i,])
    }
    arilst <- c(arilst, mclust::adjustedRandIndex(orig, mem))
    
    ## Update the log likelihood.
    #LLold <- LLnew
    LLnew <- findLL(G,Ftot,W,X,beta, Z,sigmaSq )
    if(any(c(((LLnew - LLold) < thresh), (iter > nitermax)))){
      continue <- FALSE
    }
    LLold <- LLnew
    lllst <- c(lllst, LLnew)
    iter <- iter+1
    if(o >0){
      pred <- cbind(1,Ftot) %*% t(beta)
      mse <- rbind(mse, colSums((pred - Z)^2)/dim(Z)[1])
    }
    if(k > 0){
      predLR <- apply(sigmoid(W,cbind(1,Ftot)) > 0.5,1,as.numeric)
      acc <- rep(0,k)
      for(i in 1:k){
        cm <-  confusionMatrix(factor(predLR[,i], levels = c(0,1)),factor(X[,i]))
        acc[i] <- round(cm$overall[1],2)
      }
      #cm <-  confusionMatrix(factor(as.numeric(predLR)),factor(X))
      accuracy <- rbind(accuracy, acc)
    }
   
  }
  
  ## Updating the logistic weight parameters
  return(list(Ftotold,Ftot_orig,W_orig,arilst,W,lllst,beta, Z, mse,accuracy))
}


## Creating simulation of the data
nc <- 3
## Number of binary cov
k <- 3
## Number of continuous cov
o <- 0
## Number of nodes 
N <- 100
# Probability of cluster assignment
pClust <- c(0.3,0.3,0.3)
#c(0.7,0,0,0,0.7,0,0,0,0.7)
pC <- matrix(nrow=nc, ncol =k,data = c(1,0,0,
                                       0,1,0,
                                       0,0,1))

#dist <- data.frame(matrix(nrow = nc, ncol = o, data=list(list(5,1), list(5,1),list(5,1))))

dist <- data.frame(matrix(nrow = nc, ncol = o, data=list(list(5,1), list(0,1),list(0,1),
                                                         list(0,1), list(10,1),list(0,1),
                                                         list(0,1), list(0,1),list(15,1))))
#dist <- cbind(c(5,10,15), c(1,1,1))
B <- c(0.5,0.01,0.5,0.01,0.01,0.5)

## setting alpha value
alpha <- 0.001

## setting lambda value
lambda <- 0.001

## Setting threshold for log-likelihood difference
thresh <- 0.0001

#Types of covariate values
covTypes <- c("binary")#c("continuous")#c("binary")#c("continuous")
CovNamesLin <- c()#c("cv1", "cv2","cv3") 
CovNamesLP <- c("bv1", "bv2","bv3") 
missing = 0
randomize = TRUE
G <- intergraph::asIgraph(NWSimBin(nc, k, pC, N, pClust, B, o,dist,covTypes,CovNamesLin , CovNamesLP, missing))
orig <- V(G)$Cluster

## Use the igraph for week 50-2020
# iGraph_op <- readRDS("Code/iGraph_NW.rds")
# G <- iGraph_op[["50-2020"]]
# G <- as.undirected(G, mode = "collapse")
# CovNamesLin<- c("covVal") 
# CovNamesLP <- c()
# k <- 0
# o <- 1
# nc <- 6
# alpha <- 0.001
# lambda <- 0.001
# thresh <- 0.001
# randomize = FALSE
# N <- gorder(G)
# orig <- membership(cluster_edge_betweenness(G))

# Setting iterations if log likelihood is taking too long
nitermax <- 500
Nsim <- 20
opf <- list()
ARIV <- rep(0,Nsim)
arilst <- list()
## saving the final cluster list
cLst <- c(V(G)$Cluster)
mse <- matrix(nrow = Nsim, ncol = o+k, 0)
Z <-  as.data.frame(vertex_attr(G)) %>% 
  dplyr::select(all_of(c(CovNamesLin)))
X <-  as.data.frame(vertex_attr(G)) %>% 
  dplyr::select(all_of(c(CovNamesLP)))
for(j in 1:Nsim){
  print(j)
  op <- CESNA(G, nc, k,o, N, alpha, lambda, thresh, nitermax, randomize, orig, CovNamesLin, CovNamesLP)
  
   opf[[j]] <- op
}
if(o >0){
  Ffin <- op[[1]]
  mem <- rep(NA, N)
  for(i in 1:N){
    mem[i] <- which.max(Ffin[i,])
  }
  df <- cbind(1, opf[[j]][[1]])
  bt <- opf[[j]][[7]]
  pred <- df %*% t(bt)
  pred  <- cbind(pred,Z, V(G)$Cluster, mem, op[[8]])
  #mse[j, ] <- colSums((pred - Z)^2)/dim(Z)[1]
}

if(k >0 ){
  Ffin <- op[[1]]
  mem <- rep(NA, N)
  for(i in 1:N){
    mem[i] <- which.max(Ffin[i,])
  }
  df <- cbind(1, op[[1]])
  W <- op[[5]]
  pred <- apply(sigmoid(W,df) > 0.5,1,as.numeric)
  pred <- cbind(pred,X, V(G)$Cluster, mem)
}

