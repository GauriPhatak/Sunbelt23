
library(tidyverse)
library(ggraph)


## Setting up the parameters

#3 number of networks of size 200
NWParams <- matrix(0, nrow = 0, ncol = 18)
## creating files for large networks.
##Directed graph yes? No?
dir <- "directed"

## Number of in and out binary covariates
k_out <- 5
k_in <- 0
k <- k_out + k_in

## Number of in and out continuous covariates
o_out <- 5
o_in <- 0
o <- o_in + o_out

## Types of covariates
covType <- list(c("continuous","binary"))

## Number of communities
nc_sim <- c(10,15,50,60,100,150)#rep(c(75,100,150), times = 2)#rep(c(5,6,8), times = 2)
## Beta distribution alpha and beta values
a_all <- rep(4,each = length(nc_sim))#rep(c(2,3), each =3)
b_all <- rep(9, length(a_all))

## Percent of missing data in covariates
missing <- NULL
## Setting missingnness type
missType <- rep(list(c("Random", "Random", "Random")), c(1))
MARparam <- rep(list(NULL),1)
## Number of nodes
N <- rep(c(600,6000,10000), each = 2)
count = length(N)

for(i in 1:count){
  
  nc <- nc_sim[i]
  a <- a_all[i]
  b <- b_all[i]
  # Probability of cluster assignment
  pClust <- list(rep(1/nc_sim[i],nc_sim[i]))
  
  ## trying random unequal group sizes
  min_val <- 0.05
  max_val <- 0.1
  
  # repeat {
  #   x <- runif(nc_sim[i], min_val, max_val)
  #   x <- x / sum(x)# normalize to sum to 1
  #   
  #   # after normalization values may violate bounds — check
  #   if (all(x >= min_val & x <= max_val)) break
  # }
  pClust <- list(rep(1/nc_sim[i],nc_sim[i] )) #list(sample_probabilities_scaled(nc_sim[i], 0.05, 0.40))
  
  ## Probability for overlap for clusters
  tmp <- rep(0, choose(nc_sim[i],2))
  #tmp[sample(1:length(tmp),floor(length(tmp)/2))] <- 0.05
  tmp[sample(1:length(tmp),0.3*nc_sim[i])] <- 0.05
  pClustOL <- list(c(tmp,0))
  # 
  pC1 <- matrix(
    nrow = k,
    ncol = 2,
    byrow = TRUE,
    data = c(0.9,0.1,
             0.9,0.1,
             0.9,0.1,
             0.8,0.2,
             0.8,0.2))
  
  
  
  pC <- list(pC1)#,pC2,pC3)
  ## Probability for covariate simulation for continuous covarites
  dist <- list(list(c(5,1), c(5,1),c(30,1),
                    c(20,1),c(20,1),c(65,1),
                    c(40,1),c(40,1),c(75,1),
                    c(10,1),c(10,1),c(5,1),
                    c(5,1),c(5,1),c(15,1)))
  
  ## Probability of connection between two nodes of the same community.
  pConn <- list(rep(0.5, nc_sim[i]))
  
  ## Percentage of outgoing edges in each cluster
  dirPct <- list(rep(0.5, nc_sim[i]))
  
  ## Type of network
  Type <- "Cohesive"
  epsilon <- 0.000001
  
  NWParams <- rbind(NWParams ,c(nc = nc, k_in=k_in, k_out=k_out, o_in=o_in, o_out=o_out, missing = missing,
                                N=N[i],pC = pC, pClust=pClust, pConn=pConn, dist = dist, epsilon=epsilon,
                                dirPct=dirPct, covType=covType, DirType = dir, Type = Type, 
                                pClustOL = pClustOL, a = a, b = b))
}
NWParams <- as.data.frame(NWParams)

NWConf <- function(Sim){
  
  ##Directed graph yes? No?
  dir <- Sim$DirType[[1]]
  
  ## Number of in and out binary covariates
  k_out <- Sim$k_out[[1]]
  k_in <- Sim$k_in[[1]]
  k <- k_out + k_in
  
  ## Number of in and out continuous covariates
  o_out <- Sim$o_out[[1]]
  o_in <- Sim$o_in[[1]]
  o <- o_in + o_out
  
  ## Naming the binary covariates
  CovNamesLPin <- c()
  CovNamesLPout <- c()
  CovNamesLinin <- c()
  CovNamesLinout <- c()
  
  if(k_in > 0 ){
    CovNamesLPin <- paste0("bvin", 1:k_in)
  }
  if(k_out > 0 ){
    CovNamesLPout <- paste0("bvout", 1:k_out)
  }
  
  ## Naming the continuous covariates
  if(o_in > 0 ){
    CovNamesLinin <- paste0("cvin", 1:o_in)
  }
  if(o_out > 0 ){
    CovNamesLinout <- paste0("cvout", 1:o_out)
  }
  
  
  ## Types of covariates
  covTypes <- Sim$covType[[1]]
  
  ## Number of communities
  nc_sim <- Sim$nc[[1]]
  
  ## Percent of missing data in covariates
  missing <- Sim$missing[[1]]
  ## Setting missingnness type
  missType <- Sim$mt[[1]]
  MARparam <- Sim$mar[[1]]
  ## Number of nodes
  N <- Sim$N[[1]]
  
  # Probability of cluster assignment
  pClust <- Sim$pClust[[1]]
  
  ## Probability for covariate simulation for binary covarites
  pC <- Sim$pC[[1]]
  
  ## Probability for covariate simulation for continuous covarites
  dist <- Sim$dist[[1]]
  
  ## Probability of connection between two nodes of the same community.
  pConn <- Sim$pConn[[1]]
  
  ## Percentage of outgoing edges in each cluster
  dirPct <- Sim$dirPct[[1]]
  
  ## Beta distribution alpha and beta values
  a <- Sim$a[[1]]
  b <- Sim$b[[1]]
  
  ## Type of network
  Type <- Sim$Type[[1]]
  
  ## Cluster Overlap proportion
  pClustOL <- Sim$pClustOL[[1]]
  
  seed <- sample(1:100000,1)
  ## Generating network with covariates and overlapping clusters
  NWlst <- genBipartite(N,nc_sim,pClust,k_in,k_out,o_in,o_out,pC,dist,covTypes,
                        CovNamesLPin,CovNamesLPout,CovNamesLinin,CovNamesLinout,
                        pConn,dir = "directed",dirPct,epsilon,missing, 
                        a, b, Type, pClustOL, missType, MARparam, seed)
  
  return(NWlst)
}

## number of networks for each type
NumNWs <- 10

for(i in 1:count){
  
  G_orig <- list()
  cov_orig <- list()
  
  for(j in 1:NumNWs){
    Sim <- NWParams[i,]
    NWlst <- NWConf(Sim)
    G_orig[[j]] <- NWlst[[1]]
    #delete_vertices(G_orig[[j]],degree(G_orig[[j]]) == 0 )
    print(ggraph(G_orig[[i]], layout = "fr") + 
            geom_edge_link(alpha = 0.5) + 
            geom_node_point(color = "grey25", alpha =0.5) +
            theme_void())
    cov_orig[[j]] <-  NWlst[[2]]
  }
  
  saveRDS(list(G_orig, cov_orig), paste0("NWexamples",i,".rds"))
  
}

