library(igraph)
library(aricode)  # For NMI
library(linkprediction)  # For partition density (install via devtools::install_github("arc85/LinkPrediction"))
library(stringr)
library(tidyverse)
library(gmp)
library(e1071)

make_letter_names <- function(n) {
  # Convert number to base-26 alphabetic (Excel style)
  to_letters <- function(x) {
    out <- c()
    while (x > 0) {
      x <- x - 1
      out <- c((x %% 26) + 1, out)
      x <- x %/% 26
    }
    paste0(letters[out], collapse = "")
  }
  
  sapply(1:n, to_letters)
}


## Part of log likelihood contributed by the penalty term
penLL <- function(penalty, beta){
  if(penalty =="Ridge"){
    return(sum(beta^2))
    
  }else if(penalty == "LASSO"){
    return(sum(abs(beta)))
    
  }else if(penalty == "GroupLASSO"){
    return(sum(abs(beta)))
    
  }else if(penalty == "ElasticNet"){
    alpha  <- 0.5
    return(alpha * sum(abs(beta)) +(1-alpha)*sum(beta^2))
    
  }
}

## Find correlation between values in R
findCor <- function(Fm,Hm, dir){
  if(dir  == "directed"){
    m <- cor(cbind(Fm,Hm))
    m <- m[upper.tri(m)]
  }else{
    corV <- cor(Fm)
  }
  return(m)
}

## Calculate ARI for the network outputs
ARIop <- function(Ftot,Htot,orig,nc,N){
  mem <- rep(NA, N)
  
  for (i in 1:N) {
    m <- data.frame(com = rep(make_letter_names(nc) , times = 2),
                    val = c(Ftot[i,], Htot[i,]))
    
    if(!(length(which.max(m$val)) > 0)){
      print(which.max(m$val))
    }
    mem[i] <-  m$com[which.max(m$val)]
  }
  return(mclust::adjustedRandIndex(orig, mem))
}

MSE <- function(Wtmat1,Wtmat2,Z,beta,N,dir){
  #cbind(1, Wtmat)
  if(dir=="directed"){
    pred <- cbind(1, Wtmat1, Wtmat2) %*% beta
  }else{
    pred <- cbind(1, matrix(Wtmat1, ncol = (dim(beta)[1]-1))) %*% beta
  }
  mse <- sqrt(colSums((pred - Z) ^ 2) / N)
  return(mse)
}

MSEop <- function(Ftot, Htot, covOrig, betaout,N, dir, o_in,o_out, missValsout){
  mseouttot <- matrix(nrow = 0, ncol = (o_in+o_out))
  mse_outMD <- c()#matrix(nrow = 0, ncol = (o_in + o_out))
  covOrig <- as.matrix(covOrig)
  if((o_in + o_out) > 0){
    mseout <- 0
    mseouttot <- 0
    if (o_out > 0) {
      mseouttot <- MSE(Ftot,Htot,covOrig,t(as.matrix(betaout)),N,dir)#MSE(Ftot,Htot,Z_out,betaout,N,dir)
      for(i in 1:o_out){
        if(sum(missValsout[,i]) > 0){
          mse_outMD <- c(mse_outMD, MSE(Ftot[missValsout[,i],],
                                        Htot[missValsout[,i],],
                                        covOrig[missValsout[, i],i],
                                        as.matrix(betaout[i, ]), N, dir))
        }
        else{mse_outMD <- c(mse_outMD,0)}
      }
    }
  }
  
  return(list(c(mse_outMD) , c(mseouttot)))
}

## Logistic covariates parameter updates
sigmoid <- function(W, Ftotn) {
  Q <- 1 / (1 + exp(-1 * (W %*% t(Ftotn))))
  return(Q)
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
      cm <-  caret::confusionMatrix(factor(predLR[, i], levels = c(0, 1)), factor(X[, i]))
      acc[i] <- round(cm$overall[1], 2)
    }
    accuracy <- rbind(accuracy,acc)
  }
  
  return(accuracy)
}

accuOP <- function(k_in, k_out,Wout,Ftot, Htot,covOrig,dir,missValsout){
  
  accuracy_out <- matrix(nrow = 0, ncol = k_in + k_out)
  accuracy_outMD <- matrix(nrow = 0, ncol = k_in + k_out)
  if((k_in + k_out) > 0){
    #### Calculating the accuracy for the binary covariates
    
    ## Overall Accuracy
    accuracy_out <- accuCalc(k_out,Wout,Ftot,Htot,covOrig,dir)
    ## Missing data prediction accuracy
    if(sum(missValsout) >0 ){
      accuracy_outMD <- accuCalc(k_out,Wout,
                                 Ftot[as.logical(rowSums(missValsout)),],
                                 Htot[as.logical(rowSums(missValsout)),],
                                 covOrig[as.logical(rowSums(missValsout)),],dir)
    }
  }
  return(list(c(accuracy_outMD), c(accuracy_out)))
}

getDelta <- function(N, epsilon =0 ){
  delta <- sqrt(-1*log(1-(1/N)))#sqrt(-1*log(1 - epsilon)) #
  return(delta)
}

memOverlapCalc <- function(Fm, Hm, delta, N, nc){
  ol <- as.data.frame(matrix(0, ncol = nc, nrow = N))
  
  for (i in 1:N) {
    ol[i, ] <- (as.numeric(Fm[i, ] > delta) + as.numeric(Hm[i, ] > delta)) > 0
  }
  return(ol)
}

OmegaIdx_ <- function(A_mat, B_mat,N){
  # columns = clusters of each clustering (they can be different number of clusters)
  N <- nrow(A_mat)
  P <- choose(N, 2)
  
  coA <- A_mat %*% t(A_mat)  # The observed matrix
  coA <- coA[upper.tri(coA)]
  coB <- B_mat %*% t(B_mat) # The expected matrix
  coB <- coB[upper.tri(coB)]
  
  J <- max(coA)
  K <- max(coB)
  UL <- min(J,K)
  O <- E <- 0
  for(j in 0:UL){
    N_j1 <- (coA == j) ## the number of observed values equal to j  
    N_j2 <- (coB == j) ## the number of expected values equal to j
    M <- (N_j1 == TRUE) &  (N_j2 == TRUE) ## Observed and expected values in agreement
    A_j  <- sum(M) ## Total in agreement for j
    O <- O + A_j
    E <-sum.bigz(E, mul.bigz(as.bigz(sum(N_j1)), as.bigz(sum(N_j2))))
  }
  O <- O/P
  E <- as.numeric(E/as.bigz(P)^2)
  OI <- (O-E) / (1-E)
  ##OI <- as.double(div.bigz(sub.bigz(as.bigz(O), E), sub.bigz(1, E)))
  
  return(OI)
}

omega_index_fast <- function(obs_matrix, exp_matrix) {
  # Validate inputs
  if (nrow(obs_matrix) != nrow(exp_matrix)) {
    stop("Matrices must have the same number of nodes")
  }
  
  n <- nrow(obs_matrix)
  
  # Create indicator matrices for node pairs sharing communities
  create_overlap_matrix <- function(mat) {
    # mat: binary matrix n x m
    # Returns: n x n matrix where [i,j] = number of shared communities between i and j
    tcrossprod(mat)  # mat %*% t(mat) gives shared community count
  }
  
  obs_overlap <- create_overlap_matrix(obs_matrix)
  exp_overlap <- create_overlap_matrix(exp_matrix)
  
  # Convert to binary: 1 if nodes share at least one community, 0 otherwise
  obs_binary <- (obs_overlap > 0) * 1
  exp_binary <- (exp_overlap > 0) * 1
  
  # Remove diagonal (self-pairs)
  diag(obs_binary) <- diag(exp_binary) <- 0
  
  # Calculate agreement
  agreement <- sum(obs_binary == exp_binary) / (n * (n - 1))
  
  # Calculate expected agreement
  # Probability that nodes are together in each partition
  p_obs_together <- sum(obs_binary) / (n * (n - 1))
  p_exp_together <- sum(exp_binary) / (n * (n - 1))
  
  # Expected agreement if partitions are independent
  expected_agreement <- p_obs_together * p_exp_together + 
    (1 - p_obs_together) * (1 - p_exp_together)
  
  # Calculate Omega Index
  if (abs(1 - expected_agreement) < 1e-10) {
    return(0)
  }
  
  omega <- (agreement - expected_agreement) / (1 - expected_agreement)
  
  # Clip to [0, 1] range
  return(max(0, min(1, omega)))
}

OmegaIdx <- function(G, Fm, Hm, N, delta, nc, nc_sim) {
  # A_mat, B_mat: binary membership matrices, same number of rows (items)
  A_mat <- as.matrix(memOverlapCalc(Fm, Hm, delta,N, nc))
  B_mat <- as.data.frame(vertex_attr(G)) %>% 
    dplyr::select(any_of(c(make_letter_names(nc_sim)))) %>% 
    abs() %>% 
    as.matrix()
  OI <- OmegaIdx_(A_mat,B_mat,N)#omega_index_fast(A_mat,B_mat)#OmegaIdx_(A_mat,B_mat,N)
  return(OI)
}

## Calculating the silhouette score based on covariate values
SilhouetteScore <- function(C, Z, DistMeasure = "Euclidean"){
  nc <- ncol(C)
  if(DistMeasure == "Euclidean"){
    A <- as.matrix(dist(apply(Z, 2, scales::rescale, to = c(0,1)),method = "euclidean"))
  }else if(DistMeasure == "Jaccard"){
    A <- as.matrix(dist(apply(Z, 2, scales::rescale, to = c(0,1)), method = "binary"))
  }

  NodesWoAssignment <- (rowSums(C) == 0)
  if(sum(NodesWoAssignment) > 0){
    C[,nc+1] <- as.numeric(NodesWoAssignment)
  }
  
  #Setting intra cluster distances
  D <- matrix(0, nrow = nrow(C), ncol =ncol(C))
  
  for(i in 1:ncol(C)){
    idx <- which(C[,i] == 1)
    if(length(idx) > 1){
      tempA <- A[idx,idx]
      tempA[is.infinite(tempA)] <- 0
      D[idx, i] <- rowSums(tempA, na.rm = TRUE)/(length(idx) -1)
    }
  }
  
  ## average intra cluster distance for each cluster
  a_i <- rowSums(D)/rowSums(C)
  a_i[is.nan(a_i)] <- 0
  ## calculating inter cluster distances
  ## for each node calculate distance with nodes not in its cluster
  InterDist <- matrix(0,nrow = nrow(C), ncol = ncol(C))
  Percluster <- colSums(C)
  b_i <- rep(0, nrow(C))
  for(j in 1:nrow(C)){
    nonCols <- which(!(C[j,] == 1))
    if(length(nonCols) > 0 ){
      for(nonCol in nonCols){
        idx <- which(C[,nonCol] == 1)
        tempA <- A[j, idx]
        InterDist[j,nonCol] <- sum( tempA)/Percluster[nonCol]
      }
      b_i[j] <- min( InterDist[j, nonCols], na.rm =TRUE )
    }
  }
  S <- rep(0, nrow(C))
  for(j in 1:nrow(C)){
    S[j] <- ( b_i[j] - a_i[j] ) / max(a_i[j], b_i[j]) 
  }
  
  return(mean(S, na.rm = TRUE))
}

## Community wise average triangle participation ratio
Comm_TPR <- function(G, Fm, Hm, delta, N, nc, dir){
  C <- memOverlapCalc(Fm, Hm, delta, N, nc)
  tpr <- 0
  tprVec <- c()
  
  ##Create a new community of background nodes that as unassigned
  NodesWoAssignment <- (rowSums(memOverlapCalc(Fm,Hm,delta, N, nc)) == 0)
  numNodesWoAssignment <- sum(NodesWoAssignment)
  if(sum(NodesWoAssignment) > 0){
    
 
  C[,nc+1] <- as.numeric(NodesWoAssignment)
}
  for(i in 1:dim(C)[2]){
    idx <- which(C[,i] == 1)
    iG <- induced_subgraph(G,idx)
    
    # Count triangles for each vertex
    triangle_counts <- count_triangles(iG)
    
    # Identify vertices participating in triangles
    participating_vertices <- sum(triangle_counts > 0)
    
    # Calculate the triangle participation ratio
    total_vertices <- vcount(iG)

    tprVec <- c(tprVec,(participating_vertices / total_vertices) )
  }
  ## Triangle participation ratio of the whole network
  Total_triangle_counts <- count_triangles(G)
  
  # Identify vertices participating in triangles
  Total_participating_vertices <- sum(Total_triangle_counts > 0)
  
  # Calculate the triangle participation ratio
  total_vertices <- vcount(G)
  NWTPR <-  Total_participating_vertices / total_vertices
  
  ##TPR W the background cluster
  tprW <- mean(tprVec, na.rm =TRUE)
  
  ##TPR not including the new unassigned cluster
  tprWo <- mean(tprVec[1:nc], na.rm =TRUE)
  
  ## Weighted triangle participation ratio
  ##TPR W the background cluster
  tprWW <- tprW * (N - numNodesWoAssignment)/N
  
  ##TPR not including the new unassigned cluster
  tprWoW <- tprWo * (N - numNodesWoAssignment)/N
  
return(c(tprVec,tprW,tprWo,tprWW, tprWoW))  
  
}

## Ego splitting the graph to calculate conductance
EgoSplitConductance <- function(G,Fm , Hm, dir, delta , N , nc){
  C <- memOverlapCalc(Fm, Hm, delta, N, nc)

  G <- set_vertex_attr(G, "name", value = 1:vcount(G))
  ##Create a new community of background nodes that as unassigned
  NodesWoAssignment <- (rowSums(memOverlapCalc(Fm,Hm,delta, N, nc)) == 0)
  if(sum(NodesWoAssignment) > 0){
    C[,nc+1] <- as.numeric(NodesWoAssignment)
  }
  ## Get original graph edges
  edge_orig <- as_edgelist(G)
  colnames(edge_orig) <- c("from","to")
  ## list of induced groups
  Gsub <- list()
  
  ## List of intercluster links
  dropped_edges <- matrix(0,nrow = 0,ncol =2)
  for(i in 1:vcount(G)){
    neigh <- neighbors(G, v = i, mode = "total")
    target_nodes <- neigh[which((as.matrix(C[neigh, ]) %*% t(C[i,])) == 0)]
    if(length(target_nodes) > 0){
      dropped_edges <- rbind(dropped_edges, 
                             as.data.frame(edge_orig) %>%
                               filter(from == i & to %in% target_nodes))
    }
  }
  ## find the number of communities a node is part of. set vertex color based on that
  numCom <- rowSums(C)
  VertColor <- ifelse(numCom >1, "darkblue", "yellow4")
  G <- set_vertex_attr(G, "VertColor",value = VertColor)
  for(i in 1:dim(C)[2]){
    idx <- which(C[,i] == 1)
    Gsub[[i]] <- induced_subgraph(G, vids = idx)

    Gsub[[i]] <- set_vertex_attr(Gsub[[i]], "Origname", value = vertex_attr(Gsub[[i]],"name"))
    Gsub[[i]] <- set_vertex_attr(Gsub[[i]], "Group", value = i)
    Gsub[[i]] <- set_vertex_attr(Gsub[[i]], "VertColor", value = vertex_attr(Gsub[[i]],"VertColor"))
    Gsub[[i]] <- set_vertex_attr(Gsub[[i]], "name", value = paste0(vertex_attr(Gsub[[i]],"name"),"_",i))
    }
  CombinedGraph <- disjoint_union(Gsub)
  
  if(nrow(dropped_edges) > 0 ){
  ## keep distinct edges
  dropped_edges <- dropped_edges %>% distinct()
  
  CombinedGraph <- set_edge_attr(CombinedGraph, name = "color", value = "grey28")
  
  CombinedGraph <- add_edges(CombinedGraph, as.vector(t(dropped_edges)), color = "darkred")
  
  }
  
  Vatt <- as.data.frame(vertex_attr(CombinedGraph)) %>% 
    select(c(Origname, name, Group)) %>%
    group_by(Origname) %>%
    filter(n() > 1) %>%
    tidyr::expand(item1 = name, item2 = name)
  if(nrow(Vatt) > 0 ){
    if(dir == "directed"){
      Vatt <- Vatt %>% 
        filter(item1 != item2) %>%  # Avoid duplicates and self-pairs
        ungroup() 
    }else{
      Vatt <- Vatt %>% 
        filter(item1 < item2) %>%  # Avoid duplicates and self-pairs
        ungroup()  
    }

    CombinedGraph <- add_edges(CombinedGraph, as.vector(t(Vatt[,2:3])), color = "darkred")
    
  }

  
  CombinedGraph <- igraph::simplify(CombinedGraph,remove.multiple = TRUE,remove.loops = TRUE)
  
  ConductanceVal <- rep(0,dim(C)[2])
  for (i in 1:dim(C)[2]) {
    S <- V(CombinedGraph)$name[which(as.numeric(V(CombinedGraph)$Group) == i)]
    Sc <- V(CombinedGraph)$name[which(as.numeric(V(CombinedGraph)$Group) != i)]
    
    # Calculate cut size between S and Sc
    cut_edges <- 0
    for (edge in E(CombinedGraph)) {
      ends <- ends(CombinedGraph, edge)
      if ((ends[1] %in% S && ends[2] %in% Sc) || 
          (ends[1] %in% Sc && ends[2] %in% S)) {
        cut_edges <- cut_edges + 1
      }
    }
    # Calculate volume of S
    vol_S <- sum(igraph::degree(CombinedGraph)[which(as.numeric(V(CombinedGraph)$Group) == i)])
    vol_Sc <- sum(igraph::degree(CombinedGraph)[which(as.numeric(V(CombinedGraph)$Group) != i)])
    if(min(vol_S, vol_Sc) > 0){
      ConductanceVal[i] <- cut_edges / min(vol_S, vol_Sc)
    }
  }
  ##punish internal density if nodes have no assignment
  numNodesWoAssignment <- sum( rowSums(memOverlapCalc(Fm,Hm,delta, N, nc)) == 0 )
  
  ## MEan conductance W the new background cluster
  MeanConductanceW <- mean(ConductanceVal,na.rm =T)
  
  ## Mean conductance without the new background cluster
  MeanConductanceWo <- mean(ConductanceVal[1:nc],na.rm =T )
  
  ## Punish based on unassigned nodes
  WeightedMeanConductanceW <- MeanConductanceW*(N/(N - numNodesWoAssignment ))
  WeightedMeanConductanceWo <- MeanConductanceWo*(N/(N - numNodesWoAssignment ))
  
 return(c(ConductanceVal, MeanConductanceW, MeanConductanceWo,WeightedMeanConductanceW, WeightedMeanConductanceWo))
}

##internal density of a cluster. Need to be maximized
InternalDensity <- function(G, d, epsilon, dir){
  
  Fm <- d$Ffin
  Hm <- d$Hfin
  
  N <- dim(Fm)[1]
  nc <- dim(Fm)[2]
  
  delta <- getDelta(N, epsilon)
  
  C <- memOverlapCalc(Fm, Hm, delta, N, nc)
  ##Create a new community of background nodes that as unassigned
  NodesWoAssignment <- (rowSums(memOverlapCalc(Fm,Hm,delta, N, nc)) == 0)
  if(sum(NodesWoAssignment) > 0){
    
   C[,nc+1] <- as.numeric(NodesWoAssignment)
}
  InternalDensityVal <- rep(0,dim(C)[2])
  n <- rep(0,dim(C)[2])
  for(i in 1:dim(C)[2]){
    idx <- which(C[,i] == 1)
    n[i] <- length(idx)
    Gsub <- induced_subgraph(G, vids = idx)
    if(dir == "directed"){
      possEdges <- n[i] * (n[i]-1)
    }else{
      possEdges <- n[i] * (n[i]-1)/2
    }
    InternalDensityVal[i] <- ecount(Gsub)/possEdges
  }
  ## considering the background custer
  CommInternalDensityW <- sum(InternalDensityVal, na.rm = TRUE)
  
  ## without considering the background cluster
  CommInternalDensityWo <- sum(InternalDensityVal[1:nc], na.rm = TRUE)
  
  
  ## Calculating the total internal density
  if(dir == "directed"){
    possEdges <- N * (N-1)
  }else{
    possEdges <- N * (N-1)/2
  }
  TotalInternalDensity <- ecount(G)/possEdges
  
  ##punish internal density if nodes have no assignment
  numNodesWoAssignment <- sum( rowSums(memOverlapCalc(Fm,Hm,delta, N, nc)) == 0 )
  WeightedInternalDensityW <- CommInternalDensityW* (N-numNodesWoAssignment)/N
  WeightedInternalDensityWo <- CommInternalDensityWo* (N-numNodesWoAssignment)/N
  
  return(c(InternalDensityVal, 
           CommInternalDensityW,CommInternalDensityWo,TotalInternalDensity,
           WeightedInternalDensityW,WeightedInternalDensityWo))
}

##calculate attribute cohesion
## Dissimilarity score or distance measure. Have to minimize.
## also includes non normalized variance. Does not take scales of the attributes into consideration
AverageDissimilarityScore <- function(d,epsilon){
  ##reading the covariates (cont.)
  Z <- as.data.frame(d$Zout_cov)#as.data.frame(scale(d$Zout_cov)) 
  
  ## reading the community weight matrices
  Fm <- d$Ffin
  Hm <- d$Hfin
  
  ## Number of nodes and communities
  N <- dim(Fm)[1]
  nc <- dim(Fm)[2]
  
  ##Decide the community affiliations
  delta <- getDelta(N, epsilon)
  C <- memOverlapCalc(Fm, Hm, delta, N, nc)
  
  ##Create a new community of background nodes that as unassigned
  NodesWoAssignment <- (rowSums(C) == 0)
  if(sum(NodesWoAssignment) > 0){
    C[,nc+1] <- as.numeric(NodesWoAssignment)
  }
  ## find number of nodes wihtout assignment
  numNodesWoAssignment <- sum( rowSums(memOverlapCalc(Fm,Hm,delta, N, nc)) == 0 )
  ## Set column names to communities
  ColNames <- make_letter_names(dim(C)[2])
  colnames(C) <- ColNames#letters[1:dim(C)[2]]
  ## number of covariates
  n_cov <- dim(Z)[2]
  ##Creating variance matrix
  vm <- matrix(0, nrow = dim(C)[2], ncol = n_cov)
  vmNorm <- matrix(0, nrow =  dim(C)[2], ncol = n_cov)
  
  gm <- matrix(0, nrow = 1, ncol = n_cov)
  n <- colSums(C)
  for(i in 1:n_cov){
    ##Whole covariate variance
    gm[1,i] <- var(Z[,i])#sum((Z[,i] - mean(Z[,i]))^2) / N
    for(j in 1:dim(C)[2]){
      ## for each attribute a find the variance in cluster C
      idx <- which(C[,j] == 1)
      vm[j, i] <- var(Z[idx, i])#sum((Z[idx, i] - mean(Z[idx,i]))^2) / length(idx)
      vmNorm[j,i] <- vm[j,i]/gm[1,i] ##Normalized variance per covariate per cluster. Gives dispersion score
    }
  }
  
  AvgVarianceW <- sum((rowSums(vm)/n_cov) * (n), na.rm =  TRUE)/N
  DispersionScoreW <- sum((rowSums(vmNorm)/n_cov) * (n), na.rm = TRUE)/N
  
  ## Average variance and dispersion score without taking the unassigned nodes into consideration
  if(ncol(vm) >1){
    AvgVarianceWo <- sum((rowSums(vm[1:nc, ])/n_cov) * (n[1:nc]), na.rm = TRUE)/N
    DispersionScoreWo <- sum((rowSums(vmNorm[1:nc, ])/n_cov) * (n[1:nc]), na.rm = TRUE)/N
  }else{
    AvgVarianceWo <- sum((vm[1:nc, ]/n_cov) * (n[1:nc]), na.rm = TRUE)/N
    DispersionScoreWo <- sum((vmNorm[1:nc, ]/n_cov) * (n[1:nc]), na.rm = TRUE)/N
  }

  ## Punish the value of variance and dispersion if there are nodes that have not been assigned
  WeightedAvgVarianceW <- AvgVarianceW * (N/(N - numNodesWoAssignment ))
  WeightedDispersionScoreW <- DispersionScoreW * (N/(N - numNodesWoAssignment ))
  
  return(c(AvgVarianceW, DispersionScoreW,AvgVarianceWo, DispersionScoreWo,WeightedAvgVarianceW,WeightedDispersionScoreW))
}

## Calculating trace to determine the within cluster variance
WeightedTrace <- function(d, epsilon){
  Z <-  as.data.frame(d$Zout_cov)#as.data.frame(scale(d$Zout_cov))
  
  ## reading the community weight matrices
  Fm <- d$Ffin
  Hm <- d$Hfin
  
  ## Number of nodes and communities
  N <- dim(Fm)[1]
  nc <- dim(Fm)[2]
  
  ##Decide the community affiliations
  delta <- getDelta(N, epsilon)
  C <- memOverlapCalc(Fm, Hm, delta, N, nc)
  
  ##Create a new community of background nodes that as unassigned
  NodesWoAssignment <- (rowSums(C) == 0)
  if(sum(NodesWoAssignment) > 0){
    C[,nc+1] <- as.numeric(NodesWoAssignment)
  }
  ## find number of nodes wihtout assignment
  numNodesWoAssignment <- sum( rowSums(memOverlapCalc(Fm,Hm,delta, N, nc)) == 0 )
  ## Set column names to communities
  ColNames <- make_letter_names(dim(C)[2])
  colnames(C) <- ColNames#letters[1:dim(C)[2]]
  ## number of covariates
  n_cov <- dim(Z)[2]
  
  ## within cluster trace
  tr_i <- rep(0, dim(C)[2])
  
  ## pooled within cluster trace
  W <- matrix(0, ncol = ncol(Z), nrow = ncol(Z))
  
  ##Calculate trace per cluster
  for(i in 1:dim(C)[2]){
    idx <- which(C[,i] == 1)
    nk <- length(idx)
    if(nk > 1){
      Zk <- Z[idx, ]
      cov_i <- cov(Zk)
      tr_i[i] <- sum(diag(cov_i), na.rm =T)*nk/N
      W <- W + (nk - 1) * cov_i
    }
  }
  
  W <- W / (N - dim(C)[2])
  Pooled_within_cluster_trace <- sum(diag(W))
  Trace <- sum(tr_i)
  
  ##Proportion of variance explained
  TotalTrace <- sum(diag(cov(Z)))
  R2 <- 1 - Pooled_within_cluster_trace/TotalTrace
  
  ## Punish the value of variance and dispersion if there are nodes that have not been assigned
  WeightedTraceW <- Trace * (N/(N - numNodesWoAssignment ))

  return(c(Trace, WeightedTraceW, TotalTrace, R2))
}

##Binary covariate dispersion score
##Understanding cohesiveness within and between clusters 
BinaryDispersionScore <- function(X, Fm, Hm, epsilon){
  # ##reading the covariates (cont.)
  # Z <- as.data.frame(d$Xout_cov) 
  # 
  # ## reading the community weight matrices
  # Fm <- d$Ffin
  # Hm <- d$Hfin
  
  ## Number of nodes and communities
  N <- dim(Fm)[1]
  nc <- dim(Fm)[2]
  
  ##Decide the community affiliations
  delta <- getDelta(N, epsilon)
  C <- memOverlapCalc(Fm, Hm, delta, N, nc)
  
  
  ##Create a new community of background nodes that as unassigned
  NodesWoAssignment <- (rowSums(C) == 0)
  if(sum(NodesWoAssignment) > 0){
  C[,nc+1] <- as.numeric(NodesWoAssignment)
  }
  ## find number of nodes wihtout assignment
  numNodesWoAssignment <- sum( rowSums(memOverlapCalc(Fm,Hm,delta, N, nc)) == 0 )
  
  ## Set column names to communities
  colnames(C) <- make_letter_names(dim(C)[2]) #letters[1:dim(C)[2]]
  ## number of covariates
  n_cov <- dim(X)[2]
  
  ##Number of nodes per cluster
  num_perClust <- colSums(C)
  
  ##Creating probability matrix
  pm <- matrix(0, nrow = dim(C)[2], ncol = n_cov)
  pmNorm <- matrix(0, nrow =  dim(C)[2], ncol = n_cov)
  
  gm <- matrix(0, nrow = 1, ncol = n_cov)
  n <- colSums(C)
  for(i in 1:n_cov){
    ##Whole covariate variance
    p <- sum(X[,i])/N
    gm[1,i] <- p*(1-p)   #sum((Z[,i] - mean(Z[,i]))^2) / N
    for(j in 1:dim(C)[2]){
      ## for each attribute a find the variance in cluster C
      idx <- which(C[,j] == 1)
      p_cov <- sum(X[idx,i])/length(idx)
      pm[j, i] <- p_cov*(1-p_cov) #sum((Z[idx, i] - mean(Z[idx,i]))^2) / length(idx)
      pmNorm[j,i] <- pm[j,i]/gm[1,i] ##Normalized variance per covariate per cluster. Gives dispersion score
    }
  }
  
  AvgBinarySpreadW <- sum((rowSums(pm)/n_cov) * (n), na.rm =  TRUE)/N
  BinaryDispersionScoreW <- sum((rowSums(pmNorm)/n_cov) * (n), na.rm = TRUE)/N
  
  ## Average variance and dispersion score without taking the unassigned nodes into consideration
  if(ncol(pm) >1){
    AvgBinarySpreadWo <- sum((rowSums(pm[1:nc, ])/n_cov) * (n[1:nc]))/N
    BinaryDispersionScoreWo <- sum((rowSums(pmNorm[1:nc, ])/n_cov) * (n[1:nc]))/N
  }else{
    AvgBinarySpreadWo <- sum((pm[1:nc, ]/n_cov) * (n[1:nc]))/N
    BinaryDispersionScoreWo <- sum((pmNorm[1:nc, ]/n_cov) * (n[1:nc]))/N
  }
  
  
  ## Punish the value of variance and dispersion if there are nodes that have not been assigned
  WeightedAvgBinarySpreadW <- AvgBinarySpreadW * (N/(N - numNodesWoAssignment ))
  WeightedBinaryDispersionScoreW <- BinaryDispersionScoreW * (N/(N - numNodesWoAssignment ))
  
  # ## Calculating the Within cluster dispersion for each feature for each cluster
  # ##For pure feature the dispersion value will be low.
  # disp_within <- matrix(0, ncol = n_cov, nrow = ncol(C))
  # for(i in 1:ncol(C)){
  #   idx <- which(C[,i] == 1)
  #   for(j in 1:n_cov){
  #     p_ji <-  sum(X[idx,j])/sum(X[,j])
  #     disp_within[i, j] <- p_ji*(1-p_ji)
  #   }
  # }
  # 
  # ##Average over the number of values percluster and then over the number of clusters
  # AvgDispersionWithinCluster <- mean(rowSums(disp_within)/num_perClust, na.rm = TRUE)
  # 
  # ##Weighted by the number of unassigned nodes
  # WeightedAvgDispersionWithinClusterW <- AvgDispersionWithinCluster * (N/(N - numNodesWoAssignment ))
  # 
  # ##Calculating the between cluster dispersion
  # ##Have not included this yet in assessment
  # p_k <- colSums(C)/N
  # disp_bet <- rep(0, n_cov)
  # for(i in 1:n_cov){
  #   disp_bet[i] <- sum(num_perClust * (disp_within[,i] - p_k)^2)
  # }
  
  return(c(AvgBinarySpreadW,BinaryDispersionScoreW,WeightedAvgBinarySpreadW, WeightedBinaryDispersionScoreW))
}

## similarity within each cluster
# Have to maximizehigher the better
## Need to implement Jaccard similarity for categorical variables
AverageSimilarityScore <- function(d, epsilon){
  
  Z <- as.data.frame(d$Zout_cov) 
  
  Fm <- d$Ffin
  Hm <- d$Hfin
  
  N <- dim(Fm)[1]
  nc <- dim(Fm)[2]
  
  delta <- getDelta(N, epsilon)
  C <- memOverlapCalc(Fm, Hm, delta, N, nc)
  
  ##Create a new community of background nodes that as unassigned
  NodesWoAssignment <- (rowSums(memOverlapCalc(Fm,Hm,delta, N, nc)) == 0)
  if(sum(NodesWoAssignment) > 0){
    
   C[,nc+1] <- as.numeric(NodesWoAssignment)
}
  numNodesWoAssignment <- sum( rowSums(memOverlapCalc(Fm,Hm,delta, N, nc)) == 0 )
  
  colnames(C) <- make_letter_names(dim(C)[2]) #letters[1:dim(C)[2]]
  
  n_cov <- dim(Z)[2]
  sim <- rep(0,dim(C)[2])
  numClust <- rep(0,dim(C)[2])
  Z_normalized <- apply(Z, 2, function(x) (x - min(x)) / (max(x) - min(x)))
  for(i in 1:dim(C)[2]){
    idx <- which(C[,i] == 1)
    numClust[i] <- length(idx)
    if(length(idx)>0){
      sim[i] <- AvgCosineSimilarity(Z_normalized[idx,],n_cov)
    }
  }
  
  ##Global average intracluster similarity With assigned nodes
  AvgIntraClusterSimilarityW <- sum(sim * numClust,na.rm =TRUE)/N
  
  ##Global average intracluster similarity Without assigned nodes
  AvgIntraClusterSimilarityWo <- sum(sim[1:nc] * numClust[1:nc], na.rm=TRUE)/N
  
  ## Punish the score for number of nodes not assigned
  WeightedAvgIntraClusterSimilarityW <- AvgIntraClusterSimilarityW * (N-numNodesWoAssignment)/N
  WeightedAvgIntraClusterSimilarityWo <- AvgIntraClusterSimilarityWo * (N-numNodesWoAssignment)/N
  
  return(c(AvgIntraClusterSimilarityW, AvgIntraClusterSimilarityWo, WeightedAvgIntraClusterSimilarityW, WeightedAvgIntraClusterSimilarityW))
  
}

AvgCosineSimilarity <- function(Z,n_cov){
  Z <- as.matrix(Z, ncol = n_cov)
  
  L2_Norm <- apply(Z, 1 , function(x) sqrt(sum(x^2)))
  X <- as.matrix(Z/ L2_Norm)
  CosSim <- X %*% t(X)
  AvgCosSim <- mean(CosSim[upper.tri(CosSim)], na.rm = TRUE)
  return(AvgCosSim)
}

conductance <- function(graph, communities, ground_truth) {
  #weights <- matrix(1, nrow = vcount(graph), ncol = length(communities))/ rowSums(ground_truth)
  weights <- as.matrix((ground_truth + rowSums(ground_truth))/rowSums(ground_truth))
  weights[is.nan(weights)] <- 0
  totCond <- 0
  for(k in seq_along(communities)){
    boundary <- 0
    vol_S <- 0
    comm <- unlist(communities[[k]])
    node_weights <- weights[,k]
    for (i in seq_along(comm)) {
      neighbors <- igraph::neighbors(graph, v = comm[i])
      #print(paste0(" Weight if i ",node_weights[comm[i]]))
      for (j in neighbors) {
        if (!(j %in% comm)) {
          #print(paste0("For node ", j," the weight",node_weights[j]))
          boundary <- boundary + node_weights[comm[i]] * node_weights[j]
        }
      }
      vol_S <- vol_S + igraph::degree(graph, comm[i]) * node_weights[comm[i]]
    }
    #print(vol_S)
    #print(sum(igraph::degree(graph)) - vol_S)
    vol_notS <- sum(igraph::degree(graph)) - vol_S
    #print(boundary / min(vol_S, vol_notS))
    totCond <- totCond + (boundary / min(vol_S, vol_notS))  
  }
  return(totCond)
}

## Non dominated sorting used for selecting hyperparameters
HyperParameterSelection <- function(metricsCov, cols_to_select){
  metricsCov$front <- 0
  df_list <- metricsCov %>% 
    filter(unassigned <= degree01) %>%
    group_by(bigN, OL, dir, pctMiss) %>%
    group_split() 
  fronts_list <- list()
  for(i in 1:length(df_list)){
    if(length(cols_to_select) > 1 ){
      X <-  df_list[[i]] %>%
        dplyr::select(!!!rlang::syms(cols_to_select)) %>% #WeightedMeanConductanceW, WeightedDispersionScoreW, WeightedBinaryDispersionScoreW) %>% 
        t()
      if(any(is.na(X))){
        print("Stop here")
      }else{
        fronts <- ecr::doNondominatedSorting(X)
        df_list[[i]]$front <-  fronts$ranks
      }
    }else{
      X <-  df_list[[i]] %>%
        dplyr::select(!!rlang::sym(cols_to_select)) %>% #WeightedMeanConductanceW, WeightedDispersionScoreW, WeightedBinaryDispersionScoreW) %>% 
        t() 
      fronts <- rank(X, ties.method = "first") #ecr::doNondominatedSorting(X)
      df_list[[i]]$front <-  fronts
      
    }

  }
  
  FrontTotal <- bind_rows(df_list) %>% ungroup()
  #FirstFront <- bind_rows(df_list) %>% filter(front %in% c(1)) %>% arrange(bigN, dir, OL, pctMiss)
  #FirstTwoFronts <- bind_rows(df_list) %>% filter(front %in% c(1,2)) %>% arrange(bigN, dir, OL, pctMiss, front)
  
  return(FrontTotal)
}

partition_density <- function(graph, communities) {
  densities <- sapply(communities, function(comm) {
    subg <- induced_subgraph(graph, unlist(comm))
    m <- ecount(subg)
    n <- vcount(subg)
    if (n <= 1) return(0)
    m / (n * (n - 1) / 2)
  })
  mean(densities)
}

## Newman girvan hard partition modularity
Q_HardPartition <- function(graph, alpha) {
  A <- as_adjacency_matrix(graph, sparse = FALSE)
  M <- ecount(graph)
  k <- rowSums(A)
  Q <- 0
  for (C in 1:ncol(alpha)) {
    for (i in 1:nrow(A)) {
      for (j in 1:nrow(A)) {
        if(A[i,j]!=0){
          Q <- sum(c(Q, (A[i,j] - k[i]*k[j]/(2*M)) * alpha[i,C] * alpha[j,C]), na.rm=TRUE)
        }
      }
    }
  }
  Q / (2 * M)
}

## Overlapping modularity by nicosa et al.
Q_Overlapping <- function(graph, membership_matrix, dir){
  A <- as_adjacency_matrix(graph, sparse = FALSE)
  
  if (nrow(A) != nrow(membership_matrix)) {
    stop("Adjacency and membership matrices must have the same number of nodes.")
  }
  if(dir == "U"){
    m <- sum(A) / 2  # Total edge weight (undirected)  
  }else{
    m <- sum(A)  # Total edge weight (directed)
  }
  
  if (m == 0) return(0)
  
  k <- rowSums(A)   # Node degrees
  n <- nrow(A)      # Number of nodes
  Q <- 0                     # Initialize modularity
  
  # Precompute normalization terms for each node
  norm_i <- sqrt(rowSums(membership_matrix^2))
  
  for (i in 1:n) {
    for (j in 1:n) {
      if (A[i,j] != 0) {
        # Compute the community overlap term
        overlap_term <- sum(membership_matrix[i,] * membership_matrix[j,]) / 
          (norm_i[i] * norm_i[j])
        
        # Add to modularity
        Q <- sum(c(Q, (A[i,j] - (k[i] * k[j]) / (2 * m)) * overlap_term), na.rm = TRUE)
      }
    }
  }
  
  Q <- Q / (2*m)  # Normalize
  return(Q)
}

omega_index <- function(comm1, comm2) {
  # Convert communities to pairwise memberships
  pairwise_comm1 <- crossprod(table(stack(comm1)))
  pairwise_comm2 <- crossprod(table(stack(comm2)))
  sum(pairwise_comm1 * pairwise_comm2) / sqrt(sum(pairwise_comm1^2) * sum(pairwise_comm2^2))
}

# Precision/Recall for overlaps
f1_score <- function(ground_truth, detected) {
  totF1 <- 0
  for( i in dim(ground_truth)[2]){
    true <- ground_truth[,i]
    pred <- detected[,i]
    precision <- sum(true & pred) / sum(pred)
    recall <- sum(true & pred) / sum(true)
    new <- 2 * (precision * recall) / (precision + recall)
    print(new)
    totF1 <- totF1 + 2 * (precision * recall) / (precision + recall)  
  }
  return(totF1)
}

permanence <- function(graph, comm, node){
  internal <- sum(comm == comm[graph[[node]]])
  total <- degree(graph, node)
  internal/total - (max(table(comm[graph[[node]]]))/total)
}

CQ <- function(graph, communities) {
  sum(sapply(communities, function(comm) {
    subg <- induced_subgraph(graph, unlist(comm))
    m_int <- ecount(subg)
    m_ext <- sum(igraph::degree(graph, unlist(comm))) - 2 * m_int
    (m_int - m_ext) / length(unlist(comm))
  }))
}

ONCut <- function(graph, communities) {
  sum(sapply(communities, function(comm) {
    comm <- unlist(comm)
    boundary <- sum(!(neighbors(graph, comm) %in% comm))
    vol <- sum(igraph::degree(graph, comm))
    boundary / vol
  }))
}

OP <- function(communities, g) {
  sum(sapply(1:vcount(g), function(v) {
    sum(sapply(communities, function(c) v %in% c)) > 1
  })) / vcount(g)
}

# Find clustering coef global for directed and undirected network
# High clustering coef means more clique based highly clustered Network
# Low means more sparse network
ClustCoef <- function(G, dir){
  A <- as_adjacency_matrix(G, sparse = FALSE)
  if(dir == "directed"){
    cc <- ClustF(A, type = "directed")
    return(cc$GlobaltotalCC)
  }else if(dir == "undirected"){
    cc <- ClustF(A, type = "undirected")
    return(cc$GlobalCC)
  }
  
}

## Power law fit check
PLfitcheck <- function(G){
  degD <- igraph::degree(G)
  degD <- degD[degD != 0]
  fit <- displ$new(degD)
  est <- estimate_xmin(fit)
  fit$setXmin(est)
  bs <- bootstrap_p(fit, threads = 4)  # KS test
  return(bs$p)  # p-value > 0.1 supports power law
}

# 1: Perfect homophily (nodes only connect within their group).
# 0: No homophily (random mixing).
# < 0: Heterophily (nodes prefer different groups).
assortativityF <- function(G,dir){
  return(igraph::assortativity_degree(G, directed = dir))
}

## Feature structure decomposition
Feature_struct_decomp <- function(G, nc, N, delta, noCov, FullM, covNames){
  
  ## Baseline 1 : CD using NW structure. No covariate community detection
  ## noCov
  
  ## Baseline 2: CD using NW only features (k-means)
  ## Get the covariates
  X <-  as.data.frame(vertex_attr(G)) %>% dplyr::select(covNames)
  ### using cmeans function from e1071 library
  Cov_op <- cmeans(X , nc)
  ##Find the overlap 
  cov <- memOverlapCalc(Cov_op$membership,Cov_op$membership,N, delta,nc)
  
  ## Full model: Structure + covariate community detection
  FullM <- FullM
  
  ##Compare B1 to Full model
  oi_b1 <- OmegaIdx_(noCov, FullM,nc,N)
  
  ##Compare B2 to Full model
  oi_b2 <- OmegaIdx_(cov, FullM,nc,N)
  
  ## Compare covariate only to ground truth
  OrigVal <-  as.data.frame(vertex_attr(G)) %>%
    dplyr::select(any_of(c(make_letter_names(nc) ))) %>%
    abs()
  oi_cov <- OmegaIdx_(OrigVal, cov, nc, N)
  
  return(list(oi_b1, oi_b2, oi_cov))
  
}

## null models: Randomize feature labels -> if performance drops features matter
null_models <- function(G, nc, k, o, N, alpha,lambda_lin, lambda_bin, thresh, nitermax, orig, randomize,
                        CovNamesLinin, CovNamesLinout, CovNamesLPin, CovNamesLPout, dir,
                        alphaLL, test, missing, covOrig,epsilon, impType, alphaLin, 
                        penalty, seed, covInit, specOP){
  ## for example for continuous feature Z[,i]
  ## V(G)$feature <- sample(V(G)$feature)
  ## fit the method again and compare the OI of original data and shuffled data with ground truth.
  
  ## Randomize the values of a feature vector
  ## null models: Randomize feature labels -> if performance drops features matter
  covtmp <-  as.data.frame(vertex_attr(G))
  if (sum(o) > 0) {
    Z_out <- as.matrix(covtmp %>% dplyr::select(all_of(CovNamesLinout)))
    for(i in 1:o[2]){
      Z_out[,i] <- sample(Z_out[,i])
      G <- G %>% set_vertex_attr(name = CovNamesLinout[i], value = c(Z_out[, i]))
    }
  }
  tryCatch({
    ## Algorithm with covariates + possible missing data
    start <- Sys.time()
    opf_RegcovRand <- CoDA(G, nc, k, o, N, alpha, lambda_lin, lambda_bin, thresh, nitermax, orig, randomize,
                           CovNamesLinin, CovNamesLinout, CovNamesLPin, CovNamesLPout, dir,
                           alphaLL, test = FALSE, missing, covOrig,epsilon, impType, alphaLin, 
                           penalty, seed, covInit, specOP )
    randOI <- tail(opf_RegcovRand$OmegaIndex,1)
    tme <- Sys.time() - start
    print(paste0("Total time take algo with scrambled covariates.", round(tme, 3)))
  }, error = function(e) {
    # Handle the error
    cat("An error occurred:", conditionMessage(e), "\n")
    NA})
  return(randOI)
}

## SURPRISE METRIC 
# Function to compute the Surprise metric NON OVERLAPPING ONLY!
Surprise <- function(graph, membership) {
  # Total number of nodes and edges
  n <- vcount(graph)
  m <- ecount(graph)
  
  # Total possible edges in the network (undirected)
  M <- n * (n - 1) / 2
  
  # Number of intra-community edges (m_in)
  m_in <- 0
  communities <- unique(membership)
  
  for (comm in communities) {
    nodes_in_comm <- which(membership == comm)
    if (length(nodes_in_comm) >= 2) {
      subgraph <- induced_subgraph(graph, nodes_in_comm)
      m_in <- m_in + ecount(subgraph)
    }
  }
  
  # Maximum possible intra-community edges (q)
  q <- 0
  for (comm in communities) {
    size <- sum(membership == comm)
    if (size >= 2) {
      q <- q + (size * (size - 1) / 2)
    }
  }
  
  # Compute Surprise using the hypergeometric distribution
  # Using phyper for cumulative probability (lower tail = FALSE for upper tail)
  p_value <- phyper(m_in - 1, q, M - q, m, lower.tail = FALSE)
  surprise <- -log10(p_value)
  
  return(list(
    m_in = m_in,
    m_out = m - m_in,
    q = q,
    M = M,
    p_value = p_value,
    surprise = surprise
  ))
}


