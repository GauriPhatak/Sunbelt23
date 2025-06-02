library(igraph)
library(aricode)  # For NMI
library(linkprediction)  # For partition density (install via devtools::install_github("arc85/LinkPrediction"))
library(stringr)
library(tidyverse)
## Calculate ARI for the network outputs
ARIop <- function(Ftot,Htot,orig,nc,N){
  mem <- rep(NA, N)
  
  for (i in 1:N) {
    m <- data.frame(com = rep(letters[1:nc], times = 2),
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
    pred <- cbind(1, Wtmat1) %*% beta
  }
  mse <- sqrt(colSums((pred - Z) ^ 2) / N)
  return(mse)
}

MSEop <- function(Ftot, Htot, covOrig, betaout,N, dir, o_in,o_out, missValsout){
  mseouttot <- matrix(nrow = 0, ncol = (o_in+o_out))
  mse_outMD <- c()#matrix(nrow = 0, ncol = (o_in + o_out))
  
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

getDelta <- function(N){
  delta <- sqrt(-1*log(1-(1/N)))
  return(delta)
}

memOverlapCalc <- function(Fm, Hm, delta, N, nc){
  ol <- as.data.frame(matrix(0, ncol = nc, nrow = N))
             
  for (i in 1:N) {
    ol[i, ] <- (as.numeric(Fm[i, ] > delta) + as.numeric(Hm[i, ] > delta)) > 0
  }
  return(ol)
}

OmegaIdx <- function(G, Fm, Hm, N, delta, nc) {
  if(printFlg == TRUE){
    print("In OmegaIdx Func")
  }
  memoverlap <- memOverlapCalc(Fm, Hm, delta, N, nc)
  
  
  #### Using soft thresholding for memory overlap
  #ol <- lapply(1:nc, function(i) which(memoverlap[,i] > quantile(memoverlap[,i], 0.75)))
  
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

# reading the simulated networks and the original assignments
# param <- "InitParamMiss_Coh_MCAR_LASSO_Cont_2.rds"#"InitParamMiss_Cohesive_MCAR_B15_LASSO.rds"
# param_c <- "InitParamMiss_Coh_MCAR_LASSO_Cont_2_comp.rds"#"InitParamMiss_Cohesive_MCAR_B15_LASSO_comb.rds"
# getwd()
# fileName <- list.files("Code/CoDAOP/CohesiveMCAR_2_NoDiag/")#"CoDAOP/CohesiveMCAR_B15_LASSO/")
# #k <- 0
# op_sim <- matrix(0, nrow =0 , ncol = 11)
# op_orig <- matrix(0, nrow =0 , ncol = 9)
# type <- c("Reg","Nocov","Sto","Mean")
# dir <- c("D","U")
# for(file in fileName){
#   
#   print(file)
#   idx <- str_extract_all(file, "\\d+")[[1]][2]
#   df <- readRDS(paste0("Code/CoDAOP/CohesiveMCAR_2_NoDiag/",file))
#   #df_reg <- df[[1]]
#   #df_sto <- df[[3]]
#   #df_nocov <- df[[2]]
#   grp <- 1
#   pb = txtProgressBar(min = 0, max = length(seq(from = 1, to = 2*length(df[[4]]), by = 2)), initial = 0)
#   
#   for(i in seq(from = 1, to = 2*length(df[[5]]), by = 2)){
#     setTxtProgressBar(pb,i)
#     
#     G <- df[[5]][[grp]]
#     N <- vcount(G)
#     delta <- getDelta(N)
#     nc <- 3
#     OrigVal <-  as.data.frame(vertex_attr(G)) %>%
#       dplyr::select(all_of(c(letters[1:nc])))
#     OrigVal[OrigVal == -1] <- 1
#     
#     for(j in 1:4){
#       FinOP <- list()
#       for(k in 1:2){
#         Fm <- df[[j]][[i+k-1]]$Ffin
#         Hm <- df[[j]][[i+k-1]]$Hfin
#         FinOP[[k]] <- as.data.frame(memOverlapCalc(Fm, Hm, delta, N, nc))
#         comnty <- data.frame(which(as.matrix(FinOP[[k]]) == 1, arr.ind = TRUE)) %>%
#           group_by(col) %>%
#           group_map(~.x)
#         op_sim <- rbind(op_sim, c(
#           round(conductance(graph = G, communities = comnty, FinOP[[k]]),4),
#           round(partition_density(graph = G, communities = comnty),4),
#           round(Q_HardPartition(G, FinOP[[k]]),4),
#           round(Q_Overlapping(G, FinOP[[k]], dir[k]),4),
#           #round(Suprise(G, ),4),
#           #omega_index(communities[1], communities[2]),
#           #NMI(factor(OrigVal), factor(FinOP[[k]])),
#           #f1_score(ground_truth =  OrigVal, detected =  FinOP[[k]]),
#           #permanence(g, communities[[1]], 15),
#           round(CQ(G, comnty),4),
#           round(ONCut(G, comnty),4),
#           round(OP(comnty, G),4),
#           dir[k],
#           type[j],
#           idx,
#           grp))
#       }
#     }
#     
#     comnty <- data.frame(which(as.matrix(OrigVal) == 1, arr.ind = TRUE)) %>%
#       group_by(col) %>%
#       group_map(~.x)
#     op_orig <- rbind(op_orig, c(
#       round(conductance(graph = G, communities = comnty, OrigVal),4),
#       round(partition_density(graph = G, communities = comnty),4),
#       round(Q_HardPartition(G, OrigVal),4),
#       round(Q_Overlapping(G, OrigVal, dir[1]),4),
#       #round(Suprise(G, ),4),
#       #omega_index(communities[1], communities[2]),
#       #NMI(factor(OrigVal), factor(FinOP[[k]])),
#       #f1_score(ground_truth =  OrigVal, detected =  FinOP[[k]]),
#       #permanence(g, communities[[1]], 15),
#       round(CQ(G, comnty),4),
#       round(ONCut(G, comnty),4),
#       round(OP(comnty, G),4),
#       idx,
#       grp))
#     grp <- grp + 1
#   }
#   
#   close(pb)
# }
# 
# colnames(op_sim) <- c("Conductance", "Partition Density","modHP", "Qov",  "Cover Quality", "ONCut", "Overlap Proportion", "Dir","type","idx", "grp")
# colnames(op_orig) <- c("Conductance", "Partition Density","modHP", "Qov",  "Cover Quality", "ONCut", "Overlap Proportion", "idx", "grp")
# 
# saveRDS(op_sim, "MetricsDataSim.rds")
# saveRDS(op_orig, "MetricsDataOrig.rds")

#gbg <- readRDS("MetricsDataSim.rds")

