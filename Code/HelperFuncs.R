## This file to keep all the functions used in testing.
library(igraph)
library(tidyverse)
library(ITNr)
library(data.table)
library(expm)
library(wordspace)
library(stats)

f <- function(x, n){
  #set.seed(42)
  sample(c(x,x,sample(x, n-2*length(x), replace=TRUE)))
}

## Simple SBM log likelihood
llBasicsbm <- function(g){
  
  mm <- mixing_matrix_igraph(g,"group")
  ct <- as.matrix(table(data.frame(V(g)$group)))
  llval <- sum(colSums(mm *log(mm /(ct %*% t(ct)))))
  
  return(llval)
}

## Degree corrected SBM log likelihood
llDCsbm <- function(g){
  
  d <- as.numeric(degree(g))
  d <- data.frame(cbind(deg = d,grp = as.numeric(V(g)$group)))
  d <- d %>% 
    group_by(grp) %>% 
    summarise(d = sum(deg),
              n = n())
  mm <- mixing_matrix_igraph(g,"group")
  deg <-  matrix(d$d)
  llval <- sum(colSums(mm *log(mm /(deg %*% t(deg)))))
  
  return(llval)
}


SBMrepeat <-  function(g,l,n,prevll,maxIter, thresh,FUN){
  print("inside the SBMrepeat function")
  grps <-  1:n
  iter  = 0
  gtmp <-  g
  ## Dataframe to save data for each iteration
  columns <- c("iterNum", "MaxVal", "VertChanged", "PrevGrp", "Newgrp","LTprev") # need to add neighbours here.
  opdf <-  data.frame(matrix(nrow = 0, ncol = length(columns))) 
  
  ## Get the degrees of the different nodes.
  deg <- degree(g)
  degrees <- list()
  colnames(opdf) <- columns
  verts <- V(g)
  repeat{
    Mat <- matrix(rep(-Inf, n*l), nrow = l )
    
    for(i in verts){
      for(j in grps[V(g)$group[i] != grps]){
        
        gtmp <- g
        V(gtmp)$group[i] <-  j
        
        Mat[i,j] <-  FUN(gtmp)
       
      }
    }
    
    maxVal <- max(Mat, na.rm = TRUE)
    ind <- which(Mat == maxVal, arr.ind = TRUE)
    pg <- which(Mat[ind[1,1],] == -Inf, arr.ind = TRUE)
    iter  <-  iter +1
    k <- (prevll - maxVal)
    opdf <- opdf %>% 
      add_row(iterNum = iter, 
              MaxVal= maxVal, 
              VertChanged = ind[1,1], 
              PrevGrp = pg, 
              Newgrp = ind[1,2],
              LTprev = k)
    a <- as.data.frame(cbind(degree(gtmp,as.numeric(as_ids(neighbors(gtmp,ind[1,1])))),
                             V(gtmp)$group[as.numeric(as_ids(neighbors(gtmp,ind[1,1])))]))
    a <- subset(a,a[,2] == as.numeric(pg))
    b <- as.list(a$V1)
    names(b) <-  rownames(a)
    degrees <- c(degrees, list(b))
    if( iter < maxIter){
      if((maxVal > prevll)) {
        V(g)$group[ind[1,1]] <- ind[1,2]
        prevll <-  maxVal
      }
      else if(maxVal < prevll ){
        if( (prevll - maxVal) < (prevll*thresh*(-1))){
          
          V(g)$group[ind[1,1]] <- ind[1,2]
          prevll <-  maxVal
          
        }
        else{
          
          opdf <- data.table(opdf)
          opdf$degrees <- degrees
          return(list(FinAsgn  = V(g)$group ,FinalLLMat =  Mat, Sequence = opdf))
        }
      }
    }
    else{
      #print(V(g)$group)
      #flush.console()
      opdf <- data.table(opdf)
      opdf$degrees <- degrees
      return(list(FinAsgn  = V(g)$group ,FinalLLMat =  Mat, Sequence = opdf))
      
    }
    
  }
}

##SBM function.
## Add code to give the output for \omega matrix (m_rs/2*m)
lazySBM <- function(g,n,type, thresh,maxIter){
  print("inside the lazysbm function")
  
  ## setting vriables for the group sizes and total nodes
  l <-  length(V(g))
  prevll <- NA
  k <- 0
  
  #3 checking if the number of groups is larger then the algorithm can handle.
  if(2*n > l){
    return(print("There cannot be less than two nodes per group."))
  }
  
  ## Check if the random assignment assigns groups with less than 
  ## two nodes or if there is less than one between group edge assignment.
  while( is.na(prevll) || is.infinite(prevll)){
    randMem <- f(1:n, l)
    V(g)$group <- randMem
    if(type == "simple"){
      prevll <- llBasicsbm(g)
    }
    else if(type == "DC"){
      prevll <- llDCsbm(g)
    }
  }
  origV <- V(g)$group
  if(type == "simple"){
    op <- SBMrepeat(g,l,n, prevll,maxIter,thresh ,FUN = llBasicsbm)
  }
  else if(type  == "DC"){
    op <- SBMrepeat(g,l,n, prevll,maxIter, thresh ,FUN = llDCsbm)
  }
  return(list(FromSBMRep = op, OrigRand = origV))
}

## calculate the ARI?

DataSetup <- function(v,n,pm, bs, type, thresh, maxIter, nRandStarts = 1){
  print("inside the DATASETUP function")
  
  g <- sample_sbm(v, 
                  pref.matrix=pm, 
                  block.sizes=bs, 
                  directed = FALSE)
  
  print("Generated sample graph")
  
  ## Using edge betweenness to compare the log likelihood.
  clu <- cluster_edge_betweenness(g)
  grpOrig <-  rep(1:n, times = bs)

  relations <- data.frame(get.edgelist(g))
  g <- graph_from_data_frame(relations, 
                             directed=FALSE, 
                             vertices=as.character(1:v))
  op <- list()
  for(nr in 1:nRandStarts){
  if(type == "simple"){
    op <- append(op, lazySBM(g, n, "simple", thresh, maxIter))
  }
  else if(type == "DC"){
    op <- append(op , lazySBM(g, n, "DC", thresh, maxIter))
  }
  else if(type == "both"){
    op1 <- lazySBM(g, n, "simple", thresh, maxIter)
    op2 <- lazySBM(g, n, "DC", thresh, maxIter)
    op <- append(op, list(simple = op1,dc = op2))
  }
  else {
    print("Incorrect type entered.")
  }
    
  }
  return(list(genGraph = g, FrmLazySBM = op, FrmSampleSBM = grpOrig))
  
}

RegSpectralClust <-  function(G, k, regularize = TRUE ){
  
  
## Sample graph for debugging
  # k =3
  # G <- sample_sbm(n = 90, 
  #                 pref.matrix=matrix(c(0.8,0.01,0.02,0.01,0.7,0.01,0.02,0.01,0.8),
  #                                    nrow = 3),
  #                 block.sizes=90*c(2/5,2/5,1/5), 
  #                 directed = FALSE)
  # plot(G)
  # 
  
# Step 1:
  
## Find Adjacency matrix for the given graph.
  
  A <- as_adjacency_matrix(G)
  
## find the Degree Diagonal matrix. regularized diagonal matrix here tau  =  average node degree.
  I <- diag(x = 1, nrow = length(V(G)), ncol = length(V(G)))
  tau = mean(degree(G))
  if(regularize == FALSE){tau = 0}
  Dt <- diag(degree(G)+tau, nrow = length(V(G)))
  
## find the -1/2 root of Matrix Dt

## Calculate the Regularized graph laplacian
  
  Lt <-  solve(sqrtm(Dt)) %*% A %*% solve(sqrtm(Dt))
  

# Step 2:
  ## Fnd the eigen values and vectors of the Lt matrix
  
  X <-  eigen(Lt)
  
  ## print the difference between the eigen values 
  print(X$values - lead(X$values))
  flush.console()

# Step 3:
  
  #3 Normalize each row of the matrix to have unit length
  Xt <- normalize.rows(X$vectors[,1:k]) 
  
# Step 4 
  #Run kK- means algorithm on the subset of the eigen vector matrix of size nXk
  
  op <- kmeans(Xt,k,iter.max = 30)
  
  return(op$cluster)
  
}


CovAssistedSpecClust <- function(G, X, k, Regularize = TRUE){
  
  # Step 1:
  
  ## Find Adjacency matrix for the given graph.
  
  A <- as_adjacency_matrix(G)
  
  ## find the Degree Diagonal matrix. regularized diagonal matrix here tau  =  average node degree.
  I <- diag(x = 1, nrow = length(V(G)), ncol = length(V(G)))
  tau = mean(degree(G))
  if(regularize == FALSE){tau = 0}
  Dt <- diag(degree(G)+tau, nrow = length(V(G)))
  
  ## find the -1/2 root of Matrix Dt
  
  ## Calculate the Regularized graph laplacian
  
  Lt <-  solve(sqrtm(Dt)) %*% A %*% solve(sqrtm(Dt))
  
  
  # Step 2:
  ## Fnd the eigen values and vectors of the Lt matrix
  
  X <-  eigen(Lt)
  
  ## print the difference between the eigen values 
  print(X$values - lead(X$values))
  flush.console()
  
  # Step 3:
  
  #3 Normalize each row of the matrix to have unit length
  Xt <- normalize.rows(X$vectors[,1:k]) 
  
  # Step 4 
  #Run kK- means algorithm on the subset of the eigen vector matrix of size nXk
  
  op <- kmeans(Xt,k,iter.max = 30)
  
  return(op$cluster)
}















