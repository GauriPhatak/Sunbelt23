
library(igraph)
library(tidyverse)
library(ITNr)


## Sample atleast 2 times.
## This way there are atleast 2 nodes in each group
f <- function(x, n){
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


##SBM function.
## Add code to give the output for \omega matrix (m_rs/2*m)
lazySBM <- function(g,n,type, thresh,maxIter){
  grps <-  1:n
  l <-  length(V(g))
  prevll <- NA
  k <- 0

  if(2*n > l){
    return(print("There cannot be less than two nodes per group."))
  }

  while( is.na(prevll) || is.infinite(prevll)){
    randMem <- f(grps, l)
    V(g)$group <- randMem
    if(type == 1){
    prevll <- llBasicsbm(g)
    }
    else{
      prevll <- llDCsbm(g)
    }
  }
  origV <- V(g)$group
  
  iter  = 0
  gtmp <-  g
  ## Dataframe to save data for each iteration
  columns <- c("iterNum", "MaxVal", "VertChanged", "PrevGrp", "Newgrp","LTprev")
  opdf <-  data.frame(matrix(nrow = 0, ncol = length(columns))) 
  colnames(opdf) <- columns
  verts <- V(g)
  repeat{
    Mat <- matrix(rep(-Inf, n*l), nrow = l )
    
    for(i in verts){
      for(j in grps[V(g)$group[i] != grps]){
        
        gtmp <- g
        V(gtmp)$group[i] <-  j
        
        if(type == 1){
          Mat[i,j] <- llBasicsbm(gtmp)
        }
        else{
          Mat[i,j] <- llDCsbm(gtmp)
        }
        
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
    if( iter < maxIter){
      if((maxVal > prevll)) {
        V(g)$group[ind[1,1]] <- ind[1,2]
        prevll <-  maxVal
      }
      else if(maxVal < prevll ){
        
        #if( (prevll - maxVal) < thresh){
        if( (prevll - maxVal) < (prevll*thresh)){
          
          V(g)$group[ind[1,1]] <- ind[1,2]
          prevll <-  maxVal
          
        }
        else{
          
          return(list(g , Mat, opdf, origV))
        }
      }
    }
    else{

      return(list(g , Mat, opdf, origV))
    }
   
  }
  
  return(list(g , Mat, opdf, origV))
}

## calculate the ARI?

DataSetup <- function(v,n,pm, bs, type, thresh, maxIter){
  
  g <- sample_sbm(v, 
                  pref.matrix=pm, 
                  block.sizes=bs, 
                  directed = FALSE)
  
  
  ## Using edge betweenness to compare the log likelihood.
  clu <- cluster_edge_betweenness(g)
  g$group <- clu$membership
  
  relations <- data.frame(get.edgelist(g))
  g <- graph_from_data_frame(relations, 
                             directed=FALSE, 
                             vertices=as.character(1:v))
  if(type == 1){
  op <- lazySBM(g, n, type, thresh, maxIter)
  }
  else if(type == 2){
  op <- lazySBM(g, n, type, thresh, maxIter)
  }
  else {
    print("Incorrect type entered.")
  }
  V(g)$group <- clu$membership
  
  return(list(op, clu$membership))
  
}



v <- 80 #readline(prompt = "Enter the number of vertices: ")
n  <- 3 #readline(prompt = "Enter the number of groups: ")
pm <- matrix(c(0.8, 0.01, 0.01,0.01, 0.7, 0.01,0.01,0.01,0.6), nrow = n)
bs <- c(35,25,20)
thresh = 0.05
maxIter  = 100

## Basic SBM
op <- DataSetup(v,n,pm, bs, 1, thresh, maxIter)
opDF <- op[[1]][[3]]
plot(op[[1]][[1]], vertex.color = op[[1]][[4]], vertex.size = 8, vertex.label = NA)
plot(op[[1]][[1]], vertex.color = V(op[[1]][[1]])$group, vertex.size = 8, vertex.label = NA)
plot(op[[1]][[1]], vertex.color = op[[2]], vertex.size = 8, vertex.label = NA)
llBasicsbm(op[[1]][[1]])
llBasicsbm(op[[2]])

## DC SBM
op <- DataSetup(v,n,pm, bs, 2, thresh, maxIter)
opDFDC <- op[[1]][[3]]
plot(op[[1]][[1]], vertex.color = op[[1]][[4]], vertex.size = 8, vertex.label = NA)
plot(op[[1]][[1]], vertex.color = V(op[[1]][[1]])$group, vertex.size = 8, vertex.label = NA)
plot(op[[1]][[1]], vertex.color = op[[2]], vertex.size = 8, vertex.label = NA)
llDCsbm(op[[1]][[1]])
llDCsbm(op[[2]])
