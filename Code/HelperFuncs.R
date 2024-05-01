## This file to keep all the functions used in testing.
library(igraph)
library(tidyverse)
library(ITNr)
library(data.table)
library(expm)
library(wordspace)
library(stats)
library(tidyverse)
library(ergm)
library(fs)
library(network)
library(intergraph)
#library(tidygraph)
#library(ggraph)
#library(statnet)
library(nnet)
library(caret)
#library(ggnet)
library(stringr)
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
    dplyr::summarise(d = sum(deg),
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
    
    #print("Inside the repeat loop")
    Mat <- matrix(rep(-Inf, n*l), nrow = l )
    
    for(i in verts){
      for(j in grps[V(g)$group[i] != grps]){
        
        gtmp <- g
        V(gtmp)$group[i] <-  j
        
        Mat[i,j] <-  FUN(gtmp)
       
      }
    }
    #print("test1")
    maxVal <- max(Mat, na.rm = TRUE)
    ind <- which(Mat == maxVal, arr.ind = TRUE)
    pg <- which(Mat[ind[1,1],] == -Inf, arr.ind = TRUE)
    iter  <-  iter +1
    k <- (prevll - maxVal)
    #print("test2")
    opdf <- opdf %>% 
      add_row(iterNum = iter, 
              MaxVal= maxVal, 
              VertChanged = ind[1,1], 
              PrevGrp = pg, 
              Newgrp = ind[1,2],
              LTprev = k)
    a <- as.data.frame(cbind(degree(gtmp,as.numeric(as_ids(neighbors(gtmp,ind[1,1])))),
                             V(gtmp)$group[as.numeric(as_ids(neighbors(gtmp,ind[1,1])))]))
    #print("test3")
    
    a <- subset(a,a[,2] == as.numeric(pg))
    b <- as.list(a$V1)
    names(b) <-  rownames(a)
    #print("test4")
    degrees <- c(degrees, list(b))
    if( iter < maxIter){
      #print(iter)
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
  #print(l)
  prevll <- NA
  k <- 0
  
  #3 checking if the number of groups is larger then the algorithm can handle.
  if((2*n) > l){
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
  
# Step 1:
  
## Find Adjacency matrix for the given graph.
  
  A <- as_adjacency_matrix(G)
  
## find the Degree Diagonal matrix. regularized diagonal matrix here tau  =  average node degree.
  I <- diag(x = 1, nrow = length(V(G)), ncol = length(V(G)))
  tau = mean(igraph::degree(G))
  if(regularize == FALSE){tau = 0}
  Dt <- diag(igraph::degree(G)+tau+0.5, nrow = length(V(G)))
  
## find the -1/2 root of Matrix Dt

## Calculate the Regularized graph laplacian
  
  Lt <-  solve(sqrtm(Dt)) %*% A %*% solve(sqrtm(Dt))
  

# Step 2:
  ## Fnd the eigen values and vectors of the Lt matrix
  
  X <-  eigen(Lt)
  
  ## print the difference between the eigen values 
 # print(X$values - lead(X$values))
 # flush.console()

# Step 3:
  
  #3 Normalize each row of the matrix to have unit length
  Xt <- normalize.rows(Re(X$vectors[,1:k])) 
  
# Step 4 
  #Run kK- means algorithm on the subset of the eigen vector matrix of size nXk
  
  op <- kmeans(x = Xt,centers = k ,nstart = 10, iter.max = 100)
  
  return(op$cluster)
  
}

findAlpha <- function(G, X, Regularize=TRUE){
## Find Adjacency matrix for the given graph.
  A <- as.matrix(as_adjacency_matrix(G))
  
  ## find the Degree Diagonal matrix. regularized diagonal matrix here tau  =  average node degree.
  I <- diag(x = 1, nrow = length(V(G)), ncol = length(V(G)))
  tau = mean(igraph::degree(G))
  if(Regularize == FALSE){tau = 0}
  Dt <- diag(igraph::degree(G)+tau, nrow = length(V(G)))
  
  ## find the -1/2 root of Matrix Dt
  
  ## Calculate the Regularized graph laplacian
  Lt <-  solve(sqrtm(Dt)) %*% A %*% solve(sqrtm(Dt))
  
  maxs <- apply(X, 2, max)
  mins <- apply(X, 2, min)
  X <- scale(X, center = mins, scale = maxs - mins)
  
  ## Choosing alpha based on leading eigen value of both L_t %*% L_t and X %*% X_t
  return(eigen(Lt %*% Lt)$values[1]/ eigen(X %*% t(X))$values[1])

}

CovAssistedSpecClust <- function(G, X, k, alpha, Regularize = TRUE, type = "assortative", kmeansIter = 100, retAlpha =FALSE){
  
  # Step 1:
  
  ## Find Adjacency matrix for the given graph.
  A <- as.matrix(as_adjacency_matrix(G))
  
  ## find the Degree Diagonal matrix. regularized diagonal matrix here tau  =  average node degree.
  I <- diag(x = 1, nrow = length(V(G)), ncol = length(V(G)))
  tau = mean(igraph::degree(G))
  if(Regularize == FALSE){tau = 0}
  Dt <- diag(igraph::degree(G)+tau, nrow = length(V(G)))
  
  ## find the -1/2 root of Matrix Dt
  
  ## Calculate the Regularized graph laplacian
  Lt <-  solve(sqrtm(Dt)) %*% A %*% solve(sqrtm(Dt))
  
  maxs <- apply(X, 2, max)
  mins <- apply(X, 2, min)
  X <- scale(X, center = mins, scale = maxs - mins)
  
  ## Choosing alpha based on leading eigen value of both L_t %*% L_t and X %*% X_t
  if(is.na(alpha)){
    alpha <- eigen(Lt %*% Lt)$values[1]/ eigen(X %*% t(X))$values[1]
    }
  
  #print(alpha)
  if(type == "non-assortative"){
    L_alpha <-  Lt %*% Lt + alpha *(X %*% t(X))
  }
  else if(type == "CCA"){
    L_alpha <- Lt %*% X
  }
  else if(type == "assortative"){
    L_alpha <-  Lt + alpha *(X %*% t(X))
  }
  else{
    return("No type Specified!")
  }
  
  
  
  # Step 2:
  ## Fnd the eigen values and vectors of the Lt matrix
  
  X_L <-  eigen(L_alpha)
  
  ## print the difference between the eigen values 
  #print(X$values - lead(X$values))
  #flush.console()
  
  # Step 3:
  
  #3 Normalize each row of the matrix to have unit length
  Xt <- normalize.rows(X_L$vectors[,1:k]) 
  
  # Step 4 
  #Run kK- means algorithm on the subset of the eigen vector matrix of size nXk
  
  op <- kmeans(Xt,k,iter.max = kmeansIter, nstart = 10)
  
  if(retAlpha == TRUE){
    return(list(op$cluster, alpha))
  }
  return(op$cluster)
}

findEdges <- function(geom, hwd, id, type){
  
  city_cntr <- ct %>% 
    st_as_sf() %>% 
    st_centroid() %>% 
    st_cast("POINT") %>% 
    st_transform(4326) %>% 
    st_coordinates()
  
  distMat <- distm(city_cntr)
  colnames(distMat) <- ct$city_name
  rownames(distMat) <- ct$city_name
  
  edges <- data.frame(matrix(ncol = 4, nrow=0))
  
  for(hwy_id in id) {
    dist<- c()
    sbst <- hwd[hwd$st_hwy_idx == hwy_id,]
    int <- st_intersection(sbst, geom)
    sortedV <- int %>% 
      group_by(city_name) %>%
      dplyr::summarise(m = median(beg_mp_no)) %>%
      arrange(m)
    #Subset of cities in the data
    #ct$city_name[ct$city_name %in% unique(int$city_name)]
    sortedV$city_nameConn <- c(sortedV$city_name[-1], 0)
    
    sortedCT <- as.matrix(cbind(sortedV$city_name, sortedV$city_nameConn), ncol =2 )#[1: length(sortedV$city_name)-1,]
    #print(hwy_id)
    #print(sortedCT)
    
    for (i in 1:dim(sortedV)[1]-1) {
      dist <- c(dist, distMat[sortedCT[i,1], sortedCT[i,2]])
    }
    ##Cannot use the milepost numbers as distance metric... it's different for different highways
    #dist <- c(sortedV$m[-1],0) - sortedV$m
    edges <- rbind(edges,cbind(sortedV$city_name, sortedV$city_nameConn,
                               rep(hwy_id,length(sortedV$city_name)),
                               dist)[1:length(sortedV$city_name)-1,])
  }
  edges <- cbind(edges, type)
  colnames(edges) <- c("to","from","hwy_id","distance","type")
  return(edges)
}


ShortestPathGrph <- function(ct_sbst, hwy_grph, city_name, simple = TRUE){
  
  ## keeping all the edges from above. Checking if this would add more edges to the data
  saveDirEdges <-  data.frame(matrix(ncol = 3, nrow = 0))
  
  distV <-  distances(hwy_grph, ct_sbst, 
                      ct_sbst,
                      weights = round(as.numeric(E(hwy_grph)$distance)))
  
  for (i in 1:length(ct_sbst)) {
    
    trav <- list()
    path <- shortest_paths(hwy_grph, ct_sbst[i], 
                           ct_sbst, output = "both" , 
                           weights = round(as.numeric(E(hwy_grph)$distance)))
     
    for(j in 1:length(ct_sbst)){
      curPath <- names(unlist(path$vpath[[j]]))
      subPath <- curPath[curPath %in% ct_sbst]
      if(length(subPath) > 1){
        newEdges <- cbind(subPath[1:length(subPath)-1], subPath[2:length(subPath)])
        newEdges <- cbind(newEdges, as.numeric(apply(newEdges , 1, function(x) distV[x[1], x[2]])))
        saveDirEdges <- rbind(saveDirEdges, newEdges )
      }
    }
  }
  colnames(saveDirEdges) <- c("to","from", "distance")
  saveDirEdges$distance <- as.numeric(saveDirEdges$distance)
  
  ## Creating edges from the paths that we generate using shortest path algorithm. 
  ## Coloring based on the weight of the edge i.e. the number of nodes between the two nodes of interest.
  
  DirectGrph <- graph_from_data_frame(d = saveDirEdges, vertices = city_name, directed = FALSE)
  
  if(simple == TRUE){
    DirectGrph <- simplify(DirectGrph, edge.attr.comb = "min")
  }
  
  return(DirectGrph)
}

plot_graph <- function(g, m, title){
  plot.igraph(g, vertex.size = 4, vertex.label = NA,
              vertex.color = m,
              edge.width = 1, layout = cbind(as.numeric(V(g)$lat), as.numeric(V(g)$lon)), 
              main = title,
              margin = -0.01)
}

## for all categories at random.MAR
## Ties randomly removed
removeRandEdges <- function(graph, percent){
  
  N <- gsize(graph)
  numE <- floor((N * percent)/100)
  to_be_removed <- get.edgelist(graph)[sample(1:N, numE),] #sample(1:N, numE, replace = FALSE)
  return(list("graph" =  delete_edges(graph , to_be_removed),"ids"= to_be_removed))
}

#gbg <- removeRandEdges(cowNW, 10)
#G <- gbg$graph
#ids <- gbg$ids

## For all categories based on a percentage for each categories in it. MNAR

removeCatedges <- function(graph, attr, fxn, percent, category){
  
  edges <- unlist(incident_edges(graph, which(fxn(vertex_attr(graph, attr), category))), recursive = FALSE, use.names = FALSE)
  N <- length(edges)
  numE <- ceiling((N*percent)/100)
  to_be_removed <- get.edgelist(graph, names = FALSE)[sample(edges, numE, replace = FALSE),]
  to_be_removed <- paste(to_be_removed[,1], to_be_removed[,2], sep = "|")
  print(length(to_be_removed))
  ret_grph <- delete_edges(graph , to_be_removed)
  return(list("graph" = ret_grph, "ids"= gsize(ret_grph)))
}

#gbg <- removeCatedges(cowNW, "gender", `==`,80, "2")
#G <- gbg$graph
#ids <- gbg$ids

## remove values in the attributes at random

removeRandAtt <- function(attrs, percent){
  
  N <- length(attrs)
  numE <- floor((N*percent)/100)
  to_be_removed <- sample(1:N , numE)
  attrs[to_be_removed] <- NA
  return(attrs)
}
#G <- cowNW
#V(G)$years <-  removeRandAtt(V(G)$years, 10)

removeCatAtt <- function(attr, cat, fxn, percent){
  
  L <- which(fxn(attr, cat))
  N <- length(L)
  numE <- floor((N*percent)/100)
  to_be_removed <- sample( L, numE)
  attr[to_be_removed] <- NA
  return(attr)
}

#V(G)$gender <- removeCatAtt(V(G)$gender, "2", `==`, 20 )

## Strategies for imputing the data

## Unconditional means. This includes just simple methods of imputations.
## for categorical variables impute based on the most frequent variable
## Most common value imputation

imputeCatAtt <- function(attr){
  attr[is.na(attr)] <- which.max(table(attr))#as.matrix(as.data.frame(sort(table(attr), decreasing = TRUE ))$attr[1])
  return(attr)
}

#imputeCatAtt(att$status)

imputeContAtt <- function(attr){
  attr[is.na(attr)] <- floor(mean(attr,na.rm = TRUE))
  return(attr)
} 



CatRegImp <- function(att, attr){
  
  fitDat <- subset(att, !is.na(att[[attr]]))
  pred <- subset(att, is.na(att[[attr]]))[, -which(colnames(att) %in% attr)]
  form <- as.formula(paste0(attr,"~", paste(colnames(att)[colnames(att) != attr], sep = "+", collapse = "+")))
  
  if(nlevels(att[[attr]]) > 2){
    fit  <- multinom(form, data =fitDat)
    op   <- predict(fit, "probs", newdata = pred)
    op <- ifelse(op > 0.5, 1, 0)
    
  } else {
    fit <- glm(as.formula(paste0(attr,"~", paste(colnames(att)[colnames(att) != attr], sep = "+", collapse = "+"))), 
               data = fitDat, family = binomial(link = "logit"))
    op <- predict.glm(fit, newdata = pred, type = "response")
    op <- ifelse(op > 0.5, 1, 0)
  }
  
  att[[attr]][is.na(att[[attr]])] <- op
  
  return(att)
}



## Required functions
logit2prob <- function(coef) {
  odds <- exp(coef)
  prob <- odds / (1 + odds)
  return(prob)
}

prob2logit <- function(x){
  return(log(x / (1 - x)))
}

## Create a basis network

NetworkSim <- function(N, dir=FALSE, B, C,formula, coefs){
  
  #seed <- 42
  net <- network(N, directed = dir, density= 0 )
  net %v% 'Cluster' <- C
  
  #coefs = c(prob2logit(B))
  
  if(!is.null(formula)){
    g.sim <- simulate(as.formula(paste0("net~nodemix('Cluster',levels = TRUE,levels2 = TRUE)+", formula)), 
                      coef = coefs)#,
                      #seed = seed)  
  }
  else{
    g.sim <- simulate(as.formula(paste0("net~nodemix('Cluster',levels = TRUE,levels2 = TRUE)")), 
                      coef = coefs)#,
                      #seed = seed)
  }
  
  return(g.sim)
  
}



SimNW <- function(pC,k,N,B, formula){
  C = sample(1:k,N,replace = TRUE, prob = pC)
  coefs = c(prob2logit(B))
  
  ## Simulating the basis network with clusters
  g.sim <- NetworkSim(N,FALSE,B,C,formula, coefs)
  return(g.sim)
}

AssgnCat <- function(g.sim, catName, pCat, cat,k){
  g.sim %v% catName <- NA
  j <- 1
  for (i in 1:k) {
    bk <- g.sim %v% catName
    bk[which(g.sim %v% 'Cluster' == i)] <- sample(cat,
                                                  length( which(g.sim %v% 'Cluster' == i)), 
                                                  prob = pCat[[i]],
                                                  replace = TRUE)
    g.sim %v% catName <- bk
  }
  g.sim %v% paste0(catName,"Orig") <- g.sim %v% catName
  g.sim <- asIgraph(g.sim)
  return(g.sim)
}

comdet <- function(g.sim, cat, cont = NULL){
  # running regularized spectral clustering on the original network without any covariate information.
  origmem <- RegSpectralClust(g.sim, 3,regularize = TRUE)
  
  ## Covariate assisted spectral clustering 
  ## Since we are using categorical data we have to perform one hot encoding.
  X <- as.matrix(cat)
  Xdum <- as.matrix(data.frame(predict(dummyVars("~.", data = X ), newdata <- X)))
  if(length(cont) > 0){
    Xdum <- cbind(Xdum, cont)
  }
  origmemCov <- CovAssistedSpecClust(G = g.sim, Xdum, 3,
                                     Regularize =TRUE, alpha =NA,
                                     type ="non-assortative")
  return(data.frame(origmem = origmem, origmemCov = origmemCov))
}
# 

## Check correlation between the covariates and the community generated.
# 
# table(origmemCov, V(g.sim)$LOTR)
# 
# table(origmemCov, V(g.sim)$Cluster)
# 
# table(V(g.sim)$LOTR, V(g.sim)$Cluster)

# create missingness in simulated networks

## Missingness created randomly.
Mymode <- function(codes){
  names(which.max(table(codes)))
}

Mymean <- function(contVar){
  return(mean(contVar,na.rm =TRUE ))
}

standardize <- function(x){return((x-min(x))/(max(x)-min(x)))}

IterativeImputeQual <- function(g.sim,N,coms,niter, k,covMsng, Midx, impType = "Mode"){
  #print(covMsng$name[[1]])
     ## Based on notes from 8th
  j <- 1
  g <- g.sim
  vtxat <- as.data.frame(vertex.attributes(g))
  ## start with randomly assigned covariate values
  for (i in 1:dim(covMsng)[1]) {
    ## setting random initial values to vertices 
    g <- set_vertex_attr(g, 
                         covMsng$name[i], 
                         index = c(Midx[[i]]), 
                         sample(unique(covMsng$cat[[i]]),
                                size  = length(Midx[[i]]),
                                replace = TRUE))
  }
  ## assigning random communitites to begin
  V(g)$com <- sample(1:k, size = N, replace = TRUE)
  prevComLst <- V(g)$com#coms$origmemCov

  gOrig<- g

  op <- matrix(nrow=0,ncol = 12)#data.frame(matrix(nrow=0,ncol = 12))
  covARI <- data.frame(matrix(nrow = 0,ncol =2))
  for(impT in impType){
    g <- gOrig
    for (i in 1:niter){
    ## Run community detection 
      X <- as.matrix(vertex_attr(g, covMsng$name))#as.matrix(V(g)$LOTR)
      Xdum <- as.matrix(data.frame(predict(dummyVars("~.", data = X ),newdata <- X)))
      #if(impT == "MultRegression"){
      Xdum <- cbind(Xdum, as.numeric(V(g)$Continuous))
      #}

      ## USe the eigen values from known network only to find alpha.
      #print(Midx[[1]])
      #alphaComp <- findAlpha(g, Xdum[!(1:length(vtxat[,1])) %in% Midx[[1]],], Regularize =TRUE)
      com <- CovAssistedSpecClust(G = g, Xdum, 3,
                                  Regularize =TRUE, alpha =NA, 
                                  type ="assortative",
                                  kmeansIter = 1000,
                                  retAlpha = TRUE)
      alpha <- com [[2]]
      #print(paste0("AlphaVal:", alpha))
      com <- com[[1]]
      V(g)$com <- com
    
      ## Just use neighbors within cluster for imputation. If there is a tie deal with ties randomly.
      #for(j in 1:length(Midx)){
        oldCov <- vertex_attr(g, name = covMsng$name[j])
        if(impT == "Mode"){
          for (idx in Midx[[j]]) {

            b <- as.data.frame(vertex.attributes(g, index = neighbors(g,idx)))
            ## Number of neighbours with the same community assignment as the missing node
            oldcom <- sum(b$com == V(g)$com[idx])/length(b$com)
            b <- b[b$com == V(g)$com[idx], ]
            g <- set_vertex_attr(g, covMsng$name[j], 
                                index = idx, 
                                Mymode(b[[covMsng$name[j]]]))
                                ## Stopping rule: are all the neighbors in the same community
            status <- all(vertex_attr(g,index = neighbors(g,idx))$com == V(g)$com[idx])
            op <- rbind(op, cbind(i, 
                                  idx, 
                                  oldcom,
                                  vertex_attr(g, name = covMsng$name[j],
                                              index = idx),
                                  oldCov[idx],
                                  vertex_attr(g.sim, name = covMsng$name[j],
                                              index = idx),
                                  mclust::adjustedRandIndex(coms$origmemCov, 
                                                            V(g)$com),
                                  mclust::adjustedRandIndex(prevComLst,
                                                            V(g)$com),
                                  vtxat$Cluster[as.numeric(idx)],
                                  status, 
                                  alpha,
                                  impT))
            prevComLst <-  V(g)$com
          }
        }
        if(impT == "MultRegression"){
          # Use the multinomial value as response. So categories are the response and continuous variable is the explanatory variable
          ## Subsetting the index values that are missing
          ## using both the cluster assignment and the continuous variable in multinomial regression.
          ## converting the community detection to dummy variables
          #dat <- as.matrix(as.factor(V(g)$com))
          #3dat <- as.matrix(data.frame(predict(dummyVars("~.", data = dat ),newdata <- dat)))
          #colnames(Xdum) <- as.character(1:k)
          dat <- data.frame(cbind(cat = unlist(as.character(V(g)$Categories)), 
                                  cont = standardize(vtxat$Continuous), 
                                  com = as.factor(V(g)$com)))
          #data.frame(cbind(cat = vtxat$Categories, cont = vtxat$Continuous, Xdum))
          dat$cont <- as.numeric(dat$cont)
          topred <- dat[(1:length(vtxat[,1])) %in% Midx[[j]],]
          totrain <- dat[!(1:length(vtxat[,1])) %in% Midx[[j]],]
          fit_multi <- multinom(cat ~., data =totrain)
          predClass <- predict(fit_multi, newdata = topred[, -1], type= "class")
          g <- set_vertex_attr(g, covMsng$name[j], 
                                index = unlist(Midx[[j]]), 
                                unlist(as.character(predClass)))
          ## check if all the neighbours for each missing node are in the same community
          status <- rep(NA, length(Midx[[j]]))
          comStat <- 1
          for(idx in Midx[[j]]){
          status[comStat] = all(vertex_attr(g,index = neighbors(g,idx))$com == V(g)$com[idx])
          comStat <- comStat+1
          }
          op <- rbind(op, cbind(i, 
                                  Midx[[j]], 
                                  1,
                                  vertex_attr(g, name = covMsng$name[j],
                                              index = unlist(Midx[[j]])),
                                  oldCov[Midx[[j]]],
                                  vertex_attr(g.sim, name = covMsng$name[j],
                                              index = unlist(Midx[[j]])),
                                  mclust::adjustedRandIndex(coms$origmemCov, 
                                                            V(g)$com),
                                  mclust::adjustedRandIndex(prevComLst,
                                                            V(g)$com),
                                  vtxat$Cluster[as.numeric(Midx[[j]])],
                                  status,
                                  alpha,
                                  impT))
        prevComLst <-  V(g)$com
        }

        
        v <- mclust::adjustedRandIndex(vertex_attr(g, name = covMsng$name[j],
                                                  index = Midx[[j]]),
                                      vertex_attr(g.sim, name = covMsng$name[j],
                                                  index = Midx[[j]]))
        covARI <- rbind(covARI,cbind(rep(v ,times = length(Midx[[j]]))))
      #}
    }
  }

  op <- as.data.frame(op)
  colnames(op) <- c("iteration","MissingId","PropSameComm",
                    "Imputed","prev","Assigned",
                    "ARI","PrevARI","Cluster","neighComStat","alpha","ImpType")
  op$covARI <- covARI[,1]
  op$iteration <- as.numeric(op$iteration)
  op$MissingId <- as.numeric(op$MissingId)
  #op$alpha <- alpha
  return(list(op, g))
}

IterativeImputeCont <- function(g.sim,N,coms,niter, k,covMsng, Midx, impType = "Mean"){
  # j <- 1
  g <- g.sim
  vtxat <- as.data.frame(vertex.attributes(g))
  ## start with randomly assigned covariate values
  #print("first")
  rep <- as.numeric(unlist(sample(vtxat$Continuous,#[(1:length(vtxat[,1])) %in% Midx],
                                    size  = length(Midx),
                                    replace = TRUE)))
  #print(unlist(Midx[[1]]))
  #for (i in 1:dim(covMsng)[1]) {
    ## setting random initial values to vertices 
  g <- set_vertex_attr(g, 
                       covMsng$name[1], 
                       index = as.numeric(Midx), 
                       rep)
  #}
  ## assigning random communitites to begin
  V(g)$com <- sample(1:k, size = N, replace = TRUE)
  prevComLst <- V(g)$com#coms$origmemCov

  gOrig<- g

  op <- matrix(nrow=0,ncol = 12)
  covARI <- data.frame(matrix(nrow = 0,ncol =2))
  for(impT in impType){
    #print(impT)
    g <- gOrig
    for (i in 1:niter){
    ## Run community detection 
      X <- as.matrix(vertex_attr(g, "Categories"))#as.matrix(V(g)$LOTR)
      Xdum <- as.matrix(data.frame(predict(dummyVars("~.", data = X ),newdata <- X)))
      #if(impT == "MultRegression"){
      Xdum <- cbind(Xdum, V(g)$Continuous)
      #}
      ## USe the eigen values from known network only to find alpha.
      #print(Midx[[1]])
      #alphaComp <- findAlpha(g, Xdum[!(1:length(vtxat[,1])) %in% Midx[[1]],], Regularize =TRUE)
      #print(Xdum)
      com <- CovAssistedSpecClust(G = g, Xdum, 3,
                                  Regularize =TRUE, alpha =NA, 
                                  type ="assortative",
                                  kmeansIter = 1000,
                                  retAlpha = TRUE)
      alpha <- com [[2]]
      #print(paste0("AlphaVal:", alpha))
      com <- com[[1]]
      V(g)$com <- com
     ## Just use neighbors within cluster for imputation. If there is a tie deal with ties randomly.
      #for(j in 1:length(Midx)){
        #print(length(Midx[[1]]))

        #Midx[[j]] <- as.numeric(Midx[[j]])
        oldCov <- vertex_attr(g, name = covMsng$name[1])
        if(impT == "Mean"){
          #print("mean")
          for (idx in Midx) {

            b <- as.data.frame(vertex.attributes(g, index = neighbors(g,idx)))
            ## Number of neighbours with the same community assignment as the missing node
            oldcom <- sum(b$com == V(g)$com[idx])/length(b$com)
            #print("new idx")
            #print(V(g)$com[idx])
            #print(sum(b$com == V(g)$com[idx]))
            #print(length(b$com))
            bnew <- b[b$com == V(g)$com[idx], ]
            mn <- Mymean(as.numeric(bnew$Continuous))
            #print(mn)
            if(sum(b$com == V(g)$com[idx]) == 0){
              mn <- Mymean(as.numeric(V(g)$Continuous))
            }
            #print("2nd")
            #print(Mymean(as.numeric(b$Continuous)))
            #print(b)
            g <- set_vertex_attr(g, "Continuous", 
                                index = idx, 
                                mn)
                                ## Stopping rule: are all the neighbors in the same community
            status <- all(vertex.attributes(g,index = neighbors(g,idx))$com == V(g)$com[idx])
            op <- rbind(op, cbind(i, 
                                  idx, 
                                  oldcom,
                                  vertex_attr(g, name = covMsng$name[1],
                                              index = idx),
                                  oldCov[idx],
                                  vertex_attr(g.sim, name = covMsng$name[1],
                                              index = idx),
                                  mclust::adjustedRandIndex(coms$origmemCov, 
                                                            V(g)$com),
                                  mclust::adjustedRandIndex(prevComLst,
                                                            V(g)$com),
                                  vtxat$Cluster[as.numeric(idx)],
                                  status, 
                                  alpha,
                                  impT))
            prevComLst <-  V(g)$com
          }
        }
        if(impT == "QuantRegression"){
          # Use the multinomial value as response. So categories are the response and continuous variable is the explanatory variable
          ## Subsetting the index values that are missing
          ## using both the cluster assignment and the continuous variable in multinomial regression.
          ## converting the community detection to dummy variables
          dat <- data.frame(cbind(cont = V(g)$Continuous, 
                                  cat = V(g)$Categories, 
                                  com = as.factor(V(g)$com)))
          dat$cont <- as.numeric(dat$cont)
          topred <- dat[Midx,]#topred <- dat[(1:length(vtxat[,1])) %in% Midx[[1]],]
          totrain <- dat[-c(Midx),]#totrain <- dat[!((1:length(vtxat[,1])) %in% Midx[[1]]),]
          fit_quant <- lm(cont ~., data=totrain)
          pred <- predict(fit_quant, newdata = topred[, -1])
          g <- set_vertex_attr(g, "Continuous", 
                                index = Midx, 
                                as.numeric(unlist(pred)))
          ## check if all the neighbours for each missing node are in the same community
          status <- rep(NA, length(Midx))
          comStat <- 1
          for(idx in Midx){
          status[comStat] = all(vertex.attributes(g,index = neighbors(g,idx))$com == V(g)$com[idx])
          comStat <- comStat+1
          }
          op <- rbind(op, cbind(i, 
                                  Midx, 
                                  1,
                                  vertex_attr(g, name = covMsng$name[1],
                                              index = unlist(Midx)),
                                  oldCov[Midx],
                                  vertex_attr(g.sim, name = covMsng$name[1],
                                              index = unlist(Midx)),
                                  mclust::adjustedRandIndex(coms$origmemCov, 
                                                            V(g)$com),
                                  mclust::adjustedRandIndex(prevComLst,
                                                            V(g)$com),
                                  vtxat$Cluster[as.numeric(Midx)],
                                  status,
                                  alpha,
                                  impT))
         prevComLst <-  V(g)$com
        }

        v <- mclust::adjustedRandIndex(vertex_attr(g, name = covMsng$name[1],
                                                  index = Midx),
                                       vertex_attr(g.sim, name = covMsng$name[1],
                                                  index = Midx))
        covARI <- rbind(covARI,
                        cbind(rep(v ,times = length(Midx))))
      #}
    }
  }

  op <- as.data.frame(op)
  colnames(op) <- c("iteration","MissingId","PropSameComm",
                    "Imputed","prev","Assigned",
                    "ARI","PrevARI","Cluster","neighComStat","alpha","ImpType")
  op$covARI <- covARI[,1]
  op$iteration <- as.numeric(op$iteration)
  op$MissingId <- as.numeric(op$MissingId)
  #op$alpha <- alpha
  return(list(op, g))
}


IterativeImputeSimul <- function(g.sim,N,coms,niter, k, MidxCat, MidxCon, impType = "Simul"){
  g <- g.sim
  vtxat <- as.data.frame(vertex.attributes(g))
  #print("first print")
  #print(MidxCon)
  ## start with randomly assigned covariate values
  ## setting random initial continuous values to vertices 
  rep <- as.numeric(unlist(sample(vtxat$Continuous,
                                  size  = length(MidxCon),
                                  replace = TRUE)))
  g <- set_vertex_attr(graph= g, 
                       name ="Continuous", 
                       index = as.numeric(MidxCon), 
                       value= rep)
  ##Setting random initial values for Categorical variable
  #print("2nd print")
  g <- set_vertex_attr(g, 
                       "Categories", 
                       index = c(MidxCat), 
                       sample(unique(vtxat$Categories),
                              size  = length(MidxCat),
                              replace = TRUE))

  ## assigning random communitites to begin
  V(g)$com <- sample(1:k, size = N, replace = TRUE)
  prevComLst <- V(g)$com#coms$origmemCov

  op <- matrix(nrow=0,ncol = 12)
  covARI <- data.frame(matrix(nrow = 0,ncol =2))
  gOrig<- g
  #covARI <- data.frame(matrix(nrow = 0,ncol =2))
  for(impT in impType){
    g <- gOrig
    for (i in 1:niter){
    ## Run community detection 
      X <- as.matrix(vertex_attr(g, "Categories"))
      Xdum <- as.matrix(data.frame(predict(dummyVars("~.", data = X ),newdata <- X)))
      Xdum <- cbind(Xdum, V(g)$Continuous)
      ## USe the eigen values from known network only to find alpha.
      com <- CovAssistedSpecClust(G = g, Xdum, 3,
                                  Regularize =TRUE, alpha =NA, 
                                  type ="assortative",
                                  kmeansIter = 1000,
                                  retAlpha = TRUE)
      alpha <- com [[2]]
      com <- com[[1]]
      V(g)$com <- com
        oldCov <- c(vertex_attr(g, name = "Continuous",index = unlist(MidxCon)),
                    vertex_attr(g, name = "Categories",index = unlist(MidxCat)))
        if(impT == "MeanMode"){
          
        }
        if(impT == "SimulReg"){

          ## Perform the regression to impute continuous variable.
          dat <- data.frame(cbind(cont = V(g)$Continuous, 
                                  cat = V(g)$Categories, 
                                  com = as.factor(V(g)$com)))
          dat$cont <- as.numeric(dat$cont)
          topred <- dat[MidxCon,]
          totrain <- dat[-c(MidxCon),]
          fit <- lm(cont ~., data=totrain)

          ## Predict the values of missing continuous data.
          pred <- predict(fit, newdata = topred[, -1])
          g <- set_vertex_attr(g, "Continuous", 
                                index = MidxCon, 
                                as.numeric(unlist(pred)))

          ## Perform multinomial regression to impute categorical variable.
          dat <- data.frame(cbind(cat = V(g)$Categories, 
                                  cont = V(g)$Continuous, 
                                  com = as.factor(V(g)$com)))
          dat$cont <- as.numeric(dat$cont)
          topred <- dat[MidxCat,]
          totrain <- dat[-c(MidxCat),]
          fit <- multinom(cat ~., data=totrain)

          ## Predict the values of missing categorical variable.
          pred <- predict(fit, newdata = topred[, -1], type ="class")
          g <- set_vertex_attr(g, "Categories", 
                               index = MidxCat, 
                               as.character(unlist(pred)))

          ## check if all the neighbours for each missing node are in the same community
          Midx <- c(MidxCon, MidxCat)
          status <- rep(NA, length(Midx))
          comStat <- 1
          for(idx in Midx){
          status[comStat] = all(vertex.attributes(g,index = neighbors(g,idx))$com == V(g)$com[idx])
          comStat <- comStat+1
          }
          op <- rbind(op, cbind(i, 
                                  Midx, 
                                  1,
                                  c(vertex_attr(g, name = "Continuous",
                                              index = unlist(MidxCon)),
                                    vertex_attr(g, name = "Categories",
                                              index = unlist(MidxCat))),
                                  oldCov,
                                  c(vertex_attr(g.sim, name = "Continuous",
                                              index = unlist(MidxCon)),
                                    vertex_attr(g.sim, name = "Categories",
                                              index = unlist(MidxCat))),
                                  mclust::adjustedRandIndex(coms$origmemCov, 
                                                            V(g)$com),
                                  mclust::adjustedRandIndex(prevComLst,
                                                            V(g)$com),
                                  vtxat$Cluster[as.numeric(Midx)],
                                  status,
                                  alpha,
                                  impT))
         prevComLst <-  V(g)$com
        }
        vcat <- mclust::adjustedRandIndex(vertex_attr(g, name = "Categories",
                                                  index = MidxCat),
                                          vertex_attr(g.sim, name = "Categories",
                                                  index = MidxCat))
        vcon <- mclust::adjustedRandIndex(vertex_attr(g, name = "Continuous",
                                                  index = MidxCon),
                                          vertex_attr(g.sim, name = "Continuous",
                                                  index = MidxCon))
        covARI <- rbind(covARI,
                        c(rep(vcon ,times = length(MidxCon)),
                          rep(vcat ,times = length(MidxCat))))
      #}
    }
  }

  op <- as.data.frame(op)
  colnames(op) <- c("iteration","MissingId","PropSameComm",
                    "Imputed","prev","Assigned",
                    "ARI","PrevARI","Cluster","neighComStat","alpha","ImpType")
  op$covARI <- covARI[,1]
  op$iteration <- as.numeric(op$iteration)
  op$MissingId <- as.numeric(op$MissingId)
  #op$alpha <- alpha
  return(list(op, g))
}

plot_combinations <- function(g.sim, coms, p){
  
  par(mfrow=c(2,2))
  
  plot.igraph( g.sim,
               edge.col = adjustcolor('black',alpha.f = 0.5),
               vertex.color = as.factor(V(g.sim)$LOTROrig),
               vertex.label = NA,
               vertex.size = 5,
               main = "Imputed categories")
  
  plot.igraph( g.sim,
               edge.col = adjustcolor('black',alpha.f = 0.5),
               vertex.color = as.factor(V(g.sim)$com),
               vertex.label = NA,
               vertex.size = 5,
               main = paste0("Final Community values :", as.character(round(mclust::adjustedRandIndex(V(g.sim)$Cluster,V(g.sim)$com),3))))
  
  plot.igraph( g.sim,
               edge.col = adjustcolor('black',alpha.f = 0.5),
               vertex.color = V(g.sim)$Cluster,
               vertex.label = NA,
               vertex.size = 5,
               main=paste0("Original cluster assignment"))
  plot.igraph( g.sim,
               edge.col = adjustcolor('black',alpha.f = 0.5),
               vertex.color = coms$origmem,
               vertex.label = NA,
               vertex.size = 5,
               main=paste0("Reg spectral clustering w/o cov ",
                           as.character(round(mclust::adjustedRandIndex(V(g.sim)$Cluster,coms$origmem),3))))
  mtext(paste0("Percent missing ", p))
  
}

CntEdges <- function(g) {
  ends <- as_edgelist(g)
  e1 <- V(g)[ends[,1]]$LOTR
  e2 <- V(g)[ends[,2]]$LOTR
  c1 <- V(g)[ends[,1]]$Cluster
  c2 <- V(g)[ends[,2]]$Cluster
  endcolors <- t(apply(cbind(e1, e2, c1,c2), 1, sort))
  return(table(endcolors[,1], endcolors[,2], endcolors[,3],endcolors[,4]))
}
