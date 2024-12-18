#source("Code/CoDACov.R")
source("Code/CoDACov.R")
# args=commandArgs(trailingOnly = TRUE)
# if (length(args)==0) {
#   stop("At least one argument must be supplied (input file).n", call.=FALSE)
# }
# simvec <- as.integer(gsub("[\r\n]", "",args[1]))
# print(paste0("the simulation number ", simvec))

SaveCoDASim <- function(simvec, sim, InitParamFile){
  
  S <- as.data.frame(readRDS(InitParamFile))
  
  if (sim == TRUE) {
    printFlg <<- FALSE
    Sim <- data.frame(S[simvec,])
    ## Number of simulations
    Nsim <- Sim$Nsim
    
    ##Directed graph yes? No?
    dir <- Sim$DirType
    
    ## Number of in and out binary covariates
    k_out <- Sim$k_out
    k_in <- Sim$k_in
    k <- k_out + k_in
    
    ## Number of in and out continuous covariates
    o_out <- Sim$o_out
    o_in <- Sim$o_in
    o <- o_in + o_out
    
    ## Types of covariates
    covTypes <- Sim$covType#c( "continuous")
    
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
    
    ## Number of communities
    nc <- Sim$nc
    ## Learning rate forgradient descent
    alpha <- Sim$alpha
    
    ##loglik calculation using weight for graph vs covariates
    alphaLL <- Sim$alphaLL
    
    ##penalty of logistic regression
    lambda <- Sim$lambda
    
    ## Switching to percent increase. Using 0.01% increase. 
    ## Change in loglikelihood as stopping criterion
    thresh <- Sim$thresh
    
    ## Randomize the sequence of node update or not
    randomize <- Sim$randomize
    
    ## Percent of missing data in covariates
    missing <- Sim$missing
    
    ## Number of nodes
    N <- Sim$N
    
    # Probability of cluster assignment
    pClust <- Sim$pClust[[1]]
    
    ## Probability for covariate simulation for binary covarites
    pC <- Sim$pC[[1]]
    
    ## Probability for covariate simulation for continuous covarites
    dist <- Sim$dist[[1]]
    
    ## Probability of connection between two nodes of the same community.
    pConn <- Sim$pConn[[1]]
    
    ## Baseline probability of connection
    epsilon <- Sim$epsilon
    
    ## Percentage of outgoing edges in each cluster
    dirPct <- Sim$dirPct[[1]]
    
    ## Maximum number iterations in the stopping criterion
    nitermax <- Sim$nitermax
    
    ## All test output is saved in these lists
    opf_cov <- list()
    opf_noCov <- list()
    
    ## Delta value threshold for assigning communities
    delta <- getDelta(N)
    
    ## saving the original cluster list
    mse <- matrix(nrow = Nsim, ncol = o + k, 0)
    lvl <- 1
    ncVal <- nc
    
    bigN <- 1
    Gtot<- list() 
    GTlogLikcovtot<- list() 
    GTlogLiknoCovtot <- list()
    opf_covtot <- list()
    opf_noCovtot <- list()
    TimeTaken <- list()
    
    for(m in 1:bigN){
      ## Generating network with covariates and overlapping clusters
      NWlst <- genBipartite(N,nc,pClust,k_in,k_out,o_in,o_out,pC,dist,covTypes,
                            CovNamesLPin,CovNamesLPout,CovNamesLinin,CovNamesLinout,
                            pConn,dir= "directed",dirPct,epsilon,missing)
      G_orig <- NWlst[[1]]
      # Run all the algorithms 
      # BigClam: Undirected + No Covariates
      # CoDA: Directed + No Covariates
      # CESNA: Undirected + Covariates
      # CoDACov: Directed + Covariates
      # CESNA Miss: Undirected + Missing Covariate Imputation
      # CoDA Miss: Directed + Missing Covariate Imputation
      
      for(dir in c("directed", "undirected")){
        
        if(dir == "undirected"){
          G <- convertToUndir(G_orig, nc)
          alpha <- 0.0001
        }else{
          G <- G_orig
          }
        
        gbg <- as.data.frame(vertex_attr(G))
        #plot.igraph(G, vertex.color = as.factor(V(G)$Cluster), vertex.label = NA)
        #gbg <- cbind(gbg, Ftot,Htot)#opf_cov[[1]]$Ffin, opf_cov[[1]]$Hfin, ol)
        ## find the ground truth loglikelihood
        GTlogLikcov <- GTLogLik(G,nc,pConn,NULL,
                                CovNamesLinin, CovNamesLinout,
                                CovNamesLPin,CovNamesLPout,
                                k_in,k_out,o_in,o_out, epsilon, dir)
        GTlogLiknoCov <- GTLogLik(G,nc,pConn,NULL,
                                  CovNamesLinin=c(), CovNamesLinout =c(),
                                  CovNamesLPin=c(),CovNamesLPout=c(),
                                  k_in =0 ,k_out =0 ,o_in=0,o_out=0, epsilon, dir)
        orig <<- V(G)$Cluster
        
        ## cluster assignments based on covariate assisted clustering.
        Z <-  as.data.frame(vertex_attr(G)) %>%
          dplyr::select(all_of(c(CovNamesLinin, CovNamesLinout)))
        X <-  as.data.frame(vertex_attr(G)) %>%
          dplyr::select(all_of(c(CovNamesLPin, CovNamesLPout)))
        
        if (length(Z) > 0) {
          cov <- Z
          if (missing > 0) {
            mns <- colMeans(cov, na.rm  = TRUE)
            for (i in 1:o) {
              cov[is.na(cov[, i]), i] <- mns[i]
            }
          }
          origCov <<-  CovAssistedSpecClust(G, cov, nc, alpha = 0.5)
        }
        if (length(X) > 0) {
          cov <- X
          if (missing > 0) {
            for (i in 1:k) {
              cov[is.na(cov[, i]), i] <- mode(cov[, i])
            }
          }
          origCov <<-  CovAssistedSpecClust(G, cov, nc, alpha = 0.5)
        }
        
        ## Running Nsim simulations for the above settings.
        
        for (j in 1:Nsim) {
          print(paste0("iteration ",lvl," alpha ",alpha,  " alphaLL ",alphaLL," num cmnty ",nc))
          start <- Sys.time()
          
          ## Algorithm with covariates + possible missing data
          opf_cov[[lvl]] <- CoDA(G,nc,k = c(k_in, k_out),o = c(o_in, o_out),N,alpha,
                                 lambda,thresh,nitermax,orig,randomize = TRUE,
                                 CovNamesLinin,CovNamesLinout,CovNamesLPin,CovNamesLPout, dir,
                                 alphaLL = NULL,test = TRUE)
          tme <- Sys.time() - start
          print(paste0("Total time take ", tme))
          TimeTaken[[lvl]] <- tme
          ## Algorithm without covariates + possible missing data
          opf_noCov[[lvl]] <- CoDA(G,nc,k = c(0, 0),o = c(0, 0),N,alpha,
                                   lambda,thresh ,nitermax,orig,randomize = TRUE,
                                   CovNamesLinin = c(),CovNamesLinout = c(),CovNamesLPin = c(),
                                   CovNamesLPout = c(),dir,alphaLL = NULL,test = TRUE)
          lvl <- lvl + 1
        }
        Gtot <-  list(Gtot, G)
        GTlogLikcovtot <-  list(GTlogLikcovtot, GTlogLikcov)
        GTlogLiknoCovtot <- list(GTlogLiknoCovtot, GTlogLiknoCov)
      }
    }
    opf_covtot <- opf_cov#list(opf_covtot, opf_cov)
    opf_noCovtot <- opf_noCov#list(opf_noCovtot, opf_noCov)
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
    Nsim <- 10
    #opf <- list()
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
    #EG <- expand.grid(1:Nsim, alpha, alphaLL, nc)
  }
  return(list(opf_covtot, opf_noCovtot,Gtot, GTlogLikcovtot, GTlogLiknoCovtot))
  #return(0)
}


lst <- SaveCoDASim(simvec = 1, 
                   sim = TRUE, 
                   InitParamFile = "InitParamMiss_1.rds")

saveRDS(lst, paste0(file,"_OP", simvec,".rds"))

