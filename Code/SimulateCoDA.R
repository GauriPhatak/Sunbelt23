#source("Code/CoDACov.R")
source("Code/CoDACov.R")
# args=commandArgs(trailingOnly = TRUE)
# if (length(args)==0) {
#   stop("At least one argument must be supplied (input file).n", call.=FALSE)
# }
# simvec <- as.integer(gsub("[\r\n]", "",args[1]))
# print(paste0("the simulation number ", simvec))

SaveCoDASim <- function(simvec, sim, InitParamFile){
  
  S <- as.data.frame(readRDS(paste0(getwd(),InitParamFile)))
  
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
    ## Seperate linear regression alpha
    alphaLin <- Sim$alphaLin
    
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
    missing <- Sim$missing[[1]]
    
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
    
    ## Beta distribution alpha and beta values
    a <- Sim$a
    b <- Sim$b
    
    ## Type of network
    Type <- Sim$Type
    
    ## Cluster Overlap proportion
    pClustOL <- Sim$pClustOL[[1]]
    
    ## All test output is saved in these lists
    opf_Regcov <- list()
    opf_StoRegcov <- list()
    opf_noCov <- list()
    
    ## Deciding penalty type
    penalty <- Sim$penalty
    
    ## Delta value threshold for assigning communities
    delta <- getDelta(N)
    
    ## saving the original cluster list
    mse <- matrix(nrow = Nsim, ncol = o + k, 0)
    lvl <- 1
    ncVal <- nc
    
    bigN <- Sim$bigN
    Gtot<- list() 
    GTlogLikcovtot<- list() 
    GTlogLiknoCovtot <- list()
    opf_Regcovtot <- list()
    opf_StoRegcovtot <- list()
    opf_noCovtot <- list()
    TimeTaken <- list()
    
    for(m in 1:bigN){
      ## Generating network with covariates and overlapping clusters
      NWlst <- genBipartite(N,nc,pClust,k_in,k_out,o_in,o_out,pC,dist,covTypes,
                            CovNamesLPin,CovNamesLPout,CovNamesLinin,CovNamesLinout,
                            pConn,dir= "directed",dirPct,epsilon,missing, 
                            a, b, Type, pClustOL)
      G_orig <- NWlst[[1]]
      orig_Cov <- NWlst[[2]]
      F_u <- NWlst[[3]]
      H_u <- NWlst[[4]]
      igraph::plot.igraph(G_orig,vertex.label = NA, vertex.size = 5,
                          vertex.color = as.factor(V(G_orig)$Cluster),
                          edge.arrow.size= 0.1, edge.color = "grey28")
      
      # View(cbind(V(G_orig)$bvout1, V(G_orig)$bvout2, V(G_orig)$bvout3))
      # Run all the algorithms 
      # BigClam: Undirected + No Covariates
      # CoDA: Directed + No Covariates
      # CESNA: Undirected + Covariates
      # CoDACov: Directed + Covariates
      # CESNA Miss: Undirected + Missing Covariate Imputation
      # CoDA Miss: Directed + Missing Covariate Imputation
      
        for(dir in c("directed", "undirected")){
          print(paste0("This is ", dir, " network."))
          if(dir == "undirected"){
            G <- convertToUndir(G_orig, nc)
            F_u_tot <- F_u + H_u
            F_u <- F_u_tot
            H_u <- F_u_tot
            #alpha <- 0.0001
            #thresh <- 0.0001
          }else{
            G <- G_orig
          }
          
          #gbg <- as.data.frame(vertex_attr(G))
          
          #plot.igraph(G, vertex.color = as.factor(V(G)$Cluster), vertex.label = NA)
          #gbg <- cbind(gbg, Ftot,Htot)#opf_cov[[1]]$Ffin, opf_cov[[1]]$Hfin, ol)
          ## find the ground truth loglikelihood
          GTlogLikcovLst <- GTLogLik(G,nc,pConn,NULL,
                                     CovNamesLinin, CovNamesLinout,
                                     CovNamesLPin,CovNamesLPout,
                                     k_in,k_out,o_in,o_out, epsilon, dir, 
                                     orig_Cov, Fm = F_u, Hm = H_u)
          GTlogLikcov <- GTlogLikcovLst[[1]]
          GTlogLikcovSep <- GTlogLikcovLst[[2]]
          GTlogLiknoCovLst <- GTLogLik(G,nc,pConn,NULL,
                                       CovNamesLinin=c(), CovNamesLinout =c(),
                                       CovNamesLPin=c(),CovNamesLPout=c(),
                                       k_in =0 ,k_out =0,o_in=0,o_out=0,epsilon,dir,
                                       orig_Cov, F_u, H_u)
          GTlogLiknocov <- GTlogLiknoCovLst[[1]]
          GTlogLiknocovSep <- GTlogLiknoCovLst[[2]]
          orig <<- V(G)$Cluster
          
          ## cluster assignments based on covariate assisted clustering.
          Z <-  as.data.frame(vertex_attr(G)) %>%
            dplyr::select(all_of(c(CovNamesLinin, CovNamesLinout)))
          X <-  as.data.frame(vertex_attr(G)) %>%
            dplyr::select(all_of(c(CovNamesLPin, CovNamesLPout)))
          
          if (length(Z) > 0) {
            cov <- Z
            if (!is.null(missing)) {
              mns <- colMeans(cov, na.rm  = TRUE)
              for (i in 1:o) {
                cov[is.na(cov[, i]), i] <- mns[i]
              }
            }
            origCov <<-  CovAssistedSpecClust(G, cov, nc, alpha = 0.5)
          }
          if (length(X) > 0) {
            cov <- X
            if (!is.null(missing)) {
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
            opf_Regcov[[lvl]] <- CoDA(G,nc,k = c(k_in, k_out),o = c(o_in, o_out),N,alpha,
                                   lambda,thresh,nitermax,orig,randomize = TRUE,
                                   CovNamesLinin,CovNamesLinout,CovNamesLPin,CovNamesLPout, dir,
                                   alphaLL = NULL,test = TRUE,missing = missing, covOrig = orig_Cov, 
                                   epsilon =0, impType = "Reg", alphaLin, penalty )
            tme <- Sys.time() - start
            print(paste0("Total time take algo with covariates and simple regression", tme))
            
            start <- Sys.time()
            opf_StoRegcov[[lvl]] <- CoDA(G,nc,k = c(k_in, k_out),o = c(o_in, o_out),N,alpha,
                                           lambda,thresh,nitermax,orig,randomize = TRUE,
                                           CovNamesLinin,CovNamesLinout,CovNamesLPin,CovNamesLPout, dir,
                                           alphaLL = NULL,test = TRUE,missing = missing, covOrig = orig_Cov,
                                           epsilon = 0, impType = "StochasticReg", alphaLin, penalty )
            tme <- Sys.time() - start
            print(paste0("Total time take algo with covariates and stochastic regression", tme))
            
            ## Algorithm without covariates + possible missing data
            start <- Sys.time()
            opf_noCov[[lvl]] <- 0
              # CoDA(G,nc,k = c(0, 0),o = c(0, 0),N,alpha,
              #                        lambda,thresh ,nitermax,orig,randomize = TRUE,
              #                        CovNamesLinin = c(),CovNamesLinout = c(),CovNamesLPin = c(),
              #                        CovNamesLPout = c(),dir,alphaLL = NULL,test = TRUE,
              #                        missing = missing, covOrig = orig_Cov, epsilon =0, impType = "", alphaLin, penalty)
            tme <- Sys.time() - start
            print(paste0("Total time take algo without covariates", tme))
            
            lvl <- lvl + 1
          }
          Gtot <-  append(Gtot, list(G))
          GTlogLikcovtot <-  append(GTlogLikcovtot, GTlogLikcov)
          GTlogLiknoCovtot <- append(GTlogLiknoCovtot, GTlogLiknocov)
        }
    }
    opf_Regcovtot <- opf_Regcov
    opf_StoRegcovtot <- opf_StoRegcov
    opf_noCovtot <- opf_noCov
  }
  
  return(list(opf_Regcovtot, opf_noCovtot,opf_StoRegcovtot,Gtot, GTlogLikcovtot, GTlogLiknoCovtot))
  #return(0)
}

lst2 <- SaveCoDASim(simvec = 2, 
                    sim = TRUE, 
                    InitParamFile = "/Code/InitParamMiss_Cohesive_MCAR.rds")

# lst2 <- SaveCoDASim(simvec = 1, 
#                     sim = TRUE, 
#                     InitParamFile = "/Code/InitParamMiss_Nested.rds")
# 
# lst3 <- SaveCoDASim(simvec = 1, 
#                     sim = TRUE, 
#                     InitParamFile = "/Code/InitParamMiss _2-Mode.rds")

#saveRDS(lst, paste0(file,"_OP", simvec,".rds"))
