#source("Code/CoDACov.R")
source(paste0(getwd(),"/Code/CoDACov.R"))
# args=commandArgs(trailingOnly = TRUE)
# if (length(args)==0) {
#   stop("At least one argument must be supplied (input file).n", call.=FALSE)
# }
# simvec <- as.integer(gsub("[\r\n]", "",args[1]))
# print(paste0("the simulation number ", simvec))
#predVals <<- rep(0, 5)
SaveCoDASim <- function(simvec, sim, InitParamFile){
  
  S <- as.data.frame(readRDS(paste0(getwd(),InitParamFile)))
  
  if (sim == TRUE) {
    printFlg <<- FALSE
    test <- FALSE
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
    missing <- Sim$missing[[1]]#md##c(0,0,10)#
    ## Setting missingnness type
    missType <- Sim$mt[[1]]#c("Random", "Random", "Random")#c("Random", "Random", "Random") #c("Random", "Random", "GT_MAR") #GT_MAR LT_MAR Random
    MARparam <- Sim$mar[[1]]#c(25, 80)
    #ns <- 0#newseed
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
    pClustOL <- Sim$pClustOL[[1]]#ol#Sim$pClustOL[[1]]
    
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
    opf_Regcov <- list()
    
    opf_StoRegcovtot <- list()
    opf_StoRegcov <- list()
    
    opf_noCovtot <- list()
    opf_noCov <- list()
    
    opf_MeanRegcov <- list()
    opf_MeanRegcovtot <- list()
    
    opf_ProxRegcov <- list()
    opf_ProxRegcovtot <- list()
    
    TimeTaken <- list()
    FS <- list()
    
    params <- matrix(0, nrow= 0, ncol =34)
    ## List for saving the original covariates
    orig_Cov <- list()
    for(m in 1:bigN){
      ## Generating network with covariates and overlapping clusters
      NWlst <- genBipartite(N,nc,pClust,k_in,k_out,o_in,o_out,pC,dist,covTypes,
                            CovNamesLPin,CovNamesLPout,CovNamesLinin,CovNamesLinout,
                            pConn,dir= "directed",dirPct,epsilon,missing, 
                            a, b, Type, pClustOL, missType, MARparam)
      G_orig <- NWlst[[1]]
      orig_Cov[[m]] <- NWlst[[2]]
      ## Just to compare
      covtmp <-  as.data.frame(vertex_attr(G_orig)) %>%
        dplyr::select(all_of(c(CovNamesLinin, CovNamesLinout, CovNamesLPin, CovNamesLPout)))
      #gbg <- cbind(covtmp , orig_Cov[[m]])
      F_u <- NWlst[[3]]
      H_u <- NWlst[[4]]
      
      #df <- cbind(F_u, H_u, V(G_orig)$a, V(G_orig)$b, V(G_orig)$c)
      #igraph::plot.igraph(G_orig,vertex.label = NA, vertex.size = 5,
      #                     vertex.color = as.factor(V(G_orig)$Cluster),
       #                    edge.arrow.size= 0.1, edge.color = "grey28")
      
      # View(cbind(V(G_orig)$bvout1, V(G_orig)$bvout2, V(G_orig)$bvout3))
      # Run all the algorithms 
      # BigClam: Undirected + No Covariates
      # CoDA: Directed + No Covariates
      # CESNA: Undirected + Covariates
      # CoDACov: Directed + Covariates
      # CESNA Miss: Undirected + Missing Covariate Imputation
      # CoDA Miss: Directed + Missing Covariate Imputation
      for (j in 1:Nsim) {
          seed <- Sys.time()#sample.int(100000, 1)
          for(dir in c("directed", "undirected")){
            if(dir == "undirected"){
              G <- convertToUndir(G_orig, nc)
              F_u_tot <- F_u + H_u
              F_uGT <- F_u_tot
              H_uGT <- F_u_tot
              #alpha <- 0.0001
              #thresh <- 0.0001
            }else{
              G <- G_orig
              F_uGT <- F_u
              H_uGT <- H_u
            }
            GTlogLikcovLst <- GTLogLik(G,nc,pConn,NULL,
                                       CovNamesLinin, CovNamesLinout,
                                       CovNamesLPin,CovNamesLPout,
                                       k_in,k_out,o_in,o_out, epsilon, dir, 
                                       orig_Cov[[m]], Fm = F_uGT, Hm = H_uGT)
            GTlogLikcov <- GTlogLikcovLst[[1]]
            GTlogLikcovSep <- GTlogLikcovLst[[2]]
            vifVal <- GTlogLikcovLst[[3]]
            GTlogLiknoCovLst <- GTLogLik(G,nc,pConn,NULL,
                                         CovNamesLinin=c(), CovNamesLinout =c(),
                                         CovNamesLPin=c(),CovNamesLPout=c(),
                                         k_in =0 ,k_out =0,o_in=0,o_out=0,epsilon,dir,
                                         orig_Cov[[m]], F_uGT, H_uGT)
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
            
            print(paste0("in the ", dir," iteration ",j," for network ",m," alpha ",alpha," alphaLL ",alphaLL," num cmnty ",nc))
            
            ## Algorithm with covariates + possible missing data
            start <- Sys.time()
            opf_Regcov[[lvl]] <- CoDA(G,nc,k = c(k_in, k_out),o = c(o_in, o_out),N,alpha,
                                   lambda,thresh,nitermax,orig,randomize = TRUE,
                                   CovNamesLinin,CovNamesLinout,CovNamesLPin,CovNamesLPout, dir,
                                   alphaLL, test,missing, covOrig = orig_Cov[[m]], 
                                   epsilon =0, impType = "Reg", alphaLin, penalty, seed )
            tme <- Sys.time() - start
            print(paste0("Total time take algo with covariates and simple regression ", round(tme, 3)))
            
            start <- Sys.time()
            opf_StoRegcov[[lvl]] <- CoDA(G,nc,k = c(k_in, k_out),o = c(o_in, o_out),N,alpha,
                                           lambda,thresh,nitermax,orig,randomize = TRUE,
                                           CovNamesLinin,CovNamesLinout,CovNamesLPin,CovNamesLPout, dir,
                                           alphaLL,test,missing = missing, covOrig = orig_Cov[[m]],
                                           epsilon = 0, impType = "StochasticReg", alphaLin, penalty,seed )
            tme <- Sys.time() - start
            print(paste0("Total time take algo with covariates and stochastic regression ", round(tme, 3)))
            
            ## Algorithm without covariates + possible missing data
            start <- Sys.time()
            opf_noCov[[lvl]] <- CoDA(G,nc,k = c(0, 0),o = c(0, 0),N,alpha,
                                     lambda,thresh ,nitermax,orig,randomize = TRUE,
                                     CovNamesLinin = c(),CovNamesLinout = c(),CovNamesLPin = c(),CovNamesLPout = c(),
                                     dir,alphaLL, test,missing = missing, covOrig = orig_Cov[[m]],
                                     epsilon =0, impType = "", alphaLin, penalty, seed)
            tme <- Sys.time() - start
            print(paste0("Total time take algo without covariates ", round(tme, 3)))
            
            ##Algorithm with mean imputation of missing covariates
            covtmp <-  as.data.frame(vertex_attr(G)) %>%
              dplyr::select(all_of(c(CovNamesLinin, CovNamesLinout, CovNamesLPin, CovNamesLPout))) %>% 
              mutate_all(~ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x))     
            
            G_mean <- G
            for(name in c(CovNamesLinin, CovNamesLinout, CovNamesLPin, CovNamesLPout)){
              G_mean <- set_vertex_attr(G_mean, name,index = V(G_mean), covtmp[[name]])
            }
            start <- Sys.time()
            opf_MeanRegcov[[lvl]] <- CoDA(G_mean,nc,k = c(k_in, k_out),o = c(o_in, o_out),N,alpha,
                                         lambda,thresh,nitermax,orig,randomize = TRUE,
                                         CovNamesLinin,CovNamesLinout,CovNamesLPin,CovNamesLPout, dir,
                                         alphaLL,test,missing = NULL, covOrig = orig_Cov[[m]],
                                         epsilon = 0, impType = "Reg", alphaLin, penalty,seed )
            tme <- Sys.time() - start
            print(paste0("Total time take algo with covariates and reg with mean imputation ", round(tme, 3)))
            
            # start <- Sys.time()
            # opf_ProxRegcov[[lvl]] <- CoDA(G,nc,k = c(k_in, k_out),o = c(o_in, o_out),N,alpha,
            #                              lambda,thresh,nitermax,orig,randomize = TRUE,
            #                              CovNamesLinin,CovNamesLinout,CovNamesLPin,CovNamesLPout, dir,
            #                              alphaLL,test,missing = missing, covOrig = orig_Cov[[m]],
            #                              epsilon = 0, impType = "Reg", alphaLin, "Proximal",seed )
            # tme <- Sys.time() - start
            # print(paste0("Total time take algo with covariates, regression and Proximal penalty ", round(tme, 3)))
            
            FS[[lvl]] <- cbind(F_uGT, H_uGT)
            
            lvl <- lvl + 1
            
            params <- rbind(params, append(Sim, list(bigN = m,nsim = j,dir = dir)))
          }
          GTlogLikcovtot <-  append(GTlogLikcovtot, GTlogLikcov)
          GTlogLiknoCovtot <- append(GTlogLiknoCovtot, GTlogLiknocov)
          
      }
      Gtot <-  append(Gtot, list(G))
      
    }
    opf_Regcovtot <- opf_Regcov
    opf_StoRegcovtot <- opf_StoRegcov
    opf_noCovtot <- opf_noCov
    opf_MeanRegcovtot <- opf_MeanRegcov
    opf_ProxRegcovtot <- opf_ProxRegcov
  }
  params <- as.data.frame(params)
  return(list(opf_Regcovtot, opf_noCovtot, opf_StoRegcovtot, opf_MeanRegcovtot, opf_ProxRegcovtot, Gtot, GTlogLikcovtot, GTlogLiknoCovtot, FS, orig_Cov, params))
  #return(0)
}

print(getwd())
simvec = 4
sim = TRUE
InitParamFile = "/Code/InitParamMiss_Coh_MAR_LASSO_Cont_scaled.rds"
S <- as.data.frame(readRDS(paste0(getwd(),InitParamFile)))
# c <- expand.grid(md = list(c(15,15,15),
#                            c(45,45,45),
#                            c(0,0,15),
#                            c(0,0,45)),
#                  ol = list(c(0.0, 0.0, 0.0, 0),
#                            c(0.1, 0.1, 0.1, 0),
#                            c(0.2, 0.2, 0.2, 0),
#                            c(0.3, 0.3, 0.3, 0)),
#                  initSeed = c(0,42))
#c$mar <- rep(list(NULL, c(23,60)), each = 2)
#c$mt <- rep(list(c("Random", "Random", "Random"), c("Random", "Random", "GT_MAR") ), each= 2)
#for(j in 1:32){
  # md <<- c$md[[j]]
  # ol <<- c$ol[[j]]
  # mar <<- c$mar[[j]]
  # #newseed <<- 42#c$initSeed[j]
  # mt <<- c$mt[[j]]
df_3 <- SaveCoDASim(simvec,
                    sim,
                    InitParamFile)
saveRDS(df, paste0(getwd(),"/Code/CaseStudies/CaseStudyN",(j),".rds"))
#}
