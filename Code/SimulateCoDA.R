#source("Code/CoDACov.R")
source(paste0(getwd(),"/Code/CoDACov.R"))
# args=commandArgs(trailingOnly = TRUE)
# if (length(args)==0) {
#   stop("At least one argument must be supplied (input file).n", call.=FALSE)
# }
# simvec <- as.integer(gsub("[\r\n]", "",args[1]))
# print(paste0("the simulation number ", simvec))
#predVals <<- rep(0, 5)
SaveCoDASim <- function(simvec, sim, InitParamFile,getMetrics,reg,nocov){
  
  S <- as.data.frame(readRDS(paste0(getwd(),InitParamFile)))
  
  if (sim == TRUE) {
    printFlg <<- FALSE
    test <- F
    Sim <- data.frame(S[simvec,])
    ## Number of simulations
    Nsim <- 2#Sim$Nsim
    bigN <- 7#Sim$bigN
    
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
    covTypes <- Sim$covType
    
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
    nc_sim <- Sim$nc
    
    ## Save BIc for the number of clusters
    test_nc <- Sim$test_nc
    
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
    ## Setting missingnness type
    missType <- Sim$mt[[1]]
    MARparam <- Sim$mar[[1]]
    ## Number of nodes
    N <- 200#Sim$N
    
    # Probability of cluster assignment
    pClust <- Sim$pClust[[1]]
    
    ## Probability for covariate simulation for binary covarites
    pC <- Sim$pC[[1]]
    
    ## Probability for covariate simulation for continuous covarites
    dist <- Sim$dist[[1]]
    
    ## Probability of connection between two nodes of the same community.
    pConn <- Sim$pConn[[1]]
    
    ## Baseline probability of connection
    epsilon <- 0.000001#Sim$epsilon
    
    ## Percentage of outgoing edges in each cluster
    dirPct <- Sim$dirPct[[1]]
    
    ## Maximum number iterations in the stopping criterion
    nitermax <- Sim$nitermax
    
    ## Beta distribution alpha and beta values
    a <- 5#Sim$a
    b <- 8#Sim$b
    
    ## Type of network
    Type <- Sim$Type
    
    ## Cluster Overlap proportion
    pClustOL <- Sim$pClustOL[[1]]
    
    ## Deciding penalty type
    penalty <- Sim$penalty
    
    ## Delta value threshold for assigning communities
    #delta <- getDelta(N)
    
    ## What type of inital imputation should be done on the missing data
    covInit <- Sim$covInit[[1]]
    
    ## saving the original cluster list
    mse <- matrix(nrow = Nsim, ncol = o + k, 0)
    lvl <- 1
    ncVal <- nc_sim
    
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
    
    params <- matrix(0, nrow= 0, ncol =35)
    nitermax <- 60000
    
    ## List for saving the original covariates
    orig_Cov <- list()
    
    ## Metric values
    op_sim <- matrix(0, nrow= 0, ncol = 17)
    
    for(m in 1:bigN){
      seed <- m
      ## Generating network with covariates and overlapping clusters
      NWlst <- genBipartite(N,nc_sim,pClust,k_in,k_out,o_in,o_out,pC,dist,covTypes,
                            CovNamesLPin,CovNamesLPout,CovNamesLinin,CovNamesLinout,
                            pConn,dir = "directed",dirPct,epsilon,missing, 
                            a, b, Type, pClustOL, missType, MARparam, seed)
      G_orig <- NWlst[[1]]
      orig_Cov[[m]] <- NWlst[[2]]
      ## Just to compare
      covtmp <-  cbind(as.data.frame(vertex_attr(G_orig)),
                       inDeg = igraph::degree(G_orig, mode = "in"), 
                       outDeg = igraph::degree(G_orig, mode = "out")) 
      
      F_u <- NWlst[[3]]
      H_u <- NWlst[[4]]
      print(igraph::plot.igraph(G_orig,vertex.label = NA,vertex.size = 5,
                                vertex.color = as.factor(V(G_orig)$Cluster),
                                edge.arrow.size= 0.1, edge.color = "grey28"))
      
      
      #
      # Run all the algorithms 
      # BigClam: Undirected + No Covariates
      # CoDA: Directed + No Covariates
      # CESNA: Undirected + Covariates
      # CoDACov: Directed + Covariates
      # CESNA Miss: Undirected + Missing Covariate Imputation
      # CoDA Miss: Directed + Missing Covariate Imputation
      
      for (j in 1:Nsim) {
        seed <- m + j
        for(dir in c( "undirected")){#,"undirected")){
          for(nc in c(4)){#test_nc[[1]]){
            nc <- as.numeric(nc)
          if(dir == "undirected"){
            
            G <- convertToUndir(G_orig, nc_sim)
            F_u_tot <- F_u + H_u 
            F_uGT <- F_u_tot
            H_uGT <- F_u_tot
          }else{
            
            G <- G_orig
            F_uGT <- F_u
            H_uGT <- H_u
          }
          
            # getting the epsilon value
            #E <- ecount(G)
            #epsilon <- 2*E/(N*(N-1))
            delta <- getDelta(N, epsilon)
            
          ## Block to calculate base log likelihood and clustering using spectral clustering. 
          {
             GTlogLikcovLst <-  list(0,0,0)#GTLogLik(G,nc,pConn,NULL,
            #                            CovNamesLinin, CovNamesLinout,
            #                            CovNamesLPin,CovNamesLPout,
            #                            k_in,k_out,o_in,o_out, epsilon, dir, 
            #                            orig_Cov = orig_Cov[[m]], Fm = F_uGT, Hm = H_uGT)
            GTlogLikcov <- GTlogLikcovLst[[1]]
            GTlogLikcovSep <- GTlogLikcovLst[[2]]
            vifVal <- GTlogLikcovLst[[3]]
            GTlogLiknoCovLst <- list(0,0,0)#GTLogLik(G,nc,pConn,NULL,
                                         # CovNamesLinin=c(), CovNamesLinout =c(),
                                         # CovNamesLPin=c(),CovNamesLPout=c(),
                                         # k_in =0 ,k_out =0,o_in=0,o_out=0,epsilon,dir,
                                         # orig_Cov[[m]], F_uGT, H_uGT)
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
              specOP <<-  CovAssistedSpecClust(G, cov, nc, alpha = 0.5)
            }
            if (length(X) > 0) {
              cov <- X
              if (!is.null(missing)) {
                for (i in 1:k) {
                  cov[is.na(cov[, i]), i] <- mode(cov[, i])
                }
              }
              specOP <<-  CovAssistedSpecClust(G, cov, nc, alpha = 0.5)
            }
          }
          
          print(paste0("in the ", dir," iteration ",j," for network ",m,
                       " alpha ",alpha," alphaLL ",alphaLL," alphaLin ",alphaLin, 
                       " lambda ", lambda, " initial impute type ", covInit, 
                       " Penalty, ", penalty, " nc ", nc))
          
          if(reg == TRUE){
            tryCatch({
              ## Algorithm with covariates + possible missing data
              start <- Sys.time()
              opf_Regcov[[lvl]] <- CoDA(G, nc, k = c(k_in, k_out), o = c(o_in, o_out), N, alpha,
                                        lambda, thresh, nitermax, orig, randomize = FALSE,
                                        CovNamesLinin, CovNamesLinout, CovNamesLPin, CovNamesLPout, dir,
                                        alphaLL, test, missing, covOrig = orig_Cov[[m]],
                                        epsilon, impType = "Reg", alphaLin, penalty, seed, covInit, specOP,nc_sim )
              tme <- Sys.time() - start
              print(paste0("Total time take algo with covariates and simple regression ", round(tme, 3)))
            }, error = function(e) {
              # Handle the error
              cat("An error occurred:", conditionMessage(e), "\n")
              NA})
          }
         
          if(nocov == TRUE){
            tryCatch({
              # Algorithm without covariates + possible missing data
              start <- Sys.time()
              opf_noCov[[lvl]] <- CoDA(G, nc, k = c(0, 0), o = c(0, 0), N, alpha, lambda, 
                                       thresh, nitermax, orig, randomize = TRUE, CovNamesLinin = c(), 
                                       CovNamesLinout = c(), CovNamesLPin = c(), CovNamesLPout = c(), 
                                       dir, alphaLL, test, missing = missing, covOrig = orig_Cov[[m]], 
                                       epsilon = 0, impType = "", alphaLin, penalty, seed, covInit, specOP,nc_sim)
              tme <- Sys.time() - start
              print(paste0("Total time take algo without covariates ", round(tme, 3)))
            }, error = function(e) {
              # Handle the error
              NA})              
          }
          
          FS[[lvl]] <- cbind(F_uGT, H_uGT)
          params <- rbind(params, append(Sim, list(bigN = m,nsim = j,dir = dir)))
          
          GTlogLikcovtot <-  append(GTlogLikcovtot, GTlogLikcov)
          GTlogLiknoCovtot <- append(GTlogLiknoCovtot, GTlogLiknocov)
          
          ## Code to get network metrics
          if(getMetrics == TRUE){
            Fm <- opf_Regcov[[lvl]]$Ffin
            Hm <- opf_Regcov[[lvl]]$Hfin
            Q_OL <- ifelse(dir == "directed", "D","U")
            FullM <- as.data.frame(memOverlapCalc(Fm, Hm, delta, N, nc))
            comnty <- data.frame(which(as.matrix(FullM) == 1, arr.ind = TRUE)) %>%
              group_by(col) %>%
              group_map(~.x)

            # Find clustering coef global for directed and undirected network
            # High clustering coef means more clique based highly clustered Network
            # Low means more sparse network
            CC <- ClustCoef(G, dir)

            ## Power law fit check
            PLFit <- PLfitcheck(G)

            # 1: Perfect homophily (nodes only connect within their group).
            # 0: No homophily (random mixing).
            # < 0: Heterophily (nodes prefer different groups).
            NetA <- assortativity(G)
            
            ## finc the triangle participation ratio for each community and the overall average.
            TPRvec <- Comm_TPR(G, Fm, Hm, delta, N, nc, dir)

            if(nocov == TRUE){
              ## find memoverlap for no covariates method
              Fm <- opf_noCov[[lvl]]$Ffin
              Hm <- opf_noCov[[lvl]]$Hfin
              NoCovM <- as.data.frame(memOverlapCalc(Fm, Hm, delta, N, nc))
              ## Feature structure decomposition
              fsd <- c(unlist(Feature_struct_decomp(G, nc, N, delta, NoCovM, FullM)))
            }else{
              fsd <-  c(0, 0, 0)
              }

            if(reg == TRUE){
             randOI <- null_models(G, nc, k = c(k_in, k_out), o = c(o_in, o_out), N, alpha,
                                          lambda, thresh, nitermax, orig, randomize = FALSE,
                                          CovNamesLinin, CovNamesLinout, CovNamesLPin, CovNamesLPout, dir,
                                          alphaLL, test = FALSE, missing, covOrig = orig_Cov[[m]],
                                          epsilon = 0, impType = "Reg", alphaLin, penalty, seed, covInit, specOP )

            }else{
              randOI <- 0
            }
            
            op_sim <- rbind(op_sim, c(
              round(conductance(graph = G, communities = comnty, FullM),4),
              round(partition_density(graph = G, communities = comnty),4),
              round(Q_HardPartition(G, FullM),4),
              round(Q_Overlapping(G, FullM, Q_OL),4),
              round(CQ(G, comnty),4),
              round(ONCut(G, comnty),4),
              round(OP(comnty, G),4),
              CC,
              PLFit,
              NetA,
              fsd,
              randOI,
              dir,
              m,
              j))          
          }
          
          lvl <- lvl + 1
          gc()
        }
        }
      }
      {
        opf_Regcovtot <- opf_Regcov
        opf_noCovtot <- opf_noCov
        Gtot <-  append(Gtot, list(G_orig))
        
      }
    }
    colnames(op_sim) <-c("Cond","PD","Q_HP","Q_OL","CQ","OnCut","OP","CC","PLFit","NetA","oi_b1","oi_b2", "oi_cov","randOI","dir","bigN","nsim")
    params <- as.data.frame(params)
    print(params[1,])
    
    return(list(opf_Regcovtot, opf_noCovtot, Gtot, GTlogLikcovtot, GTlogLiknoCovtot, FS,orig_Cov, params, op_sim))
    #return(0)
  }
}

simvec = 3
sim = TRUE
InitParamFile = "/Code/4clustersims_new.rds"
getMetrics = FALSE
reg = TRUE
nocov = FALSE
df_cc <- list()
i <- 1
for(simvec in c(20,23,26,38,41,44)){
  df_cc[[i]] <- SaveCoDASim(simvec,
                     sim,
                     InitParamFile,
                     getMetrics, 
                     reg, 
                     nocov)
  i <- i+1

  }

saveRDS(df_cc, "CrossCheck.rds")
