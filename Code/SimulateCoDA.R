#source("Code/CoDACov.R")
source(paste0(getwd(),"/Code/CoDACov.R"))
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
      
      df <- cbind(F_u, H_u, V(G_orig)$a, V(G_orig)$b, V(G_orig)$c)
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
      for (j in 1:Nsim) {
          seed <- Sys.time()#sample.int(100000, 1)
          for(dir in c("directed", "undirected")){
            print(paste0("This is ", dir, " network."))
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
                                       orig_Cov, Fm = F_uGT, Hm = H_uGT)
            GTlogLikcov <- GTlogLikcovLst[[1]]
            GTlogLikcovSep <- GTlogLikcovLst[[2]]
            vifVal <- GTlogLikcovLst[[3]]
            GTlogLiknoCovLst <- GTLogLik(G,nc,pConn,NULL,
                                         CovNamesLinin=c(), CovNamesLinout =c(),
                                         CovNamesLPin=c(),CovNamesLPout=c(),
                                         k_in =0 ,k_out =0,o_in=0,o_out=0,epsilon,dir,
                                         orig_Cov, F_uGT, H_uGT)
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
            
            print(paste0("iteration ",lvl," alpha ",alpha,  " alphaLL ",alphaLL," num cmnty ",nc))
            
            ## Algorithm with covariates + possible missing data
            start <- Sys.time()
            opf_Regcov[[lvl]] <- CoDA(G,nc,k = c(k_in, k_out),o = c(o_in, o_out),N,alpha,
                                   lambda,thresh,nitermax,orig,randomize = TRUE,
                                   CovNamesLinin,CovNamesLinout,CovNamesLPin,CovNamesLPout, dir,
                                   alphaLL = NULL, test,missing = missing, covOrig = orig_Cov, 
                                   epsilon =0, impType = "Reg", alphaLin, penalty, seed )
            tme <- Sys.time() - start
            print(paste0("Total time take algo with covariates and simple regression ", round(tme, 3)))
            
            start <- Sys.time()
            opf_StoRegcov[[lvl]] <- CoDA(G,nc,k = c(k_in, k_out),o = c(o_in, o_out),N,alpha,
                                           lambda,thresh,nitermax,orig,randomize = TRUE,
                                           CovNamesLinin,CovNamesLinout,CovNamesLPin,CovNamesLPout, dir,
                                           alphaLL = NULL,test,missing = missing, covOrig = orig_Cov,
                                           epsilon = 0, impType = "StochasticReg", alphaLin, penalty,seed )
            tme <- Sys.time() - start
            print(paste0("Total time take algo with covariates and stochastic regression ", round(tme, 3)))
            
            ## Algorithm without covariates + possible missing data
            start <- Sys.time()
            opf_noCov[[lvl]] <- CoDA(G,nc,k = c(0, 0),o = c(0, 0),N,alpha,
                                     lambda,thresh ,nitermax,orig,randomize = TRUE,
                                     CovNamesLinin = c(),CovNamesLinout = c(),CovNamesLPin = c(),CovNamesLPout = c(),
                                     dir,alphaLL = NULL, test,missing = missing, covOrig = orig_Cov,
                                     epsilon =0, impType = "", alphaLin, penalty, seed)
            tme <- Sys.time() - start
            print(paste0("Total time take algo without covariates ", round(tme, 3)))
            
            ##Algorithm with mean imputation of missing covariates
            covtmp <-  as.data.frame(vertex_attr(G)) %>%
              dplyr::select(all_of(c(CovNamesLinin, CovNamesLinout, CovNamesLPin, CovNamesLPout))) %>% 
              mutate_all(~ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x))          

            for(name in c(CovNamesLinin, CovNamesLinout, CovNamesLPin, CovNamesLPout)){
              G_mean <- G %>% 
                set_vertex_attr(name,index = V(G), covtmp[[name]])
            }
            start <- Sys.time()
            opf_MeanRegcov[[lvl]] <- CoDA(G_mean,nc,k = c(k_in, k_out),o = c(o_in, o_out),N,alpha,
                                         lambda,thresh,nitermax,orig,randomize = TRUE,
                                         CovNamesLinin,CovNamesLinout,CovNamesLPin,CovNamesLPout, dir,
                                         alphaLL = NULL,test,missing = NULL, covOrig = orig_Cov,
                                         epsilon = 0, impType = "Reg", alphaLin, penalty,seed )
            tme <- Sys.time() - start
            print(paste0("Total time take algo with covariates and reg with mean imputation ", round(tme, 3)))
            
            start <- Sys.time()
            opf_ProxRegcov[[lvl]] <- CoDA(G,nc,k = c(k_in, k_out),o = c(o_in, o_out),N,alpha,
                                         lambda,thresh,nitermax,orig,randomize = TRUE,
                                         CovNamesLinin,CovNamesLinout,CovNamesLPin,CovNamesLPout, dir,
                                         alphaLL = NULL,test,missing = missing, covOrig = orig_Cov,
                                         epsilon = 0, impType = "Reg", alphaLin, "Proximal",seed )
            tme <- Sys.time() - start
            print(paste0("Total time take algo with covariates, regression and Proximal penalty ", round(tme, 3)))
            
            FS[[lvl]] <- cbind(F_uGT, H_uGT)
            
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
    opf_MeanRegcovtot <- opf_MeanRegcov
    opf_ProxRegcovtot <- opf_ProxRegcov
  }
  
  return(list(opf_Regcovtot, opf_noCovtot, opf_StoRegcovtot, opf_MeanRegcovtot, opf_ProxRegcovtot, Gtot, GTlogLikcovtot, GTlogLiknoCovtot, FS))
  #return(0)
}


simvec = 194
sim = TRUE
InitParamFile = "/Code/InitParamMiss_Coh_MCAR_LASSO_Cont_2.rds"
S <- as.data.frame(readRDS(paste0(getwd(),InitParamFile)))

lstNst194_1 <- SaveCoDASim(simvec ,
                       sim ,
                       InitParamFile)
#"/Code/InitParamMiss_Cohesive_Nested_MCAR_B15_LASSO.rds")
df <- lstNst194_1

ends <- matrix(0, nrow =0, ncol = 18)
for(i in 1:12){
  #G <- df[[4]][[i]]
  F_u <- df[[1]][[i]]$Ffin
  H_u <- df[[1]][[i]]$Hfin
  #df <- cbind(F_u, H_u, V(G)$a, V(G)$b, V(G)$c)
  # igraph::plot.igraph(G,vertex.label = NA, vertex.size = 5,
  #                     vertex.color = as.factor(V(G)$Cluster),
  #                     edge.arrow.size = 0.1, edge.color = "grey28")
  ends <- rbind(ends, c(tail(df[[1]][[i]]$OmegaIndex,1), tail(df[[3]][[i]]$OmegaIndex,1), 
                        tail(df[[2]][[i]]$OmegaIndex,1), 
                        tail(df[[4]][[i]]$OmegaIndex,1), tail(df[[5]][[i]]$OmegaIndex,1),
                        tail(df[[1]][[i]]$Loglik,1)[1], tail(df[[3]][[i]]$Loglik,1)[1], 
                        tail(df[[2]][[i]]$Loglik,1)[1],
                        tail(df[[4]][[i]]$Loglik,1)[1],tail(df[[5]][[i]]$Loglik,1)[1],
                        mean(tail(df[[1]][[i]]$MSEMD,1)), mean(tail(df[[3]][[i]]$MSEMD,1)), 
                        mean(tail(df[[5]][[i]]$MSEMD,1)),
                        length(df[[1]][[i]]$Loglik)/4,length(df[[2]][[i]]$Loglik)/4,
                        length(df[[3]][[i]]$Loglik)/4,
                        length(df[[4]][[i]]$Loglik)/4,length(df[[5]][[i]]$Loglik)/4))
  #tail(df[[2]][[i]]$OmegaIndex)
  #tail(df[[1]][[i]]$Acc)
  #tail(df[[3]][[i]]$Acc)

  # p <- data.frame(seq = 1:length(df[[1]][[i]]$OmegaIndex) ,
  #            OIReg = df[[1]][[i]]$OmegaIndex)  %>%
  #   ggplot() +
  #   geom_line(aes(x = seq, y = OIReg), color = "green") +
  #   geom_line(data  = data.frame(seq = 1:length(df[[3]][[i]]$OmegaIndex) ,
  #                                OISto = df[[3]][[i]]$OmegaIndex),
  #             aes(x = seq, y = OISto), color = "blue") +
  #   geom_line(data  =data.frame(seq = 1:length(df[[2]][[i]]$OmegaIndex) ,
  #                             OINoCov = df[[2]][[i]]$OmegaIndex),
  #             aes(x = seq, y = OINoCov))
  # print(p)
  # print(data.frame(seq = 1:length(df[[1]][[i]]$Loglik[,1]) ,
  #                  OIReg = df[[1]][[i]]$Loglik[,1])  %>%
  #         ggplot() +
  #         geom_line(aes(x = seq, y = OIReg), color = "green") +
  #         geom_line(data  = data.frame(seq = 1:length(df[[3]][[i]]$Loglik[,1]) ,
  #                                      OISto = df[[3]][[i]]$Loglik[,1]),
  #                   aes(x = seq, y = OISto), color = "blue") +
  #         geom_line(data  =data.frame(seq = 1:length(df[[2]][[i]]$Loglik[,1]) ,
  #                                     OINoCov = df[[2]][[i]]$Loglik[,1]),
  #                   aes(x = seq, y = OINoCov)))

  ## plotting the final value of community weights
  # cw <- as.data.frame(cbind(1:100, df[[1]][[i]]$Ffin, df[[3]][[i]]$Ffin,
  #                           df[[1]][[i]]$Hfin, df[[3]][[i]]$Hfin,
  #                           #df[[1]][[i]]$Forig, df[[1]][[i]]$Horig,
  #                           df[[7]][[i]])
  #                     )
  # colnames(cw) <- c("n","f-r-1","f-r-2","f-r-3","f-s-1","f-s-2","f-s-3",
  #                   "h-r-1","h-r-2","h-r-3","h-s-1","h-s-2","h-s-3",
  #                   #"f-n-1","f-n-2","f-n-3","h-n-1","h-n-2","h-n-3",
  #                   "f-o-1","f-o-2","f-o-3","h-o-1","h-o-2","h-o-3")
  # print( cw %>% pivot_longer(!n, values_to = "cw", names_to = "type")%>%
  #          separate(type,c ("dir","RegType", "col") ,sep = "-") %>%
  #          ggplot()+
  #          geom_point(aes(x = n, y = cw , color = RegType), alpha =0.5)+
  #          facet_wrap(~dir+col))
}

colnames(ends) <- c("OI_R","OI_S","OI_N", "OI_M","OI_PR","LL_R","LL_S","LL_N","LL_M","LL_PR","MSE_R","MSE_S","MSE_PR","Len_R","Len_N","Len_S","Len_M","Len_PR")
ends <- data.frame(ends)
ends$type <- c(rep(c("D","U"), 6))
ends$NW <- rep(c(1:2), each =6)

ends <- data.frame(ends) %>%
  pivot_longer(-c(type,NW))#, values_to = "Values", names_to = "type")

ends %>%
  filter(name %in% c("OI_R","OI_S","OI_N", "OI_M","OI_PR")) %>%
  group_by(NW,type) %>%
  aggregate(value~ . , FUN =mean) %>%
  ggplot() +
  geom_boxplot(aes(y = value, x = name)) +
  geom_point(aes(x = name, y = value, color = factor(NW))) +
  facet_wrap(~type)

ends %>%
  filter(name %in% c("OI_R","OI_S","OI_N", "OI_M","OI_PR")) %>%
  group_by(NW,type) %>%
  #aggregate(value~ . , FUN =mean) %>%
  ggplot() +
  geom_boxplot(aes(y = value, x = name)) +
  geom_point(aes(x = name, y = value, color = factor(NW))) +
  facet_wrap(~type)


ends %>%
  filter(name %in% c("MSE_R","MSE_S","MSE_PR")) %>%
  group_by(NW,type) %>%
  aggregate(value~ . , FUN =mean) %>%
  ggplot()+
  geom_boxplot(aes(y = value, x = name))+
  geom_point(aes(x = name, y = value, color = factor(NW)))+
  facet_wrap(~type)

ends %>%
  filter(name %in% c("Len_R","Len_N","Len_S","Len_M","Len_PR")) %>%
  group_by(NW,type) %>%
  #aggregate(value~ . , FUN =mean) %>%
  ggplot() +
  geom_boxplot(aes(y = value, x = name)) +
  geom_point(aes(x = name, y = value, color = factor(NW))) +
  facet_wrap(~type)

# lst2 <- SaveCoDASim(simvec = 2, 
#                     sim = TRUE, 
#                     InitParamFile = "/Code/InitParamMiss_Cohesive_MCAR.rds")

# lst2 <- SaveCoDASim(simvec = 1, 
#                     sim = TRUE, 
#                     InitParamFile = "/Code/InitParamMiss_Nested.rds")
# 
# lst3 <- SaveCoDASim(simvec = 1, 
#                     sim = TRUE, 
#                     InitParamFile = "/Code/InitParamMiss _2-Mode.rds")

#saveRDS(lst, paste0(file,"_OP", simvec,".rds"))
