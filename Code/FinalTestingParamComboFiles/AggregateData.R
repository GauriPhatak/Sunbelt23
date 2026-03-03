library(tidyverse)
library(dplyr)
library(igraph)
library(png)
library(RColorBrewer)
library(gridExtra)
library(ggpubr)
library(ecr)
library(nsga2R)
source("C:/Users/gauph/Documents/StatisticsMS_PhD/Wastewater-Surveillance-OSU/Sunbelt23/Code/CoDACov.R")

## Reading combination of parameters
filePath <- "C:/Users/gauph/Documents/StatisticsMS_PhD/Wastewater-Surveillance-OSU/Sunbelt23/Code/FinalTestingParamComboFiles/4clustersims_comb_MultipleStarts.rds"
S <- readRDS(filePath)
S$penalty <- as.character(S$penalty)
## Reading names of files
dataFolder <- "C:/Users/gauph/Documents/StatisticsMS_PhD/Wastewater-Surveillance-OSU/Sunbelt23/Code/FinalTestingParamComboFiles/ResultsFolder/MultiInitializationScaledFullOP/"
files <- list.files(dataFolder)
N <- S$N[1]
nsimV <- S$Nsim[1]
epsilon <- 0.00001
delta <- getDelta(N,epsilon)

params <- matrix(0, nrow = 0, ncol = 9)
#zerocoef <- matrix(0,nrow = 0, ncol = 3)
#numNodesPerCluster <- matrix(0, nrow = 0, ncol = 8)
conduct <- matrix(0, nrow = 0, ncol = 8)

## setting up output for covariate simulations
cn <- c("MinMeanConductanceW","MaxMeanConductanceW","SdMeanConductanceW", "MeanConductanceW","MeanConductanceWo","WeightedMeanConductanceW","WeightedMeanConductanceWo","ConductanceWDensityPen",
        "tprW","tprWo","tprWW","tprWoW",
        "AvgVarianceW","DispersionScoreW","AvgVarianceWo","DispersionScoreWo","WeightedAvgVarianceW","WeightedDispersionScoreW",
        "AvgIntraClusterSimilarityW","AvgIntraClusterSimilarityWo","WeightedAvgIntraClusterSimilarityW","WeightedAvgIntraClusterSimilarityWo",
        "CommInternalDensityW","CommInternalDensityWo","TotalInternalDensity","WeightedInternalDensityW","WeightedInternalDensityWo",
        "SilhouetteLin","SilhouetteBin","OI","AvgBinarySpreadW","BinaryDispersionScoreW","WeightedAvgBinarySpreadW",
        "WeightedBinaryDispersionScoreW","Trace","WeightedTrace","TotalTrace", "R2","sIdx","bigN","nsim","i","effective_nc",
        "MSE1","MSE2","MSE3","Acc1","Acc2","Acc3","NumIters", "unassigned","degree0","degree1","BIC","ICL","LL","nc2","Nnodes","Nedges",
        "dir","groupSize", "Failure","Filename")
metricsCov <- matrix(0, nrow =0,ncol = length(cn))
cov_opCols <- c("Dir","effective_nc","bigN","nsim","BIC","ICL","LL","nc2","Nnodes","Nedges",
                "OI","MSE1","MSE2","MSE3", "Acc1","Acc2","Acc3","NumIters", "sIdx", "newOI","unassigned","degree0","degree1", "i")
cov_op <- matrix(0, nrow = 0, ncol = length(cov_opCols))

## setting up output for no covariate values
nocovCols <- c("Dir","effective_nc","bigN","nsim","BIC","ICL","LL","nc2", "Nnodes","Nedges",
               "NumIters","OI", "sIdx", "newOI","unassigned","degree0","degree1","i")
nocov_op <- matrix(0, nrow = 0, ncol = length(nocovCols))

cn_nc <- c("MinMeanConductanceW","MaxMeanConductanceW","SdMeanConductanceW", "MeanConductanceW","MeanConductanceWo","WeightedMeanConductanceW","WeightedMeanConductanceWo"," ConductanceWDensityPen",
           "tprW","tprWo","tprWW","tprWoW",
           "CommInternalDensityW","CommInternalDensityWo","TotalInternalDensity","WeightedInternalDensityW","WeightedInternalDensityWo",
           "OI","sIdx","bigN","nsim","i","effective_nc","NumIters","unassigned","degree0","degree1","BIC","ICL","LL","nc2", "Nnodes","Nedges", "dir","groupSize","Failure")
metricsNoCov <- matrix(0,nrow =0, ncol = length(cn_nc)) 
## find the number of nodes without assignmment

numNodesWoAssignment <- c()
nc_sim <- 4

##replacing overlap column with:
OL_combo <- as.data.frame(cbind(prop = c(0,0.05,0.1,0.15,0.2,0.25), OL = c("NoOL","5Pct","10Pct","15Pct","20Pct","25Pct")))
OL_combo$prop <- as.numeric(OL_combo$prop)

S$sIdx <- 1:nrow(S)
print(paste0("The number of files: ", length(files)))
total <- length(files)
range_size = 10
seqID <- 1 #as.numeric(args[1])
starts <- seq(1, total, by = range_size)
ends <- pmin(starts + range_size - 1, total)
print(paste0("Start: ",starts[seqID], " ends: ", ends[seqID]))
itr <- 0

for(fn in files[91:102]){#files[starts[seqID]:ends[seqID]]){
  itr <- itr+1
  print(paste0(fn, " seq #: ", itr))
  df_rds <- readRDS(paste0(dataFolder,fn))
  sIdx <- as.numeric(str_extract_all(fn, "\\d+")[[1]])[1]
  
  alphaV <- df_rds[[8]]$V2[1][[1]]
  lambdaV_bin <- df_rds[[8]]$V5[1][[1]]
  lambdaV_lin <- df_rds[[8]]$V6[1][[1]]
  pClustOLV <- df_rds[[8]]$V26[[1]]
  penaltyV <- as.character(df_rds[[8]]$V30[1][[1]])
  pctMiss <- ifelse(is.null(unlist(df_rds[[8]]$V11[1])[1]) , 0 , unlist(df_rds[[8]]$V11[1])[1])#S$missing[sIdx][[1]]#
  params <- rbind(params, c(sIdx, alphaV,lambdaV_bin, lambdaV_lin,pClustOLV,penaltyV,
                            pctMiss,length(df_rds[[1]]),length(df_rds[[2]])))
  i <- 1
  
  for(bigN in 1:S$bigN[1]){
    G_orig <-  df_rds[[3]][[bigN]]
    ## Original binary covariates
    #Xout_cov <- as.data.frame(vertex_attr(G_orig))[,5:7]#df_rds[[7]][[bigN]][[1]]
    for(nsim in 1:nsimV){
      for(dir in c("directed","undirected")){
        for(nc in S$test_nc[[1]]){
          
          if(dir == "directed"){
            G <- G_orig
          }else{
            G <- convertToUndir(G_orig,nc_sim)
          }
          
          ## find the number of nodes with degree 0 or 1
          degree0 <- sum(igraph::degree(G) == 0 )
          degree1 <- sum(igraph::degree(G) == 1 )
          degree01 <- degree0 + degree1
          
          ## Pct missing data
          if(length(df_rds[[1]]) > 0){
            if(length(df_rds[[1]][[i]]) > 0){
              if(!is.null(df_rds[[1]][[i]])){
                d <- df_rds[[1]][[i]]
                bicOP <- d$bic
                numIter <- dim(d$Loglik)[1]
                Fm <- d$Ffin
                Hm <- d$Hfin
                FAILURE <- d[[22]]
                if(FAILURE == FALSE){
                  
                  ## find the number of nodes without assigned cluster
                  C <- memOverlapCalc(Fm,Hm,delta, N, nc)
                  Z <- d$Zout_cov
                  X <- d$Xout_cov
                  
                  
                  ## if a community consists of more than 80% of the nodes drop it
                  idx <- which(colSums(C)/N > 0.8 )
                  if(length(idx) > 0 ){
                    C <- as.matrix(C[, -idx], ncol = nc - length(idx), nrow = N)
                    nc <- nc - length(idx)
                    #print("Done dropping large communities")
                  }
                  ## drop tiny communities
                  if(nc > 0){
                    drop_cols <- DropTinyCommunities(degree01, C)
                    if(length(drop_cols)  == 2){
                      C <- drop_cols[[1]]
                      nc <- drop_cols[[2]]
                      #print("Done dropping tiny communities")
                      #print(drop_cols)
                    }
                  }
                  
                  ## drop highly overlapping or contained clusters
                  maxOL <- 0.70
                  if(nc > 0){
                    drop_cols <- Percent_overlap(C, maxOL) #DropTinyCommunities(degree01, C)
                    if(length(drop_cols)  == 2){
                      C <- drop_cols[[1]]
                      nc <- drop_cols[[2]]
                      #print("done merging similar communities")
                    }
                  }
                  effective_nc <- nc
                  if(nc > 0){
                    ## Calculate new OI
                    newOI   <- OmegaIdx(G,C, N, nc, nc_sim)
                    
                    ## function to create community of background nodes
                    
                    ##Create a new community of background nodes that as unassigned
                    NodesWoAssignment <- (rowSums(C) == 0)
                    ##punish internal density if nodes have no assignment
                    numNodesWoAssignment <- sum(NodesWoAssignment)
                    
                    if(sum(NodesWoAssignment) > 0){
                      C <- cbind(C, new_col = as.numeric(NodesWoAssignment))
                      nc <- nc+1
                    }
                    
                    if(ncol(C) > 1){
                      ## Calculating the metrics
                      Silhouette_lincov <- SilhouetteScore(C, Z,"Euclidean")
                      Silhouette_bincov <- SilhouetteScore(C, X, "Jaccard")
                      conduct <- EgoSplitConductance(G , C, dir, degree01, N, effective_nc, numNodesWoAssignment)
                      TPR     <- Comm_TPR(G, C , N, effective_nc, dir)
                      ADS     <- AverageDissimilarityScore(Z, C, degree01, numNodesWoAssignment)
                      AS      <- AverageSimilarityScore(Z, C, degree01, numNodesWoAssignment)
                      ID      <- InternalDensity(G, d, epsilon, dir)
                      
                      
                      #newOI <- newOI[1]
                      BinaryDispersion <- BinaryDispersionScore(X, C, degree01, numNodesWoAssignment)
                      traceVal <- WeightedTrace(d, epsilon)
                      
                      ##find the size of clusters
                      gs <- groupSize(C) 
                      
                      ## finding min, max and sd of mean conductance
                      #indices_to_remove <- (length(conduct) - 3):length(conduct)
                      conduct_sbst <- conduct[1:effective_nc]
                      MinMeanConductanceW <- min(conduct_sbst)
                      MaxMeanConductanceW <- max(conduct_sbst)
                      SdMeanConductanceW <- sd(conduct_sbst)
                      
                      metricsCov <- rbind(metricsCov, c(MinMeanConductanceW,MaxMeanConductanceW,SdMeanConductanceW,
                                                        tail(conduct,5),tail(TPR,4),ADS,AS,tail(ID,5), 
                                                        Silhouette_lincov,Silhouette_bincov,newOI, 
                                                        BinaryDispersion,traceVal, sIdx, bigN, nsim, i, 
                                                        effective_nc,
                                                        d$MSE,
                                                        d$Acc,
                                                        numIter,
                                                        numNodesWoAssignment,
                                                        degree0,
                                                        degree1, bicOP,dir, gs, FAILURE, fn))
                      
                      cov_op <- rbind(cov_op, unlist(c(dir, effective_nc, bigN, nsim,
                                                       bicOP,
                                                       d$OmegaIndex,
                                                       d$MSE,
                                                       d$Acc,
                                                       numIter,
                                                       sIdx,
                                                       newOI,
                                                       numNodesWoAssignment,
                                                       degree0,
                                                       degree1,
                                                       i )))
                      
                    }
                  }
                  
                }else{print(paste0("Failed at ", sIdx, " at ",i))}
              }
            }
          }
          
          if(length(df_rds[[2]]) > 0){
            if(!is.null(df_rds[[2]][[i]])){
              dn <- df_rds[[2]][[i]]
              Fm <- dn$Ffin
              Hm <- dn$Hfin
              FAILURE <- dn[[22]]
              numIter <- dim(dn$Loglik)[1]
              C <- memOverlapCalc(Fm,Hm,delta, N, nc)
              
              if(FAILURE == FALSE){
                ## if a community consists of more than 90% of the nodes drop it
                idx <- which(colSums(C)/N > 0.8 )
                if(length(idx) > 0 ){
                  C <- C[, -idx]
                  nc <- nc - length(idx)
                }
                
                ## drop tiny communities
                if(nc > 1){
                  drop_cols <- DropTinyCommunities(degree01, C)
                  if(length(drop_cols)  == 2){
                    C <- drop_cols[[1]]
                    nc <- drop_cols[[2]]
                  }
                }
                
                ## drop highly overlapping or contained clusters
                maxOL <- 0.70
                if(nc > 0){
                  drop_cols <- Percent_overlap(C, maxOL) #DropTinyCommunities(degree01, C)
                  if(length(drop_cols)  == 2){
                    C <- drop_cols[[1]]
                    nc <- drop_cols[[2]]
                  }
                }
                
                effective_nc <- nc
                if(nc > 0){
                  ## Calculate new OI
                  newOI   <- OmegaIdx(G,C, N, nc, nc_sim)
                  
                  ## function to create community of background nodes
                  
                  ##Create a new community of background nodes that as unassigned
                  NodesWoAssignment <- (rowSums(C) == 0)
                  ##punish internal density if nodes have no assignment
                  numNodesWoAssignment <- sum(NodesWoAssignment)
                  
                  if(sum(NodesWoAssignment) > 0){
                    C <- cbind(C, new_col = as.numeric(NodesWoAssignment))
                    nc <- nc+1
                  }
                  
                  ## Calculating the metrics
                  conduct <- EgoSplitConductance(G , C, dir, degree01, N, effective_nc, numNodesWoAssignment)
                  
                  TPR     <- Comm_TPR(G, C , N, effective_nc, dir)
                  ID      <- InternalDensity(G, dn, epsilon, dir)
                  
                  gs <- groupSize(C)  
                  
                  ## finding min, max and sd of mean conductance
                  #indices_to_remove <- (length(conduct) - 3):length(conduct)
                  conduct_sbst <- conduct[1:nc]
                  MinMeanConductanceW <- min(conduct_sbst)
                  MaxMeanConductanceW <- max(conduct_sbst)
                  SdMeanConductanceW <- sd(conduct_sbst)
                  
                  metricsNoCov <- rbind(metricsNoCov, c(MinMeanConductanceW,MaxMeanConductanceW,SdMeanConductanceW,
                                                        tail(conduct,5),tail(TPR,4),tail(ID,5),newOI, sIdx, bigN, nsim,
                                                         i,effective_nc, 
                                                        numIter,
                                                        numNodesWoAssignment,
                                                        degree0,
                                                        degree1,
                                                        dn$bic,dir,
                                                        gs,FAILURE))
                  
                  
                  nocov_op <- rbind(nocov_op, unlist(c(dir, effective_nc, bigN, nsim, dn$bic,numIter,
                                                       dn$OmegaIndex,sIdx,newOI,numNodesWoAssignment,
                                                       degree0,degree1,i)))
                }
              }
            }
          }
          
          i <- i + 1
        }
      }
    }
  }
}

params <- t(as.data.frame(apply(params,1,unlist)))[,c(1:5,12:15)]
colnames(params) <- c("sIdx", "alpha","lambda_bin","lambda_lin","prop","penalty","PctMiss","lenCov","lenNoCov")
params <- as.data.frame(cbind(apply(params[,c(1:5,8:9)],2,as.numeric), params[,c(6,7)]))
params$prop <- as.numeric(params$prop)
params$sIdx <- as.numeric(params$sIdx)
params <- left_join(params, OL_combo )
params <- params %>% select(-prop)
params <- params %>% arrange(sIdx)
#params$penalty <- rep(c("LASSO","Ridge"), each = 72)
#params <- as.data.frame(cbind(params, ifelse(params[,5] == 0,"NoOL","2GrpOL1")))[,c(1:4,6:10)]
if(dim(cov_op)[1] > 0){
  
  cov_op <- as.data.frame(cov_op)
  colnames(cov_op) <- cov_opCols
  cov_op[,c(2:24)] <- as.data.frame(apply(cov_op[,c(2:24)],2,as.numeric))
  cov_op$Dir <- unlist(cov_op$Dir)
  cov_op <- left_join(cov_op, params)
  
  ## Recalculating the BIc cause the formula chnaged?
  cov_op <-  cov_op %>%
    dplyr::mutate(MSEavg = rowMeans(select(., MSE1, MSE2, MSE3), na.rm = TRUE),
                  Accavg = rowMeans(select(., Acc1, Acc2, Acc3), na.rm = TRUE))
  
  ## adding column names to metrics 
  metricsCov <- as.data.frame(metricsCov)
  colnames(metricsCov) <- cn
  
  metricsCov[,1:59] <- as.data.frame(apply(metricsCov[,1:59],2,as.numeric))
  metricsCov <- left_join(metricsCov, params)
  
  #metricsCov$OL <- cov_op$OL
  #metricsCov$alpha <- cov_op$alpha
  #metricsCov$penalty <- cov_op$penalty
  #metricsCov$MSEavg <- cov_op$MSEavg
  #metricsCov$Accavg <- cov_op$Accavg
  #metricsCov$lambda_bin <- cov_op$lambda_bin
  #metricsCov$lambda_lin <- cov_op$lambda_lin
  #metricsCov$degree0 <- cov_op$degree0
  #metricsCov$degree1 <- cov_op$degree1
  #metricsCov$degree01 <- metricsCov$degree0+metricsCov$degree1
  #metricsCov$unassigned <- cov_op$unassigned
  #metricsCov$pctMiss <- cov_op$PctMiss
  #metricsCov$NumIters <- cov_op$NumIters
  #metricsCov$LL <- cov_op$LL
  metricsCov <- metricsCov %>% arrange(sIdx)
  
  ## Saving the files
  file <- "C:/Users/gauph/Documents/StatisticsMS_PhD/Wastewater-Surveillance-OSU/Sunbelt23/Code/FinalTestingParamComboFiles/ResultsFolder/SubsetResults/"
  ## saving the data in files
  saveRDS(cov_op, paste0(file,"MultInitializationScale_covOP",as.character(seqID),".rds"))
  saveRDS(metricsCov, paste0(file,"MultInitializationScale_metricsCov",as.character(seqID),".rds"))
  
}

if(dim(nocov_op)[1] > 0){
  cov_op <- nocov_op
  metricsCov <- metricsNoCov
  cov_op <- as.data.frame(cov_op)
  colnames(cov_op) <- nocovCols#cov_opCols
  cov_op[,c(2:18)] <- as.data.frame(apply(cov_op[,c(2:18)],2,as.numeric))
  cov_op$Dir <- unlist(cov_op$Dir)
  cov_op <- left_join(cov_op, params)
  
  ## adding column names to metrics 
  metricsCov <- as.data.frame(metricsCov)
  colnames(metricsCov) <- cn_nc
  
  metricsCov[,1:33] <- as.data.frame(apply(metricsCov[,1:33],2,as.numeric))
  metricsCov <- left_join(metricsCov, params)
  metricsCov$degree01 <- metricsCov$degree0+metricsCov$degree1
  
  # metricsCov$OL <- cov_op$OL
  # metricsCov$alpha <- cov_op$alpha
  # metricsCov$penalty <- cov_op$penalty
  # metricsCov$degree0 <- cov_op$degree0
  # metricsCov$degree1 <- cov_op$degree1
  # metricsCov$degree01 <- metricsCov$degree0+metricsCov$degree1
  # metricsCov$unassigned <- cov_op$unassigned
  # metricsCov$NumIters <- cov_op$NumIters
  # metricsCov$LL <- cov_op$LL
  metricsCov <- metricsCov %>% arrange(sIdx)
  
  ## Saving the files
  file <- "C:/Users/gauph/Documents/StatisticsMS_PhD/Wastewater-Surveillance-OSU/Sunbelt23/Code/FinalTestingParamComboFiles/ResultsFolder/SubsetResults/"
  ## saving the data in files
  saveRDS(cov_op, paste0(file,"MultInitializationScale_nocovOP",as.character(seqID),".rds"))
  saveRDS(metricsCov, paste0(file,"MultInitializationScale_metricsnoCov",as.character(seqID),".rds"))
  
}

