
library(tidyverse)
library(dplyr)
library(igraph)

##This reads the command line arguments
# args=commandArgs(trailingOnly = TRUE)
# if (length(args)==0) {
#   stop("At least one argument must be supplied (input file).n", call.=FALSE)
# }
# 
getwd()
source(paste0(getwd(),"/Code/CoDACov.R"))

## Read full graph
df <- readRDS(paste0(getwd(),"/Code/FinalTestingParamComboFiles/LargeSimulatedData/NWexamples200_MAR1.rds"))

## Read full graph
#df <- readRDS(paste0(args[2],"NWexamples200_MAR",args[1],".rds"))
k_out <- 3
k_in <- 0
k <- k_in+k_out

## Number of in and out continuous covariates
o_out <- 3
o_in <- 0
o <- o_in+o_out

alphaLL <- 1

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
test_nc <- c(4)
#c(4,6,8,10,13,15)
thresh  <- 0.00005
epsilon <- 0.000001
nitermax <- 30000
alphaV <- c(0.00001, 0.00005,0.0001, 0.0005,0.001,0.005,0.01)
test = FALSE
randomize = TRUE
missingV <- list(c(15,0,0,0,0,0),
                 c(15,15,0,0,0,0),
                 c(25,0,0,0,0,0),
                 c(15,0,0,0,0,0),
                 c(15,15,0,0,0,0),
                 c(25,0,0,0,0,0))
missing <- missingV[[as.numeric(1)]]
alphaLin <- 1#0.001
penaltyV <- c("LASSO")
seed <- sample(1:100000, 1)
lambda_linV  <- c(0.00001 ,0.0001, 0.001, 0.01)#,0.0001, 0.001, 0.01)
lambda_binV  <- c(0.00001 ,0.0001, 0.001, 0.01)#c(0.00001,0.0001, 0.001,0.01)
lambda_grphV <- c(0.0001)#, 0.0001)#,0.0001, 0.001,0.01, 0.1)
covInitOpt <- c("Nmean")
printFlg <- FALSE

## setting flags for running different comm det and nw metrics
nocov <- FALSE
reg <- TRUE

## setting up metrics dataframe
cn <- c("MeanConductanceW","MeanConductanceWo","WeightedMeanConductanceW","WeightedMeanConductanceWo",
        "tprW","tprWo","tprWW","tprWoW",
        "AvgVarianceW","DispersionScoreW","AvgVarianceWo","DispersionScoreWo","WeightedAvgVarianceW","WeightedDispersionScoreW",
        "AvgIntraClusterSimilarityW","AvgIntraClusterSimilarityWo","WeightedAvgIntraClusterSimilarityW","WeightedAvgIntraClusterSimilarityWo",
        "CommInternalDensityW","CommInternalDensityWo","TotalInternalDensity","WeightedInternalDensityW","WeightedInternalDensityWo",
        "Silhouette", "numNodesWoAssignment" ,"nc","dir", "alpha", "lambda","penalty","covInit","MissingCity")
metricsCov <- matrix(0, nrow =0,ncol = length(cn))

## Save output ##MSE, MSEmd, lambda, week, number of available, penalty, initializing of covariate, alpha 
finOP <- matrix(0, nrow =0, ncol = 12)

## Metric values
op_sim <- matrix(0, nrow= 0, ncol = 14)

## Saving the imputed values for the missing covid data
imputedV <- list()

## set up list for final community assignments
finMCov <- list()
finMNoCov <- list()

opf_Regcov <- list()
opf_noCov <- list()

covInit<- "Nmedian"
#nc_sim <- params$nc[1][[1]]
N <- 200
lvl <- 1
#IDs <- as.numeric(args[[1]])
NumNWs <- 1:15

for(ID in NumNWs){#:length(test_nc)){
  print(paste0("The network ID number: ", ID))
  NWlst <- df[[ID]]
  G_orig <- NWlst[[1]]
  orig_Cov <- NWlst[[2]]
  nc <- test_nc#[ID]
  nc_sim <- test_nc#[ID]
  ## Just to compare
  covtmp <-  cbind(as.data.frame(vertex_attr(G_orig)),
                   inDeg = igraph::degree(G_orig, mode = "in"), 
                   outDeg = igraph::degree(G_orig, mode = "out")) 
  
  F_u <- NWlst[[3]]
  H_u <- NWlst[[4]]
  
  for(dir in c("directed","undirected")){
    for(lambda_bin in lambda_binV){
      for(lambda_lin in lambda_linV){
        for(alpha in alphaV){
          for(lambda_grph in lambda_grphV){
            for(penalty in penaltyV){
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
              
              if(reg == TRUE){
                print(paste0("in the ", dir,
                             " alpha ",alpha," initial impute type ", covInit, 
                             " Penalty, ", penalty, " nc ", nc))
                tryCatch({
                  ## Algorithm with covariates + possible missing data
                  start <- Sys.time()
                  opf_Regcov[[lvl]] <- CoDA(G, nc, k = c(k_in, k_out), o = c(o_in, o_out), N, alpha,
                                            lambda_lin, lambda_bin, thresh, nitermax, orig, randomize,
                                            CovNamesLinin, CovNamesLinout, CovNamesLPin, CovNamesLPout, dir,
                                            alphaLL, test, missing, covOrig=orig_Cov,epsilon, 
                                            impType = "Reg", alphaLin, penalty, seed, covInit, specOP,nc_sim,lambda_grph )
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
                  opf_noCov[[lvl]] <- CoDA(G, nc, k = c(0, 0), o = c(0, 0), N, alpha, lambda_lin, lambda_bin, 
                                           thresh, nitermax, orig, randomize, CovNamesLinin = c(), 
                                           CovNamesLinout = c(), CovNamesLPin = c(), CovNamesLPout = c(), 
                                           dir, alphaLL, test, missing = missing, covOrig=orig_Cov, 
                                           epsilon = 0, impType = "", alphaLin, penalty, seed, covInit, specOP,nc_sim,lambda_grph)
                  tme <- Sys.time() - start
                  print(paste0("Total time take algo without covariates ", round(tme, 3)))
                }, error = function(e) {
                  # Handle the error
                  NA})              
              }
              
              #params <- rbind(params, append(Sim, list(bigN = m,nsim = j,dir = dir)))
              
              lvl <- lvl + 1
              gc()
            }
          }
        }
      }
    }
  }
}
saveRDS(opf_Regcov, paste0("/home/phatakg/novus/Output_rds/MARSimulation_",args[1],".rds"))  
