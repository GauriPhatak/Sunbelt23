
source("/nfs/stak/users/phatakg/ResearchCode/Sunbelt23/Code/HelperFuncs.R")
args=commandArgs(trailingOnly = TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}
print(args)

Impute = TRUE
w =gsub("[\r\n]", "", as.character(args[2]))#"Simul"
print(paste0("The simulation type ",w))
## Number of clusters
k <- 3
## Define the number of nodes 
N <- 200
## Define the number of categories we are interested in
ncat = 5
cat <- c("Dachshund","Saluki","Labrador","FlatCoat","LhasaApso") 
#"BorderCollie","Aussie","Pitbull","JapaneseSpitz",
# "Hound","Hound","Retriever","Retriever","Terrier","Herding","Herding","Terrier","Spitz",
origcat <- cat

clusters <- list(c(0.3,0.3,0.3), c(0.5,0.3,0.2))#c("equal", "unequal")
clustEdgeP <- list(c(0.5,0.05,0.5,0.05,0.05,0.5), c(0.2,0.05,0.2,0.05,0.05,0.2))#c("high", "low")
catProb <- list(list(c(rep(0.6,2),rep(0.02,3)), # Cluster 1
                     c(rep(0.02,2),rep(0.6, 2),rep(0.02,1)), # Cluster 2
                     c(rep(0.02,4),rep(0.6,1)))) # Cluster 3#c("structured")
#          list(list(c(rep(0.6,1),rep(0.02,2)), # Cluster 1
#                    c(rep(0.02,1),rep(0.6, 1),rep(0.02,1)), # Cluster 2
#                    c(rep(0.02,2),rep(0.6,1)))) # Cluster 3#c("structured")
           
catEdgeP <- list(c(0.5,0.02),c(0.1,0.02))#c("high","low")
normCoef <- 0.51#c("low")
normCenters <- list(list(mean = c(15,10,5), var = c(1,1,1)), list(mean = c(7,6,5), var = c(1,1,1)))#c("close","far")

## Setting values for cluster dependent missingness
#Missing type
mt <- "Cluster"#"Random" #
propmsng <- list(c(0.3,0.3,0.3),
                 c(0.25,0.25,0.2,0.2,0.1),#c(0.4,0.3,0.3), 
                 c(0.2,0.2,0.6))

simgrid <- expand.grid(clusters = clusters, ClustEdgeP = clustEdgeP, 
                       catProb = catProb, catEdgeP = catEdgeP, 
                       normCoef = normCoef, normCenters = normCenters)

simvec <- as.integer(gsub("[\r\n]", "",args[1]))
print(paste0("the simulation number ", simvec))
for(sim in simvec){
  simnum <- simgrid[sim,]

  ## Setting distributions for the continous covariate value
  p = unlist(simnum[[1]])
  ## Define the probabilty of connection i.e. edge between nodes belonging to different/same cluster.
  B = unlist(simnum[[2]])
  ## Assigning probability of a categorical level for a  particular cluster being assigned to the node. 
  pC <- simnum[[3]][[1]]
  ## Probability of connecting based on category assignment. Keeping it simple for now. If the nodes belong to the same category then the probability is higher. If not it is lower.
  CtProb <- unlist(simnum[[4]])
  ## Need to assign coefficient value to the continuous variable
  normCoef = unlist(simnum[[5]])
  ## Setting distributions for the continuous covariate value
  pCont <- data.frame(cbind(mean = c(simnum[[6]][[1]]$mean), var = c(simnum[[6]][[1]]$var)))

  ## Creating the dataframe with category combinations
  combi <- data.frame(rbind(t(combn(sort(cat, decreasing =FALSE),2, simplify = TRUE)), cbind(cat,cat)))

  combi$prob <- CtProb[2]
  for(i in 1:length(combi$cat)){
    if(combi[i,1] == combi[i,2]){
      combi[i,3]= CtProb[1]
    }
  }
  colnames(combi) <- c("cat1","cat2","prob")
  combi <- combi[with(combi, order(cat2, cat1)),]
  pCat <- data.frame(cat =cat,clust = rep(c("1","2","3"), each = length(cat)), pCat  = unlist(pC)) 
  pCat <- pCat[with(pCat, order(clust, cat)),]


  superiter <- 100
  niter <- 7
  pctMsng <- c(15,25,50,70)
  op <- data.frame(matrix(nrow=0,ncol = 16))
  comDetVals <- data.frame(matrix(nrow = 0, ncol = 5))
  CovARI <- data.frame(matrix(nrow = 0, ncol = 3))
  edgeDF <- data.frame(matrix(nrow= 0 ,ncol = 6))
  degVal <- data.frame(matrix(nrow = 0, ncol = 7))
  IC <- data.frame(matrix(nrow = 0, ncol = 7))
  clust9 <- data.frame(matrix(nrow = 0, ncol = 3))

  for(val in 1:superiter){
      
    print(paste0(sim,"-",val))
    ## Create random assignments for assigning clusters with probability vector p
    C <- 0
    while(length(table(C)) < k){
        C = sample(1:k,N,replace = TRUE, prob = p)
    }

    ## Create an empty network i.e. no cluster assignments or category assignments or edges
    net <- network(N, directed = FALSE, density= 0 )
    
    ## Assign clusters to the nodes 
    net %v% 'Cluster' <- C
    
    ## Based on the probability matrix pC assign cats to nodes
    cats <- NA
    normVal <- NA
    for (i in 1:k) {
      cats[which(as.numeric(net %v% 'Cluster') == i)] <- sample(x = origcat,
                                                              size = length(which(as.numeric(net %v% 'Cluster') == i)), 
                                                              prob = pC[[i]],
                                                              replace = TRUE)
      normVal[which(as.numeric(net %v% 'Cluster') == i)]  <- rnorm(n = length(which(as.numeric(net %v% 'Cluster') == i)), 
                                                              mean = pCont[i,1], 
                                                              sd = pCont[i,2])                                                  
    }
    
    net %v% "Categories" <- cats
    net %v% "Continuous" <- normVal

    ## change name to g.sim for some nonsense reason!
    g.sim <- net
    
    ## We would like to simulate networks using different ergm terms. We can compare combination of both to just one of them.
    ## Use the probability of connection values defined earlier for both cluster and category for fitting the data. 
    g.sim_C_B <- simulate(g.sim ~ nodemix("Cluster", levels = TRUE, levels2 = TRUE) + 
                                  nodemix(c("Categories"), levels=TRUE, levels2 = TRUE) + 
                                  nodecov("Continuous"),
                          nsim = 1,
                          coef = prob2logit(c(B, combi$prob, normCoef)),
                          control=control.simulate( MCMC.burnin=10000, MCMC.interval=1000))

    ## This function return a table of number of edges between the different books
    g.sim_C_B <- asIgraph(g.sim_C_B)

    ## List all the network characteristics 
    ## create a dataframe of all the node characterictics.
    ## Cluster assignment, book assignment. Add groupby.
    g <- g.sim_C_B
    ends <- as_edgelist(g)
    e1 <- V(g)[ends[,1]]$Categories
    e2 <- V(g)[ends[,2]]$Categories
    c1 <- V(g)[ends[,1]]$Cluster
    c2 <- V(g)[ends[,2]]$Cluster
    endcolors <- t(apply(cbind(e1, e2, c1,c2), 1, sort))
    CntEdgesop <- table(endcolors[,1], endcolors[,2], endcolors[,3],endcolors[,4])

    M <- as.data.frame(CntEdgesop)
    edgeDF <- rbind(edgeDF,cbind(M[M$Freq !=0,],val))
    
    ## This function find the empirical mean and variance of the degrees for different titles.
    degVal <- rbind(degVal,
                    c(data.frame(degree = igraph::degree(g.sim_C_B), 
                                cat = V(g.sim_C_B)$Categories,
                                normVal = V(g.sim_C_B)$Continuous,
                                cluster = V(g.sim_C_B)$Cluster) %>%
                        group_by(cat,cluster) %>%
                        summarise(meanDeg = mean(degree),
                                  minDeg = min(degree),
                                  maxDeg = max(degree),
                                  varDeg = var(degree),
                                  meanN = mean(normVal),
                                  minN = min(normVal),
                                  maxN = max(normVal),
                                  varN = var(normVal),
                                  count = n(),
                                  val)))
    
    g <- asNetwork(g.sim_C_B)
    fit_C_B_N <- ergm(g ~ edges() + nodemix(c("Cluster"),levels = TRUE,levels2 = TRUE) +
                      nodemix(c("Categories"),levels = TRUE,levels2 = TRUE)+
                      nodecov("Continuous"))
    fit_C <- ergm(g ~ edges()+nodemix(c("Cluster"),levels = TRUE,levels2 = TRUE))
    fit_B <- ergm(g ~ edges()+nodemix(c("Categories"),levels = TRUE,levels2 = TRUE))
    fit_N <- ergm(g ~ edges()+nodecov(c("Continuous")))
    
    ## creating dataframe for finding AIC / BIC values
    
    IC <- rbind(IC,data.frame(AIC_B = AIC(fit_B),
                              AIC_C = AIC(fit_C),
                              AIC_C_B_N = AIC(fit_C_B_N),
                              AIC_N = AIC(fit_N),
                              BIC_B = BIC(fit_B),
                              BIC_C = BIC(fit_C),
                              BIC_C_B_N = BIC(fit_C_B_N),
                              BIC_N = BIC(fit_N),
                              iter = val))

    ### Collecting data from Covariate assisted regularized spectral clustering and corresponding alpha values
    X <- as.matrix(V(g.sim_C_B)$Categories)
    Xdum <- as.matrix(data.frame(predict(dummyVars("~.", data = X ), newdata <- X)))
    ## attach the continuous variable to the categories. 
    Xdum <- cbind(Xdum, as.numeric(V(g.sim_C_B)$Continuous))
    origmemCov <- CovAssistedSpecClust(G = g.sim_C_B, Xdum, 3,
                                      Regularize =TRUE, alpha =NA,
                                      type ="non-assortative",
                                      retAlpha = TRUE)
    CovARI <- rbind(CovARI, c(mclust::adjustedRandIndex(V(g.sim_C_B)$Cluster, origmemCov[[1]]),
                              origmemCov[[2]]))
    if(Impute ==TRUE){
      g.imp <- list()
      for(pt in pctMsng){
        if(w == "Qual"){
          covMsng <- data.frame(name = c("Categories") ,
                            p = pt,
                            cat = I(list(V(g.sim_C_B)$Categories)))
          Midx <- list()
          for (i in 1:dim(covMsng)[1]) {
            Midx[[i]] <- sample(1:N,ceiling(covMsng$p[i] *N/100))
            }
          coms <- comdet(g.sim_C_B, V(g.sim_C_B)$Categories)
          imptype = c("Mode","MultRegression")
          df <- IterativeImputeQual(g.sim_C_B, N, coms,niter,k,covMsng,Midx, imptype)
          g.imp[[i]] <- df[[2]] ###******MAYBE NEED TO ADD REGRESSION OUTPUT AS WELL******###
          # keeping track of community detection values for percentage missing
          comDetVals <- rbind(comDetVals,
                              c(val,pt,round(mclust::adjustedRandIndex(V(g.imp[[i]])$Cluster,V(g.imp[[i]])$com),3), 
                              round(mclust::adjustedRandIndex(V(g.imp[[i]])$Cluster,coms$origmem),3)))
          op <- rbind(op,cbind(df[[1]], pt,val,w))
        }

        if(w == "Cont"){
        ## Continuous variable imputation
          covMsng <- data.frame(name = c("Continuous"),
                                p = pt,
                                cat = unique(floor(V(g.sim_C_B)$Continuous )))
          a <- 1:N
          numSamp <- ceiling(pt*N/100)
          car <- caret::createDataPartition(as.factor(V(g.sim_C_B)$Categories), 
                                              p = 1-(pt/N), 
                                              list = FALSE)
          Midx <- a[-car]
          coms <- comdet(g.sim_C_B, V(g.sim_C_B)$Categories, as.numeric(V(g.sim_C_B)$Continuous))
          imptype = c("Mean","QuantRegression")
          df <- IterativeImputeCont(g.sim_C_B, N, coms,niter,k,covMsng,Midx, imptype)
          g.imp[[i]] <- df[[2]] ###******MAYBE NEED TO ADD REGRESSION OUTPUT AS WELL******###
          # keeping track of community detection values for percentage missing
          comDetVals <- rbind(comDetVals,
                              c(val,pt,round(mclust::adjustedRandIndex(V(g.imp[[i]])$Cluster,V(g.imp[[i]])$com),3), 
                              round(mclust::adjustedRandIndex(V(g.imp[[i]])$Cluster,coms$origmem),3)))
          op <- rbind(op,cbind(df[[1]], pt,val,w))
        }

        if(w == "Simul"){
          #print("In SimulImpute")
          # a <- 1:N
          # ## Setting missing data for categorical variable
          # car <- caret::createDataPartition(as.factor(V(g.sim_C_B)$Categories), 
          #                                     p = round(1-(pt/(2*N)),2), 
          #                                     list = FALSE)
          # MidxCon <- a[-car]
          # ## Setting missing data for Continuous variable. Here we are sampling from indicies that are not Categorical variable
          # MidxCat <- sample(car, ceiling((pt/(2*N))*100), replace = FALSE)
          if(mt == "Random"){
            Midx <- RandomMissingness(N,V(g.sim_C_B)$Categories,pt)
          } else if (mt == "Cluster") {
            Midx <- ClusterMissingness(propmsng[[2]], N, V(g.sim_C_B)$Categories, pt, C)
            #print(Midx)
          }
          MidxCon <- unlist(Midx[[1]])
          MidxCat <- unlist(Midx[[2]])
          coms <- comdet(g.sim_C_B, V(g.sim_C_B)$Categories, as.numeric(V(g.sim_C_B)$Continuous))
          imptype = c("SimulReg")
          df <- IterativeImputeSimul(g.sim_C_B, N, coms,niter,k,MidxCat,MidxCon, imptype)
          g.imp[[i]] <- df[[2]]
          ## keeping track of community detection values for percentage missing
          comDetVals <- rbind(comDetVals,
                              c(val,pt,round(mclust::adjustedRandIndex(V(g.imp[[i]])$Cluster,V(g.imp[[i]])$com),3), 
                              round(mclust::adjustedRandIndex(V(g.imp[[i]])$Cluster,coms$origmem),3)))
          op <- rbind(op,cbind(df[[1]], pt,val,w))
        }
        i <- i+1
      }
    }
    clust9 <- rbind(clust9 , cbind(RegSpectralClust(g.sim_C_B,9), V(g.sim_C_B)$Cluster, rep(val, N)))
  }

saveRDS(list(op,comDetVals,CovARI,edgeDF,degVal,IC),paste0(w,"U",sim,".rds"))
}