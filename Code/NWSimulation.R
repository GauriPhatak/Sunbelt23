source("HelperFuncs.R")


## Number of clusters
k <- 3

## Probability of a node belonging to a cluster. This defines the size of the clusters.
p = c(0.3,0.3,0.3)

## Define the number of nodes 
N <- 100

## Define the probabilty of connection i.e. edge between nodes belonging to different/same cluster.
B = c(0.2,0.05,0.2,0.05,0.05,0.2)

## Setting up a category for LOTR book universe
cat <- c("FOTR","TT","ROTK","Hob","Sil", "UT","CoH","BL","FoG")
## Assigning color to each book
catcol  <- c("#B22222","#FF7F50","#FF8C00",
             "#6B8E23","#7FFF00","#32CD32",
             "#6495ED","#00008B","#87CEEB")
## Assigning book to color
col <- data.frame(as.matrix(cbind(catcol, cat)))

## Assigning probability of a book from a particular cluster being assigned to the node. 
pC <- list(c(rep(0.5,3),rep(0.02,6)), # Cluster 1
           c(rep(0.02,3),rep(0.5, 3),rep(0.02,3)), # Cluster 2
           c(rep(0.02,6),rep(0.5,3))) # Cluster 3

## Creating book group. Books that have higer probability of connection.
c1 <- c("FOTR|TT|ROTK")
c2 <- c("Hob|Sil|UT")
c3 <- c("CoH|BL|FoG")
clusts <- c(c1,c2,c3)

## Here we assign probability of a connection between two nodes with particular book values. Higher probability of nodes within the same group. Lower for intergroup connection. 
## Define dataframe for mapping probability to book combinations
cat <- sort(cat, decreasing = TRUE)
combi <- data.frame(rbind(t(combn(cat,2, simplify = TRUE)), cbind(cat,cat)))
combi$prob <- 0
for(i in 1:length(combi$cat)){
  if( which(str_detect(combi[i,1], clusts) == TRUE,arr.ind = TRUE) == which(str_detect(combi[i,2], clusts) == TRUE,arr.ind = TRUE)){
    if(combi[i,1] == combi[i,2]){
      combi$prob[i] = 0.3
    }
    else{
      combi$prob[i] = 0.15
    }
  }
  else{
    combi$prob[i] = 0.02
  }
}
colnames(combi) <- c("cat1","cat2","prob")
combi <- combi[with(combi, order(cat1, cat2)),]

pCat <- data.frame(cat =cat,clust = rep(c("1","2","3"), each = length(cat)), pCat  = unlist(pC)) 
pCat <- pCat[with(pCat, order(clust, cat)),]

superiter <- 500
pctMsng <- c(5,7,10,20,30,40,50)
op <- data.frame(matrix(nrow=0,ncol = 9))

opLD <- data.frame(matrix(nrow =0 , ncol = 15))
comDetVals <- data.frame(matrix(nrow = 0, ncol = 5))
CovARI <- data.frame(matrix(nrow = 0, ncol = 3))
edgeDF <- data.frame(matrix(nrow= 0 ,ncol = 5))
degVal <- data.frame(matrix(nrow = 0, ncol = 6))
IC <- data.frame(matrix(nrow = 0, ncol = 7))

for(val in 1:superiter){
  print(val)
  ## Create random assignments for assigning clusters with probability vector p
  C = sample(1:k,N,replace = TRUE, prob = p)
  
  ## Create an empty network i.e. no cluster assignments or category assignments or edges
  net <- network(N, directed = FALSE, density= 0 )
  
  ## Assign clusters to the nodes 
  net %v% 'Cluster' <- C
  
  ## Based on the probability matrix pC assign books to nodes
  bk <- NA
  for (i in 1:k) {
    #bk <- net %v% "LOTR"
    bk[which(as.numeric(net %v% 'Cluster') == i)] <- sample(x = cat,
                                                            size = length(which(as.numeric(net %v% 'Cluster') == i)), 
                                                            prob = pC[[i]],
                                                            replace = TRUE)
  }
  
  net %v% "LOTR" <- bk
  
  ## Assigning book group 
  BookGrp <- NA
  for (i in 1:N) {
    BookGrp[i] <- which(str_detect((net %v% "LOTR")[i], clusts) == TRUE,arr.ind = TRUE)
  }
  net %v% "BookGrp" <- BookGrp
  
  
  ## change name to g.sim for some nonsense reason!
  g.sim <- net
  
  ## We would like to simulate networks using different ergm terms. We can compare combination of both to just one of them.
  ## Use the probability of connection values defined earlier for both cluster and category for fitting the data. 
  g.sim_C_B <- simulate(g.sim ~ nodemix("Cluster", levels = TRUE, levels2 = TRUE) +  nodemix(c("LOTR"), levels=TRUE, levels2 = TRUE),
                        nsim = 1,
                        coef = prob2logit(c(B, combi$prob)),
                        control=control.simulate( MCMC.burnin=10000, MCMC.interval=1000))
  fit_C_B <- ergm(g.sim_C_B ~ nodemix(c("Cluster"),levels = TRUE,levels2 = TRUE)+
                    nodemix(c("LOTR"),levels = TRUE,levels2 = TRUE))
  
  #theta <- logit2prob(coef(fit_C_B))
  
  
  ## This function return a table of number of edges between the different books
  g.sim_C_B <- asIgraph(g.sim_C_B)
  M <- CntEdges(g.sim_C_B)
  edgeDF <- rbind(edgeDF,
                  data.frame(NumEdges = as.numeric(M[!lower.tri(M)][order(row(M)[!lower.tri(row(M))])]),
                             Cat1 = rep(rownames(M), times = nrow(M))[!lower.tri(M)][order(row(M)[!lower.tri(row(M))])],
                             Cat2 = rep(rownames(M), each = nrow(M))[!lower.tri(M)][order(row(M)[!lower.tri(row(M))])],
                             val))
  
  
  ## This function find the empirical mean and variance of the degrees for different titles.
  degVal <- rbind(degVal,
                  c(data.frame(degree = igraph::degree(g.sim_C_B), 
                               cat = V(g.sim_C_B)$LOTR) %>%
                      group_by(cat) %>%
                      summarise(meanDeg = mean(degree),
                                min = min(degree),
                                max = max(degree),
                                var = var(degree),
                                val)))
  
  g <- asNetwork(g.sim_C_B)
  fit_C_B <- ergm(g ~ nodemix(c("Cluster"),levels = TRUE,levels2 = TRUE) +
                    nodemix(c("LOTR"),levels = TRUE,levels2 = TRUE))
  fit_C <- ergm(g ~ nodemix(c("Cluster"),levels = TRUE,levels2 = TRUE))
  fit_B <- ergm(g ~ nodemix(c("LOTR"),levels = TRUE,levels2 = TRUE))
  
  ## creating dataframe for finding AIC / BIC values
  
  IC <- rbind(IC,data.frame(AIC_B = AIC(fit_B),
                            AIC_C = AIC(fit_C),
                            AIC_C_B = AIC(fit_C_B),
                            BIC_B = BIC(fit_B),
                            BIC_C = BIC(fit_C),
                            BIC_C_B = BIC(fit_C_B),
                            iter = val))
  
  ### Collecting data from Covariate assisted regularized spectral clustering and corresponding alpha values
  X <- as.matrix(V(g.sim_C_B)$LOTR)
  Xdum <- as.matrix(data.frame(predict(dummyVars("~.", data = X ), newdata <- X)))
  origmemCov <- CovAssistedSpecClust(G = g.sim_C_B, Xdum, 3,
                                     Regularize =TRUE, alpha =NA,
                                     type ="non-assortative",
                                     retAlpha = TRUE)
  CovARI <- rbind(CovARI, c(mclust::adjustedRandIndex(V(g.sim_C_B)$Cluster, origmemCov[[1]]),
                            origmemCov[[2]]))
  
  for(pt in pctMsng){
    covMsng <- data.frame(name = c("LOTR") ,
                          p = pt,
                          cat = I(list(V(g.sim_C_B)$LOTR)))
    coms <- comdet(g.sim_C_B, V(g.sim_C_B)$LOTR)
    df <- IterativeImpute(g.sim_C_B, N, coms,val,k,covMsng)
    #g.imp[[i]] <- df[[2]]
    
    ## keeping track of community detection values for percentage missing
    # comDetVals <- rbind(comDetVals,
    #                     c(val,pt,round(mclust::adjustedRandIndex(V(g.imp[[i]])$Cluster,V(g.imp[[i]])$com),3), 
    #                       round(mclust::adjustedRandIndex(V(g.imp[[i]])$Cluster,coms$origmem),3)))
    # 
    op <- rbind(op,df[[1]])
    #i <- i+1
  }
  
  
  ## Use the probability of connection values defined earlier for only cluster for fitting the data. 
  g.sim_C <- simulate(g.sim ~ nodemix("Cluster", levels = TRUE, levels2 = TRUE),
                      nsim = 1,
                      coef = prob2logit(c(B)),
                      control=control.simulate( MCMC.burnin=10000, MCMC.interval=1000))
  
  
  ## Use the probability of connection values defined earlier for only category for fitting the data. 
  g.sim_B <- simulate(g.sim ~ nodemix(c("LOTR"), levels=TRUE, levels2 = TRUE),
                      nsim = 1,
                      coef = prob2logit(c(combi$prob)),
                      control=control.simulate( MCMC.burnin=10000, MCMC.interval=1000))
  
}

#summary(fit_C_B)
#summary(fit_C)
#summary(fit_B)
