## Iterative sampling Base code Script

source("HelperFuncs.R")

### Categorical variable correlated to the assigned clusters and different number of levels.(binary variable where one is correlated to two clusters and 2nd is correlated to the last cluster. Or more than one level correlated with cluster)
### Categorical variable correlated to the assigned clusters and same number of levels.
## Using LOTR Universe books.
niter= 5
superiter <- 1000
k <- 3
N <- 100
## Almost equally defined clusters
pC = c(0.3,0.4,0.3)
## Highly defined clusters
B = c(0.3,0.1,0.3,0.05,0.1,0.2)
## setting color codes
catcol  <- c("#B22222","#FF7F50","#FF8C00",
             "#6B8E23","#7FFF00","#32CD32",
             "#6495ED","#00008B","#87CEEB")
cat <- c("FOTR","TT","ROTK","Hob","Sil", "UT","CoH","BL","FoG")
col <- data.frame(as.matrix(cbind(catcol, cat)))
pCat <- list(c(0.3,0.4,0.4,0.02,0.01,0.02,0.02,0.01,0.01), # Cluster 1
             c(0.01,0.02,0.01,0.2,0.3,0.4,0.02,0.01,0.01), # Cluster 2
             c(0.01,0.01,0.02,0.01,0.2,0.02,0.3,0.4,0.2)) # Cluster 3

c1 <- c("FOTR|TT|ROTK")
c2 <- c("Hob|Sil|UT")
c3 <- c("CoH|BL|FoG")

opLD <- data.frame(matrix(nrow =0 , ncol = 14))

for(val in 1:superiter){
  ## Simulate network
  g.sim <- SimNW(pC,k,N,B)
  
  
  g.sim <- AssgnCat(g.sim, "LOTR",pCat, cat)
  
  g <- data.frame(as.matrix(V(g.sim)$LOTR))
  colnames(g) <- "cat"
  vtx <- left_join(g,col, by = join_by(cat== cat))
  #g.sim <- asIgraph(g.sim)
  V(g.sim)$col <- vtx$catcol
  
  ## Regularized and cova assisted spectral clustering
  coms <- comdet(g.sim, V(g.sim)$LOTR)
  
  pctMsng <- c(5,7,10,20,30,40,50)
  op <- data.frame(matrix(nrow=0,ncol = 9))
  g.imp <- list()
  i <-  1
  for(p in pctMsng){
    covMsng <- data.frame(name = c("LOTR") ,
                          p = p,
                          cat = I(list(cat)))
    df <- IterativeImpute(g.sim, N, coms,niter,k,covMsng)
    g.imp[[i]] <- df[[2]]
    op <- rbind(op,df[[1]])
    i <- i+1
  }
  
  op$pcnt <- rep(pctMsng, times = niter*pctMsng)
  op$sameClust <- NA
  for(i in 1:dim(op)[1]){
    str <- c(op[i,4], op[i,6])
    op$sameClust[i] <- all(sapply(str, grepl, c1)) | all(sapply(str, grepl, c2)) | all(sapply(str, grepl, c3))
    
  }
  
  op$superiter <- val
  opLD <- rbind(opLD,op)
  print(val)
  
}
propCorClust <- opLD %>%
  group_by(iteration, pcnt,superiter) %>%
  dplyr::summarise(prop = sum(sameClust)/n())
# for(i in 1:length(pctMsng)){
#   plot_combinations(g.imp[[i]], coms, pctMsng[i])
# }

#op <- op[with(op,order(op[,2],op[,1],op[,11])),]
