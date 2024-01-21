library(igraph)
library(tidyverse)
library(mclust)
library(utils)
library("RColorBrewer")
library(ggplot2)
library(gg)
library(wesanderson)
library(gplots)
library(kernlab)
library(gridExtra)
library(mltools)
library("RColorBrewer")
source("Code/HelperFuncs.R")


## Lazega lawyer network data.

#Coworker network data
cow <- read.table("../../Research/ResearchCode/Data/lawyers/LazegaLawyers/ELwork.dat")

#Friendship network data
adv <- read.table("../../Research/ResearchCode/Data/lawyers/LazegaLawyers/ELadv.dat")

#Advice network data
frn <- read.table("../../Research/ResearchCode/Data/lawyers/LazegaLawyers/ELfriend.dat")

#Attributes data 

# The three networks refer to co-worker, friendship, and advice. The first 36 respondents are the partners in the firm. 
# The attribute variables in the file ELattr.dat are:
# seniority
# status (1=partner; 2=associate)
# gender (1=man; 2=woman)
# office (1=Boston; 2=Hartford; 3=Providence)
# years with the firm
# age
# practice (1=litigation; 2=corporate)
# law school (1: harvard, yale; 2: ucon; 3: other)

attri <- read.table("../../Research/ResearchCode/Data/lawyers/LazegaLawyers/ELattr.dat")
colnames(attri) <- c('name','status','gender', 'office', 'years','age','practice','school')
att <- attri
att$status <- as.factor(att$status)
att$gender <- as.factor(att$gender)
att$office <- as.factor(att$office)
att$practice <- as.factor(att$practice)
att$school <- as.factor(att$school)
att$years <- scale(att$years)
att$age <- scale(att$age)
## Tis is the only undirected network in this dataset
cowNW <- graph_from_adjacency_matrix(as.matrix(cow), mode = "undirected")
plot.igraph(cowNW, margin=0, frame=F, vertex.size = 3,vertex.frame = NA, 
            vertex.label=NA, edge.color="cornsilk4", 
            edge.arrow.size = 0.01)

advNW <- graph_from_adjacency_matrix(as.matrix(adv))
plot.igraph(advNW, margin=0, frame=F, vertex.size = 3,vertex.frame = NA, 
            vertex.label=NA, edge.color="cornsilk4",
            edge.arrow.size = 0.01)

frnNW <- graph_from_adjacency_matrix(as.matrix(frn))
plot.igraph(frnNW, margin=0, frame=F, vertex.size = 3,vertex.frame = NA, 
            vertex.label=NA, edge.color="cornsilk4",
            edge.arrow.size = 0.01)

vertex_attr(cowNW) <- att

ggplot(att) + geom_point(aes(x = name, y = age,color = as.factor(gender)))

## Strategies of generating missingness in data
## Missing Ties

## definding the percentage of missing ties

pcntMsgn <-c(50) #c(0, 2, 5, 10, 20, 30, 40, 50, 60, 70, 75, 80, 85, 90) 


## Functions created
### unit non-response : all outgoing ties for a node are missing

### item non-response : occasionally missing ties of nodes

## Lets start with choosing a particular variables and remove the edges:



#imputeContAtt(att$years)
#getwd()
## Running the imputation method N times each and check accuracy of simple imputations
N <- 100
origmem <- RegSpectralClust(cowNW, 3,regularize = TRUE)
#origmemadv <- RegSpectralClust(advNW, 3, regularize = TRUE)
#adjustedRandIndex(origmem, origmemadv)


# origmemCov <- CovAssistedSpecClust(cowNW, 
#                                    as.matrix(att),
#                                    3,1, Regularize = TRUE,
#                                    "assortative", 10000)

#adjustedRandIndex(origmem, origmemCov)

plot.igraph(cowNW, margin=0, frame=F, vertex.size = 3,vertex.frame = NA, 
            vertex.label=NA, edge.color="cornsilk4", vertex.color = origmem,
            edge.arrow.size = 0.01)

plot.igraph(cowNW, margin=0, frame=F, vertex.size = 3,vertex.frame = NA, 
            vertex.label=NA, edge.color="cornsilk4", vertex.color = origmemadv,
            edge.arrow.size = 0.01)

plot.igraph(advNW, margin=0, frame=F, vertex.size = 3,vertex.frame = NA, 
            vertex.label=NA, edge.color="cornsilk4", vertex.color = origmemadv,
            edge.arrow.size = 0.01)

ARIcomparison <- as.data.frame(pcntMsgn)
ARIcomparison$rmRandEdges <- 0
ARIcomparison$rmsts1 <- 0 ##status 1
ARIcomparison$rmsts2 <- 0
ARIcomparison$rmgndr1 <- 0
ARIcomparison$rmgndr2 <- 0
ARIcomparison$rmofc1 <-0
ARIcomparison$rmofc2 <-0
ARIcomparison$rmofc3 <-0 ## office 3 has very few people
ARIcomparison$rmyrLT5 <- 0
ARIcomparison$rmyrGT5 <- 0
ARIcomparison$rmageLT40 <- 0
ARIcomparison$rmageGT40 <- 0
ARIcomparison$rmprc1 <- 0
ARIcomparison$rmprc2 <- 0
ARIcomparison$rmsch1 <- 0 
ARIcomparison$rmsch2 <- 0
ARIcomparison$rmsch3 <- 0

## Adding columns to save ARI after imputing attributes. 
## These values were removed at random.
# 
# ARIcomparison$rmAttSts <- 0
# ARIcomparison$rmAttGndr <- 0
# ARIcomparison$rmAttOfc <- 0
# ARIcomparison$rmAttYr <- 0
# ARIcomparison$rmAttAge <- 0
# ARIcomparison$rmAttPrc <- 0
# ARIcomparison$rmAttSch <- 0

## Adding columns to save ARI after imputing attributes. 
## These values were removed at random for SPECIFIC attributes.


j <- 0
G <- cowNW
for(percent in pcntMsgn){
  j <- j +1
  ariVal <- 0
  print(percent)
  #pb <- txtProgressBar(min = 0, max = N, style = 3)
  for(i in 1:N){
    #Randomly remove edges and run the above for 100 runs. compare with Spectral clustering and find ARI
    op <- removeRandEdges(G, percent)
    specmem <- RegSpectralClust(op$graph, 3, regularize = TRUE)
    v <- adjustedRandIndex(origmem,specmem)
    #print(v)
    ARIcomparison$rmRandEdges[j] <-  ARIcomparison$rmRandEdges[j] + v/N
    #print(1)
    op <- removeCatedges(G,"status", `==`,percent, "1" )
    specmem <- RegSpectralClust(op$graph, 3, regularize = TRUE)
    ARIcomparison$rmsts1[j] <-  ARIcomparison$rmsts1[j] + adjustedRandIndex(origmem, specmem)/N ##status 1

    op <- removeCatedges(G,"status", `==`,percent, "2" )
    specmem <- RegSpectralClust(op$graph, 3, regularize = TRUE)
    ARIcomparison$rmsts2[j] <- ARIcomparison$rmsts2[j]+ adjustedRandIndex(origmem, specmem)/N

    op <- removeCatedges(G,"gender", `==`,percent, "1" )
    specmem <- RegSpectralClust(op$graph, 3, regularize = TRUE)
    ARIcomparison$rmgndr1[j] <- ARIcomparison$rmgndr1[j]+ adjustedRandIndex(origmem, specmem)/N

    op <- removeCatedges(G,"gender", `==`,percent, "2" )
    specmem <- RegSpectralClust(op$graph, 3, regularize = TRUE)
    ARIcomparison$rmgndr2[j] <- ARIcomparison$rmgndr2[j]+ adjustedRandIndex(origmem, specmem)/N

    op <- removeCatedges(G,"office", `==`,percent, "1" )
    specmem <- RegSpectralClust(op$graph, 3, regularize = TRUE)
    ARIcomparison$rmofc1[j] <- ARIcomparison$rmofc1[j]+ adjustedRandIndex(origmem, specmem)/N

    op <- removeCatedges(G,"office", `==`,percent, "2" )
    specmem <- RegSpectralClust(op$graph, 3, regularize = TRUE)
    ARIcomparison$rmofc2[j] <- ARIcomparison$rmofc2[j]+ adjustedRandIndex(origmem, specmem)/N

    # op <- removeCatedges(G,"office", `==`,percent, "3" )
    # specmem <- RegSpectralClust(op$graph, 3, regularize = TRUE)
    # ARIcomparison$rmofc3[j] <-ARIcomparison$rmofc3[j]+ adjustedRandIndex(origmem, specmem)/N ## office 3 has very few people
    #print(2)
    op <- removeCatedges(G,"years", `<`,percent, median(att$years) )
    specmem <- RegSpectralClust(op$graph, 3, regularize = TRUE)
    ARIcomparison$rmyrLT5[j] <- ARIcomparison$rmyrLT5[j]+ adjustedRandIndex(origmem, specmem)/N

    op <- removeCatedges(G,"years", `>`,percent, median(att$years) )
    specmem <- RegSpectralClust(op$graph, 3, regularize = TRUE)
    ARIcomparison$rmyrGT5[j] <- ARIcomparison$rmyrGT5[j]+ adjustedRandIndex(origmem, specmem)/N
    #print(3)
    op <- removeCatedges(G,"age", `<`,percent, median(att$age) )
    specmem <- RegSpectralClust(op$graph, 3, regularize = TRUE)
    ARIcomparison$rmageLT40[j] <- ARIcomparison$rmageLT40[j]+ adjustedRandIndex(origmem, specmem)/N

    op <- removeCatedges(G,"age", `>`,percent, median(att$age) )
    specmem <- RegSpectralClust(op$graph, 3, regularize = TRUE)
    ARIcomparison$rmageGT40[j] <- ARIcomparison$rmageGT40[j]+ adjustedRandIndex(origmem, specmem)/N

    op <- removeCatedges(G,"practice", `==`,percent, "1" )
    specmem <- RegSpectralClust(op$graph, 3, regularize = TRUE)
    ARIcomparison$rmprc1[j] <- ARIcomparison$rmprc1[j]+ adjustedRandIndex(origmem, specmem)/N

    op <- removeCatedges(G,"practice", `==`,percent, "2" )
    specmem <- RegSpectralClust(op$graph, 3, regularize = TRUE)
    ARIcomparison$rmprc2[j] <- ARIcomparison$rmprc2[j]+ adjustedRandIndex(origmem, specmem)/N

    op <- removeCatedges(G,"school", `==`,percent, "1" )
    specmem <- RegSpectralClust(op$graph, 3, regularize = TRUE)
    ARIcomparison$rmsch1[j] <- ARIcomparison$rmsch1[j]+ adjustedRandIndex(origmem, specmem)/N

    op <- removeCatedges(G,"school", `==`,percent, "2" )
    specmem <- RegSpectralClust(op$graph, 3, regularize = TRUE)
    ARIcomparison$rmsch2[j] <- ARIcomparison$rmsch2[j]+ adjustedRandIndex(origmem, specmem)/N
    #print(4)
    op <- removeCatedges(G,"school", `==`,percent, "3" )
    specmem <- RegSpectralClust(op$graph, 3, regularize = TRUE)
    ARIcomparison$rmsch3[j] <- ARIcomparison$rmsch3[j]+ adjustedRandIndex(origmem, specmem)/N

    #setTxtProgressBar(pb, i)
  }
}

## cluster detection using Kmeans directly on the covariates.
kmem<- matrix(nrow = N, ncol = 0)
#for(i in 1:N){
kmem <- kmeans(att[,-c(1)], 3, nstart = 1000)
#}
adjustedRandIndex(kmem$cluster , origmem)
plot.igraph(cowNW, margin=0, frame=F, vertex.size = 3,vertex.frame = NA, 
            vertex.label=NA, edge.color="cornsilk4", vertex.color = kmem$cluster,
            edge.arrow.size = 0.01)

## performing spectral clustering on just the covariate data.
Dat <- one_hot(as.data.table(att[,-c(1)]))
Dat <- as.data.frame(Dat)
spec <- specc(as.matrix(Dat),centers =  3)
adjustedRandIndex(spec@.Data , origmem)
plot.igraph(cowNW, margin=0, frame=F, vertex.size = 3,vertex.frame = NA, 
            vertex.label=NA, edge.color="cornsilk4", vertex.color = spec@.Data,
            edge.arrow.size = 0.01)

adjustedRandIndex(kmem$cluster, spec@.Data)


### Creating a heatmap of the ari values for the ARI comparison dataframe.
rownames(ARIcomparison) <- ARIcomparison$pcntMsgn

col <- colorRampPalette(brewer.pal(10, "RdYlBu"))(256)
heatmap.2(as.matrix(ARIcomparison[, -c(1,2)]), scale = "none", Rowv = NA, 
          Colv = NA, trace ="none" ,col = bluered(300), density.info = "none")

heatmap(as.matrix(ARIcomparison[, -c(1,2)]), Rowv = NA, 
          Colv = NA,scale = "none", col = bluered(200))
## What if we consider the cowNW to be ground truth and 

### USe a combinations of all the different covariates and compare with the spectral without covariates
att <- att[,-1]
names <- colnames(att)
combi <- do.call("c",lapply(seq_along(names),function(i) utils::combn(names,i,FUN=base::list)))

categoricalCols <- 'status|gender|office|practice|school'

continousCols <- 'years|age'

combiOPCov <- list()
combiOPNoCov <- list()

for (com in combi[1:7]) {
  print(com)
  subDat <- subset(att, select = com)
  #subDat <- one_hot(as.data.table(subDat))
  #subDat <- as.data.frame(subDat)
  
  subDatopNoCov <- as.data.frame(matrix(0,ncol = ncol(subDat), nrow = length(pcntMsgn)))
  colnames(subDatopNoCov) <- colnames(subDat)
  
  subDatopCov <- as.data.frame(matrix(0,ncol = ncol(subDat), nrow = length(pcntMsgn)))
  colnames(subDatopCov) <- colnames(subDat)
  
  j <- 0
  
  catCols <- which(str_detect(colnames(subDat) , 
                              categoricalCols) == TRUE)
  contCols <- which(str_detect(colnames(subDat) , 
                               continousCols) == TRUE)
  
  op <- cowNW
  ## find the communities for original covariates assisted graphs
  CovMem <- CovAssistedSpecClust(op, as.matrix(as.data.frame(one_hot(as.data.table(subDat)))), 3,1)
  
  for(percent in pcntMsgn){
    j <- j + 1
    #ariVal <- 0
    print(percent)
    
    
    op <- cowNW
    ## assigning the vertex attributes to graph
    #vertex_attr(op) <- subDat
    
    
    for(i in 1:N){
      if(!is_empty(catCols)){
        for (col in catCols) {
          covDat <- subDat
          covDat[,col] <-  removeRandAtt(covDat[,col], percent)
          a <- att
          a[, colnames(covDat)[col]] <- covDat[,col]
          CatRegImp(a, names[col])
          #if(is.null(covDat[,col])){print(covDat)}
          covDat[,col] <-  imputeCatAtt(covDat[,col])
          mem <- CovAssistedSpecClust(op, as.matrix(as.data.frame(one_hot(as.data.table(covDat)))), 3,1)
          subDatopNoCov[j,col] <- subDatopNoCov[j,col] + adjustedRandIndex(origmem, mem)/N
          subDatopCov[j,col] <- subDatopCov[j,col] + adjustedRandIndex(CovMem, mem)/N
        }
      }
      #print("done with categorical cols") # nolint
      if(!is_empty(contCols)){
        for (col in contCols) {
          covDat <- subDat
          covDat[,col] <-  removeRandAtt(covDat[,col], percent)
          covDat[,col] <-  as.integer(imputeContAtt(covDat[,col]))
          mem <- CovAssistedSpecClust(op, as.matrix(as.data.frame(one_hot(as.data.table(covDat)))), 3,1)
          subDatopNoCov[j,col] <- subDatopNoCov[j,col] + adjustedRandIndex(origmem, mem)/N
          subDatopCov[j,col] <- subDatopCov[j,col] + adjustedRandIndex(CovMem, mem)/N
        }
      }
    }
  }
  combiOPCov <- append(combiOPCov, cbind(pcntMsgn, subDatopCov))
  combiOPNoCov <- append(combiOPNoCov, cbind(pcntMsgn, subDatopNoCov))
}



## saving the combination in rds file

s <- cbind(seq(from = 1, to = 127, by = 3), seq(from = 3, to = 127, by = 3))

for(i in 1:dim(s)[1]){
  
  saveRDS(combi[s[i,1]:s[i,2]], paste0("Input/Number",s[i,1], ".rds"))
}


saveRDS(combi[1:2],"Number1_2.rds")
saveRDS(combi[1:28],"Number1_28.rds")
saveRDS(combi[29:46],"Number29_46.rds")
saveRDS(combi[47:63],"Number47_63.rds")
saveRDS(combi[64:74],"Number64_74.rds")
saveRDS(combi[75:84],"Number75_84.rds")
saveRDS(combi[85:98],"Number85_98.rds")
saveRDS(combi[99:110],"Number99_110.rds")
saveRDS(combi[111:120],"Number111_120.rds")
saveRDS(combi[121:127],"Number121_127.rds")

########## Creating plots using the output from Cluster ###########


ARIcomp <- readRDS("output/ARIcomp2.rds")
ARIcomp$rmRandEdges <- ARIcomp$rmRandEdges/100
### Creating a heatmap of the ari values for the ARI comparison dataframe.
rownames(ARIcomp) <- ARIcomp$pcntMsgn

col <- wes_palette("Zissou1", 500, type = "continuous")#colorRampPalette(brewer.pal(10, "RdYlBu"))(256)
heatmap.2(as.matrix(subset(ARIcomp, select = -c(pcntMsgn, rmofc3))), scale = "none", 
          Rowv = NA, Colv = NA, trace ="none" ,col = col, density.info = "none")

Nedges <- readRDS("output/Nedges2.rds")

col <- wes_palette("Zissou1", 500, type = "continuous")#colorRampPalette(brewer.pal(10, "RdYlBu"))(256)
rownames(Nedges) <- Nedges[,1]
heatmap.2(as.matrix(subset(Nedges, select = -c(pcntMsgn, rmofc3))), scale = "none", 
          Rowv = NA, Colv = NA, trace ="none" ,col = col, density.info = "none")


#### Plotting the remove/impute category data.


files <- list.files(path = "output/outputN10", pattern = "^Number")
files <- as.data.frame(cbind(files, str_extract(files, "[0-9]+")))
files$V2 <- as.integer(files$V2)
files <- arrange(files, V2)
covcomp <- list()
Noncovcomp <- list()
for(file in files$files){
  print(file)
  op <- readRDS(paste0("output/outputN10/",file))
  covcomp <- append(covcomp,op[[1]])
  Noncovcomp <- append(Noncovcomp, op[[2]])
}

## Reading the singles
createHeatmap <- function(data, palette = "Royal1", title){
  singles <- as.data.frame(do.call(cbind, data)) %>%
    dplyr::select(-starts_with("pcnt"))
  rownames(singles) <- c(2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90) 
  col <- wes_palette(palette, 500, type = "continuous")#colorRampPalette(brewer.pal(10, "RdYlBu"))(256)
  hm <-heatmap.2(as.matrix(singles), scale = "none", 
            Rowv = NA, Colv = NA, trace ="none" ,col = col, density.info = "none",main = title)  
  return(hm)
  }


hmcov <- createHeatmap(covcomp[1:14], "Darjeeling1", "Comparing to CD with covariates")

hmnocov <- createHeatmap(Noncovcomp[1:14], "FantasticFox1", "Comparing to Cd without covariates")

singles <- as.data.frame(do.call(cbind, covcomp[1:14])) %>%
  dplyr::select(-starts_with("pcnt"))
singles$pcnt <- c(2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90) 

singles <- singles %>% pivot_longer(cols = !pcnt)

ggplot(singles)+ geom_point(aes(y = value, x = pcnt, color =name))+
  geom_line(aes(y = value, x = pcnt, color =name))

## Doubles
doubles <- as.data.frame(do.call(cbind,covcomp[15:71]))
createHeatmap <- function(data, palette, title, by_val){
  gbg <- as.data.frame(matrix(ncol = by_val+1, nrow = 0))
  for(i in seq(1,length(colnames(data)), by=by_val)){
    gbg <- rbind(gbg, cbind(pivot_longer(data[,i:(i+by_val-1)], !pcntMsgn), i))
    
  }
  col <- wes_palette(palette, 500, type = "continuous")#colorRampPalette(brewer.pal(10, "RdYlBu"))(256)
  p <- gbg %>%  ggplot() +
    geom_tile(aes(x = name, y = factor(pcntMsgn), fill = value)) +
    scale_fill_gradientn(colours = col) +
    facet_wrap(~ i, scales = "free") +
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank()
    )
  return(p)
}

createHeatmap(doubles, "FantasticFox1", "W/ covariates comparison", 3)

##Triples
triples <- as.data.frame(do.call(cbind,covcomp[72:217]))
createHeatmap(triples, "FantasticFox1", "W/ covariates comparison", 4)

## Create a map for each variable in the data## example. Compare status in presence of all the other variables.##


coldoubles <- as.data.frame(matrix(colnames(doubles)[-grep("pcnt", colnames(doubles))], ncol = 2, byrow = TRUE))
colnames(doubles) <- make.unique(colnames(doubles))
colnames(coldoubles) <- c("a","b")
coldoubles$pair <- 1:dim(coldoubles)[1]
coldoubles <- pivot_longer(coldoubles, cols = c("a","b") )


doubleHeatmaps <- function(var){
  Dat <- doubles %>% dplyr::select(starts_with(var))
  colnames(Dat) <- coldoubles$value %>% 
    subset(coldoubles$pair %in% coldoubles$pair[grep(var,coldoubles$value)] & !grepl(var, coldoubles$value), ) 
 return(heatmap.2(as.matrix(Dat), scale = "none",  Rowv = NA, Colv = NA, trace ="none" , breaks = seq(0,1, length = 500), key=FALSE,
            col =  wes_palette("Zissou1", 499, type = "continuous"), density.info = "none", main = var))
}

doubleHeatmaps("status")
doubleHeatmaps("gender")
doubleHeatmaps("office")
doubleHeatmaps("practice")
doubleHeatmaps("school")
doubleHeatmaps("age")
doubleHeatmaps("years")



misng <- as.data.frame(1- Nedges[,-c(1)]/gsize(cowNW)) %>% 
  pivot_longer(everything()) %>% subset(name != "rmofc3")

ariVal <- ARIcomp[, -c(1)] %>% pivot_longer(everything()) %>% subset(name != "rmofc3")

dat <- cbind(misng, ariVal$value)

colnames(dat) <-c("Attribute", "pcntMissing","ARI")

ggplot(dat) +
  geom_point(aes(x = pcntMissing, y = ARI, color = Attribute))+
  geom_line(aes(x = pcntMissing, y = ARI, color = Attribute))

###### Node importance measures for original network

closeness.cent <- closeness(cowNW, mode = "all")

degree.cent<- degree(cowNW, mode = "all")

bet <- betweenness(cowNW)

neighbourhood <- ego_size(cowNW, order = 2, mode = "all")

ARILOO <- list()
### Remove one node at atime and perform community detection. See if there is a change.

for(v in V(cowNW)){
  nw <- delete_vertices(cowNW, v)
  ## spectral clustering for the removed graph
  mem <- RegSpectralClust(nw, 3)
  aricomp <- adjustedRandIndex(mem, origmem[-which(V(cowNW) == v)])
  ARILOO <- append(ARILOO, aricomp)
  
}

nodecomp <- cbind(ARI = unlist(ARILOO), closeness.cent, degree.cent, bet, neighbourhood)

## creating a correlation matrix to compare node importance measures and the ARIs

cor(nodecomp, use ="pairwise.complete.obs")

## Remove based on node importance measures. Multiple at a time. and check if there is any difference. 

####################### data imputation using regression

files <- list.files(path = "output/RegressionImputationAlpha5", pattern = "^Number")
files <- as.data.frame(cbind(files, str_extract(files, "[0-9]+")))
files$V2 <- as.integer(files$V2)
files <- arrange(files, V2)
covcomp <- list()
Noncovcomp <- list()
RegComp <- list()
diff <- list()
for(file in files$files){
  print(file)
  op <- readRDS(paste0("output/RegressionImputationAlpha5/",file))
  covcomp <- append(covcomp,op[[1]])
  Noncovcomp <- append(Noncovcomp, op[[2]])
  RegComp <-  append(RegComp, op[[3]])
  diff <- append(diff, op[[4]])
}

as.data.frame(do.call(cbind, RegComp[1:14])) %>%
  dplyr::select(-starts_with("pcnt")) %>%
  mutate(pcnt=c(2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90)) %>%
  pivot_longer(cols = !pcnt)%>%
  ggplot() + 
  geom_point(aes(y = value, x = pcnt, color = name)) +
  geom_line(aes(y = value, x = pcnt, color = name)) +ggtitle("Using regression imputation")

as.data.frame(do.call(cbind, covcomp[1:14])) %>%
  dplyr::select(-starts_with("pcnt")) %>%
  mutate(pcnt=c(2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90)) %>%
  pivot_longer(cols = !pcnt)%>%
  ggplot() + 
  geom_point(aes(y = value, x = pcnt, color = name)) +
  geom_line(aes(y = value, x = pcnt, color = name))+ggtitle("Using simple imputation")

# as.data.frame(do.call(cbind, Noncovcomp[1:14])) %>%
#   dplyr::select(-starts_with("pcnt")) %>%
#   mutate(pcnt=c(2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90)) %>%
#   pivot_longer(cols = !pcnt)%>%
#   ggplot() + 
#   geom_point(aes(y = value, x = pcnt, color = name)) +
#   geom_line(aes(y = value, x = pcnt, color = name))+ggtitle("Using simple imputation Comparing with non covariance CD")
# 

regrI <- as.data.frame(do.call(cbind, RegComp[1:14])) %>%
  dplyr::select(-starts_with("pcnt")) %>%
  mutate(pcnt=c(2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90)) %>%
  pivot_longer(cols = !pcnt) %>%
  mutate(type = "Regression")

sImp <- as.data.frame(do.call(cbind, covcomp[1:14])) %>%
  dplyr::select(-starts_with("pcnt")) %>%
  mutate(pcnt=c(2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90)) %>%
  pivot_longer(cols = !pcnt) %>%
  mutate(type = "Simple Imputation")

df <- rbind(regrI, sImp)

df %>%
  ggplot() + 
  geom_point(aes(y = value, x = as.factor(pcnt), color = type, group = value)) +
  geom_line(aes(y = value, x = as.factor(pcnt), color = type, group = type)) +
  facet_wrap(vars(name), ncol = 1)+theme_minimal() + ggtitle("compare regression and simple imputation")

Sdiff <- as.data.frame(do.call(cbind, diff[1:14])) %>%
  dplyr::select(-starts_with("pcnt")) %>%
  mutate(pcnt=c(2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90)) %>%
  pivot_longer(cols = !pcnt)

colnames(regrI) <- c("pcnt","name","Regression","type")
compare <- regrI
compare$SimpleImpute <- sImp$value
compare$Valdiff <- compare$Regression - compare$SimpleImpute
compare$numdiff <- Sdiff$value

compare %>% 
  ggplot() + 
  geom_point(aes(y = Valdiff, x = numdiff, color = name)) +
  geom_line(aes(y = Valdiff, x = numdiff, color = name))

compare %>% 
  ggplot() + 
  geom_point(aes(y = Valdiff, x = log(numdiff), color = name)) +
  geom_line(aes(y = Valdiff, x = log(numdiff), color = name))

compare %>% 
  ggplot() + 
  geom_point(aes(y = Valdiff, x = pcnt, color = name)) +
  geom_line(aes(y = Valdiff, x = pcnt, color = name))


compare %>% 
  ggplot() + 
  geom_point(aes(y = numdiff, x = pcnt, color = name)) +
  geom_line(aes(y = numdiff, x = pcnt, color = name))





