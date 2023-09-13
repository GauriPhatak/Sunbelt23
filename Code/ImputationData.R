library(igraph)
library(tidyverse)


### Reading and cleaning data 
## Reading EU core email data. But this dataset does not contain any node attributes
#e <- as.matrix(read.table(gzfile("../../Research/ResearchCode/Data/EmailEUCore/.txt.gz")))
# e[,1] <- as.numeric(e[,1])+1
# e[,2] <- as.numeric(e[,2])+1
# c <- read.table(gzfile("../../Research/ResearchCode/Data/EmailEUCore/email-Eu-core-department-labels.txt.gz")) 
# 
# ## Twitter Dataset obtained from : https://github.com/he-tiantian/Attributed-Graph-Data/tree/master
# ## Reading edge list
# e <- as.matrix(read.table("../../Research/ResearchCode/Data/TwitterWAtt/twitter/edgelist"))
# e[,1] <- as.numeric(e[,1])+1
# e[,2] <- as.numeric(e[,2])+1
# 
# ## Vertex class info
# #c <- read.table("../../Research/ResearchCode/Data/TwitterWAtt/twitter/Class_info",
# #               quote = "\"'",)
# c <- matrix(readLines(con <- file("../../Research/ResearchCode/Data/TwitterWAtt/twitter/Class_info"))[-1], 
#             ncol = 2, byrow = TRUE)
# close(con)
# 
# ## Attributes
# a <- read.table("../../Research/ResearchCode/Data/TwitterWAtt/twitter/Attribute",
#                 sep = "\t", comment.char= "", quote = "")
# 
# at <- read.table("../../Research/ResearchCode/Data/TwitterWAtt/twitter/vertex2aid")
# 
# ## Reading Flickr data
# ##In this network, the ground-truth communities are defined as user-created interest-based groups that have more than five members.
# #e <- as.matrix(read.table(gzfile("../../Research/ResearchCode/Data/Flickr/flickrEdges.txt.gz")))
# 
# 
# ## Reading facebook dataset
# ## reading the combined data file. 
# e <- as.matrix(read.table(gzfile("../../Research/ResearchCode/Data/FaceboookCirc/facebook_combined.txt.gz")))
# e[,1] <- as.numeric(e[,1])+1
# e[,2] <- as.numeric(e[,2])+1
# 
# # Reading the 
# 
# ## Reading the ground truth community of the ego networks
# # Create list of text files
# txt_files_ls = list.files(path="../../Research/ResearchCode/Data/FaceboookCirc/facebook/", pattern="*.circles") 
# # Read the files in, assuming comma separator
# txt_files_df <- lapply(txt_files_ls, function(x) {read.table(file = paste0("../../Research/ResearchCode/Data/FaceboookCirc/facebook/", x), sep =" ")})
# # Combine them
# combined_df <- do.call("rbind", lapply(txt_files_df, as.data.frame)) 
# combined_df <- separate(combined_df, V1, sep = "\\t", into = c("circle","nodes"),extra = "merge")
# combined_df$com <- 1:length(combined_df$circle)
# combined_df <- combined_df %>% separate_longer_delim(nodes, delim = "\t")
# 
# ## Reading the edges form the egonetwork data
# 
# ## Reading the ground truth community of the ego networks
# # Create list of text files
# txt_files_ls = list.files(path="../../Research/ResearchCode/Data/FaceboookCirc/facebook/", pattern="*.edges")
# # Read the files in, assuming comma separator
# txt_files_df <- lapply(txt_files_ls, 
#                        function(x) { m <- read.table(file = paste0("../../Research/ResearchCode/Data/FaceboookCirc/facebook/", x), sep =" ")
#                        rbind(m, cbind(as.numeric(gsub("\\D","",x)), unique(unlist(m))))})
# # Combine them
# combined_e <- do.call("rbind", lapply(txt_files_df, as.data.frame)) 
# combined_e[,1] <- combined_e[,1]+1
# combined_e[,2] <- combined_e[,2]+1
# #combined_df <- separate(combined_df, V1, sep = "\\s", into = c("circle","nodes"),extra = "merge")
# #combined_df$com <- 1:length(combined_df$circle)
# #combined_df <- combined_df %>% separate_longer_delim(nodes, delim = "\t")
# 
# 
# #a <- read.table("../../Research/ResearchCode/Data/TwitterWAtt/twitter/Attribute",
# #                sep = "\t", comment.char= "", quote = "")
# 

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

att <- read.table("../../Research/ResearchCode/Data/lawyers/LazegaLawyers/ELattr.dat")

colnames(att) <- c('name','status','gender', 'office', 'years','age','practice','school')

## Tis is the only undirected network in this dataset
cowNW <- graph_from_adjacency_matrix(as.matrix(cow), mode = "undirected")
plot.igraph(cowNW, margin=0, frame=F, vertex.size = 3,vertex.frame = NA, 
            vertex.label=NA, edge.color="cornsilk4", vertex.color =  V(cowNW)$status,
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

ggplot(att) + geom_point(aes(x = seniority, y = age,color = as.factor(gender)))

## Strategies of generating missingness in data
## Missing Ties

## definding the percentage of missing ties

pcntMsgn <- c(0, 2, 5, 10, 20, 30, 40) 

### unit non-response : all outgoing ties for a node are missing

### item non-response : occasionally missing ties of nodes

## Lets start with choosing a particular variables and remove the edges:

## for all categories at random.MAR
## Ties randomly removed
removeRandEdges <- function(graph, percent){
  
  N <- gsize(graph)
  numE <- floor((N * percent)/100)
  to_be_removed <- sample(1:N, numE, replace = FALSE)
  return(delete_edges(graph, to_be_removed))
}

cowNV10 <- removeRandEdges(cowNW, 10)


## For all categories based on a percentage for each categories in it. MNAR

removeCatedges <- function(graph, attr, category, fxn, percent){
  
  N <- gsize(graph)
  numE <- floor((N*percent)/100)
  to_be_removed <-  sample(unlist(incident_edges(graph, 
                                                 which(fxn(vertex_attr(graph, attr), category))), 
                                  recursive = FALSE), 
                           numE)
  
  return(delete_edges(graph, to_be_removed))
}

gbg <- removeCatedges(cowNW, "gender", "2", `==`,10)

## remove values in the attributes at random

removeRandAtt <- function(){
  
}

removeCatAtt <- function(){
  
}

## Strategies for imputing the data

## Unconditional means. This includes just simple methods of imputations.
## for categorical variables impute based on the most frequent variable
## Most common value imputation

imputeCatAtt <- function(attr){
  return(as.matrix(as.data.frame(sort(table(attr), decreasing = TRUE ))$attr[1]))
}

imputeCatAtt(att$status)

imputeContAtt <- function(attr){
  return(floor(mean(attr,na.rm = TRUE)))
}

imputeContAtt(att$years)

## Unconditional distributions

## Conditional means

## Conditional distributions


##  1st test
##    1. Impute
##    2. Community Detection
