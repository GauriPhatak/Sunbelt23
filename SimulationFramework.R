library(tidyverse)
library(igraph)
library(ergm)
library(fs)
library(network)
library(intergraph)
library(tidygraph)
library(ggraph)
library(statnet)

## Terms under consideration nodecov, star, triangle, transtive, dyad

###################################################3
##Synthetic population parameters
## Going to use only undirected network
## Setting the required coefficient values are the required argument values
## Arguments
k <- 3
N <- 100
pC = c(0.3,0.4,0.3)
B = c(0.7,0.02,0.7,0.03,0.02,0.7)
C = sample(1:k,N,replace = TRUE, prob = pC)
coefs = c(prob2logit(B))


# params <- list(
#   ## number of nodes
#   N =  100,
#   ## number of simulated networks
#   nsim = 1,
#   ## Number of communities 
#   k = 3,
#   ## proportion of nodes in each community
#   pC = c(0.3,0.4,0.3),
#   ## Assigning cluster to each node in the network
#   C = sample(1:k,N,replace = TRUE, prob = pC),
#   ## Mixture probability matrix
#   B = c(0.7,0.02,0.7,0.03,0.02,0.7),
#   ## coefficient values for the desired networks models
#   coefs = c(prob2logit(B), -1),
#   seed = 42,
#   edgeC = 0,
#   mutualC = 0,
#   TrianglesC = 0,
# )

## Need to create a nodematch coefficent vector for differential homophily. 
## Check nodematch function in ergm terms

## Required functions
logit2prob <- function(coef) {
  odds <- exp(coef)
  prob <- odds / (1 + odds)
  return(prob)
}

prob2logit <- function(x){
  return(log(x / (1 - x)))
}

## Create a basis network

NetworkSim <- function(N, dir=FALSE, B, C,formula, coefs){
  
  seed <- 42
  net <- network(N, directed = dir, density= 0 )
  net %v% 'Cluster' <- C
 
  #form <- as.formula(paste0("net~"))
  if(!is.null(formula)){
    g.sim <- simulate(as.formula(paste0("net~nodemix('Cluster',levels = TRUE,levels2 = TRUE)+", formula)), 
                      coef = coefs,
                      seed = seed)  
  }
  else{
    g.sim <- simulate(as.formula(paste0("net~nodemix('Cluster',levels = TRUE,levels2 = TRUE)")), 
                      coef = coefs,
                      seed = seed)
  }
  
  # g.sim <- simulate(net~nodemix('Cluster',
  #                               levels = TRUE,
  #                               levels2 = TRUE)+kstar(3), 
  #                   coef = c(prob2logit(B), -1),
  #                   seed = seed)
  
  return(g.sim)
  
}


g.sim <- NetworkSim(N,FALSE,B,C, NULL, coefs)

#Simulate a network with only 
fit <- ergm(g.sim~nodemix('Cluster', levels = TRUE, levels2 = TRUE))
summary(fit)

plot( g.sim,
      edge.col = adjustcolor('black',alpha.f = 0.5),
      #main=paste('homophily =', round(logit2prob(Hseq[i]), digits = 2)),
      vertex.col = g.sim %v% 'Cluster')

## 1.2 2.2 1.3 2.3 3.3
logit2prob(coef(fit))

Hseq <- seq(-4,4, length.out = 9)

par(mfrow=c(3,3), mar=c(2,2,2,2))

for(i in seq_along(Hseq)){
  g.sim <- simulate(~edges+nodematch('Cluster', diff=TRUE), 
                    coef = c(-4,Hseq[i]), 
                    seed = seed,
                    basis= net)
  coef(g.sim)
  #summary(g.sim)
  plot( g.sim,
        edge.col = adjustcolor('black',alpha.f = 0.5),
        #main=paste('homophily =', round(logit2prob(Hseq[i]), digits = 2)),
        vertex.col = g.sim %v% 'Cluster')
}

### Use degree term? Probably not

############# EXPLANATION OF UNIFORM AND DIFFERENTIAL HOMOPHILY
#There are two broad flavors of homophily: uniform and differential. 
#Consider the case of assortative mating based on ethnicity. 
#Uniform homophily refers to the tendency for people to form sexual unions with people of the 
#same ethnicity and is the same regardless of which ethnicity is considered. 
#Differential homophily accounts for the fact that the assortative tendencies of 
#sexual partnerships are different depending on the ethnicities of the individuals involved. 
#In the United States, homophilous (or “ethnically-concordant”) partnerships 
#are much more likely among African Americans than among white couples.

## https://rstudio-pubs-static.s3.amazonaws.com/471073_d45a4acd780b4987932dc8fc47c46dd5.html
## https://eehh-stanford.github.io/SNA-workshop/ergm-intro.html#:~:text=ERGMs%20are%20analogous%20to%20logistic,network%20to%20Exponential%20Random%20Graphs.
## https://zalmquist.github.io/ERGM_Lab/ergm-terms.html
## https://seng.netlify.app/meetings/introduction-to-network-simulation-with-statnet/
## https://statnet.org/workshop-advanced-ergm/advanced_ergm_tutorial.pdf

#plot(g.sim)
## Use network.initialize toc create an empty network. 
## Assign values to the ERGM terms to generate networks for various combinations of interest.

### Types of communities in networks: guided by Yang and Leskovich 2012

# 1: Dense within communities. SBM or DCSBM type

# 2: High star values Look at 3,4 and 5 star

# 3: High triangle values. 





## Multiple scenario ground truth

## Relationship between covariates and communities. 

## Terms in ERGM to consider while modelling.

## Dyads

## Triads

## Stars

## Covariates/Attributes


## Likelihood function for Combination of attributes and above terms.


## Bayesian model for network development