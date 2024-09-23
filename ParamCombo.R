k_out <- 3
k_in <- 3
k <- k_out + k_in
## Number of continuous cov
o_out <- 0
o_in <- 0
o <- o_in + o_out
covTypes <- c("binary")
CovNamesLinin <- c()
CovNamesLPin <- paste0("bvin", 1:k_in)
CovNamesLinout <- c()
CovNamesLPout <- paste0("bvout", 1:k_out)
## setting alpha value
nc <- c(3)
alpha <- c(0.0002,0.0005,0.0008,0.001)
alphaLL <- c(1)#c(0.3,0.5, 0.7,0.9,1)
lambda <- 0.001

## Switching to percent increase. Using 0.01% increase
thresh <- c(0.1, 0.01, 0.001, 0.0001)
randomize = c(TRUE)
missing = 0

## Number of nodes
N <- 100
# Probability of cluster assignment and cluster overlap.
pClust <- list(c(0.3, 0.3, 0.3, 0.0, 0.0, 0.0, 0),
                c(0.3, 0.3, 0.3, 0.1, 0.1, 0.1, 0),
                c(0.3, 0.3, 0.3, 0.0, 0.0, 0.0, 0.1),
                c(0.3, 0.3, 0.3, 0.1, 0.1, 0.1, 0.1))

pC1 <- matrix(
  nrow = k,
  ncol = k,
  byrow = TRUE,
  data = c(
    1,0,0,0,0,0,
    0,1,0,0,0,0,
    0,0,1,0,0,0,
    0,0,0,1,0,0,
    0,0,0,0,1,0,
    0,0,0,0,0,1)
)
pC2 <- matrix(
  nrow = k,
  ncol = k,
  byrow = TRUE,
  data = c(
    0.7,0.01,0.01,0.01,0.01,0.01,
    0.01,0.7,0.01,0.01,0.01,0.01,
    0.01,0.01,0.7,0.01,0.01,0.01,
    0.01,0.01,0.01,0.7,0.01,0.01,
    0.01,0.01,0.01,0.01,0.7,0.01,
    0.01,0.01,0.01,0.01,0.01,0.7)
)
pC <- list(pC1)#, pC2)
## this is Pc i.e. the probability of connection between two nodes of the same community.
pConn <- list(c(0.5, 0.5, 0.5))
              #c(0.3,0.4,0.5))
dirPct <-  list(c(0.5, 0.5, 0.5))
               # c(0.75, 0.75, 0.75))
epsilon <- 0.001
nitermax <- 700
Nsim <- 100

EG <- expand.grid(Nsim=Nsim,alpha=alpha,alphaLL=alphaLL,nc = nc,
                  k_in=k_in,k_out=k_out,lambda=lambda,missing=missing,
                  randomize=randomize,N=N,pClust=pClust,pC=pC,pConn=pConn,
                  epsilon=epsilon,dirPct=dirPct,nitermax=nitermax, thresh = thresh)

saveRDS(EG, "InitParam2.rds")

EG <- readRDS("InitParam.rds")
#simvec <- EG[1,]
#simvec$Nsim

## 10 sim run base case
#saveRDS(opf_cov,"BCcov_10.rds")
#saveRDS(opf_noCov,"BCnocov_10.rds")

