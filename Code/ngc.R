#' Plot DAG of network or ring graph showing Granger causality
#' @param fit object of class ngc
#' @param ngc.type type of plot to create
plot.ngc <- 
  function(
    fit, #object of class ngc
    ngc.type = "dag", #"dag" or "granger"
    sparsify = FALSE #whether to remove covariates with no edges in the high-dimensional case
  ){
    if (class(fit) != "ngc")
    {
      stop("Class of argument must be ngc")
    }
    p <- fit$p
    d <- fit$d
    if (ngc.type == "dag" & p > 20)
    {
      cat("Warning: plot.ngc is not designed for plotting networks of more than 20 covariates.")
      if (sparsify)
      {
        cat("plot.ngc will remove covariates with no edges in order to plot the network.")
      }
    }
    covNames <- fit$covNames
    group <- fit$group
    if (ngc.type == "granger")
    {
      plot_ring(fit$ring, p, d)
    }
    else
    {
      plot_network(fit$dag, p, d, group, fit$fitMethod == "lasso", sparsify)
    }
    if (!is.null(covNames))
    {
      legend(d+1.5, p, paste(1:p, covNames, sep = " - "), cex = labelCex, ncol = p%/%10 + 1, title = "Legend")
    }
  }

plot_ring <- function(g, p, d)
{
  if (is.null(E(g)$weight))
  {
    edgeThickness = 0
  }
  else
  {
    edgeThickness = E(g)$weight^2/mean(E(g)$weight^2)
  }
  #control maximum and minimum thickness
  edgeThickness <- ifelse(edgeThickness > 0.2, edgeThickness, 0.2)
  edgeThickness <- ifelse(edgeThickness < 5, edgeThickness, 5)
  labelCex <- max(min(10/p, 1), 0.3)
  arrowSize <- 0.5*labelCex
  plot(g, layout = layout_in_circle(g), vertex.label.cex = labelCex,
       edge.arrow.size = arrowSize, vertex.shape = "none", edge.width = edgeThickness)
}



plot_network <- function(g, p, d, group =  NULL, signed = TRUE, sparsify = FALSE, title = NULL)
{
  edgeTails <- tail_of(g, E(g))
  edgeHeads <- head_of(g, E(g))
  if (sparsify)
  {
    nodes <- as.numeric(unique(c(edgeHeads%%p, edgeTails%%p)))
    nodes[nodes==0] <- p
    nodes <- sort(nodes)
    toRemove <- (1:p)[-nodes] + rep(seq(0, p*d, p), each = p - length(nodes))
    g <- delete_vertices(g, toRemove)
    p <- length(nodes)
    vertexLabels <- rep(nodes, d+1)
  }
  else
  {
    vertexLabels <- rep(1:p, d+1)
  }
  xcoords = rep(1:(d+1), each=p)
  ycoords = rep(p:1, d+1) 
  layout_matrix = matrix(c(xcoords, ycoords), ncol=2)
  groupList <- NULL
  if (!is.null(group))
  {
    groupList <- lapply(unique(group),function(x){which(group==x)})
  }
  par(mar=c(2.5, 2.5, 2.5, 2.5))
  if (signed)
  {
    edgeColor = ifelse(E(g)$weight > 0, "blue", "red")
  }
  else
  {
    edgeColor = NULL
  }
  if (is.null(E(g)$weight))
  {
    edgeThickness = 0
  }
  else
  {
    edgeThickness <- E(g)$weight^2/mean(E(g)$weight^2)
  }
  #control maximum and minimum thickness
  edgeThickness <- ifelse(edgeThickness > 0.2, edgeThickness, 0.2)
  edgeThickness <- ifelse(edgeThickness < 5, edgeThickness, 5)
  labelCex <- max(min(10/p, 1), 0.3)
  arrowSize <- 0.5*labelCex
  #curve edges that are more than 1 lag
  edgeCurvature <- (edgeTails <= p*(d-1))*0.25
  edgeCurvature <- edgeCurvature*(-1)^((head_of(g, E(g)) %% p) < (edgeTails %% p))
  aRatio <- ((d+3)/p)/2
  plot(g, asp = aRatio, layout = layout_matrix, main = title,
       mark.groups = groupList, mark.border = NA,
       vertex.label.cex = labelCex,
       vertex.label = vertexLabels, vertex.shape = "none",
       edge.color = edgeColor, edge.width = edgeThickness,
       edge.arrow.size = arrowSize, edge.curved = edgeCurvature,
       rescale = FALSE, xlim = c(0, d+2), ylim = c(0, p))
  text(0, -0.25*labelCex, "Lag", cex = labelCex)
  lagStep <- ifelse(d < 10, 1, 5)
  for (i in seq(lagStep, d, lagStep))
  {
    text(i, -0.25*labelCex, d-i+1, cex = labelCex)
  }
}


plot_ring <- function(g, p, d)
{
  if (is.null(E(g)$weight))
  {
    edgeThickness = 0
  }
  else
  {
    edgeThickness = E(g)$weight^2/mean(E(g)$weight^2)
  }
  #control maximum and minimum thickness
  edgeThickness <- ifelse(edgeThickness > 0.2, edgeThickness, 0.2)
  edgeThickness <- ifelse(edgeThickness < 5, edgeThickness, 5)
  labelCex <- max(min(10/p, 1), 0.3)
  arrowSize <- 0.5*labelCex
  plot(g, layout = layout_in_circle(g), vertex.label.cex = labelCex,
       edge.arrow.size = arrowSize, vertex.shape = "none", edge.width = edgeThickness)
}


############### trying out ngc function #############################

ngc <-
  function(
    X, #input array dim=(n,p,T) (longitudinal), or (p,T) (time series); last time=Y
    d = NULL, #number of time lags to consider
    fitMethod = "lasso", #"lasso" for Granger lasso fit; "sam" or "hierbasis" for sparse additive modeling
    estMethod = "regular", #method to use. options are "regular", "truncate", and "threshold"
    group = NULL, #vector of group indices of length p; if null, no group structure
    groupByTime = FALSE, #whether to group covariates across time points
    typeIerr = NULL, #acceptable type I error rate; if provided, error-based lasso is fitted
    typeIIerr = 0.1, #acceptable type II error rate
    weights = NULL, #wmatrix of weights for Alasso. If no weights are provided, use regular lasso.
    thresholdConstant = NULL, #constant used for calculating threshold value
    refit = FALSE, #whether to refit a linear regression after initial thresholding
    covNames = NULL #covariate names
  ){
    ####### START OF FN #######
    if (fitMethod != "lasso" & fitMethod != "sam" & fitMethod != "hierbasis")
    {
      stop("Invalid fit method specified")
    }
    
    if (estMethod != "regular" & estMethod != "truncate" & estMethod != "threshold")
    {
      stop("Invalid estimation method specified")
    }
    
    if (fitMethod != "lasso" & estMethod != "regular")
    {
      stop("This estimation method is not supported for additive modeling")
    }
    
    if (!is.array(X) | length(dim(X)) <2 | length(dim(X)) > 3)
    {
      stop("Invalid X")
    }
    
    if (length(dim(X))==2)
    {
      n <- 1
      p <- dim(X)[1]
      len <- dim(X)[2]
    }
    else
    {
      n <- dim(X)[1]
      p <- dim(X)[2]
      len <- dim(X)[3]
    }
    
    if (is.null(d))
    {
      d <- len-2
    }
    
    if (d >= len)
    {
      stop("Number of time lags to consider cannot exceed number of time points")
    }
    
    #Set up replicates for the time series case
    #Put X into array format
    #The transformed matrix has (len-d) replicates over d+1 time points
    if (n == 1)
    {
      if (d >= len-1)
      {
        stop("Number of time lags to consider must be restricted in order to fit time series")
      }
      cat("Warning: stationarity assumption is required for time series data")
      xMat <- X
      n <- len-d
      len <- d+1
      X <- array(0, c(n,p,len))
      for (i in 1:n)
      {
        X[i,,] <- xMat[,i:(d+i)]
      }
    }
    
    if (!is.null(covNames))
    {
      if (length(covNames) != p)
      {
        stop("Number of covariate names must match number of covariates")
      }
    }
    
    if (!is.null(group))
    {
      if (length(group)!=p)
      {
        stop("Invalid group specification")
      }
      if (!is.numeric(group))
      {
        stop("Groups must be specified with consecutive integers")
      }
      if (!all.equal(order(group), 1:p))
      {
        stop("Groups must be specified with consecutive integers")
      }
      # apply groups across time points
      if (groupByTime)
      {
        group <- rep(group, d)
      }
      else
      {
        ngrp = length(unique(group))
        group <- group + rep(seq(0, (d-1)*ngrp, by = ngrp), each = p)
      }
    }
    
    if (fitMethod == "lasso")
    {
      if (estMethod == "regular")
      {
        fit <- grangerLasso(X, d = d, group = group, typeIerr = typeIerr,
                            weights = weights) 
      }
      else if (estMethod == "truncate")
      {
        fit <- grangerTlasso(X, d = d, group = group, typeIerr = typeIerr,
                             typeIIerr = typeIIerr, weights = weights)
      }
      
      else #threshold
      {
        fit <- grangerThrLasso(X, d = d, group = group, typeIerr = typeIerr,
                               typeIIerr = typeIIerr, weights = weights,
                               thresholdConstant = thresholdConstant,
                               refit = refit)
      }
    }
    else if (fitMethod == "sam")
    {
      fit <- grangerSam(X, d = d)
    }
    else if (fitMethod == "hierbasis")
    {
      fit <- grangerHierBasis(X, d = d)
    }
    
    dagMat <- Matrix(0, nrow=p*(d+1), ncol=p*(d+1), sparse = TRUE)
    ringMat <- Matrix(0, nrow=p, ncol=p)
    edgeIx <- which(fit$estMat != 0, arr.ind = T)
    edgeCount <- dim(edgeIx)[1]
    if (is.null(fit$tsOrder))
    {
      tsOrder <- ifelse(edgeCount > 0, max(d-edgeIx[,3]+1), 0)
      fit$tsOrder <- ifelse(!is.null(tsOrder), tsOrder, 0)
    }
    if (edgeCount > 0)
    {
      for (i in 1:edgeCount)
      {
        edge <- edgeIx[i,]
        pStart <- edge[2]
        pEnd <- edge[1]
        lag <- edge[3]
        dagMat[((lag-1)*p + pStart),(d*p + pEnd)] <- fit$estMat[pEnd, pStart, lag]
        ringMat[pStart, pEnd] <- ringMat[pStart, pEnd] + fit$estMat[pEnd, pStart, lag]^2
      } 
    }
    fit$dag <- graph_from_adjacency_matrix(dagMat, mode = "directed", weighted = TRUE)
    fit$ring <- graph_from_adjacency_matrix(sqrt(ringMat), mode = "directed", weighted = TRUE)
    fit$fitMethod <- fitMethod
    fit$estMethod <- estMethod
    fit$n <- n
    fit$p <- p
    fit$len <- len
    fit$d <- d
    fit$group <- group
    fit$X <- X
    fit$covNames <- covNames
    class(fit) <- "ngc"
    return(fit)
  }

grangerLasso <-
  function(
    X, 				  #input array dim=(n,p,T) (longitudinal), or (p,T) (time series); last time=Y
    d = NULL, 	  #number of time lags to consider
    group = NULL, #group indices
    penCoefMethod = 'errBased', #choose between errBased and direct
    typeIerr = 0.10, 		  #sig level for lambda (...Method=errBased)
    lambda = NULL,		  #value of lambda (...Method=direct)
    weights = NULL,			  #matrix of weights for Alasso
    Alasso.power = 1,		  #power for Adaptive lasso
    eps = 1e-8			  #used for setting small values to zero
  ){
    ####### START OF FN #######
    n <- dim(X)[1]
    p <- dim(X)[2]
    tp <- dim(X)[3]
    
    useAlasso <- !is.null(weights)
    Onep = matrix(1,p,(p*d))
    estMat <- array( 0, c(p, p, d) )
    
    #scale the X matrix
    for (i in 1:(tp-1))
    {
      X[,,i] <- scale( X[,,i] )*sqrt(n/(n-1))
    }
    
    #first put all X matrices into one big matrix
    XX <- array2mat(X[,,(tp-d):(tp-1)])
    YY <- X[,,tp]
    
    if (!useAlasso)
    {
      temp = pldag.set(XX, YY, group = group, sigLevel=typeIerr, wantScale=TRUE)
    }
    else
    {
      W = abs( array2mat(weights[,,1:d]) )
      W = ( W + eps*Onep ) ^ (-Alasso.power)
      if (sum(W < 1) > 0)
      {
        W[(W < 1)] = 1
      }
      
      temp = pldag.set(XX, YY, group = group, sigLevel=typeIerr,
                       useWghts=TRUE, wghts=W, wantScale=TRUE)
    }
    AA <- as.matrix(temp$AA)
    lambda <- temp$lambda
    sigma <- temp$sigma
    intercepts <- temp$intercepts
    ##Put the matrix output of pldag.set into an array to make it
    ##compatible with other parts of the code
    for(i in 1:p)
    {
      estMat[i,,] = AA[i,]
    }
    rm(temp)
    
    return(list(estMat = estMat, lambda = lambda, sigma = sigma, intercepts = intercepts))
  }


pldag.set <-
  function(
    X1,				##predictor set, nxp1 matrix
    X2, 				##outcome set, nxp2 matrix (X2~X1)
    group = NULL, #group indices
    sigLevel = NULL, 		##siglevel for MB & BdE penalties
    useWghts = FALSE,		##use wghts for the penalty?
    ##if TRUE, wghts should be provided
    wghts=NULL, 		##wghts for penalty, e.g. Alasso
    excld=NULL,			##which variabels to exclude (this
    ##is equiv to setting wght=Inf)
    wantScale = FALSE,	##whether to scale
    useNaive=TRUE,		##see below!
    useTLASSO=FALSE,		##is truncating lasso used?
    d = NULL			##number of time lags to cosider
  ){
    ####### START OF FN #######
    method <- "naive"		#method used in glmnet (see glmnet help)
    if (!useNaive)
    {
      method <- "covariance"
    }
    
    if (dim(X1)[1] != dim(X2)[1])
    {
      stop("Number of observations does not  match")
    }
    n <- dim(X2)[1]
    
    p1 <- dim(X1)[2]
    p2 <- dim(X2)[2]
    prdctrIndx = 1:p1
    
    if (is.null(wghts))
    {
      if (useWghts)
      {
        cat("WARNING: No weights provided, using weights=1", "\n")
      }
      if (is.null(group))
      {
        wghts <- matrix(1,p2,p1)
      }
      else
      {
        wghts <- matrix(rep(as.vector(sqrt(table(group))), p2), nrow=p2)
      }
    }
    
    ##calculate penalty coefficient (lambda)
    nvar <- p2		#no of variables, to use in calculation of lambda
    ncov <- p1		#no of covariates, to use in calculation of lambda
    
    #in case of TLASSO, one set of values of X is given at each time, but
    #ncov needs to be adjusted based on the truncating pealty
    if (useTLASSO)
    {
      if (is.null(d))
      {
        stop('Number of effective time lags needed for TLASSO')
      }
      else
      {
        ncov <- d*p1
      }
    }
    
    AA <- matrix(0, p2, p1)
    
    if ( is.null(excld) )
    {
      excld <- matrix(FALSE,p2,p1)
    }
    
    if ( (dim(excld)[1]!= dim(AA)[1]) || (dim(excld)[2]!= dim(AA)[2]) )
    {
      stop("Wrong dimension for variables to exclude")
    }
    
    lambda <- NULL
    if(!is.null(sigLevel))
    {
      lambda <- (1/sqrt(n))*qnorm(1-sigLevel/(2*ncov*nvar))
    }
    
    lambdas <- rep(NA, p2)
    sigmas <- rep(NA, p2)
    sigmas2 <- rep(NA, p2)
    intercepts <- rep(0, p2)
    ##main estimation loop
    for (i in 1:p2)
    {
      y <- X2[ , i]
      ww <- wghts[i, ]
      
      temp <- excld[i, ]
      excldIndx <-prdctrIndx[temp]
      
      if (length(excldIndx) < p1-1)
      {
        #estimate sigma with deviance
        if (!is.null(lambda))
        {
          if (is.null(group))
          {
            fit1 <- glmnet(X1, y, lambda = lambda*sd(y)*sqrt((n-1)/n), penalty.factor = ww,
                           standardize = wantScale, exclude = excldIndx)
          }
          else
          {
            fit1 <- gglasso(X1, y, group = sort(group), loss="ls", 
                            lambda = lambda)
          }
          betas <- coef(fit1)
          intercepts[i] <- betas[1]
          betas <- betas[-1]
          dev <- deviance(fit1)
          if (!is.null(dev))
          {
            residDf <- max(n-1-sum(betas!=0), 1)
            sigmas[i] <- sqrt(dev/residDf)
          }
          else
          {
            sigmas[i] = 1
          }
        }
        else #estimate sigma with cvm
        {
          if (is.null(group))
          {
            fit1 <- cv.glmnet(X1, y, penalty.factor = ww,
                              standardize = wantScale, exclude = excldIndx)
          }
          else
          {
            fit1 <- cv.gglasso(X1, y, group = sort(group), pred.loss = "L1", pf  = ww)
          }
          lambdas[i] <- fit1$lambda.1se
          betas <- coef(fit1, s="lambda.1se")
          intercepts[i] <- betas[1]
          betas <- betas[-1]
          sigmas[i] <- mean(sqrt(fit1$cvm))
        }
        if (length(betas) > 0)
        {
          AA[i, ] <- betas
        }
        rm(fit1)
        rm(betas)
      }
    }
    if (is.null(lambda))
    {
      lambda <- lambdas
    }
    return(list(AA = AA, lambda = lambda, sigma = mean(sigmas, na.rm = TRUE), intercepts = intercepts))
  }


