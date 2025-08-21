library(terra)
library(sf)
library(tidyverse)
library(dplyr)
library(mapview)
library(sp)
library(usmap)
library(raster)
library(igraph)
library(GGally)
library(ggnet)
library(geosphere)
library(lubridate)
library(gganimate)
library(gifski)
library(png)
library(legendry)
library(RColorBrewer)
library(ggiraph)
library(plotly)
library(network)
library(intergraph)
library(ergm)
library(gridExtra)
library(ggpubr)
library(kableExtra)
library(ggside)
library(grid)
library(ggnewscale)

##This reads the command line arguments
# args=commandArgs(trailingOnly = TRUE)
# if (length(args)==0) {
#   stop("At least one argument must be supplied (input file).n", call.=FALSE)
# }
# 
week <- 102#as.integer(gsub("[\r\n]", "",args[1]))

source(paste0(getwd(),"/Code/CoDACov.R"))

## Read full graph
G_FullGrph <- readRDS(paste0(getwd(),"/Code/G_FullGrph.rds"))

## Aggregate data from 7/11/25
aggr <- readRDS(paste0(getwd(),"/Code/aggr_7_11.rds"))
unique_count <- as.data.frame(table(aggr$weekno))

## geolocations of cities in oregon. Includes column for included or not included wastewater locations
O_attr <- readRDS(paste0(getwd(),"/Code/O_attr.rds"))

test_nc <- c(2,3,4)#,5,6,7,8) ## Can change
N <- vcount(G_FullGrph)
thresh  <- 0.00005
nitermax <- 30000
dir <- "undirected"
alphaLL <- 0.001
alphaOpt <- c(0.00001)#, 0.0001, 0.001, 0.01)
test = FALSE
missing <- NULL
alphaLin <- 0.001
penaltyOpt <- c("Ridge")#,"LASSO","ElasticNet")
seed <- 5
lambdaOpt <- c(0.00001)#, 0.0001, 0.001, 0.01)
covInitOpt <- c("Nmedian")#,"Nmode","Nmean")
printFlg <- FALSE
delta <- getDelta(N)

## setting colors
colors <-  c("cornflowerblue", "coral2","chocolate3", "chartreuse4","blue4",
             "darkorchid", "orangered","gold1","cyan2")

## setting flags for running different comm det and nw metrics
nocovF <- TRUE
covF <- TRUE
nwmetricF <- TRUE
imputedF <- TRUE
plotF <- TRUE
## Read full graph
G_nocov <- G_FullGrph
fg_attr <- O_attr[O_attr$included == "orangered3",]
fg_attr <- fg_attr[fg_attr$name %in% unique(aggr$Location), ]
G_nocov <- set_vertex_attr(G_nocov, name = "X", value = as.matrix(fg_attr$X ))
G_nocov <- set_vertex_attr(G_nocov, name = "Y", value = as.matrix(fg_attr$Y ))
G_FullGrph <- G_nocov


## Save output ##MSE, MSEmd, lambda, week, number of available, penalty, initializing of covariate, alpha 
finOP <- matrix(0, nrow =0, ncol = 12)

## Metric values
op_sim <- matrix(0, nrow= 0, ncol = 14)

## Number of nodes. Saving the number of communities each conde belongs to
numCom <- matrix(0, nrow = vcount(G_FullGrph), ncol = 0)
numCom_nocov <- matrix(0, nrow = vcount(G_FullGrph), ncol = 0)

## Saving the imputed values for the missing covid data
imputedV <- list()
if(plotF == TRUE){
  pdf(paste0("COVID_data_plots",week,".pdf"), width = 10, height = 6)  # Adjust dimensions as needed
  
}
for(nc in test_nc){
  
  specOP <- RegSpectralClust(G_nocov, nc, regularize = TRUE )
  C <- letters[specOP]
  op <- model.matrix( ~ C - 1)
  colnames(op) <- letters[1:nc]
  orig <- specOP
 
  for (i in 1:nc){
    G_nocov <- G_nocov %>% 
      set_vertex_attr(name = letters[i],value = c(op[, i])) 
  }
  
  for(alpha in alphaOpt){
    
    if(nocovF == TRUE){
      op_noCov <- CoDA(G_nocov, nc, k = c(0, 0), o = c(0, 0), N, alpha, lambda,thresh, nitermax, 
                       orig, randomize = TRUE, CovNamesLinin = c(),CovNamesLinout = c(), 
                       CovNamesLPin = c(), CovNamesLPout = c(),dir, alphaLL, test, missing = missing,
                       covOrig ,epsilon = 0, impType = "", alphaLin, penalty, seed, covInit, specOP)
      Fin_mem <- op_noCov$Ffin > delta
      colnames(Fin_mem) <- letters[1:nc]
      
      Fm <- op_noCov$Ffin
      Hm <- op_noCov$Hfin
      NoCovM <- as.data.frame(memOverlapCalc(Fm, Hm, delta, N, nc))
      ## saving number of communities each node belongs to.
      numCom_nocov <- cbind(numCom_nocov, rowSums(NoCovM))
    }
    
    for(covInit in covInitOpt){
      for(penalty in penaltyOpt){
        
        p <- list()
        j <- 1
        
        for(lambda in lambdaOpt){
          
          print(paste0("The week number ", week," avail ", unique_count$Freq[week],
                       " cov init ", covInit, " penalty ", penalty," alpha ", alpha,
                       " lambda ",lambda, " # Clusters ", nc))
          sub_aggr <- aggr %>% filter(weekno == week) %>% ungroup()
          
          ## joining the data and creating na's for where there is no data
          fg_attr <- O_attr[O_attr$included == "orangered3",]
          fg_attr <- fg_attr[fg_attr$name %in% unique(aggr$Location), ]
          fg_attr <- left_join(fg_attr,sub_aggr %>% 
                                 select(CopiesPerul, Location), 
                               by = join_by(name == Location)) %>% distinct()
          
          G_FullGrph <- set_vertex_attr(G_FullGrph,name = "CopiesPerul", value = fg_attr$CopiesPerul )
          cov <- fg_attr$CopiesPerul
          cov[is.na(cov)] <- mean(cov, na.rm=T)
          covOrig <- as.data.frame(cov)
          
          specOP <- CovAssistedSpecClust(G_FullGrph, covOrig , nc , 0.5 )
          C <- letters[specOP]
          op <- model.matrix( ~ C - 1)
          colnames(op) <- letters[1:nc]
          orig <- specOP
          for (i in 1:nc){
            G_FullGrph <- G_FullGrph %>% 
              set_vertex_attr(name = letters[i],value = c(op[, i])) 
          }
          
          G <- G_FullGrph
          N <- nrow(fg_attr)
          if(covF == TRUE){
            op_cov <- CoDA(G, nc, k = c(0, 0), o = c(0, 1), N, alpha,lambda, thresh, nitermax,orig, randomize = FALSE,
                           CovNamesLinin = c(),CovNamesLinout= c("CopiesPerul"),CovNamesLPin = c(),CovNamesLPout= c(), 
                           dir,alphaLL, test, missing = NULL, covOrig,epsilon = 0, impType = "Reg",alphaLin, penalty, 
                           seed, covInit, specOP )
            
          }
          
          if(nwmetricF == TRUE){
            
            ## Code to get network metrics
            Fm <- op_cov$Ffin
            Hm <- op_cov$Hfin
            Q_OL <- ifelse(dir == "directed", "D","U")
            FullM <- as.data.frame(memOverlapCalc(Fm, Hm, delta, N, nc))
            comnty <- data.frame(which(as.matrix(FullM) == 1, arr.ind = TRUE)) %>%
              group_by(col) %>%
              group_map(~.x)
            
            # Find clustering coef global for directed and undirected network
            # High clustering coef means more clique based highly clustered Network
            # Low means more sparse network
            CC <- ClustCoef(G, dir)
            
            ## Power law fit check
            PLFit <- PLfitcheck(G)
            
            # 1: Perfect homophily (nodes only connect within their group).
            # 0: No homophily (random mixing).
            # < 0: Heterophily (nodes prefer different groups).
            NetA <- assortativityF(G,dir)
            
            
            ## Feature structure decomposition
            G <- G %>% set_vertex_attr(name = "cov",value = cov) 
            fsd <- c(unlist(Feature_struct_decomp(G, nc, N, delta, NoCovM, FullM, c("cov") ) ) )
            
            randOI <- null_models(G, nc, k = c(0,0), o = c(0,1), N, alpha,lambda, thresh, nitermax, orig, randomize = FALSE,
                                  CovNamesLinin = c(),CovNamesLinout= c("CopiesPerul"),CovNamesLPin = c(),CovNamesLPout= c(), 
                                  dir,alphaLL, test = FALSE, missing = NULL, covOrig,epsilon = 0,impType = "Reg", alphaLin, 
                                  penalty, seed, covInit, specOP )
            
            op_sim <- rbind(op_sim, c(
              round(conductance(graph = G, communities = comnty, FullM),4),
              round(partition_density(graph = G, communities = comnty),4),
              round(Q_HardPartition(G, FullM),4),
              round(Q_Overlapping(G, FullM, Q_OL),4),
              round(CQ(G, comnty),4),
              round(ONCut(G, comnty),4),
              round(OP(comnty, G),4),
              CC,
              PLFit,
              NetA,
              fsd,
              randOI))          
            
          }
          
          
          delta <- getDelta(N)
          Fin_mem <- op_cov$Ffin > delta
          colnames(Fin_mem) <- letters[1:nc]
          ## saving number of communities each node belongs to.
          numCom <- cbind(numCom, rowSums(Fin_mem))
          if(plotF == TRUE){
            attr <- cbind(as.data.frame(vertex_attr(G)) %>% select(!letters[1:nc]), Fin_mem)
            eList <- as_edgelist(G)
            eList <- as.data.frame(cbind(eList, pair = c(1:nrow(eList)))) %>% 
              pivot_longer(!pair, values_to = "name",names_to = "gbg") %>%
              left_join(attr %>% select(name, X, Y))
            attr$Available <- !is.na(attr$CopiesPerul)
            q <- ggplot()+ 
              geom_line(data = eList, aes(x = X, y = Y, group = pair), color = "grey80") +
              geom_point(data = attr, aes(x = X, y = Y, color = Available)) +
              scale_color_manual(values = c("#661100", "#CC6677")) +
              new_scale_color()
            for(k in 1:nc){
              q <- q + ggforce::geom_mark_hull(data = attr[attr[,letters[k]] == TRUE,],
                                      aes(x = X, y = Y,fill = .data[[letters[k]]] ), 
                                      concavity = 3, expand = unit(2, "mm"),
                                      alpha = 0.1, fill = colors[k], color = colors[k])
            }
              q <- q+xlab("Latitude") + ylab("Longitude")+
              theme(legend.position = "bottom", panel.background = element_blank()) 
              p[[j]] <- q
            
            j <- j + 1
          }
          
          
          ## filling the final output
          #MSE, MSEmd, lambda, week, number of available, penalty, initializing of covariate, alpha 
          finOP <- rbind(finOP, c(c(tail(op_cov$Loglik,1)) ,tail(op_cov$MSE,1),tail(op_cov$MSEMD,1),
                                  week,unique_count$Freq[week],penalty,covInit,alpha,lambda))
          
          if(imputedF == TRUE){
            imputedV <- append(imputedV, list(op_cov$Zout_cov))
          }
          
        }
        if(plotF == TRUE){
          pt <- ggarrange(plotlist = p, nrow = 2, ncol = 2)
          fig <- annotate_figure(pt, top = text_grob(paste0("Week ", week, 
                                                            " avail ", unique_count$Freq[week], 
                                                            " penalty ", penalty, 
                                                            " alpha ", alpha,
                                                            " nc ",nc)))
          # Convert to grob and draw
          grob <- ggplotGrob(fig)
          grid.newpage()
          grid.draw(grob)
        }
      }
    }
  }
  
}
if(plotF == TRUE){
  dev.off()
}
colnames(finOP) <- c("LLTot", "LLGrph","LLBCov","LLCCov","MSE","MSEMD","Week","Freq","Penalty","covInit","Alpha","Lambda")
colnames(op_sim) <- c("Cond","PD","Q_HP","Q_OL","CQ","OnCut","OP","CC","PLFit","NetA","oi_b1","oi_b2", "oi_cov","randOI")

saveRDS(finOP, paste0("realCOVID_DataOP", week,".rds"))
saveRDS(numCom, paste0("realCOVID_NumberOfCommunities",week,".rds"))
saveRDS(numCom_nocov, paste0("realCOVID_NumberOfCommunitiesNoCov",week,".rds"))
saveRDS(op_sim, paste0("realCOVID_NetworkMetricsOutput",week,".rds"))
saveRDS(imputedV, paste0("ImputedValues",week,".rds"))

# readRDS(paste0("realCOVID_DataOP", week,".rds"))
# readRDS(paste0("realCOVID_NumberOfCommunities",week,".rds"))
# readRDS(paste0("realCOVID_NumberOfCommunitiesNoCov",week,".rds"))
# readRDS(paste0("realCOVID_NetworkMetricsOutput",week,".rds"))
# readRDS(paste0("ImputedValues",week,".rds"))

