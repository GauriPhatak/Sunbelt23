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
week <- 34#as.integer(gsub("[\r\n]", "",args[1]))
getwd()
source(paste0(getwd(),"/Code/CoDACov.R"))

## Read full graph
G_FullGrph <- readRDS(paste0(getwd(),"/Code/MapNetworks/G_FullGrph45Km_WS_SR_RC.rds"))

## Cities to be set missing
MissCities <- c("None")#readRDS(paste0(getwd(),"/Code/CitiesofImport.rds"))

## Aggregate data from 7/11/25
aggr <- readRDS(paste0(getwd(),"/Code/MapNetworks/aggr_04-10.rds"))
unique_count <- as.data.frame(table(aggr$weekno))

## geolocations of cities in oregon. Includes column for included or not included wastewater locations
O_attr <- readRDS(paste0(getwd(),"/Code/MapNetworks/O_attr_WS_SR_RC.rds"))

test_nc <- c(2)#,3,4,5,6)#c(2,3,4,5,6) ## Can change
N <- vcount(G_FullGrph)
thresh  <- 0.00005
nitermax <- 30000
dir <- "undirected"
alphaLL <- 0.001
alphaOpt <- c(0.00001)#,0.0001, 0.001, 0.01, 0.1)#, 0.0001)
test = FALSE
missing <- NULL
alphaLin <- 0.001
penaltyOpt <- c("LASSO")
seed <- 5 #sample(1:100000, 1)
lambdaOpt <- c(0.00001)#, 0.0001, 0.001, 0.01, 0.1)
covInitOpt <- c("Nmean")
printFlg <- FALSE
delta <- getDelta(N)

## setting colors
colors <-  c("cornflowerblue", "coral2","chocolate3", "chartreuse4","blue4",
             "darkorchid", "orangered","gold1","cyan2")

## setting flags for running different comm det and nw metrics
nocovF <- FALSE
covF <- TRUE
nwmetricF <- FALSE
imputedF <- TRUE
plotF <- TRUE
## Read full graph
G_nocov <- G_FullGrph
fg_attr <- O_attr[O_attr$included == "orangered3",]
fg_attr <- fg_attr[fg_attr$name %in% unique(aggr$Location), ]
G_nocov <- set_vertex_attr(G_nocov, name = "X", value = as.matrix(fg_attr$X ))
G_nocov <- set_vertex_attr(G_nocov, name = "Y", value = as.matrix(fg_attr$Y ))
G_FullGrph <- G_nocov

## setting up metrics dataframe
cn <- c("MeanConductanceW","MeanConductanceWo","WeightedMeanConductanceW","WeightedMeanConductanceWo",
        "AvgVarianceW","DispersionScoreW","AvgVarianceWo","DispersionScoreWo","WeightedAvgVarianceW","WeightedDispersionScoreW",
         "numNodesWoAssignment" ,"nc","dir", "alpha", "lambda","penalty","covInit","MissingCity")
metricsCov <- matrix(0, nrow =0,ncol = length(cn))

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
  
  pdf(paste0("COVID_data_plotsmeanLogCopiesPerL",week,".pdf"), width = 10, height = 6)  # Adjust dimensions as needed
  
}
getwd()
## set up list for final community assignments
finMCov <- list()
finMNoCov <- list()


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
      op_noCov <- CoDA(G_nocov, nc, k = c(0, 0), o = c(0,0), N, alpha,lambda_lin=0, lambda_bin=0, thresh, nitermax, orig, randomize =TRUE ,
                       CovNamesLinin= c(), CovNamesLinout= c(), CovNamesLPin= c(), CovNamesLPout= c(), 
                       dir,alphaLL, test, missing =NULL, covOrig =NULL, epsilon =0, impType = "Reg", alphaLin, penalty=NULL, 
                       seed, covInit=NULL, specOP,nc_sim = nc,lambda_grph=0 )
      
      Fin_mem <- op_noCov$Ffin > delta
      colnames(Fin_mem) <- letters[1:nc]
      ## saving nocov output
      finMNoCov[[length(finMNoCov)+1 ]] <- Fin_mem
      
      Fm <- op_noCov$Ffin
      Hm <- op_noCov$Hfin
      NoCovM <- as.data.frame(memOverlapCalc(Fm, Hm, delta, N, nc))
      ## saving number of communities each node belongs to.
      numCom_nocov <- cbind(numCom_nocov, rowSums(NoCovM))
      MissCity <- "None"
    }
    
    for(covInit in covInitOpt){
      for(penalty in penaltyOpt){
        
        p <- list()
        j <- 1
        
        for(lambda_lin in lambdaOpt){
          
          print(paste0("The week number ", week," avail ", unique_count$Freq[week],
                       " cov init ", covInit, " penalty ", penalty," alpha ", alpha,
                       " lambda ",lambda_lin, " # Clusters ", nc))
          sub_aggr <- aggr %>% filter(weekno == week) %>% ungroup()
          
          ## joining the data and creating na's for where there is no data
          fg_attr <- O_attr[O_attr$included == "orangered3",]
          fg_attr <- fg_attr[fg_attr$name %in% unique(aggr$Location), ]
          fg_attr <- left_join(fg_attr,sub_aggr %>% 
                                 select(meanLogCopiesPerL, Location), 
                               by = join_by(name == Location)) %>% distinct()
          #fg_attr$meanLogCopiesPerL <- NA
          #fg_attr$meanLogCopiesPerL[fg_attr$name %in% sub_aggr$Location] <- sub_aggr$meanLogCopiesPerL
          
          for(MissCity in MissCities){
            
            ## setting a particular city missing.
            attr_tmp <- as.data.frame(vertex_attr(G_FullGrph)) %>% select(name, X, Y) #fg_attr$meanLogCopiesPerL
            attr_tmp <- attr_tmp %>% left_join(fg_attr %>% select(name, meanLogCopiesPerL), by = join_by(name == name))
            attr_tmp$meanLogCopiesPerL[attr_tmp$name == MissCity] <- NA
            #fg_attr$meanLogCopiesPerL[fg_attr$name %in% sub_aggr$Location] <- sub_aggr$meanLogCopiesPerL
            #fg_attr$meanLogCopiesPerL[fg_attr$name == MissCity] <- NA
            
            G_FullGrph <- set_vertex_attr(G_FullGrph, name = "meanLogCopiesPerL", value = attr_tmp$meanLogCopiesPerL )
            cov <- fg_attr$meanLogCopiesPerL
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
            
            ## Origignal covariate list first binary then covariate
            covOrigLst <- list(NULL, covOrig)
            if(covF == TRUE){
              op_cov <- CoDA(G, nc, k = c(0, 0), o = c(0,1), N, alpha,lambda_lin, lambda_bin=0, thresh, nitermax, orig, randomize =TRUE ,
                   CovNamesLinin= c(), CovNamesLinout= c("meanLogCopiesPerL"), CovNamesLPin= c(), CovNamesLPout= c(), 
                   dir,alphaLL, test, missing =NULL, covOrigLst, epsilon =0, impType = "Reg", alphaLin, penalty, 
                   seed, covInit, specOP,nc_sim=nc,lambda_grph=0 )
            }
            
            if(nwmetricF == TRUE){  }
            
            
            delta <- getDelta(N)
            Fin_mem <- op_cov$Ffin > delta
            colnames(Fin_mem) <- letters[1:nc]
            ## saving number of communities each node belongs to.
            numCom <- cbind(numCom, rowSums(Fin_mem))
            
            ## Saving the community assignments
            finMCov[[length(finMCov) + 1]] <- Fin_mem 
            
            if(plotF == TRUE){
              attr <- cbind(as.data.frame(vertex_attr(G)) %>% select(!letters[1:nc]), Fin_mem)
              eList <- as_edgelist(G)
              eList <- as.data.frame(cbind(eList, pair = c(1:nrow(eList)))) %>% 
                pivot_longer(!pair, values_to = "name",names_to = "gbg") %>%
                left_join(attr %>% select(name, X, Y))
              attr$Available <- !is.na(attr$meanLogCopiesPerL)
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
                                    week,unique_count$Freq[week],penalty,covInit,alpha,lambda_lin))
            
            if(imputedF == TRUE){
              imputedV <- append(imputedV, list(op_cov$Zout_cov))
            }
            
            ## Calculating the metrics
            d <- op_cov
            G <- G_FullGrph
            Fm <- op_cov$Ffin
            Hm <- op_cov$Hfin
            C <- memOverlapCalc(Fm,Hm,delta, N, nc)
            epsilon <- 0.000001
            ##Create a new community of background nodes that as unassigned
            NodesWoAssignment <- (rowSums(memOverlapCalc(Fm,Hm,delta, N, nc)) == 0)
            numNodesWoAssignment <- sum(NodesWoAssignment)
            
            ## find the number of nodes with degree 0 or 1
            degree0 <- sum(igraph::degree(G) == 0 )
            degree1 <- sum(igraph::degree(G) == 1 )
            degree01 <- degree0 + degree1
            
            conduct <- EgoSplitConductance(G , as.matrix(C), dir, degree01, N, nc, numNodesWoAssignment)
            ADS     <- AverageDissimilarityScore(covOrig, C, degree01, numNodesWoAssignment)
            
            newOI   <- OmegaIdx(G,as.matrix(C),N,nc,nc)
            metricsCov <- rbind(metricsCov, c(tail(conduct,4),ADS,numNodesWoAssignment, nc, dir, alpha, lambda_lin,penalty,covInit, MissCity))
          }
          
        }
        if(plotF == TRUE){
          pt <- ggarrange(plotlist =  p,  ncol = 2, nrow = 2)
         # pt <- cowplot::plot_grid(p, ncol =2)
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


# colnames(finOP) <- c("LLTot", "LLGrph","LLBCov","LLCCov","MSE","MSEMD","Week","Freq","Penalty","covInit","Alpha","Lambda")
# #colnames(op_sim) <- c("Cond","PD","Q_HP","Q_OL","CQ","OnCut","OP","CC","PLFit","NetA","oi_b1","oi_b2", "oi_cov","randOI")
# colnames(metricsCov) <- cn
# 
# saveRDS(metricsCov,paste0(getwd(),"/Code/COVIDNWOutput/MetricsPO", week,".rds") )
# saveRDS(finOP, paste0(getwd(),"/Code/COVIDNWOutput/realCOVID_DataOPmeanLogCopiesPerL", week,".rds"))
# saveRDS(numCom, paste0(getwd(),"/Code/COVIDNWOutput/realCOVID_NumberOfCommunitiesmeanLogCopiesPerL",week,".rds"))
# saveRDS(numCom_nocov, paste0(getwd(),"/Code/COVIDNWOutput/realCOVID_NumberOfCommunitiesNoCovmeanLogCopiesPerL",week,".rds"))
# #saveRDS(op_sim, paste0(getwd(),"/Code/COVIDNWOutput/realCOVID_NetworkMetricsOutputmeanLogCopiesPerL",week,".rds"))
# saveRDS(imputedV, paste0(getwd(),"/Code/COVIDNWOutput/ImputedValuesmeanLogCopiesPerL",week,".rds"))
# saveRDS(finMCov, paste0(getwd(),"/Code/COVIDNWOutput/Fin_mem_Cov",week,".rds"))
# saveRDS(finMNoCov, paste0(getwd(),"/Code/COVIDNWOutput/Fin_mem_NoCov",week,".rds"))
# 
# 
