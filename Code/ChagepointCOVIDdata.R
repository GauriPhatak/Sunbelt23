library(tidyverse)
library(mcp)
library(changepoint)
library(strucchange)
library(readxl)
library(RColorBrewer)
library(imputeTS)
library(fs)
library(glmnet)
#library(NGC)
library(igraph)
library(influential)
library(simts)
library(tseries)


nodeImpfunc <- function(nw){
  
  colnames <-  c("Section","Indegree","Outdegree","Betweenness",
                 "NeighConnect_in","NeighConnect_out", 
                 "H_Index_in","H_Index_out","Coll_inf_in","Coll_ing_out",
                 "ivi_in","ivi_out","ivi_all","Closeness","Eigen",
                 "Strength_in","Strength_out","City") #"Closeness","EigenCentrality",
  Node_imp <- data.frame(matrix(ncol =length(colnames), nrow = 0))
  
  
  for (i in 1:(Nbrk+1)) {
    nw[[i]] <- igraph::simplify(nw[[i]],  remove.multiple = FALSE )
    plot.igraph(igraph::simplify(nw[[i]]),directed = T,
                main = paste("subsection: ", i),
                edge.arrow.size = .2,
                vertex.label = V(nw[[i]])$names,
                vertex.size = 5,
                vertex.color = c("lightblue"),
                vertex.frame.color = "blue",
                vertex.label.size=0.001)
    
    ## gathering some ndoe importance measures from all the nodes
    Node_imp <- rbind(Node_imp, 
                      cbind(i,
                            igraph::degree(nw[[i]], mode = "in"),
                            igraph::degree(nw[[i]], mode = "out"),
                            round(igraph::betweenness(nw[[i]],normalized = TRUE ),2),
                            round(neighborhood.connectivity(nw[[i]], mode = "in"),2),
                            round(neighborhood.connectivity(nw[[i]], mode = "out"),2),
                            h_index(nw[[i]], mode= "in"),
                            h_index(nw[[i]], mode= "out"),
                            collective.influence(nw[[i]],mode = "in"),
                            collective.influence(nw[[i]],mode = "out"),
                            
                            ### Calculating the integrated value of influence from a graph
                            round(ivi(nw[[i]], directed = TRUE, mode = "in"),2),
                            round(ivi(nw[[i]], directed = TRUE, mode = "out"),2),
                            round(ivi(nw[[i]], directed = TRUE, mode = "all"),2),
                            igraph::closeness(nw[[i]], mode = "out"),
                            igraph::eigen_centrality(nw[[i]], directed = TRUE)$vector,
                            igraph::strength(nw[[i]], vids = V(nw[[i]]), mode = "in", weights = E(nw[[i]])$weight),
                            igraph::strength(nw[[i]], vids = V(nw[[i]]), mode = "out", weights = E(nw[[i]])$weight),
                            V(nw[[i]])$names
                      )
    )
    
  }
  colnames(Node_imp) <- colnames
  return(Node_imp)
}

LassoVar <- function(df_ip, d, Nbrk, retFit = FALSE){
  #edgeIx <- as.data.frame(matrix(nrow = 0, ncol = 4))
  if(retFit ==TRUE){
    grphs <- list()
    fits <- list()
  }
  for (m in 1:(Nbrk+1)) {
    print(paste0("Subsection number ,", m))
    k <- df_ip[[m]] 
    cityNames <- colnames(k)
    k <- t(k)
    if(any(is.na(k))){
      print("k has missing values")
    }
    fit1 = NGC::ngc(k, d)
    
    if(retFit == TRUE){
      V(fit1$ring)$names <- cityNames
      grphs[[m]] <- fit1$ring
      fits[[m]] <- fit1
      
    }
    
  }
 
  return(list(grphs, fits))
}

dat <- read_excel(paste0("C:/Users/gauph/Documents/StatisticsMS_PhD/Wastewater-Surveillance-OSU/Sunbelt23/Code/Combined_aggregated_OHA_normalized_2026-04-24.xlsx"),
                  sheet = "COVID",
                  guess_max = 10000)

df_n <- dat  %>% 
  subset(select = c(Sample_Date,LogCopiesPerL,LogCopiesPerDayPerPerson_NormtoFlowPopRec, Location, County, Site))

weekly_summary <- df_n %>%
  mutate(
    Sample_Date = as.Date(Sample_Date),
    week_start = floor_date(Sample_Date, "week", week_start = 1)
  ) %>%
  group_by(Location, week_start) %>%
  summarise(
    mean_logcopiesperL = mean(LogCopiesPerL, na.rm = TRUE),
    .groups = "drop"
  )
df <- weekly_summary %>%
  pivot_wider( names_from = Location, values_from = c(mean_logcopiesperL, -week_start) )

## nas by date
nas_date <- data.frame(week = 1:(dim(df)),
                       val = as.numeric(rowSums(!is.na(df[,-c(1)]))))
nas_date %>% ggplot()+geom_point(aes(x = week, y = val)) +geom_line(aes(x = week, y = val), color = "gray40")

df <- df[11:290,]

svfldr <- "C:/Users/gauph/Box/FinalTestingParamComboFiles/FinalPlots_TS/"

## Find locations that have more than 50% missing values.

## Find the counties for each lcoation
df_county <- na.omit(distinct(data.frame(loc =make.names(dat$Location), county = dat$County)))

nas <- data.frame(names= colnames(df),
                  val = as.numeric(colSums(is.na(df))/dim(df)[1])*100)
summary(nas$val[4:length(nas$val)])
hist(nas$val[4:length(nas$val)], breaks = 70)

## nas by date
nas_date <- data.frame(week = 1:(dim(df)),
                       val = as.numeric(rowSums(!is.na(df[,-c(1)]))))
nas_date %>% ggplot()+geom_point(aes(x = week, y = val))+
  geom_line(aes(x = week, y = val), color = "gray40")+ geom_vline(aes(xintercept = 19 ))+geom_vline(aes(xintercept = 282 ))
## removing locations that have more than 50% missing values
pctThrs <- 40
df <- df %>% subset(select = which(nas[,2] < pctThrs))
colnames(df) <- make.names(colnames(df))
## Remvoing st.Helens , Siletz and Siletz tribe ,Ontario, Ontario Prison, Redmond

df <- df %>% subset(select = -c(St..Helens,
                                Silverton,
                                Siletz
                                
))

nas$names <- make.names(nas$names)
## Find locations that have more than 50% missing values.
Nval <- dim(df)[2]-1#27
Nbrk <- 7
nas <- data.frame(names= colnames(df),
                  val = as.numeric(colSums(is.na(df))/dim(df)[1])*100)
nas <- nas[-c(1),]
df_samp <- df[, -c(1)]
ww_samp <- df[, nas$names[sort(nas$val, index.return=TRUE)$ix[1:Nval]]]
ww_samp <- as.matrix(ww_samp)
op_bp <- list()
col <-  data.frame(loc = colnames(ww_samp))
col <- left_join(col, df_county,by = join_by(loc == loc ))

bp_nat <- data.frame(matrix(ncol = 2, nrow = 0, 0))

## Getting average value measured over the state
avgVal <- rowMeans(ww_samp, na.rm = TRUE)

ww_ip <- data.frame(matrix(ncol = ncol(ww_samp), nrow =nrow(ww_samp), 0))

for(i in 1:dim(ww_samp)[2]){
  ### OLS CUSUM changepoint det with missing data
  
  ww_ip[,i] <- na_kalman(ww_samp[,i],model = "auto.arima" )
 
  bpIp <- data.frame(y = ww_ip[,i],
                     t = 1:(length(ww_ip[,i])),
                     y_lag1 = c(lag(ww_ip[,i])),#%>%
                     y_lag2 = c(lag(ww_ip[,i], 2)))# %>%
  
  op_efp <- breakpoints(y ~ t,h =dim(ww_samp)[2], data = bpIp )
  plot(op_efp)
  
  ## Finding breakpoints without any lag considerations
  op_bp[[i]] <- breakpoints(breakpoints(y ~ t, #y_lag1 + y_lag2, 
                                        h = dim(ww_samp)[2], 
                                        data = bpIp), 
                            breaks = Nbrk)
  #print(confint(breakpoints(y ~ 1, h =dim(ww_samp)[2], data = bpIp, breaks = Nbrk)))
  ## finding breakpoints with lag consideration
  bp_nat <- rbind(bp_nat, cbind(col[i,1],col[i,2], c(1,op_bp[[i]]$breakpoints, dim(ww_samp)[1])))
  
}


colnames(ww_ip) <- colnames(ww_samp)
colnames(bp_nat) <- c("Location", "County","Breakpoint")
bp_nat$Breakpoint <- as.numeric(bp_nat$Breakpoint)


## Adding information about the breakpoint numbering
bp_nat$brkpt <- rep(0:(Nbrk+1), times = Nval)

highest_rows <- bp_nat %>%
  filter(!(Breakpoint == 1) & !(Breakpoint ==dim(df)[1])) %>%
  group_by(brkpt) %>%
  slice_max(Breakpoint) 
# Select lowest row within each group
lowest_rows <- bp_nat %>%
  filter(!(Breakpoint == 1) & !(Breakpoint ==dim(df)[1])) %>%
  group_by(brkpt) %>%
  slice_min(Breakpoint) 
# To combine results, you can use `bind_rows`:
combined_results <- bind_rows(highest_rows, lowest_rows)

bp_nat$Location <- factor(bp_nat$Location, levels = bp_nat$Location[bp_nat$brkpt == 2])

bp_nat %>% ggplot() +
  geom_point(aes(y = interaction(Location, County), x = Breakpoint, group = Location)) +
  geom_line(aes(y = interaction(Location, County), x = Breakpoint, group = Location)) +
  geom_vline(data=highest_rows, aes(xintercept = Breakpoint), color = "darkblue", linetype = 2)+
  geom_vline(data=lowest_rows, aes(xintercept = Breakpoint), color = "darkred", linetype = 2)+
  scale_x_continuous(guide = guide_axis(angle = 45), limits = c(0,dim(ww_samp)[1])) + 
  ggtitle("Breakpoint distribution for locations") + 
  ylab("Locations") + xlab("Breakpoints") +
  theme_minimal()

bp_nat %>% ggplot() +
  geom_point(data =  bp_nat ,aes(y = interaction(Location), x = Breakpoint, group = Location))+
  geom_line(aes(y = interaction(Location), x = Breakpoint, group = Location), alpha =1) +
  scale_x_continuous(guide = guide_axis(angle = 45), limits = c(0,dim(ww_samp)[1])) + 
  ggtitle("Breakpoint distribution for locations") + labs(color = "Breakpoint number") +   guides(color = guide_legend(nrow = 1))+
  ylab("Locations") + xlab("Breakpoints") + theme(legend.position = "bottom")
theme_minimal()
#ggsave("brkpts.pdf")

ggsave(paste0(svfldr,"brkPts.png"))

colors <-  c("brown4","blue4","chartreuse4", "darkorchid","cornflowerblue", "coral2",
             "salmon4","cyan2","slateblue3","gold1",
             "darkolivegreen1", "aquamarine","firebrick2", "maroon2","dodgerblue" )
d <- 2
bp <- data.frame(matrix(nrow = 0, ncol =3))
for(i in 1:length(op_bp)){
  tmp <- c(1,op_bp[[i]]$breakpoints,dim(ww_ip)[1])
  bp <- rbind(bp, cbind(rep(i, length(tmp)),tmp, rep(col$loc[i], length(tmp))))
}

colnames(bp) <- c("TS","BkPtS", "Location")
bp$BkPtS <- as.numeric(bp$BkPtS)

bp$Date <- df$week_start[bp$BkPtS]

bp <- bp %>% 
  group_by(TS) %>% 
  mutate(BkPtE = as.numeric(lead(BkPtS, 1, default = NA))) %>%
  filter(!is.na(BkPtE)) %>%
  mutate(range = paste0(BkPtS,"-",BkPtE)) %>%
  ungroup()%>%
  mutate(xRangeL =round(rep(seq(from = min(rowMeans(ww_ip, na.rm = TRUE)), 
                                to = max(rowMeans(ww_ip, na.rm = TRUE)), 
                                length.out= length(unique(bp$TS))), 
                            each =Nbrk+1), 2) ) %>%
  mutate(xRangeH = rep(c(unique(xRangeL)[-1], 5.5) , each = Nbrk+1))%>%
  mutate(id = rep(1:(Nbrk+1), times = Nval)) %>%
  ungroup()

## find the minimum and maximum value for each changepoint window. 
ggplot()+
  geom_rect(aes(xmin = as.numeric(bp$BkPtS), xmax = as.numeric(bp$BkPtE), 
                ymin = bp$xRangeL, ymax = bp$xRangeH,
                fill = as.factor(bp$id)), alpha = 1) +
  geom_line(aes(y = rowMeans(ww_ip, na.rm = TRUE), x = as.numeric(1:dim(ww_ip)[1])))+
  geom_text(aes(x = -45, y =as.numeric(unique(bp$xRangeL)), label = col$loc), 
            size = 3, vjust = 0, hjust = 0, color = "blue4")+
  scale_fill_manual(values = colors)+
  xlab("Sample week")+ ylab("Average mean log copies.")+
  theme_minimal() + 
  theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Mean value of all cities")
ggsave(paste0(svfldr,"CPperCity.png"))
## Separate out all the time windows based on their start and end change points
## Imputing data and creating equal length time windows depending on stationarity of the data.

rng <- bp %>% group_by(id) %>%
  summarise(minV = min(BkPtS),
            maxV = max(BkPtE))%>%
  mutate(minDate  = df$Date[minV],
         maxDate  = df$Date[maxV])

TSsec <- list()

for (i in 1:Nval) {
  g <- list()
  for(j in 1:(Nbrk+1)){
    bp_ss <- bp[bp$TS ==i,]
    strt <- bp_ss[bp_ss$id==j,]$BkPtS
    end <- bp_ss[bp_ss$id==j,]$BkPtE
    t = ww_ip[strt:end,i]
    rngstrt <- rng[j,]$minV   
    rngend <- rng[j,]$maxV
    ## Front padding
    fp <- strt - rngstrt
    if(fp > 0){
      t <- c(rep(NA, fp), t) 
    }
    ## End padding for the data
    ep <- rngend - end
    if(ep >0){
      t <- c(t, rep(NA,ep))
    }
    # Example 5:  Perform imputation with KalmanSmooth and user created model. This imputes the missing data
    ## eventually need to replace this with multiple imputation
    g[[j]] <- na_kalman(t) 
  }
  TSsec[i] <- list(g)
}

df_split <- list()
for(i in 1:length(TSsec[[i]])){
  df_split[[i]] <- as.matrix(sapply(TSsec,"[[",i)) 
  colnames(df_split[[i]]) <- col$loc
}
StatDat <- list()

## Data frame for saving output from dickey fuller test
adfOp <- data.frame(matrix(ncol = 1+1+2+2, nrow = 0, 0))


## check if individual sections are stationary
for(i in 1:length(df_split)){
  print(paste0("Subsection,", i))
  dd <- df_split[[i]]
  cn <- colnames(dd)
  diffedDat <- matrix(ncol = ncol(dd), nrow = nrow(dd)-1,0)
  for(j in 1:dim(dd)[2]){
    diffedDat[,j] <- diff(dd[,j], lag =1)
    ## Running dickey fuller test on the original data 
    origTst <- adf.test(dd[,j])
    ## Running dickey fuller on differenced data
    diffTst <- adf.test(diffedDat[,j])
    
    adfOp <- rbind(adfOp, c(i,cn[j],origTst$statistic[[1]], origTst$p.value, diffTst$statistic[[1]], diffTst$p.value ))
  }
  colnames(diffedDat) <- cn
  StatDat[[i]] <- diffedDat
}
colnames(adfOp) <- c("Subsection","Timeseries", "origStat","origPval","diffStat","diffpVal")
## If not then take difference and check pacf to see lag.
## If the above is OK, run VAR using LASSO
op <- LassoVar(df_ip = StatDat,d = d,Nbrk = Nbrk, retFit = TRUE)

## Create maps using the LASSO penalty
G <- op[[1]]
nodeImp <- nodeImpfunc(G)

## find the max values for different node importance measures
maxImp <- nodeImp %>%
  group_by(Section) %>%
  summarise(maxOutdegree = max(Outdegree),
            maxIndegree = max(Indegree),
            maxBetweenness = max(Betweenness),
            maxBetCity = City[which(Betweenness == max(Betweenness))]
  )

## creating network based on lat lon of the locations
grphDat <- readRDS("GraphData.rds")
g <- grphDat$WWgrph
l <- data.frame(names = make.names(V(g)$name),
                lat = as.numeric(V(g)$lat),
                lon = as.numeric(V(g)$lon))
## adding missing cities
l <- data.frame(rbind(as.matrix(l),
                      rbind(c("Rock.Creek",-122.8808,45.5554),
                            c("Dallas",-123.3170,44.9193),
                            c("Sunriver",-121.4334,43.8694))))
nameVal <- data.frame(names = make.names(V(G[[2]])$names))
layout <- left_join(nameVal, l, join_by(names == names))
## Setting labels for the nodes that have highest outdegree and indegree with different color code

lvl <- 3
var <- c("Strength_out") # "Outdegree","Indegree","Closeness","Eigen", "Strength_in",
#pdf("SegmentPlots.pdf")

for(v in var){
  gbg <- nodeImp %>% 
    group_by(Section) %>%
    subset(select = c(Section, get(v), City))%>%
    mutate(min_rank = dense_rank(as.numeric(get(v))), 
           max_rank = dense_rank(-as.numeric(get(v)))) %>%
    mutate(label1 = ifelse( max_rank <= lvl, City, "")) %>%
    mutate(label2 = ifelse( min_rank <= lvl, City, "")) %>%#min_rank <= 1 |
    mutate(color1 = case_when(max_rank <= lvl ~ "green4")) %>%
    mutate(color2 = case_when(min_rank <= lvl ~ "red"))
  
  for(i in 1:(Nbrk+1)){
    V(G[[i]])$label1 = gbg$label1[gbg$Section == i]
    V(G[[i]])$color1 = gbg$color1[gbg$Section == i]
    
    g <- simplify(G[[i]])
    Edgwt <- E(g)$weight * 10
    par(mar=c(0.8,0,0.8,0))
    png(paste0(svfldr,v,i,".png"), width = 800, height = 800)
    plot.igraph(g,directed = T,
                edge.arrow.size = 0.3,
                edge.width = Edgwt,
                edge.color = "gray60",
                edge.alpha =0.5,
                layout = cbind(as.numeric(layout$lat),
                               as.numeric(layout$lon)),
                vertex.label = V(g)$label1,
                vertex.size = 5,
                vertex.color =V(G[[i]])$color1,
                vertex.label.dist=1,
                vertex.frame.color = NULL,
                vertex.label.font = 2
    )
    dev.off()
  }
}
