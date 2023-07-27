library(tidyverse)
library(grpreg)
library(HDGCvar)
library(epitools)
library(mlVAR)
library(lmtest)
library(igraph)
##Reading the combined all ddpcr data generated using the old code.
dat <- read_csv("C:/Users/gauph/Box/Preliminary Results Coronavirus Sewer Surveillance/ddPCR results/Covid/R output data/7.18.23/Combined_all_sample_data_aggregated_bysample2023-07-18.csv")
dat <- dat[dat$SampleType == "Influent",]
dat <- subset(dat, select = c(Location, County, LogCopiesPerL, Date, PCRplate))
dat$Date <- parse_date_time(dat$Date, orders = "mdy")
## Reading highway graph data
grphDat <- readRDS("Code/GraphData.rds")

## Case data and population data
df <- read.csv("Data/Graph_Full Data_data.csv")


## Assumptions of VAR models:
#The error term has a conditional mean of zero.
#The variables in the model are stationary.
#Large outliers are unlikely.
#No perfect multicollinearity.

#In general, for a VAR(p) model, 
#the first p lags of each variable in the system would be used as regression predictors for each variable

## Check Assumpotion of stationatiry, non stationariy or unit root

## Using the HDGCvar model to do pairwise granger causality.
# p is the upper bound for lag.


## Need to pivot data to get data for each Location at the given date
dat$epiweek <- paste(epiweek(dat$Date), as.character(year(dat$Date)), sep= "-")
ts <- dat %>% 
  subset(select=c(epiweek, Location, LogCopiesPerL)) %>% 
  group_by(Location, epiweek) %>%
  dplyr::summarise(meanLC = mean(LogCopiesPerL, na.rm = TRUE)) %>%
  ungroup() %>%
  pivot_wider(values_from = meanLC, names_from = Location, values_fill = NA) %>%
  separate_wider_delim(epiweek, delim = "-", names = c("week","year"))


## Subsetting just the 2021 data

ts$week <- as.numeric(ts$week)

ts <- ts[ts$year == "2021" & ts$week <= 46 & ts$week >= 26,]

ts <- ts[, colSums(!is.na(ts)) > 8]

ts <- ts %>% subset(select = -c(week,year))

op <- cbind(expand_grid(x = 1:ncol(ts), y=1:ncol(ts)), 
            expand_grid(x = colnames(ts) , y = colnames(ts)))
colnames(op) <- c("a","b","c1","c2")

op <- op[op[,1] != op[,2],]
op$pValO1 <- unlist(apply(op, 1, function(x) grangertest(ts[,as.numeric(x[1])], 
                                                         ts[,as.numeric(x[2])], 
                                                         order = 1)$Pr[2]))

op$SigO1 <- op$pValO1 <= 0.05

## subsetting only the significant nodes

opSig <- op[op$SigO1,]

# Creating a graph using the output from granger test. 
g <- grphDat$WWgrph
l <- cbind(V(g)$name ,as.numeric(V(g)$lat), as.numeric(V(g)$lon))
G <- graph_from_edgelist(as.matrix(opSig[, c(3,4)]), directed = TRUE)
G <- simplify(G)
lsig <- cbind(l[match(V(G)$name, l[,1]),],V(G)$name)
lsig <- lsig[!is.na(lsig[,1]),]
G2 <- delete.vertices(G, c("Warm Springs", "Dallas","Rock Creek"))

plot.igraph(G2, layout= cbind(as.numeric(lsig[,2]), as.numeric(lsig[,3])), 
            vertex.colo= "darkblue", vertex.frame= NA, edge.color= "cornsilk4",
            vertex.label= NA, edge.width= 1, edge.arrow.size= 0.01, 
            vertex.size= 4)








