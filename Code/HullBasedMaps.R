library(tidyverse)
library(dplyr)
library(sf)
library(mapview)
library(sp)
library(usmap)
require(rgdal)
require(ggplot2)
library(raster)
library(igraph)
library(GGally)
library(ggnet)
library(geosphere)
library(png)

source(paste0(getwd(),"/Code/HelperFuncs.R"))

aggr <- readRDS(paste0(getwd(),"/Code/MapNetworks/aggr_10_7.rds"))

##Loading city limits data 
hwpath = "C:/Users/gauph/Documents/StatisticsMS_PhD/Wastewater-Surveillance-OSU/RScripts/OSDLDataExtract/data.gdb"
ct <- st_read(hwpath, layer = "City_Limits__2016_")
ct <- st_transform(ct, 2994)

sunriver <- st_read("C:/Users/gauph/Documents/StatisticsMS_PhD/Wastewater-Surveillance-OSU/Sunbelt23/Data/SunriverShapefile/Sunriver.shp")
sunriver<-  st_transform(sunriver, crs = 2994)
sunriver <- rename(sunriver, city_name = id)
sunriver$city_name[1] <- "Sunriver"
sunriver <- rename(sunriver, shape = geometry)
sunriver <- st_cast(sunriver, "MULTIPOLYGON")

## adding rock creek shapefile
Rockcreek <- st_read("C:/Users/gauph/Documents/StatisticsMS_PhD/Wastewater-Surveillance-OSU/Sunbelt23/Data/RockcreekShapefile/Rockcreek.shp")
Rockcreek <-  st_transform(Rockcreek, crs = 2994)
Rockcreek <- rename(Rockcreek, city_name = id)
Rockcreek$city_name[1] <- "Rock Creek"
Rockcreek <- rename(Rockcreek, shape = geometry)
Rockcreek <- st_cast(Rockcreek, "MULTIPOLYGON")

## Adding warmsprings shapefile
WarmSprings <- st_read("C:/Users/gauph/Documents/StatisticsMS_PhD/Wastewater-Surveillance-OSU/Sunbelt23/Data/WarmSpringsShapefile/WarmSprings.shp")
WarmSprings <-  st_transform(WarmSprings, crs = 2994)
WarmSprings <- rename(WarmSprings, city_name = id)
WarmSprings$city_name[1] <- "Warm Springs"
WarmSprings <- rename(WarmSprings, shape = geometry)
WarmSprings <- st_cast(WarmSprings, "MULTIPOLYGON")

ct <- bind_rows(ct, sunriver,Rockcreek, WarmSprings)
city_cntr <- ct %>% 
  st_as_sf() %>% 
  st_centroid() %>% 
  st_cast("POINT") %>% 
  st_transform(4326) %>% 
  st_coordinates()

## cities present in ct and wastewater data
ct_sbst <- ct$city_name[ct$city_name %in% unique(aggr$Location)]

## Loading highway data
hwpath = "C:/Users/gauph/Documents/StatisticsMS_PhD/Wastewater-Surveillance-OSU/RScripts/HighwayDat/data.gdb"
hwd <- st_read(hwpath, layer="Highway_Network__2015_")
hwd <- st_transform(hwd, 2994)
hwd <- st_zm(hwd)
id <- unique(hwd$st_hwy_idx)

ct_hull <- st_convex_hull(ct)
ct_hull_buff <- list()
dist <- c(10,15,20,25,30,35,40,45,50)
for(i in 1:length(dist)){
  ct_hull_buff[[i]] <- st_buffer(ct_hull, dist  = dist[i]*1000)
}

#Using the find edges function to find the edges between the different convex hulls
edges_hull<- list()
for(i in 1:length(dist)){
  edges_hull[[i]] <- findEdges(ct_hull_buff[[i]], hwd, id, as.character(dist[i]))
}

for (i in 1:length(dist)) {
  edgesfull <- bind_rows(edges_hull[1:i])
  hwy_grph <- graph_from_data_frame(d = edgesfull ,vertices = ct$city_name, directed = FALSE)
  
  mults <- simplify_and_colorize(hwy_grph)
  pdf(file= paste0(getwd(),"/Code/MapNetworks/HighWayGraphHull",dist[i],"km_WS_SR_RC.pdf"))
  par(mar=c(0,0.2,0,0.2)+.55)
  plot.igraph(mults, vertex.label = NA, vertex.size = 1, 
              edge.width = 1, layout = city_cntr)
  dev.off()
  
  incl <- V(hwy_grph)$name %in% ct_sbst
  V(hwy_grph)$included <- ifelse(incl == TRUE, "orangered3", "grey37")
  E(hwy_grph)$weight <- edgesfull$distance
  
  ## keeping all the edges from above. Checking if this would add more edges to the data
  saveDirEdges <-  data.frame(matrix(ncol = 2, nrow = 0))
  colnames(saveDirEdges) <- c("to","from")
  
  
  for (k in 1:length(ct_sbst)) {
    
    trav <- list()
    path <- shortest_paths(hwy_grph, ct_sbst[k], 
                           ct_sbst, output = "both" , 
                           weights = round(as.numeric(E(hwy_grph)$distance)))
    
    for(j in 1:length(ct_sbst)){
      curPath <- names(unlist(path$vpath[[j]]))
      subPath <- curPath[curPath %in% ct_sbst]
      if(length(subPath) > 1){
        saveDirEdges <- rbind(saveDirEdges, cbind(subPath[1:length(subPath)-1], subPath[2:length(subPath)]))
      }
    }
  }
  ## Creating edges from the paths that we generate using shortest path algorithm. 
  ## Coloring based on the weight of the edge i.e. the number of nodes between the two nodes of interest.
  
  DirectGrph <- graph_from_data_frame(d = saveDirEdges, vertices = ct$city_name, directed = FALSE)
  V(DirectGrph)$included <- ifelse(incl == TRUE, "orangered3", "grey37")
  DirectGrph <- simplify(DirectGrph)
  ## Save the graph for use
  saveRDS(DirectGrph, paste0(getwd(), "/Code/MapNetworks/ShortestPathGrph",dist[i],"Km_WS_SR_RC.rds"))
  
  G_FullGrph <- delete_vertices(DirectGrph, V(DirectGrph)[included == "grey37"])
  O_attr <- as.data.frame(cbind(as.data.frame(vertex_attr(DirectGrph)), city_cntr))
  fg_attr <- O_attr[O_attr$included == "orangered3",]
  fg_attr <- fg_attr[fg_attr$name %in% unique(aggr$Location), ]
  G_FullGrph <- induced_subgraph(G_FullGrph,vids = fg_attr$name)
  
  
  G_FullGrph <- set_vertex_attr(G_FullGrph, name = "X", value = as.matrix(fg_attr$X ))
  G_FullGrph <- set_vertex_attr(G_FullGrph, name = "Y", value = as.matrix(fg_attr$Y ))
  saveRDS(G_FullGrph, paste0(getwd(),"/Code/MapNetworks/G_HullGrph",dist[i],"Km_WS_SR_RC.rds"))
  
  pdf(file= paste0(getwd(),"/Code/MapNetworks/G_HullGraph",dist[i],"km_WS_SR_RC.pdf"))
  par(mar=c(0,0.2,0,0.2)+.55)
  plot.igraph(G_FullGrph, vertex.label = NA, vertex.size = 1,
              vertex.color = V(G_FullGrph)$included, 
              vertex.frame.color = V(G_FullGrph)$included, 
              edge.color = "darkgrey", edge.width = 1, 
              layout = city_cntr[O_attr$name %in% fg_attr$name, ])
  dev.off()

}

saveRDS(O_attr,paste0(getwd(),"/Code/MapNetworks/O_attr_WS_SR_RC.rds"))

