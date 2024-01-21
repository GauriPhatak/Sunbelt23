## This script is to try out code using random dump
library(igraph)
library(geosphere)
library(ape)

##Loading city data
#hwpath = "C:/Users/gauph/Documents/StatisticsMS_PhD/Wastewater-Surveillance-OSU/RScripts/OSDLDataExtract/data.gdb"
#rgdal::ogrListLayers(hwpath)
#ct <- st_read(hwpath, layer = "City_Limits__2016_")


## loading the census block data to see if we can discrene granularity in it
cncsPth = "C:/Users/gauph/Documents/StatisticsMS_PhD/Wastewater-Surveillance-OSU/Sunbelt23/Data/Reference_Census_Block_Groups__2010_/data.gdb"
rgdal::ogrListLayers(cncsPth)
cncs <- st_read(cncsPth, layer = "Census_Block_Groups__2010_")


data(afcon, package="spData")
oid <- order(afcon$id)
resI <- localmoran(afcon$totcon, nb2listw(paper.nb))
printCoefmat(data.frame(resI[oid,], row.names=afcon$name[oid]),
             check.names=FALSE)
hist(resI[,5])
mean(resI[,1])
sum(resI[,1])/Szero(nb2listw(paper.nb))
moran.test(afcon$totcon, nb2listw(paper.nb))
# note equality for mean() only when the sum of weights equals
# the number of observations (thanks to Juergen Symanzik)
resI <- localmoran(afcon$totcon, nb2listw(paper.nb))
printCoefmat(data.frame(resI[oid,], row.names=afcon$name[oid]),
             check.names=FALSE)
hist(p.adjust(resI[,5], method="bonferroni"))
totcon <-afcon$totcon
is.na(totcon) <- sample(1:length(totcon), 5)
totcon
resI.na <- localmoran(totcon, nb2listw(paper.nb), na.action=na.exclude,
                      zero.policy=TRUE)
if (class(attr(resI.na, "na.action")) == "exclude") {
  print(data.frame(resI.na[oid,], row.names=afcon$name[oid]), digits=2)
} else print(resI.na, digits=2)
resG <- localG(afcon$totcon, nb2listw(include.self(paper.nb)))
print(data.frame(resG[oid], row.names=afcon$name[oid]), digits=2)
set.seed(1)
resI_p <- localmoran_perm(afcon$totcon, nb2listw(paper.nb))
printCoefmat(data.frame(resI_p[oid,], row.names=afcon$name[oid]),
             check.names=FALSE)