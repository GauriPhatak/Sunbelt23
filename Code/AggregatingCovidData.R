library(tidyverse)
library(dplyr)
require(ggplot2)

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


WD <- paste0(getwd(),"/Code/Combined_aggregated_OHA_normalized_2025-10-07.xlsx")
#"/Box/Preliminary Results Coronavirus Sewer Surveillance/ddPCR results/Covid/R output data/7.11.25/Combined_all_data_2025-07-11.xlsx")

dat <- readxl::read_excel(WD,sheet = "COVID", guess_max = 10000)

## removing KF S Suburban, Portland Tryon Creek, Ontario prison, Warm Spring SSL,
df <- dat %>% filter(!(Location %in%c("KF S Suburban", "Portland Tryon Creek", "Ontario-Prison", "Warm Springs-SSL", "Siletz Tribe") ))

## cities present in ct and wastewater data
ct_sbst <- ct$city_name[ct$city_name %in% unique(df$Location)]

unique_loc <- unique(df$Location)
unique_week <- table(df$week)
df <- df[!(is.na(df$Location)),]

gbg <- data.frame(table(df$Location, df$County))
gbg <- gbg[gbg$Freq > 0,]

## Removing gilliam county
df <- df[!(df$County == "Gilliam County"),]

## remove columns that we wont use
df <- df %>%
  dplyr::mutate(Date = lubridate::ymd(Sample_Date )) %>%
  dplyr::mutate(year = lubridate::year(Date), 
                month = lubridate::month(Date), 
                day = lubridate::day(Date))

df <- df %>%
  dplyr::mutate(Yr_Mth = format(as.Date(Date),"%Y-%m")) %>% 
  dplyr::mutate(weekno = as.numeric(Date - min(Date, na.rm = TRUE)),
                weekno = (weekno %/% 7) + 1) %>%
  group_by(weekno) %>% 
  arrange(desc(LogCopiesPerL)) %>%
  dplyr::mutate(
    maxByVal = if_else(row_number() == 1, "1", "0")
  ) %>%
  arrange(weekno)%>%
  ungroup()

df <- df %>% 
  dplyr::select(c(LogCopiesPerL ,Date,Location,MasterID,County,weekno,Yr_Mth,
                  maxByVal,year,month,day )) %>%
  distinct()


aggr <- df %>%
  group_by(weekno,County, Location) %>%
  dplyr::summarise(meanLogCopiesPerL = mean(LogCopiesPerL, na.rm = TRUE),
                   across(everything() ,first)) %>% 
  dplyr::select(-c(LogCopiesPerL))

saveRDS(aggr,paste0(getwd(), "/Code/MapNetworks/aggr_10_7.rds"))

