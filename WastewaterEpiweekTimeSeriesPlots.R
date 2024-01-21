library(tidyverse)
library(epitools)
library(readxl)
library(fs)
library(ggiraph)
library(gridExtra)

dat <- read_excel(paste0(path_home(),"/Box/Preliminary Results Coronavirus Sewer Surveillance/ddPCR results/Covid/R output data/11.7.23 (1)/Combined_all_data_2023-11-08.xlsx"), sheet = "COVID",guess_max = 10000)
rank <- read_csv(paste0(path_home(),"/Box/Preliminary Results Coronavirus Sewer Surveillance/Data files/covid data files/RScripts/RequiredFiles/HealthRegions.csv"))
df <- dat  %>% 
  filter(Study == "OHA") %>% 
  subset(select = c(Date, CalcCopies, Location, County, Target, Site))

epiweek <- epitools::as.week(df$Date, format  = "%m/%d/%y")
df$epiweek <- as.numeric(epiweek$week)
df$epimonth <- month(as.Date(df$Date, format = "%m/%d/%y"))
df$year <- year(as.Date(df$Date, format = "%m/%d/%y"))

df_HR <- merge(df, rank, by = "Location", all.x = TRUE)

df_HR <- df_HR %>%
  group_by(Health.Region, year, epiweek) %>%
  summarise(meanlogcopies =  mean(logCopies, na.rm = TRUE)) %>%
  ungroup()
df_HR$size <- ifelse(df_HR$year == 2023, 0.5,0.1)
pal <- c("#F0E442", "#E69F00", "#D55E00", "#000000",  "#56B4E9", "#009E73" , "#0072B2", "#CC79A7")

############## Creating plots by health region.

for(var in c(1,2,3,5,6,7,9)){
  
  print( df_HR  %>% filter(complete.cases(.)) %>%
           filter(Health.Region == var) %>%
           ggplot()+ #df[df$Health.Region == var & !is.na(df$epiweek),])+
           geom_point(aes(y = meanlogcopies, 
                          x = factor(epiweek),
                          color = factor(year), size = factor(size), group = 1)) +
           scale_size_discrete(range = c(1, 2))+
           guides(size = FALSE, color = guide_legend(title = "Year"))+
           #scale_size(range = c(0,0.5))+
           geom_line(aes(y = meanlogcopies, 
                         x = epiweek,
                         color = factor(year))) +
           scale_colour_manual(values = pal)+
           ggtitle(paste("Comparison of SARS-CoV-2 trend by Epiweek between years for Health Region", var, sep = " "))+ 
           xlab("Epidemiological Week")+
           ylab(paste("Average log copies/L for HR",var, sep = " "))+
           theme_minimal() + theme( axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                                    panel.grid.major=element_line(colour="gray92"),
                                    panel.grid.minor=element_line(colour="gray92")))
}


### creating approx 40 graphs by location.

df_loc <- merge(df, rank, by = "Location", all.x = TRUE)

df_loc <- df_loc %>%
  subset(!(Location %in% c("Broadman","Condon","Gold Beach"))) %>%
  group_by(Location, year, epiweek) %>%
  summarise(meancopies =  mean(CalcCopies, na.rm = TRUE)) %>%
  ungroup()
#df_HR$size <- ifelse(df_HR$year == 2023, 0.5,0.1)
pal <- c( "#009E73" , "#0072B2", "#CC79A7", "#F0E442",  "#000000", "#E69F00", "#D55E00",  "#56B4E9", "#009E73" , "#0072B2", "#CC79A7")

############## Creating plots by health region.
plot.list <- list()
i <- 1
for(var in unique(df_loc$Location)){
  print(var)
  plot.list[[i]] <- df_loc  %>% filter(complete.cases(.)) %>%
           filter(Location == var) %>%
           ggplot()+ #df[df$Health.Region == var & !is.na(df$epiweek),])+
           geom_point(aes(y = meancopies, 
                          x = factor(epiweek),
                          color = factor(year),  group = 1)) +
           guides(color = guide_legend(title = "Year"))+
           #scale_size(range = c(0,0.5))+
           geom_line(aes(y = meancopies, 
                         x = factor(epiweek),
                         color = factor(year), group = year)) +
           scale_colour_manual(values = pal)+
           ggtitle(paste("Comparison of SARS-CoV-2 trend by Epiweek between years for ", var, sep = " "))+ 
           xlab("Epidemiological Week")+
           ylab(paste("Average log copies/L for ",var, sep = " "))+
           theme_minimal() + theme( axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                                    panel.grid.major=element_line(colour="gray92"),
                                    panel.grid.minor=element_line(colour="gray92"))
  
  i <- i+1
}
#for(p in plot.list){
#  print(p)
#}
#plots <- do.call(marrangeGrob, c(plot.list, list(ncol = 1, nrow = 2)))

# To save to file, here on A4 paper
#ggsave(filename = "multipage_plot.pdf", plot = plots, width = 21, height = 29.7, units = "cm")
pdf("PerLocation.pdf", height = 7, width = 8)
for(p in plot.list){
  plot(p)
}
dev.off()


## plots by county


df_county <- merge(df, rank, by = "Location", all.x = TRUE) 
df_county <- df_county %>%
  group_by(Health.Region, County, epiweek, year) %>%
  summarise(meanlogcopies =  mean(logCopies, na.rm = TRUE)) %>%
  ungroup()

pal <- c("#F0E442", "#E69F00", "#D55E00", "#000000",  "#56B4E9", "#009E73" , "#0072B2", "#CC79A7")
for(var in c(1,2,3,5,6,7,9)){
  print(df_county %>% 
          filter(complete.cases(.)) %>%
          filter(year == 2023, Health.Region == var) %>%
          ggplot()+
          geom_point(aes(y = meanlogcopies, 
                         x = factor(epiweek),
                         color = factor(County))) +
          guides(size = FALSE, color = guide_legend(title = "County"))+
          #scale_size(range = c(0,0.5))+
          geom_line(aes(y = meanlogcopies, 
                        x = factor(epiweek),
                        color = factor(County), group = County)) +
          scale_colour_manual(values = pal)+
          ggtitle(paste("Comparison of SARS-CoV-2 trend by Epiweek between counties for Health Region", var, sep = " "))+ 
          xlab("Epidemiological Week")+
          ylab(paste("Average log copies/L for HR",var, sep = " "))+
          theme_minimal() + theme( axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), panel.grid.major=element_line(colour="gray92"),panel.grid.minor=element_line(colour="gray92")))
}

############# Create a plpots by location

df_loc <- merge(df, rank, by = "Location", all.x = TRUE) 
df_loc <- df_loc %>%
  group_by(County,Location, epiweek, year) %>%
  summarise(meanlogcopies =  mean(logCopies, na.rm = TRUE)) %>%
  ungroup()

pal <- c("#F0E442", "#E69F00", "#D55E00", "#000000",  "#56B4E9", "#009E73" , "#0072B2", "#CC79A7")
for(c in unique(df$County)){
  print(df_loc %>% 
          filter(complete.cases(.)) %>%
          filter(year == 2023, County == c) %>%
          ggplot()+
          geom_point(aes(y = meanlogcopies, 
                         x = factor(epiweek),
                         color = factor(Location))) +
          guides(size = FALSE, color = guide_legend(title = "Location"))+
          #scale_size(range = c(0,0.5))+
          geom_line(aes(y = meanlogcopies, 
                        x = factor(epiweek),
                        color = factor(Location), group = Location)) +
          scale_colour_manual(values = pal)+
          ggtitle(paste("Comparison of SARS-CoV-2 trend by Epiweek between locations in", c,"County", sep = " "))+ 
          xlab("Epidemiological Week")+
          ylab(paste("Average log copies/L for",c,"County", sep = " "))+
          theme_minimal() + theme( axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), panel.grid.major=element_line(colour="gray92"),panel.grid.minor=element_line(colour="gray92")))
}



########################### Sample plots to show the team ###############################

var = 9
df_HR  %>% filter(complete.cases(.)) %>%
  filter(Health.Region == var) %>%
  ggplot()+ #df[df$Health.Region == var & !is.na(df$epiweek),])+
  geom_point(aes(y = meanlogcopies, 
                 x = factor(epiweek),
                 group = 1)) +
  scale_size_discrete(range = c(1, 2))+
  guides(size = FALSE, color = guide_legend(title = "Year"))+
  #scale_size(range = c(0,0.5))+
  geom_line(aes(y = meanlogcopies, 
                x = epiweek)) +
  scale_colour_manual(values = pal)+
  ggtitle(paste("Comparison of SARS-CoV-2 trend by Epiweek between years for Health Region", var, sep = " "))+ 
  xlab("Epidemiological Week")+
  ylab(paste("Average log copies/L for HR",var, sep = " "))+
  theme_minimal() + 
  theme( axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
         panel.grid.major=element_line(colour="gray92"),
         panel.grid.minor=element_line(colour="gray92"))+
  facet_wrap(vars(year), ncol = 1)

pal <- c("#b2df8a","#a6cee3", "#33a02c",  "#1f78b4")

if( requireNamespace("dplyr", quietly = TRUE)){
  gg <- df_HR  %>% filter(complete.cases(.)) %>%
    filter(Health.Region == var) %>%
    ggplot(aes(hover_css = "fill:none;")) +
    geom_point_interactive(aes(y = meanlogcopies, 
                   x = factor(epiweek),
                   tooltip = factor(year),
                   data_id = factor(year),
                   color = factor(year),
                   group = 1)) +
    geom_line_interactive(aes(y = meanlogcopies, 
                              x = epiweek,
                              tooltip = factor(year),
                              data_id = factor(year),
                              color = factor(year)), size = .75) +
    scale_colour_manual(values = pal) +
    guides(color = guide_legend(title = "Year")) +
    ggtitle(paste("Comparison of SARS-CoV-2 trend by Epiweek between years for Health Region", var, sep = " "))+ 
    xlab("Epidemiological Week")+
    ylab(paste("Average log copies/L for HR",var, sep = " "))+
    theme_minimal() + theme( axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                             panel.grid.major=element_line(colour="gray"),
                             panel.grid.minor=element_line(colour="gray"))
  x <- girafe(ggobj = gg, width_svg = 10, height_svg = 7)
  x <- girafe_options(x = x, opts_hover(css = "stroke:darkred;fill:darkred"), opts_hover_inv(css = "stroke:grey"))
  if( interactive() ) print(x)
}

htmlwidgets::saveWidget(file = paste0(getwd(),"/WWplotYearsHR.html"), x)

pal <- c("aquamarine", "aquamarine2","aquamarine3","aquamarine4")
if( requireNamespace("dplyr", quietly = TRUE)){
  gg <- df_HR  %>% filter(complete.cases(.)) %>%
    filter(Health.Region == var) %>%
    ggplot(aes(hover_css = "fill:none;")) +
    geom_point_interactive(aes(y = meanlogcopies, 
                               x = factor(epiweek),
                               tooltip = factor(year),
                               data_id = factor(year),
                               color = factor(year),
                               group = 1)) +
    geom_line_interactive(aes(y = meanlogcopies, 
                              x = epiweek,
                              tooltip = factor(year),
                              data_id = factor(year),
                              color = factor(year))) +
    scale_colour_manual(values = pal) +
    guides(color = guide_legend(title = "Year")) +
    ggtitle(paste("Comparison of SARS-CoV-2 trend by Epiweek between years for Health Region", var, sep = " "))+ 
    xlab("Epidemiological Week")+
    ylab(paste("Average log copies/L for HR",var, sep = " "))+
    theme_minimal() + theme( axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                             panel.grid.major=element_line(colour="gray"),
                             panel.grid.minor=element_line(colour="gray"))
  x <- girafe(ggobj = gg, width_svg = 10, height_svg = 7)
  x <- girafe_options(x = x, opts_hover(css = "stroke:darkred;fill:darkred"), opts_hover_inv(css = "stroke:grey"))
  if( interactive() ) print(x)
}

htmlwidgets::saveWidget(file = paste0(getwd(),"/WWplotYearsHR.html"), x)


plot_countiesCOVID <- function(data){
  #6. Influents are subset out from fulloha
  samples.aggr <- SampleAggr(data)
  fulloha <- samples.aggr[samples.aggr$Study == "OHA",]
  inf <- subset(fulloha, SampleType == "Influent")
  
  #tail(inf)
  oha.county <- dcast(inf, County ~ ., fun.aggregate = mean, value.var = c("LogCopiesPerL"), na.rm = TRUE)
  #head (oha.county)
  
  counties <- sf::st_as_sf(maps::map("county", plot = FALSE, fill = TRUE))
  counties <- subset(counties, grepl("oregon,", counties$ID))
  counties$ohaAvg <- NA
  for (i in 1:nrow(oha.county)) {
    m <- agrep(oha.county$County[i], counties$ID, ignore.case = T)
    counties$ohaAvg[m] <- oha.county[i, 2]
  }
  #head(counties)
  
  
  #7. Plot for Concentration by County is generated
  p <- ggplot() +  
    geom_sf(data = counties, aes(fill = ohaAvg), color = gray(0.9)) + 
    scale_fill_viridis_c(name = "Gene Copies Per L", labels = unlog) + 
    ggtitle("Geometric Mean WW SARS-CoV-2 Concentrations by County") + 
    xlab("Longitude") + 
    ylab("Latitude")
  return(p)
  
}

plot_countiesCOVID(dat)
#getwd()
#path_home()

#library(webshot)
#webshot::install_phantomjs()
#webshot("paste_your_html_here" , "output.png", delay = 0.2 , cliprect = c(440, 0, 1000, 10))
# 
# dat <- data.frame(cbind(num = 1:40 ,a = rnorm(40), b= rnorm(40, 2,1), c= rnorm(40, 1, 1), d = rnorm(40, 1,2)))
# dat <- pivot_longer(data = dat, !num)
# if( requireNamespace("dplyr", quietly = TRUE)){
#   gg <- dat %>%
#     ggplot(aes(hover_css = "fill:none;")) +
#     geom_line_interactive(aes(y = value, 
#                               x = num,
#                               tooltip = name,
#                               data_id = name,
#                               color = name), size = .75) +
#     scale_colour_manual(values = pal) +
#     guides(color = guide_legend(title = "Year")) +
#     ggtitle("Dummy plot")+ 
#     xlab("sample number")+
#     ylab("Values")+
#     theme_minimal() + theme( axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
#                              panel.grid.major=element_line(colour="gray"),
#                              panel.grid.minor=element_line(colour="gray"))
#   x <- girafe(ggobj = gg, width_svg = 10, height_svg = 7)
#   x <- girafe_options(x = x, opts_hover(css = "stroke:darkred;fill:darkred"))
#   if( interactive() ) print(x)
# }
# htmlwidgets::saveWidget(file = paste0(getwd(),"/dummy.html"), x)
