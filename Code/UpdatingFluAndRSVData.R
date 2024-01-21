library(tidyverse)
#library(xlsx)
library(openxlsx)
boxPath <- "C:/Users/gauph/Box/"
WD <- paste0(boxPath,"Preliminary Results Coronavirus Sewer Surveillance/Data files/rsv data files/processed rsv data files (R input)") 
allFiles <- list.files(WD)
#file <- allFiles[1]

for(file in allFiles){
  
  sheet1 <- read_excel(paste(WD,"/",file,sep = ""),sheet=1)
  dat <- read_excel(paste(WD,"/",file,sep = ""),sheet=2)
  dat$Target[dat$Target == "FLUA"] <- "InfA"
  dat$Target[dat$Target == "FLUB"] <- "InfB"
  #dat$Target[dat$Target == "RPP30"] <- "RP"
  
  xl_lst <- list("from QSAP" = sheet1, "for R" = dat)
  write.xlsx(xl_lst, 
             file = paste0(boxPath,"Preliminary Results Coronavirus Sewer Surveillance/Data files/rsv data files/RSV_Data_RInput/",file))
  
}
