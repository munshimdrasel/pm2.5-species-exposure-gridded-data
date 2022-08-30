#PART1: Getting data
library(downloader)

library(tidyverse)
library(sf)
library(scales)
library(fst)
library(plyr)
library(RCurl)
library(stringi)

setwd ("/Volumes/GoogleDrive/My Drive/R/pm2.5-species-exposure-gridded-data")

# setwd ("/Volumes/GoogleDrive/My Drive/R/pm2.5-species-exposure-gridded-data")

#downloading gridded PM2.5 species data from 
# https://www.ciesin.columbia.edu/data/aqdh/pm25component-EC-NH4-NO3-OC-SO4-2000-2019/data/

species <- c ("ec", "nh4", "no3", "oc", "so4")
year <- 2001:2002
location <- c ("urban", "non-urban")

url1 <- "https://www.ciesin.columbia.edu/data/aqdh/pm25component-EC-NH4-NO3-OC-SO4-2000-2019/data/aqdh-pm25component-"
url2 <- "-"
url3 <- "-"
url4 <- "-rds.zip"
res <- as.vector(do.call(paste0,expand.grid(url1, species, url2, year, url3, location, url4)))

  
  
#non-urban data
# url.non.urban = paste0("https://www.ciesin.columbia.edu/data/aqdh/pm25component-EC-NH4-NO3-OC-SO4-2000-2019/data/aqdh-pm25component-ec-",
#              year, "-non-urban-rds.zip")

for (i in 1:length(res)) {
  
  download.file(res[i], dest= paste0("./data/", sub(".*/data/", "", res[i])), mode="wb") 
}
download(url.non.urban, dest= paste0("./data/", "aqdh-pm25component-ec-",year, "-non-urban-rds.zip"), mode="wb") 


# getting all the zip files
zipF <- list.files(path = "./data",   pattern = "*.zip", full.names = TRUE)

# unzip all your files

#non-urban
ldply(.data = zipF, .fun = unzip, exdir = "./data/non-urban")




# get the rds files
rds_files <- list.files(path = "./data/non-urban", pattern = "*.rds",  full.names = TRUE)

# read the csv files
my_data <- ldply(.data = rds_files, .fun = readRDS)

write.fst(my_data, "/Volumes/GoogleDrive/My Drive/R/get_data_ampd/emissions-raw.fst")