rm(list = ls())

#PART1: Getting data
library(downloader)
library(tidyverse)
library(fst)
library(stringi)
library(plyr)
library( sf)
library( raster)
library( data.table)
library(tidycensus, quietly = TRUE)
library(tigris, quietly = TRUE)
library( fasterize)
library( USAboundaries)
library( magrittr)
library( ncdf4)
options(tigris_use_cache = TRUE)
options(tigris_class = "sf")

# setwd ("/Volumes/GoogleDrive/My Drive/R/pm2.5-species-exposure-gridded-data")

setwd ("/projects/HAQ_LAB/mrasel/R/pm2.5-species-exposure-gridded-data")

#downloading gridded PM2.5 species data from 
# https://www.ciesin.columbia.edu/data/aqdh/pm25component-EC-NH4-NO3-OC-SO4-2000-2019/data/

# selected.species <- c ("so4","ec", "nh4", "no3", "oc")
selected.species <- "oc"
years <- 2000:2019
us_states <- USAboundaries::us_states()  %>%  filter(!stusps %in% c("PR", "HI", "AK") )
states.vec <- as.vector(us_states$stusps)
states.name <- as.vector(us_states$name)
location <- c ("urban", "non-urban") #change this

# url1 <- "https://www.ciesin.columbia.edu/data/aqdh/pm25component-EC-NH4-NO3-OC-SO4-2000-2019/data/aqdh-pm25component-"
# url2 <- "-"
# url3 <- "-"
# url4 <- "-rds.zip"
# res <- as.vector(do.call(paste0,expand.grid(url1, species, url2, year, url3, location, url4)))

options(timeout=600)


#download data
# for (i in 1:length(res)) {
# 
#   download.file(res[i], dest= paste0("./data/", sub(".*/data/", "", res[i])), mode="wb")
# }


# getting all the zip files
# zipF <- list.files(path = "./data",   pattern = "*so4*", full.names = TRUE)

# unzip all your files

# ldply(.data = zipF, .fun = unzip, exdir = paste0("./data/", location, "/", selected.species, "/"))


# get the rds files
rds_files_urban <- as.vector(list.files(path = paste0("./data/", location [1], "/", selected.species), 
                                        pattern = "*.rds",  full.names = TRUE))
rds_files_non_urban <- as.vector(list.files(path = paste0("./data/", location [2], "/", selected.species), 
                                            pattern = "*",  full.names = TRUE))

# tx.block.groups <- block_groups("Texas")



# we want to use an equal area projection, here's one I often use:
p4s <- "+proj=lcc +lat_1=33 +lat_2=45 +lat_0=40 +lon_0=-97 +a=6370000 +b=6370000"

ras_dom<-raster(xmn=-126, xmx=-66, ymn=23, ymx=51,
                resolution=c(0.1,0.1))

all.list <- list()
df.list <- list()

# for (i in 1:length(states.vec)) {
#   for (i in 1:length(selected.species)){
#     for (i in 1:length(years)) {
for (i in 1:length(rds_files_non_urban)) {
  
  dt.non.urban <- ldply(.data = rds_files_non_urban[i], .fun = readRDS)
  dt.urban <- ldply(.data = rds_files_urban[i], .fun = readRDS)
  dt <- rbind(dt.non.urban, dt.urban)
  dt <- dt %>% dplyr::select(lon,lat,final.predicted.oc)
  dt.raster<- rasterize(x=dt[, 1:2], y=ras_dom,field=dt[, 3], fun = mean)
  df.sf <- rasterToPolygons( dt.raster ) %>%st_as_sf()
  
  df.sf.transform <- df.sf %>% st_transform(crs = p4s)
  
  for (j in 1:length(states.vec)) {
    block.groups <- block_groups(states.vec[j])
    block.groups.transform <- block.groups %>% st_transform(crs=p4s)
    bg.selected <- block.groups.transform %>% dplyr::select(STATEFP, COUNTYFP, TRACTCE, BLKGRPCE, geometry)
    
    conc.bg <- st_interpolate_aw( df.sf.transform["layer"], block.groups.transform,   extensive = F)
    vec <- as.vector(intersect(bg.selected$geometry, conc.bg$geometry))
    bg.new <- bg.selected %>%  filter(geometry %in% vec)
    conc.bg$species <- sub("^([^-]+-){3}([^-]+).*", "\\2", rds_files_non_urban[i]) #FOR urban change 3 to 2
    conc.bg$year <- as.numeric(sub("^([^-]+-){4}([^-]+).*", "\\2", rds_files_non_urban[i]))
    conc.bg$state <- states.vec[j]
    conc<- cbind (bg.new,conc.bg) %>% dplyr::select(-geometry.1)
    # ggplot( conc) + geom_sf( aes( fill = layer,  geometry = geometry),
    #                        color = NA) +scale_fill_viridis_c( option="plasma")
    # 
    # cbind(bg.selected, conc.bg)
    df.list[[j]] <- conc
    
  }  
  
  save(df.list, file=paste0("./data/",selected.species, "-", sub("^([^-]+-){4}([^-]+).*", "\\2", rds_files_non_urban[i]),".RData"))
  
}


# load (paste0("./data/",selected.species, "-", sub("^([^-]+-){4}([^-]+).*", "\\2", rds_files_non_urban[i]),".RData"))
#   
# test <- df.list
#   xx <- do.call (rbind, test)


#   
#   st_write(xx, paste0("./data/",selected.species, "-", sub("^([^-]+-){4}([^-]+).*", "\\2", rds_files_non_urban[i]),".shp")) 
#   
#   all.list[[i]] <- xx
#   
# }  

