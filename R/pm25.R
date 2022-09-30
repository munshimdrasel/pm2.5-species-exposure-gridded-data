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

#pm data from Randal Martins group
# https://wustl.app.box.com/v/ACAG-V5GL02-GWRPM25/folder/148055284876

RM.dir <- "./data/randal_martin_pm2.5/"


years <- 1998:2020
us_states <- USAboundaries::us_states()  %>%  filter(!stusps %in% c("PR", "HI", "AK") )
states.vec <- as.vector(us_states$stusps)
states.name <- as.vector(us_states$name)


options(timeout=600)


# we want to use an equal area projection, here's one I often use:
p4s <- "+proj=lcc +lat_1=33 +lat_2=45 +lat_0=40 +lon_0=-97 +a=6370000 +b=6370000"

df.list <- list()

#Using 2019 gridded data as 2020 also since 2020 data pm2.5 not available

for (i in 1:length(years)) {
  print( years[i])
  file_in<- paste0( RM.dir, "V5GL02.HybridPM25.NorthAmerica.", years[i],"01-",years[i], "12.nc")
  
  print(file_in)
  
  # see ?raster::brick, scroll down to details for info on reading NetCDF
  ncin_raster <- raster::brick( file_in, varname = "GWRPM25",
                                lvar = 1, level = 5)

  
  #zooming over CONUS US
  box_use <- c( xmin = -126,
                xmax = -66,
                ymin = 23,
                ymax = 51)
  
  new_raster <-crop(x = ncin_raster, y = box_use)
  
  new_raster <- aggregate( new_raster, fact=10) #increasing resolution
  
  # dt.raster<- rasterize(x=dt[, 1:2], y=ras_dom, field=dt[, 3], fun = mean)
  
  spatpoly <- rasterToPolygons( new_raster)
  spatpoly <- spTransform(spatpoly, "+proj=longlat +datum=WGS84")
  
  # convert to simple features (best integration with ggplot)
  spatpoly_sf <- as(spatpoly,'sf')
  
  spatpoly_sf <- st_transform( spatpoly_sf, crs = p4s)

  for (j in length(states.vec)) {
    block.groups <- block_groups(states.vec[j])
    block.groups.transform <- block.groups %>% st_transform(crs=p4s)
    bg.selected <- block.groups.transform %>% dplyr::select(STATEFP, COUNTYFP, TRACTCE, BLKGRPCE, geometry)
    
    conc.bg <- st_interpolate_aw( spatpoly_sf["layer"], block.groups.transform,   extensive = F)
    vec <- as.vector(intersect(bg.selected$geometry, conc.bg$geometry))
    bg.new <- bg.selected %>%  filter(geometry %in% vec)
    conc.bg$species <- "pm2.5"
    conc.bg$year <- years[i]
    conc.bg$state <- states.vec[j]
    conc<- cbind (bg.new,conc.bg) %>% dplyr::select(-geometry.1)
    # ggplot( conc) + geom_sf( aes( fill = layer,  geometry = geometry),
    #                        color = NA) +scale_fill_viridis_c( option="plasma")
    # 
    # cbind(bg.selected, conc.bg)
    df.list[[j]] <- conc
    
  }  
  save(df.list, file=paste0("./data/pm2.5-",years[i],".RData"))
  
}


# 
# 
#   }
#   
#   }
#   
# pm_data <- lapply( 2015:2020,
#                    function( yr){
#                      print( yr)
#                      file_in<- paste0( RM.dir, "V5GL02.HybridPM25.NorthAmerica.", yr,"01-",yr, "12.nc")
#                      
#                      print(file_in)
#                      
#                      # see ?raster::brick, scroll down to details for info on reading NetCDF
#                      ncin_raster <- raster::brick( file_in, varname = "GWRPM25",
#                                                    lvar = 1, level = 5)
#                      
#                     
#                      
#                      
#                      #zooming over Eagle ford shale area
#                      box_use <- c( xmin = -126,
#                                    xmax = -66,
#                                    ymin = 23,
#                                    ymax = 51)
#                      
#                      new_raster <-crop(x = ncin_raster, y = box_use)
#                      
#                      # dt.raster<- rasterize(x=dt[, 1:2], y=ras_dom, field=dt[, 3], fun = mean)
# 
#                      spatpoly <- rasterToPolygons( new_raster)
#                      spatpoly <- spTransform(spatpoly, "+proj=longlat +datum=WGS84")
#                      
#                      # convert to simple features (best integration with ggplot)
#                      spatpoly_sf <- as(spatpoly,'sf')
#                      
#                      spatpoly_sf <- st_transform( spatpoly_sf, crs = p4s)
#                      
#                      spatpoly_sf$year <- yr
#                      
#                      interp_col <- c( 'layer')
#                      pm_poly <- st_interpolate_aw( spatpoly_sf  [spatpoly_sf$year == yr, interp_col], exp_inverse_dist_all.sf, 
#                                                    extensive = F)
#                      
#                      pm_efs <- pm_poly  %>%  filter(geometry %in% unique.geometry)
#                      
#                      
#                      pm_efs$year <- yr
#                      return(pm_efs)
#                      
#                    }) 
# 
# pm_poly <- do.call(rbind, pm_data)
# 
# pm_poly <- pm_poly  %>%  filter(geometry %in% pop_race$geometry)

