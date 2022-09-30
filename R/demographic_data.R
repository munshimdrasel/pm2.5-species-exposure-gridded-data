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

#variables list
# https://api.census.gov/data/2020/acs/acs5/variables.html

years <- 2016:2019
us_states <- USAboundaries::us_states()  %>%  filter(!stusps %in% c("PR", "HI", "AK") )
states.vec <- as.vector(us_states$stusps)
states.name <- as.vector(us_states$name)


options(timeout=600)


# we want to use an equal area projection, here's one I often use:
p4s <- "+proj=lcc +lat_1=33 +lat_2=45 +lat_0=40 +lon_0=-97 +a=6370000 +b=6370000"

ras_dom<-raster(xmn=-126, xmx=-66, ymn=23, ymx=51,
                resolution=c(0.1,0.1))

all.list <- list()
df.list <- list()

race_vars <- c(
  White = "B03002_003",
  Black = "B03002_004",
  Native = "B03002_005",
  Asian = "B03002_006",
  HIPI = "B03002_007",
  Hispanic = "B03002_012"
)

med_household_income <- c( med_household_income = "B19013_001",
                           white_med_household_income= "B19013A_001E",
                           black_med_household_income = "B19013B_001E",
                           native_med_household_income= "B19013C_001E",
                           asian_med_household_income= "B19013D_001E",
                           hipi_med_household_income="B19013E_001E",
                           hispanic_med_household_income="B19013I_001E")


df.list <- list()

for (i in 1:length(years)) {
  for (j in 1:length(states.vec)) {
    
    #population
    race_block_groups <- get_acs(
      geography = "block group",
      state = states.vec[j],
      variables = race_vars,
      summary_var = "B03002_001",
      year =  years[i],
      geometry=TRUE,
      survey = "acs5",
      output="wide" )  %>% st_transform(crs = p4s)
    
    
    income_block_groups <- get_acs(geography = "block group",
                                   state = states.vec[j],
                                   geometry = "TRUE",
                                   year = years[i],
                                   survey = "acs5",
                                   variables = med_household_income,
                                   output="wide")  %>% st_transform(crs = p4s) %>% 
      as.data.frame() %>% dplyr::select(-NAME, -geometry) 
    
    
    pop.income <- merge(race_block_groups,income_block_groups, by ="GEOID") 
    df.list[[j]] <- pop.income
  }
  
  save(df.list, file=paste0("./data/pop.income.blk.grp-", years[i],".RData"))
  
}
