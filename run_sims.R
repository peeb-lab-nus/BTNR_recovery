### sim_hpc.R
# Run initialisation and 
# This script was designed to be run on a high performance computing cluster (HPC)

## PACKAGES AND DIRECTORIES ===============
rm(list = ls())
library(plyr)
library(raster)
library(sp)
library(doParallel)
library(foreach)
library(dplyr)
library(sf)
library(stringr)
library(doRNG)

#cluster <- makeCluster(detectCores())
#registerDoParallel(cluster)
registerDoParallel(detectCores())

source("utils.R")

mort_rate <- readRDS("mortality_rate_para.rds")
mort_thresh <- readRDS("mortality_thresh_para.rds")
grow_rate <- readRDS("growth_rate_para.rds")
recr_rate <- readRDS("recruit_rate_para.rds")

#recr_dist <- readRDS("results/recruit_dist_para.rds") # change this to modify dispersal distance
recr_dist <- readRDS("recruit_dist_far_para.rds")
recr_thresh <- readRDS("recruit_thresh_para.rds")

btnr <- readRDS(file = "BT_polygon.rds")
btnr <- sf::st_as_sf(btnr)

bigtree_init <- readRDS("bigtree_init.rds") # initial conditions: just the big trees

bigtree_init$id <- 1:nrow(bigtree_init) # assign IDs to the original big trees so we can track them
bigtree_init <- bigtree_init %>% dplyr::select(c("sp", "dbh", "gx", "gy", "id"))

target_sp<- unique(bigtree_init$sp)

grid_btnr <- readRDS("grid_btnr_1_2.rds") # using 1.2 m grid cells

nstart=1
nend=1000

for ( m in 1:length(target_sp)){
  init_alt2 <- subset(bigtree_init, sp == target_sp[m])
  
  # Need to assign trees to gridIDs
  init_alt2$gridID <- init_alt2 %>% 
    st_as_sf(coords = c("gx", "gy")) %>%
    st_within(., grid_btnr) %>%
    unlist()
  
  
 foreach(q = nstart:nend,.packages = c("plyr","dplyr","stringr","sf")) %dorng% {
    
   initialized_pop<-initialize_BDG_model(init_pop = init_alt2,
                         growth_rate = grow_rate,
                         grow_exceptions = init_alt2$id, # og big trees don't grow
                         mort_rate_type = "size",
                         mort_rate = mort_rate,
                         mort_thresh = mort_thresh,
                         mort_exceptions = init_alt2$id, # og big trees don't die
                         recruit_thresh = recr_thresh,
                         recruit_rate = recr_rate,
                         recruit_dist = recr_dist,
                         boundary = btnr,
                         grid = grid_btnr,
                         grid_limit = TRUE,
                         output_dir = "output",
                         output_name = paste0(target_sp[m], "_", q),
                         save_ts = FALSE)
    
    
    run_BDG_model(init_pop = initialized_pop,
                           growth_rate = grow_rate,
                           mort_rate = mort_rate,
                           mort_thresh = mort_thresh,
                           mort_rate_type = "size",
                           recruit_rate = recr_rate,
                           recruit_thresh = recr_thresh,
                           recruit_dist = recr_dist,
                           boundary = btnr,
                           print_output = TRUE,
                           grid = grid_btnr,
                           grid_limit = TRUE,
                           output_name = paste0(target_sp[m],"_",q),
                           output_dir = "output",
                           nyear = 1000,
                           print_freq = 50)
    
    
  }
}