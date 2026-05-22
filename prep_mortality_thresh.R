## prep_mortality_thresh.R
# Set bar on how large trees can get

rm(list = ls())

## PACKAGES ==================
library(plyr)
library(dplyr)
library(reshape2)
library(readxl)
library(magrittr)

## IMPORT DATA ==================
bigtree_init <- readRDS("bigtree_init.rds") # only big trees within reserve boundaries
bigtree_sp <- unique(bigtree_init$sp)

range(bigtree_init$dbh) # 300 to 1547 


## GENERATE PARAMETER FILE =========================
mort_size_thresh <- list()
for(i in 1:length(bigtree_sp)){
  mort_size_thresh[[bigtree_sp[i]]] <- 1600 
}

saveRDS(mort_size_thresh, file = "mortality_thresh_st.rds")

