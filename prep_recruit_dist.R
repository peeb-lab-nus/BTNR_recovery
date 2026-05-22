## prep_recruit_dist.R:
# Assigning dispersal distances to each species based on dispersal syndrome

## PACKAGES & DIRECTORIES ==================
rm(list = ls())

## IMPORT DATA ==================
# Import big tree dataset
bigtree <- readRDS("bigtree_init.rds")

btbt_dispersal <- read.csv("species_trait.csv")

## PREPARE RECRUITMENT DISTANCE ==================
table(btbt_dispersal$dispersal1)
table(btbt_dispersal$dispersal2)

recr_dist <- list()
for(i in unique(bigtree$sp)){
  if(btbt_dispersal$dispersal2[btbt_dispersal$sp == i] %in% c("animal_large", "animal_small", "wind", "wind/water")){
    recr_dist[[i]] <- 25
  } else {
    recr_dist[[i]] <- 5
  }
}
saveRDS(recr_dist, "recruit_dist_para.rds")

recr_dist_far <- list()
for(i in unique(bigtree$sp)){
  if(btbt_dispersal$dispersal2[btbt_dispersal$sp == i] %in% c("animal_large", "animal_small", "wind", "wind/water")){
    recr_dist_far[[i]] <- 100
  } else {
    recr_dist_far[[i]] <- 20
  }
}
saveRDS(recr_dist_far, "recruit_dist_far_para.rds")