## prep_recruit_thresh.R:
# Recruitment rates are annualized rates for the whole population and not at the individual-level

## PACKAGES & DIRECTORIES ==================
rm(list = ls())

## IMPORT DATA ==================
bigtree <- readRDS("bigtree_init.rds") # only big trees within reserve boundaries
target_sp <- unique(bigtree$sp)

recruit_thresh_para <- list()
for(i in 1:length(target_sp)){
  recruit_thresh_para[[target_sp[i]]] <- 300
}

saveRDS(recruit_thresh_para, file = "recruit_thresh_para.rds")