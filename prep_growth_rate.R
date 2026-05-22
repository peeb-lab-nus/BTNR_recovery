## prep_growth_rate_pooled.R: Generate size-dependent growth model with pooled data from pri and sec plots

## PACKAGES & DIRECTORIES ==================
rm(list = ls())
library(plyr)
library(reshape2)
library(readxl)
library(dplyr)

## IMPORT DATA ==================
bigtree_init <- readRDS("bigtree_init.rds")
btbt_eco <- read.csv("species_trait.csv")
target_sp <- unique(bigtree_init$sp)

pri_plot_data_fnames <- list.files(pattern = "bukittimah_primary", full.names = T)
for (i in pri_plot_data_fnames) {
  load(i)
}

# fix known errors for primary forest
bukittimah_primary.stem6$dbh[bukittimah_primary.stem6$tag == 10564] <- 121 # 12.1
bukittimah_primary.stem6$dbh[bukittimah_primary.stem6$tag == 4873] <- 512 # 51.2
bukittimah_primary.stem6$dbh[bukittimah_primary.stem6$tag == 4890] <- 145 # 14.5
bukittimah_primary.stem6$dbh[bukittimah_primary.stem6$tag == 10874] <- 165 # 16.5
bukittimah_primary.stem6$dbh[bukittimah_primary.stem6$tag == 10789] <- 384 # 38.4
bukittimah_primary.stem6$dbh[bukittimah_primary.stem6$tag == 1437] <- 1000 # 100

pri_plot_data_list <- list(
  bukittimah_primary.stem1,
  bukittimah_primary.stem2,
  bukittimah_primary.stem4,
  bukittimah_primary.stem5,
  bukittimah_primary.stem6,
  bukittimah_primary.stem7
)

sec_plot_data_fnames <- list.files(pattern = "bukittimah_secondary", full.names = T)
for (i in sec_plot_data_fnames) {
  load(i)
}

sec_plot_data_list <- list(
  bukittimah_secondary.stem1,
  bukittimah_secondary.stem2,
  bukittimah_secondary.stem3,
  bukittimah_secondary.stem4
)

# Clean up census data ===============
cleanCensusData <- function(x) {
  x$dbh <- as.numeric(ifelse(grepl(x$dbh, pattern = "[A-z]"), NA, x$dbh))
  x$date <- as.numeric(ifelse(grepl(x$date, pattern = "[A-z]"), NA, x$date))
  return(x)
}

pri_plot_data_clean <- pri_plot_data_list %>%
  # Convert NULLs to numeric
  lapply(FUN = cleanCensusData) %>%
  lapply(FUN = function(x) {
    left_join(x, btbt_eco, by = "sp")
  }) %>%
  lapply(FUN = function(x) {
    x <- x %>%
      mutate(dbh = case_when(status != "A" | DFstatus != "alive" ~ NA,
        .default = dbh
      )) ## all non-alive trees have NA for dbh
  })

sec_plot_data_clean <- sec_plot_data_list %>%
  # Convert NULLs to numeric
  lapply(FUN = cleanCensusData) %>%
  lapply(FUN = function(x) {
    left_join(x, btbt_eco, by = "sp")
  }) %>%
  lapply(FUN = function(x) {
    x <- x %>%
      mutate(dbh = case_when(status != "A" | DFstatus != "alive" ~ NA,
        .default = dbh
      )) ## all non-alive trees have NA for dbh
  })


# CALCULATE GROWTH RATES ===============
# Combine data
primary_forest_census_list <- list()
for (i in 1:(length(pri_plot_data_clean) - 1)) {
  dbh_pre <- pri_plot_data_clean[[i]]$dbh
  dbh_post <- pri_plot_data_clean[[i + 1]]$dbh
  date_pre <- pri_plot_data_clean[[i]]$date
  date_post <- pri_plot_data_clean[[i + 1]]$date
  hom_pre <- as.numeric(pri_plot_data_clean[[i]]$hom)
  hom_post <- as.numeric(pri_plot_data_clean[[i + 1]]$hom)
  hom_change <- hom_post - hom_pre
  time <- (date_post - date_pre) / 365.25
  sp <- pri_plot_data_clean[[i]]$sp
  stemID <- pri_plot_data_clean[[i]]$stemID
  plot <- "primary"
  primary_forest_census_list[[i]] <- data.frame(
    dbh_pre,
    dbh_post,
    date_pre,
    date_post,
    sp,
    time,
    stemID,
    plot,
    hom_change
  )
}
primary_forest_census_df <- do.call("rbind", primary_forest_census_list)

secondary_forest_census_list <- list()
for (i in 1:(length(sec_plot_data_clean) - 1)) {
  dbh_pre <- sec_plot_data_clean[[i]]$dbh
  dbh_post <- sec_plot_data_clean[[i + 1]]$dbh
  date_pre <- sec_plot_data_clean[[i]]$date
  date_post <- sec_plot_data_clean[[i + 1]]$date
  time <- (date_post - date_pre) / 365.25
  hom_pre <- as.numeric(sec_plot_data_clean[[i]]$hom)
  hom_post <- as.numeric(sec_plot_data_clean[[i + 1]]$hom)
  hom_change <- hom_post - hom_pre
  sp <- sec_plot_data_clean[[i]]$sp
  stemID <- sec_plot_data_clean[[i]]$stemID
  plot <- "secondary"
  secondary_forest_census_list[[i]] <- data.frame(
    dbh_pre,
    dbh_post,
    date_pre,
    date_post,
    sp,
    time,
    stemID,
    plot,
    hom_change
  )
}

secondary_forest_census_df <- do.call("rbind", secondary_forest_census_list)

combined_df <- rbind(primary_forest_census_df, secondary_forest_census_df) %>%
  left_join(btbt_eco[c("sp", "ecology")], by = "sp") %>%
  subset(ecology == "shade-tolerant")

### removing outliers and trees with changes in hom
### following Condit (2017), remove all trees >75 mm growth or shrank by >4s
# s = 0.006214dbh + 0.9036

combined_df2 <- combined_df %>%
  subset(dbh_pre >= 10 & dbh_post >= 10) %>%
  .[complete.cases(.), ] %>%
  filter(hom_change == "0") %>%
  mutate(growth = (dbh_post - dbh_pre) / time) %>%
  filter(growth < 75) %>%
  filter(dbh_post > dbh_pre - 4 * (0.006214 * dbh_pre + 0.9036)) %>%
  mutate(size_bin = cut(dbh_pre,
    breaks = c(10, 20, 50, 100, 300, 10000),
    right = FALSE,
    include.lowest = TRUE
  ))


nrow(combined_df2) # 24452 observations
table(combined_df2$size_bin)

# Average rates
all_species_growth_rate <- tapply(combined_df2$growth, INDEX = combined_df2$size_bin, FUN = mean)
all_species_growth_rate_df <- data.frame(
  "avg_growth_rate" = all_species_growth_rate,
  "size_bin" = factor(levels(combined_df2$size_bin),
    levels = levels(combined_df2$size_bin)
  ),
  "species" = "All shade-tolerant species"
)

### export parameter files
growth_rates <- all_species_growth_rate_df$avg_growth_rate
names(growth_rates) <- all_species_growth_rate_df$size_bin

sp_rates_list <- rep(list(growth_rates), length(target_sp))
names(sp_rates_list) <- target_sp

saveRDS(sp_rates_list, "growth_rate_para.rds")
