## prep_recruit_rate.R
# Calculate recruitment rates from censuses

## PACKAGES & DIRECTORIES ==================
rm(list = ls())
library(plyr)
library(dplyr)
library(reshape2)
library(readxl)
library(magrittr)

source(file.path("CTFSRPackage", "startCTFS", "startCTFS.r"))
startCTFS(folder = "CTFSRPackage")

## IMPORT DATA ==================
# Import big tree data -----------
bigtree <- readRDS("bigtree_init.rds") # only big trees within reserve boundaries
length(unique(bigtree$sp)) # 408 species

# Import groupings -----------
btbt_eco <- read.csv("species_trait.csv")
splist <- unique(bigtree$sp)

## CALCULATE RECRUITMENT RATES FROM CENSUS DATA ==================
# Import census data
pri_plot_data_fnames <- list.files(pattern = "bukittimah_primary", full.names = T)

for(i in pri_plot_data_fnames){load(i)}

# Correct some errors
bukittimah_primary.stem6$dbh[bukittimah_primary.stem6$tag == 10564] <- 121 # 12.1
bukittimah_primary.stem6$dbh[bukittimah_primary.stem6$tag == 4873] <- 512 # 51.2
bukittimah_primary.stem6$dbh[bukittimah_primary.stem6$tag == 4890] <- 145 # 14.5
bukittimah_primary.stem6$dbh[bukittimah_primary.stem6$tag == 10874] <- 165 # 16.5
bukittimah_primary.stem6$dbh[bukittimah_primary.stem6$tag == 10789] <- 384 # 38.4
bukittimah_primary.stem6$dbh[bukittimah_primary.stem6$tag == 1437] <- 1000 # 100

pri_census_list <- list(bukittimah_primary.stem1,
                        bukittimah_primary.stem2,
                        bukittimah_primary.stem4,
                        bukittimah_primary.stem5,
                        bukittimah_primary.stem6,
                        bukittimah_primary.stem7)


sec_plot_data_fnames <- list.files(pattern = "bukittimah_secondary", full.names = T)
for(i in sec_plot_data_fnames){load(i)}

sec_census_list <- list(bukittimah_secondary.stem1,
                        bukittimah_secondary.stem2,
                        bukittimah_secondary.stem3,
                        bukittimah_secondary.stem4)

# Define size thresholds
size_threshold_df <- data.frame(ecology = c("shade-intolerant",
                                            "shade-tolerant",
                                            "sub-canopy",
                                            "understorey"),
                                size_threshold = c(150,300,150,30))

# Bulk merge and clean
cleanCensusData <- function(x){
  x$dbh <- as.numeric(ifelse(grepl(x$dbh, pattern = "[A-z]"), NA, x$dbh))
  x$date <- as.numeric(ifelse(grepl(x$date, pattern = "[A-z]"), NA, x$date))
  return(x)
}

pri_census_list_clean <- pri_census_list %>% 
  # Convert NULLs to numeric
  lapply(FUN = cleanCensusData) %>%
  lapply(FUN = function(x) {left_join(x, btbt_eco, by = "sp")}) %>%
  lapply(FUN = function(x) {left_join(x, size_threshold_df, by = "ecology")})

sec_census_list_clean <- sec_census_list %>% 
  # Convert NULLs to numeric
  lapply(FUN = cleanCensusData) %>%
  lapply(FUN = function(x) {left_join(x, btbt_eco, by = "sp")}) %>%
  lapply(FUN = function(x) {left_join(x, size_threshold_df, by = "ecology")})

# Calculate recruitment rates
calcRecruitment <- function(x, y, census_list){
  # N2 = no. of individuals at time 2
  # R = no. of recruits
  # time = average time elapsed between censuses
  
  nrepro <- nrow(subset(census_list[[x]], status == "A" & dbh >= size_threshold))
  ntotal <- nrow(subset(census_list[[x]], status == "A"))
  recr <- recruitment(census1 = census_list[[x]],
                      census2 = census_list[[y]],
                      mindbh = 10)
  census <- x
  nrecr <- recr$R
  N2 <- recr$N2
  time <- recr$time
  
  data.frame(census, nrepro, ntotal, nrecr, N2, time)
}

# wrapper function
summarizeRecruitment <- function(census_list){
  # iterate across census intervals
  res <- mapply(x = 1:(length(census_list)-1),
         y = 2:length(census_list),
         FUN = calcRecruitment,
         MoreArgs = list("census_list" = census_list),
         SIMPLIFY = FALSE)
  res_df <- do.call("rbind",res)
  # calculate annual recruitment
  res_df$nrecr_ann <- res_df$nrecr / res_df$time
  # calculate per capita annual recruitment
  res_df$per_cap_nrecr_ann <- res_df$nrecr_ann / res_df$nrepro
  return(res_df)
}


pri_shadetol_rate <- summarizeRecruitment(lapply(pri_census_list_clean, FUN = function(x) subset(x, ecology == "shade-tolerant") ))
pri_shadetol_rate$ecology <- "shade-tolerant"

# Secondary forest
sec_shadetol_rate <- summarizeRecruitment(lapply(sec_census_list_clean, FUN = function(x) subset(x, ecology == "shade-tolerant") ))
sec_shadetol_rate$ecology <- "shade-tolerant"

# Average rate
plot_repro_rate <- rbind(pri_shadetol_rate, sec_shadetol_rate) %>% 
  mutate("rate" = per_cap_nrecr_ann)

avg_repro_rate <- plot_repro_rate %>%
  ddply(.variables = .(ecology),
        .fun = summarize,
        rate = mean(per_cap_nrecr_ann))


recr_rates_df <- subset(btbt_eco, sp %in% splist) %>%
  left_join(avg_repro_rate, by = "ecology")

recr_sp_rates_list <- dlply(.data = recr_rates_df,
                          .variables = .(sp),
                          .fun = function(x) unlist(x[["rate"]]))

saveRDS(recr_sp_rates_list, file = "recruit_rate_para.rds")
