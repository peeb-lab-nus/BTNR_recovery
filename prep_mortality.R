## prep_mortality.R
# fit mortality model
# assign them to species in our simulation

## DIRECTORY ==================
rm(list = ls())

## PACKAGES ==================
library(RColorBrewer)
library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)
library(readxl)
library(magrittr)
library(stringr)
library(plyr)
library(dplyr)
library(stringr)
library(greta)
library(bayesplot)
library(MCMCvis)
library(tidybayes)

source(file.path(src.dir, "utils.R"))

## IMPORT DATA ==================
# Import big tree data -----------
bigtree_init <- readRDS("bigtree_init.rds") # only big trees within reserve boundaries
bigtree_sp <- unique(bigtree_init$sp)

# Import functional group data -----------
btbt_grouping <- read.csv("species_trait.csv")

# Import plot data -----------
pri_plot_data_fnames <- list.files(pattern = "bukittimah_primary", full.names = T)
for(i in pri_plot_data_fnames){load(i)}

#fix known errors for primary forest
bukittimah_primary.stem6$dbh[bukittimah_primary.stem6$tag == 10564] <- 121 # 12.1
bukittimah_primary.stem6$dbh[bukittimah_primary.stem6$tag == 4873] <- 512 # 51.2
bukittimah_primary.stem6$dbh[bukittimah_primary.stem6$tag == 4890] <- 145 # 14.5
bukittimah_primary.stem6$dbh[bukittimah_primary.stem6$tag == 10874] <- 165 # 16.5
bukittimah_primary.stem6$dbh[bukittimah_primary.stem6$tag == 10789] <- 384 # 38.4
bukittimah_primary.stem6$dbh[bukittimah_primary.stem6$tag == 1437] <- 1000 # 100

pri_plot_data_list <- list(bukittimah_primary.stem1, # 1993
                           bukittimah_primary.stem2, # 1995
                           bukittimah_primary.stem4, # 2003
                           bukittimah_primary.stem5, # 2007
                           bukittimah_primary.stem6, # 2012
                           bukittimah_primary.stem7) # 2018

sec_plot_data_fnames <- list.files(pattern = "bukittimah_secondary", full.names = T)
for(i in sec_plot_data_fnames){load(i)}

sec_plot_data_list <- list(bukittimah_secondary.stem1,
                           bukittimah_secondary.stem2,
                           bukittimah_secondary.stem3,
                           bukittimah_secondary.stem4)


# Generate mortality observations
target.col <- c("sp","stemID", "CensusID", "DFstatus", "ExactDate", "dbh")
pri_mort_data_list <- list()
for (i in 1:(length(pri_plot_data_list)-1)){
  pri_mort_data_list[[i]] <- left_join(pri_plot_data_list[[i]][,target.col],
                                       pri_plot_data_list[[i+1]][,target.col],
                                       by = c("stemID", "sp"),
                                       suffix = c(".before", ".after"))
}
pri_mort_df <- do.call("rbind", pri_mort_data_list)

sec_mort_data_list <- list()
for (i in 1:(length(sec_plot_data_list)-1)){
  sec_mort_data_list[[i]] <- left_join(sec_plot_data_list[[i]][,target.col],
                                       sec_plot_data_list[[i+1]][,target.col],
                                       by = c("stemID", "sp"),
                                       suffix = c(".before", ".after"))
}
sec_mort_df <- do.call("rbind", sec_mort_data_list)

# Clean up mortality observations
cleanDBHdata <- function(x){
  return(as.numeric(ifelse(grepl(x, pattern = "[A-z]"), NA, x)))
}

pri_mort_df$dbh.before <- cleanDBHdata(pri_mort_df$dbh.before)
pri_mort_df$dbh.after <- cleanDBHdata(pri_mort_df$dbh.after)
sec_mort_df$dbh.before <- cleanDBHdata(sec_mort_df$dbh.before)
sec_mort_df$dbh.after <- cleanDBHdata(sec_mort_df$dbh.after)

pri_mort_df_clean <- pri_mort_df %>% 
  subset(dbh.before >= 10) %>%
  subset(DFstatus.before == "alive") %>%
  subset(DFstatus.after %in% c("alive","dead")) %>%
  subset(!is.na(ExactDate.before)) %>%
  subset(!is.na(ExactDate.after)) %>%
  subset(!is.na(dbh.before)) %>%
  left_join(y = btbt_grouping[,c("sp","ecology")], by = "sp") %>%
  mutate(Thresh = dbh.before >= median(dbh.before)) %>%
  mutate(Time = as.numeric(ExactDate.after - ExactDate.before) / 365.25) %>% 
  mutate(Died = ifelse(DFstatus.before == "alive" &
                         DFstatus.after == "dead", 1, 0))

sec_mort_df_clean <- sec_mort_df %>% 
  subset(dbh.before >= 10) %>%
  subset(DFstatus.before == "alive") %>%
  subset(DFstatus.after %in% c("alive","dead")) %>%
  subset(!is.na(ExactDate.before)) %>%
  subset(!is.na(ExactDate.after)) %>%
  subset(!is.na(dbh.before)) %>%
  left_join(y = btbt_grouping[,c("sp","ecology")], by = "sp") %>%
  mutate(Thresh = dbh.before >= median(dbh.before)) %>%
  mutate(Time = as.numeric(ExactDate.after - ExactDate.before) / 365.25) %>% 
  mutate(Died = ifelse(DFstatus.before == "alive" &
                         DFstatus.after == "dead", 1, 0))

## FIT MORTALITY MODEL ==================
CENSUS_survival_df<-rbind(pri_mort_df_clean,sec_mort_df_clean)
CENSUS_survival_df <- CENSUS_survival_df %>% 
  filter(ecology=="shade-tolerant")
t<-as_data(CENSUS_survival_df$Time)
dbh<-as_data(CENSUS_survival_df$dbh.before/100)

mortality_stat <- as_data(CENSUS_survival_df$Died) ##1 means dead

# define priors and relationships between covariates
b <- normal(0, 5)
c <- normal(0, 5)
a <- normal(0, 5)

# Define mortality model
h <- exp((a/dbh) + (b*dbh) + c)
mort <- 1 - exp(-h*t) #this is between 0 to 1 so it doesn't need a ilogit link

distribution(mortality_stat) <- bernoulli(mort)

model_mort <- model(a,b,c)

draws_mort <- mcmc(model_mort, n_samples = 4000,warmup=4000,chains=4)

#model_mort <- model

####diagnostics####

mcmc_trace(draws_mort, pars = c("a","b","c"))

mcmc_trace(draws_mort, pars = c("a_int","a_sd","b_int",
                                "b_sd","c_int","c_sd")) #trace plot for convergence
mcmc_trace(draws_mort, pars = c("a_z[1,1]","c_z[1,1]","b_z[1,1]")) #trace plot for convergence


#we want rhat<1.1, neff>2000
draws.sum<-MCMCsummary(draws_mort)
print(length(which(draws.sum$Rhat>1.1)))
non.converged<-print(draws.sum[which(draws.sum$Rhat>1.1),])

params_ecology<-draws_mort %>%
  spread_draws(a,b,c) %>%
  median_hdci() %>% 
  select(a,b,c) %>% 
  mutate(ecology="shade-tolerant")#get only median values for now. Alternatively, we can draw from the posterior distribution.

#save ecology only parameters
write.csv(params_ecology,file="mort_rates_st_combined.csv")

## 
mort_rate_st <- read.csv("mort_rates_st_combined.csv")

mort_rate_grp<- subset(btbt_grouping, sp %in% bigtree_sp)  %>%
  filter(ecology=="shade-tolerant") %>% 
  left_join(mort_rate_st, by = "ecology") %>% 
  dlply(.variables = .(sp), .fun = summarise,
        a = a,
        b = b,
        c = c)

saveRDS(mort_rate_grp, file =  "mortality_rate_para.rds")
