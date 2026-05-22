# utils_alt.R: Utility functions for the simulation model
# Functions to run simulation models without area-specific rates

grow_things_categ <- function(x, growth_rate, exceptions = NULL){
  # Deterministic growth function, categorical relationship with size
  # Function is iterated over species using ddply
  #
  # Arguments:
  #     x = data.frame, population data
  #     growth_rate = named list of vectors, species-specific growth rates for each size-category
  #
  # Returns:
  #     data.frame, population data
  
  sp <- x$sp[1]
  size_class <- cut(x$dbh, breaks = c(10, 20, 50, 100, 300, 10000),
                           right = FALSE, include.lowest = TRUE)
  if(is.null(exceptions)){
    # if no exceptions provided, then each stem grows by the model provided
    x$dbh <- x$dbh + growth_rate[[sp]][as.numeric(size_class)]
  } else {
    # if exceptions provided, only stems that not excepted grow
    x$dbh[!x$id %in% exceptions] <- x$dbh[!x$id %in% exceptions] + growth_rate[[sp]][as.numeric(size_class)[!x$id %in% exceptions]]
  }
  
  return(x)
}

calc_mort_rate <- function(a, b, c, x, t){
  # Calculates the annual mortality probability, with some modifications
  # given a set of parameters for the Bohlman & Pacala (2010)
  # mortality model:
  # logit(P) = a * dbh-1 + b * dbh + c
  # 
  # Arguments:
  #    a = numeric, coefficient for reciprocal term
  #    b = numeric, coefficient for linear term
  #    c = numeric, coefficient for intercept
  #    x = numeric, dbh of individual
  #    s = numeric, duration of observations (years)
  #    dbh_scale = numeric, scaling factor for dbh
  #    t_scale = numeric, scaling factor for time
  
  q<- x /100
  h <- exp((a / q) + b*q + c)
  mu = 1-exp(-h*t) # annual probability, t = 1
  return(mu)
}

kill_things_equal_rate <- function(x, 
                                   mort_rate,
                                   exceptions = NULL, 
                                   mort_thresh = NULL){
  # Stochastic death function using continuous diameter
  # Exceptions has priority (i.e., if mort_rate set to 1 and/or individual is above size-threshold
  # but individual is excepted, then mort probability = 0)
  # Assume mortality rate is not size dependent
  #
  # Arguments:
  #     x = data.frame, starting population
  #     mort_rate = list, contain model parameters for mortality model
  #     exceptions = numeric, vector of ids to be immortal
  #
  # Returns:
  #     data.frame, population after mortality applied
  sp <- x$sp[1]

  mort_rate_vector <- rep(mort_rate, length(x$sp))

  # If individual is above size threshold, probability of mortality = 1
  if(!is.null(mort_thresh)){
    mort_rate_vector[x$dbh >= mort_thresh[[sp]] ] <- 1
  }

  # If individual is an exception, probability of mortality = 0
  if(!is.null(exceptions)){
    mort_rate_vector[x$id %in% exceptions] <- 0
  }
  
  # Generate random uniform numbers
  prob <- runif(nrow(x), min = 0, max = 1)
  #deaths <- prob <= mort_rate_vector

  # If random uniform smaller than probability of mortality, then mark for death
  x$deaths <- prob <= mort_rate_vector
  
  return(x)

}

kill_things <- function(x,
                        mort_rate,
                        exceptions = NULL,
                        mort_thresh = NULL){
  # Stochastic death function using continuous diameter
  # Assume mortality rate is size dependent
  # Applied to each species
  #
  # Arguments:
  #     x = data.frame, starting population
  #     mort_rate = list, contain model parameters for mortality model
  #     exceptions = numeric, vector of ids to be immortal
  #
  # Returns:
  #     data.frame, population after mortality applied
  
  
  sp <- x$sp[1]
  
  mort_rate_vector <- calc_mort_rate(a =mort_rate[[sp]]$a,
                                     b =mort_rate[[sp]]$b,
                                     c =mort_rate[[sp]]$c, 
                                     x = x$dbh,
                                     t = 1)
  
  # If individual is above size threshold, probability of mortality = 1
  if(!is.null(mort_thresh)){
    mort_rate_vector[x$dbh >= mort_thresh[[sp]] ] <- 1
  }

  #If individual is an exception, probability of mortality = 0
  if(!is.null(exceptions)){
    mort_rate_vector[x$id %in% exceptions] <- 0
  }
  
  # Generate random uniform numbers
  prob <- runif(nrow(x), min = 0, max = 1)
  #deaths <- prob <= mort_rate_vector

  # If random uniform smaller than probability of mortality, then mark for death
  x$deaths <- prob <= mort_rate_vector
  

  return(x)
}



recruit_things <- function(x,
                            recruit_thresh,
                            recruit_rate,
                            recruit_dist,
                            boundary){
  # Stochastic recruitment function
  # 
  # Arguments:
  #     x: population
  #     recruit_thresh: named list, with each element containing the reproductive size threshold for each species
  #     recruit_rate: named list, with each element containing the per individual annual recruitment rate for each species
  #     recruit_dist: named list, with each element containing the average dispersal distance for each species
  #     boundary: SpatialPolygons, boundaries where dispersal is limited. Same coordinate systems as coordinates of individuals in population
  # 
  # Returns:
  #     data.frame, population after recruitment applied
  # x = subset(total_init, sp == "CRATCO")
  # recruit_thresh = bigtree_recruit_thresh
  # recruit_rate = bigtree_recruit_rate
  # recruit_dist = bigtree_recruit_dist
  # bounded  = FALSE
  # boundary = btnr_sf
  
  sp <- x$sp[1]

  nindiv <- nrow(x) # number of live individuals
  x_repro <- x[x$dbh >= recruit_thresh[[ sp ]],] 
  
  nrepro <- nrow(x_repro) # number of reproductive individuals
    
  # If there are reproductive individuals
  if(nrepro > 0){
      
    # Generate the number of recruits for each reproductive individual
    nrecruit <- rpois(n = nrepro, lambda = recruit_rate[[ sp ]] )
    recruits_sp <- list()
      
    # For each reproductive individual
    for(j in 1:length(nrecruit)){
      
      # If the individual recruits 
      if(nrecruit[j] > 0){

          # For each new recruit, generate new coordinates. 
        
        new_coord <- recruit_coordinates(n = nrecruit[j],
                                         parent_x = x_repro$gx[j],
                                         parent_y = x_repro$gy[j],
                                         recruit_dist = recruit_dist[[ sp ]] )
        
        # Store information on recruits for each reproductive individual
        recruits_sp[[j]] <- data.frame(sp = sp,
                                       gx = new_coord$x,
                                       gy = new_coord$y,
                                       dbh = 10,
                                       id = NA)
      }
        
      }
      recruits  <- do.call("rbind", recruits_sp)
      return(recruits)
  } else {
    # If there are no reproductive individuals, return NULL
    return(NULL)
  }
}


recruit_coordinates <- function(n, parent_x, parent_y, recruit_dist){
  # Basic isotropic dispersal function based on log-normal distribution
  # probability distribution
  # 
  # Arguments:
  #   parent_x: numeric, x-coordinates of parent individual
  #   parent_y: numeric, y-coordinates of parent individual
  #   recruit_dist: numeric, average dispersal distance (equivalent to lambda)
  #
  # Returns:  
  #   list, containing x and y coordinates of recruit
  disp_angle <- runif(min = 0, max = 360, n = n)
  #disp_dist <- rexp(n = n, rate = 1 / recruit_dist ) # exponential function
  disp_dist <- rlnorm(n = n, meanlog = log(recruit_dist), sdlog = 0.5)
  new_x = parent_x + (disp_dist  * cos(disp_angle*pi/180 ))
  new_y = parent_y + (disp_dist * sin(disp_angle*pi/180 ))
  return(data.frame("x" = new_x, "y" = new_y))
}


initialize_BDG_model <- function(init_pop,
                          growth_rate,
                          grow_exceptions = NULL,
                          mort_rate_type = "size",
                          mort_rate,
                          mort_thresh = NULL,
                          mort_exceptions = NULL,
                          recruit_rate,
                          recruit_thresh,
                          recruit_dist,
                          grid_limit = FALSE,
                          grid,
                          boundary,
                          save_ts = FALSE, ##whether to save entire time series or to save only the final time step
                          output_name,  
                          output_dir,
                          print_freq = 10){
  # Run birth-death-recruitment Markov chain simulations
  # 
  # Arguments:
  #     init_pop = data.frame, starting population must contain the following columns: "sp", "dbh", "gx", "gy", "habitat", "id"
  #     growth_rate = list, species-specific, size-dependent growth rates
  #     mort_rate_type = character, either "constant" or "size"
  #     mort_rate = list, species-specific mortality rates
  #     mort_thresh = list, functional-group specific mortality thresholds
  #     recruit_rate = list, species-specific recruitment rates
  #     recruit_thresh = list, species-specific size-thresholds for reproductive maturity
  #     recruit_dist = list, species-specific average dispersal distance
  #     const_pop = logical, whether population number should be constant
  #     pop_cap = numeric, number of stems that simulation is limited to. Value set must be higher than the starting population. When the cap is first exceeded, model culls recruits until population in the next time step is equal to the population cap. In subsequent generations, number of recruits will be forced to equal the number of deaths
  #     boundary = SpatialPolygon object, geographic limit of simulations. Recruitment outside boundary is prohibted
  #     output_dir = character, directory where output is printed
  #     print_freq = numeric, how frequently is output saved
  #
  # Returns:
  #     either list of data.frame objects, where each data.frame contains population at each time step.
  #     OR the final data.frame when the first recruited individual has reached above 300 mm in dbh
  
  pop <- init_pop
  
  if (save_ts==TRUE){
    pop_timeseries <- list()
  }

  i = 1 #counter
  
  recruit_max_dbh=0
  while(recruit_max_dbh<300){
    print(paste0("Begin time step = ", i))
    # let's assume yearly recruitment
    nlivetree <- nrow(pop)
    print(paste0("Number of live trees = ", nlivetree ))
    
  
    # RECRUIT THINGS =====
    print("Simulating recruitment ...")
    
    # Recruitment
    recruits <- ddply(.data = pop,
                      .variables = .(sp),
                      .fun = recruit_things,
                      recruit_thresh = recruit_thresh,
                      recruit_rate = recruit_rate,
                      recruit_dist = recruit_dist,
                      .progress = "text")
    
    # Recruits outside bounds are effectively killed
    if(nrow(recruits) > 0){
      
      recruit_bound_logical <- recruits %>% 
        st_as_sf(coords = c("gx", "gy")) %>%
        st_within(., boundary, sparse = FALSE)

      recruits <- recruits[recruit_bound_logical,]

      nrecruits <- nrow(recruits)

    } else {  
      
      nrecruits <- 0

    }

    #print(paste0("Number of recruits = ", nrecruits))
    
    # KILL THINGS =====
    print("Simulating mortality ...")
    if(mort_rate_type == "size"){
      temp <- ddply(.data = pop,
                    .variables = .(sp),
                    .fun = kill_things, 
                    mort_rate = mort_rate,
                    mort_thresh = mort_thresh,
                    exceptions = mort_exceptions,
                    .progress = "text")
    } else if (mort_rate_type == "constant"){
      temp <- ddply(.data = pop,
                    .variables = .(sp),
                    .fun = kill_things_equal_rate,
                    mort_rate = mort_rate)
      
    }
    survivors <- temp[temp$deaths=="FALSE",]
    survivors$deaths<-NULL
    
    
    ndeaths <- sum(temp$deaths)
    print(paste0("Number of deaths = ", ndeaths ))
    
    
    # GROW THINGS ====================
    print("Simulating growth ...")
    survivors_grown <- ddply(.data = survivors,
                             .variables = .(sp),
                             .fun = grow_things_categ, 
                             growth_rate = growth_rate,
                             exceptions = grow_exceptions,
                             .progress = "text")
    
    # GRID CULLING ====================
    if(grid_limit){

      recruits$gridID <- recruits %>% 
        st_as_sf(coords = c("gx", "gy")) %>%
        st_within(., grid) %>%
        unlist()
      
      # Recruits cannot be in already occupied grid cells
      recruits <- subset(recruits, ! gridID %in% unique(survivors_grown$gridID))      

      pop <- rbind(survivors_grown, recruits)

    } else {

      recruits$gridID<-"NA" #assign NA to gridID otherwise
      pop <- rbind(survivors_grown, recruits)
    }

    print(paste0("Number of recruits = ", nrow(recruits) ))


    ###simulation stops only when a recruit reaches 300 mm  
    # Update recruit_max_dbh
    recruit_max_dbh=max(pop[is.na(pop$id),]$dbh)
    
    if (save_ts==TRUE){
      pop_timeseries[[i]]<-pop
      
    }
    
    i = i + 1 # update counter
  }
  

  if(save_ts == TRUE){
    # Return all time points as a list
    
          saveRDS(pop_timeseries,
                  file = file.path(output_dir, paste0(output_name,"_ts", ".rds" )))
    rm(pop) # remove from memory once finished
    print("Simulation complete")
    
    return(pop)
    
  } else {
    # Return only the last time point
    res <- list("pop" = pop, "time" = i)
    saveRDS(res,
            file = file.path(output_dir, paste0(output_name, "_final",".rds" )))
    print("Simulation complete")
    
    return(pop)
    
  }
}


run_BDG_model <- function(init_pop,
                            growth_rate,
                            grow_exceptions = NULL,
                            mort_rate_type = "size",
                            mort_rate,
                            mort_thresh = NULL,
                            mort_exceptions = NULL,
                            recruit_rate,
                            recruit_thresh,
                            recruit_dist,
                            boundary,
                            grid_limit = FALSE,
                            grid,
                            nyear,
                            print_output = FALSE,
                            output_name,  
                            output_dir,
                            print_freq = NULL){
  # Run birth-death-recruitment Markov chain simulations
  # Alt2 = assumes that information on grid IDs is stored.
  # 
  # Arguments:
  #     init_pop = data.frame, starting population must contain the following columns: "sp", "dbh", "gx", "gy", "id"
  #     growth_rate = list, species-specific, size-dependent growth rates
  #     mort_rate_type = character, either "constant" or "size"
  #     mort_rate = list, species-specific mortality rates
  #     mort_thresh = list, functional-group specific mortality thresholds
  #     recruit_rate = list, species-specific recruitment rates
  #     recruit_thresh = list, species-specific size-thresholds for reproductive maturity
  #     recruit_dist = list, species-specific average dispersal distance
  #     const_pop = logical, whether population number should be constant
  #     pop_cap = numeric, number of stems that simulation is limited to. Value set must be higher than the starting population. When the cap is first exceeded, model culls recruits until population in the next time step is equal to the population cap. In subsequent generations, number of recruits will be forced to equal the number of deaths
  #     boundary = SpatialPolygon object, geographic limit of simulations. Recruitment outside boundary is prohibted
  #     nyear = numeric, number of cycles (years)
  #     print_output = logical, whether output should be saved
  #     output_dir = character, directory where output is printed
  #     print_freq = numeric, how frequently is output saved
  #
  # Returns:
  #     list of data.frame objects. Each data.frame contains population at each time step.
  
  pop_timeseries <- list()
  pop <- init_pop
  
  for(i in 1:nyear){
    print(paste0("Begin time step = ", i))
    # let's assume yearly recruitment
    nlivetree <- nrow(pop)
    print(paste0("Number of live trees = ", nlivetree ))
    
    # RECRUIT THINGS ====================
    print("Simulating recruitment ...")
    
    # Recruitment
    recruits <- ddply(.data = pop,
                      .variables = .(sp),
                      .fun = recruit_things,
                      recruit_thresh = recruit_thresh,
                      recruit_rate = recruit_rate,
                      recruit_dist = recruit_dist,
                      .progress = "text")
    
    # Recruits outside bounds are effectively killed
    if(nrow(recruits) > 0){
      
      recruit_bound_logical <- recruits %>% 
        st_as_sf(coords = c("gx", "gy")) %>%
        st_within(., boundary, sparse = FALSE)

      recruits <- recruits[recruit_bound_logical,]
      nrecruits <- nrow(recruits)
    } else {  
      nrecruits <- 0
    } 

    #print(paste0("Number of recruits = ", nrecruits))

    # KILL THINGS ====================
    print("Simulating mortality ...")
    if(mort_rate_type == "size"){
        ## NEEDS TO BE CHANGED ##
        temp <- ddply(.data = pop,
                    .variables = .(sp),
                    .fun = kill_things, 
                    mort_rate = mort_rate,
                    mort_thresh = mort_thresh,
                    exceptions = mort_exceptions,
                    .progress = "text")      
    } else if (mort_rate_type == "constant"){
      temp <- ddply(.data = pop,
                    .variables = .(sp),
                    .fun = kill_things_equal_rate,
                    mort_rate = mort_rate)
    }
    survivors <- temp[temp$deaths=="FALSE",]
    survivors$deaths<-NULL
    
    
    ndeaths <- sum(temp$deaths)
    print(paste0("Number of deaths = ", ndeaths ))
    
      
    # GROW THINGS ====================
    print("Simulating growth ...")
    survivors_grown <- ddply(.data = survivors,
                             .variables = .(sp),
                             .fun = grow_things_categ, 
                             growth_rate = growth_rate,
                                    exceptions = grow_exceptions,
                                   .progress = "text")

    # This threshold is for initialisation purposes
    # (prevent the ones that are greater than 30 cm dbh from recruiting in the next time step)
    #survivors_grown <- subset(survivors_grown, (dbh <= 300 & is.na(id)) | !is.na(id))

    
    # GRID CULLING ====================
    if(grid_limit){

      recruits$gridID <- recruits %>% 
        st_as_sf(coords = c("gx", "gy")) %>%
        st_within(., grid) %>%
        unlist()
      
      # Recruits cannot be in already occupied grid cells
      recruits <- subset(recruits, ! gridID %in% unique(survivors_grown$gridID))

      pop <- rbind(survivors_grown, recruits)

    } else {
      recruits$gridID<-"NA" #assign NA to gridID otherwise
      
      pop <- rbind(survivors_grown, recruits)
    }

    # SAVING AND EXPORTING RESULTS ====================
    if(print_output == TRUE){
      if( i %% print_freq == 0){
        #print(paste("Saving generation ", i))
        #print(paste("Memory used: ", mem_used())
        saveRDS(pop,
                file = file.path(output_dir, paste0(output_name, "_t", i, ".rds" )))  
      }
    } else {
      if( i %% print_freq == 0){
        pop_timeseries[[ length(pop_timeseries) + 1 ]] <- pop
      } 
      if( is.null(print_freq) ){
        # If output is not printed, save anyway
        pop_timeseries[[i]]  <- pop
      }
      
    }
  }
  
  return(pop_timeseries) 
}
