
# Extract growth from model

gg_df <- function(model){
  
  param <- model$param
  tmax <- model$param$tmax
  tmax <- tmax-1
  # time_cutoff <- tmax-500
  time_cutoff <- tmax/2

  # Extract growth matrix, average over last 50% of save time steps
  df_growth <- model$gg
  gg_ave <- apply(df_growth[time_cutoff:tmax,,], c(2,3), FUN = mean, na.rm = TRUE)
  
  df_gg <- melt(gg_ave)
  
  names(df_gg) <- c("Predator_Species", "Predator_Size_Class", "Growth")
  
  ## make predators and prey factor variables
  df_gg$Predator_Species <- as.factor(df_gg$Predator_Species)
  
  df_size_classes <- model$w # Create new object with list of size classes in log10 g
  df_size_classes <- as.data.frame(df_size_classes) # Convert to a dataframe
  names(df_size_classes) <- c("Predator_Size") # Create new name Size for the variable
  df_size_classes$Predator_Size_Class <- 1:length(model$w) # Create Size_class variable to join with phytodiet by
  
  df_species_predator <- model$param$groups$species # Create new object with species names from model
  df_species_predator <- as.data.frame((df_species_predator)) # Convert to a dataframe
  names(df_species_predator) <- c("Predator") # Rename variable
  df_species_predator$Predator_Species <- 1:length(model$param$groups$species) # Create new integer variable in ascending order with 'Species'
  df_species_predator$Predator_Species <- as.factor(df_species_predator$Predator_Species) # Convert to factor
  
  df_combined <- left_join(df_gg, df_size_classes,
                           by="Predator_Size_Class")
  
  df_combined <- df_combined %>%
    left_join(df_species_predator,
              by="Predator_Species")
  
  df_plot <- df_combined %>%
    mutate(dwdt = Growth*Predator_Size) %>%
    mutate(RGR = Growth/Predator_Size) %>%
    # filter(Time >= 500) %>%
    # group_by(Predator,Predator_Size) %>%
    # summarise(tavg_dwdt = mean(dwdt)) %>%
    filter(dwdt > 0)
  
  return(df_plot)
  
}

# Extract Production:Biomass from models

PB_df <- function(projection){
  
  param <- projection$param
  tmax <- projection$param$tmax
  tmax <- tmax-1
  # time_cutoff <- tmax-500
  time_cutoff <- tmax/2
  w <- projection$w
  dx <- projection$param$dx
  
  # Extract mortality matrix, average over last 50% of save time steps
  df_Z <- projection$Z
  Z_ave <- apply(df_Z[time_cutoff:tmax,,], c(2,3), FUN = mean, na.rm = TRUE)
  
  # Extract abundance matrix, average over last 50% of save time steps
  df_N <- projection$N
  N_ave <- apply(df_N[time_cutoff:tmax,,], c(2,3), FUN = mean, na.rm = TRUE)
  
  Biom_ave <- sweep(N_ave, 2, w, '*') # Calculate biomass in each size class
  Biom_int <- rowSums(Biom_ave*dx) # Biomass integral for each group
  
  # Calculate production integral for each group
  Prod_ave <- rowSums(Z_ave*Biom_ave*dx)
  
  # Production to Biomass ratio
  PB_ave <- Prod_ave/Biom_int
  names(PB_ave) <- projection$param$groups$species
  
  return(PB_ave)
  
}

ppmr_biomass_weighted <-function(logbase, df_RTP, df_B){
  
  df_B_sum <- df_B %>%
    filter(Time >=500) %>%
    group_by(Predator) %>%
    summarise(Biomass = mean(Biomass_sum))
  
  df_B_sum_total <- df_B_sum %>%
    summarise(Total_biomass = sum(Biomass))
  
  df_B_active <- df_B_sum %>%
    filter(Predator == "Active Cephs")
  
  df_B_inactive <- df_B_sum %>%
    filter(Predator == "Inactive Cephs")
  
  df_B_fish <- df_B_sum %>%
    filter(Predator == "Fish")
  
  df_B_zoo <- df_B_sum %>%
    filter(Predator == "Zooplankton")
  
  df_RTP$log10_SC <- log10(df_RTP$Predator_Size)
  
  df_RTP_tidy <- df_RTP %>%
    filter(!is.na(RTP))
  
  df_RTP_Zoo <- df_RTP_tidy %>%
    filter(Predator == "Zooplankton")
  
  df_RTP_Fish <- df_RTP_tidy %>%
    filter(Predator == "Fish")
  
  df_RTP_Active <- df_RTP_tidy %>%
    filter(Predator == "Active Cephs")
  
  df_RTP_Inactive <- df_RTP_tidy %>%
    filter(Predator == "Inactive Cephs")
  
  lm_RTP_Zoo <- lm(RTP~log10_SC, data=df_RTP_Zoo)
  lm_RTP_Fish <- lm(RTP~log10_SC, data=df_RTP_Fish)
  lm_RTP_Active <- lm(RTP~log10_SC, data=df_RTP_Active)
  lm_RTP_Inactive <- lm(RTP~log10_SC, data=df_RTP_Inactive)
  
  b_active <- lm_RTP_Active$coefficients[2]
  b_inactive <- lm_RTP_Inactive$coefficients[2]
  b_fish <- lm_RTP_Fish$coefficients[2]
  b_zoo <- lm_RTP_Zoo$coefficients[2]
  
  ppmr_active <-logbase^(1/b_active)
  ppmr_inactive <-logbase^(1/b_inactive)
  ppmr_fish <-logbase^(1/b_fish)
  ppmr_zoo <-logbase^(1/b_zoo)
  
  weighted_slope <- (ppmr_active*df_B_active$Biomass)+(ppmr_inactive*df_B_inactive$Biomass)+(ppmr_fish*df_B_fish$Biomass)+(ppmr_zoo*df_B_zoo$Biomass)
  
  ppmr_B_weighted <- weighted_slope/df_B_sum_total$Total_biomass
  
  names(ppmr_B_weighted) <- "Biomass_weighted_ppmr"
  
  return(ppmr_B_weighted)
  
}



Diet_df <- function(projection, file_name){
  
  df_diet <- colMeans(projection$diet[(ceiling(0.5*dim(projection$diet)[1])):(dim(projection$diet)[1]),,,,], dim = 1)
  
  # Create long form diet dataframe
  df_diet <- melt(df_diet)
  
  # Assign names to all variables in df_diet
  names(df_diet) <- c("Predator_Species", "Predator_Size_Class", "Prey_Species", "Prey_Size_Class", "Prey_Biomass")
  
  # make predators and prey factor variables
  df_diet$Predator_Species <- as.factor(df_diet$Predator_Species)
  df_diet$Prey_Species <- as.factor(df_diet$Prey_Species)
  
  df_size_classes <- projection$w # Create new object with list of size classes in log10 g
  df_size_classes <- as.data.frame(df_size_classes) # Convert to a dataframe
  names(df_size_classes) <- c("Prey_Size") # Create new name Size for the variable 
  df_size_classes$Prey_Size_Class <- 1:length(projection$w)
  
  df_size_classes_1 <- projection$w # Create new object with list of size classes in log10 g
  df_size_classes_1 <- as.data.frame(df_size_classes_1) # Convert to a dataframe
  names(df_size_classes_1) <- c("Predator_Size") # Create new name Size for the variable 
  df_size_classes_1$Predator_Size_Class <- 1:length(projection$w)
  
  df_species_prey <- projection$param$groups$species # Create new object with species names from model
  df_species_prey <- as.data.frame((df_species_prey)) # Convert to a dataframe
  names(df_species_prey) <- c("Prey") # Rename variable
  df_species_prey$Prey_Species <- 1:length(param$groups$species) # Create new integer variable in ascending order with 'Species'
  df_species_prey$Prey_Species <- as.factor(df_species_prey$Prey_Species) # Convert to factor
  
  df_species_predator <- projection$param$groups$species # Create new object with species names from model
  df_species_predator <- as.data.frame((df_species_predator)) # Convert to a dataframe
  names(df_species_predator) <- c("Predator") # Rename variable
  df_species_predator$Predator_Species <- 1:length(param$groups$species) # Create new integer variable in ascending order with 'Species'
  df_species_predator$Predator_Species <- as.factor(df_species_predator$Predator_Species) # Convert to factor
  
  ## join df_N with df_size_class data in order to create log10 size classes rather then the size     class integer 1:184 that is in df_N
  df_combined <- left_join(df_diet, df_size_classes,
                           by="Prey_Size_Class") %>%
    # then join df_species to df_combined so that there is a variable with actual species names rather than a numeric variable 1:n (in this case 1:7)
    left_join(df_size_classes_1,
              by="Predator_Size_Class") %>%
    
    left_join(df_species_prey,
              by="Prey_Species") %>%
    
    left_join(df_species_predator,
              by="Predator_Species") %>%
    
    filter(Prey_Biomass > 0)
  
  
  
  write_csv(df_combined, paste(file_name, "diet.csv"))
  
}



N_to_Biomass_df <- function(projection){
  
  param <- projection$param  
  tmax <- projection$param$tmax
  #tmax <- 3000
  # cutoff_time <- 2000
  
  df_N <- melt(projection$N) # Extract adundance from model projection and use melt to convert 3D                                  array into long form data
  names(df_N) <- c("Time", "Species_Names", "Size_Class", "N") # Rename variable names
  df_N$Species_Names <- factor(df_N$Species_Names) # Change Species_Names to a factor
  
  df_size_classes <- projection$w # Create new object with list of size classes in log10 g
  df_size_classes <- as.data.frame(df_size_classes) # Convert to a dataframe
  names(df_size_classes) <- c("Size") # Create new name Size for the variable 
  df_size_classes$Size_Class <- 1:length(projection$w) # Create new variable 1:184 in ascending order with log10 g                                             size classes, which will be used to left join later
  bin_width <- diff(df_size_classes$Size)
  start_bin <- 5.011872e-13 - 3.981072e-13
  bin_width_total <- prepend(bin_width, start_bin)
  df_size_classes$bin_width <- bin_width_total
  df_size_classes <- df_size_classes %>%
    mutate(half_bin_width = bin_width_total/2) %>%
    mutate(size_class_midpoint = Size + half_bin_width)
  
  
  df_species <- projection$param$groups$species # Create new object with species names from model
  df_species <- as.data.frame((df_species)) # Convert to a dataframe
  names(df_species) <- c("Species") # Rename variable
  df_species$Species_Names <- 1:length(param$groups$species) # Create new integer variable in ascending order with 'Species'
  df_species$Species_Names <- as.factor(df_species$Species_Names) # Convert to factor
  
  ## join df_N with df_size_class data in order to create log10 size classes rather then the size     class integer 1:184 that is in df_N
  df_combined <- left_join(
    select(df_N, Time, Species_Names, Size_Class, N),
    select(df_size_classes, Size, Size_Class, bin_width, size_class_midpoint),
    by="Size_Class") %>%
    # then join df_species to df_combined so that there is a variable with actual species names rather than a numeric variable 1:n (in this case 1:7)
    left_join(df_species,
              by="Species_Names")
  
  
  df_norm_biomass_tavg  <- df_combined %>%
    group_by(Time,Species, Size) %>%
    #group_by(Species, Size, Time) %>%
    summarise(Mean_N = mean(N)) %>%
    mutate(Biomass = Mean_N*Size) %>%
    filter(Biomass > 0)
  
  df_Species_Biomass <- df_norm_biomass_tavg %>%
    group_by(Species, Time) %>%
    summarise(Annual_Biomass = sum(Biomass))
  
  df_Community_Biomass <- df_norm_biomass_tavg %>%
    group_by(Time) %>%
    summarise(Community_Annual_Biomass = sum(Biomass))
  
  df_Species_Biomass$Species_Order <- factor(df_Species_Biomass$Species, 
                                             levels = c('Active Cephs',
                                                        'Inactive Cephs',
                                                        'Zooplankton',
                                                        'Fish'),
                                             ordered = TRUE)
  
  return(df_Species_Biomass)
  
}


Abundance_df <- function(projection){
  
  model <- projection
  
  param <- projection$param  
  tmax <- projection$param$tmax
  #tmax <- 3000
  # cutoff_time <- 2000
  
  tmax <- model$param$tmax
  tmax <- tmax-1
  time_cutoff <- tmax/2
  w <- model$w
  
  ## Extract abundance matrix, average over last 50% of save time steps
  df_N <- model$N
  N_ave <- apply(df_N[time_cutoff:tmax,,], c(2,3), FUN = mean, na.rm = TRUE)
  
  N_df <- melt(N_ave)
  
  names(N_df) <- c("Predator_Species", "Predator_Size_Class", "Predator_Abundance")
  
  N_df$Predator_Species <- as.factor(N_df$Predator_Species)
  
  df_size_classes <- model$w # Create new object with list of size classes in log10 g
  df_size_classes <- as.data.frame(df_size_classes) # Convert to a dataframe
  names(df_size_classes) <- c("Predator_Size") # Create new name Size for the variable
  df_size_classes$Predator_Size_Class <- 1:length(model$w)
  
  df_species_predator <- model$param$groups$species # Create new object with species names from model
  df_species_predator <- as.data.frame((df_species_predator)) # Convert to a dataframe
  names(df_species_predator) <- c("Predator") # Rename variable
  df_species_predator$Predator_Species <- 1:length(model$param$groups$species) # Create new integer variable in ascending order with 'Species'
  df_species_predator$Predator_Species <- as.factor(df_species_predator$Predator_Species) # Convert to factor
  
  ## join df_N with df_size_class data in order to create log10 size classes rather then the size     class integer 1:184 that is in df_N
  df_combined <- N_df %>%
    left_join(df_species_predator,
              by="Predator_Species")
  
  df_combined <- left_join(df_combined, df_size_classes,
                           by="Predator_Size_Class")
  
  
  df_combined$Species_Order <- factor(df_combined$Predator, 
                                               levels = c('Zooplankton',
                                                          'Active Cephs',
                                                          'Inactive Cephs',
                                                          'Fish'),
                                               ordered = TRUE)
  
  return(df_combined)
  
}

Biomass_df <- function(projection){
  
  model <- projection
  
  dx <- model$dx
  
  param <- projection$param  
  tmax <- projection$param$tmax
  #tmax <- 3000
  # cutoff_time <- 2000
  
  tmax <- model$param$tmax
  tmax <- tmax-1
  time_cutoff <- tmax/2
  w <- model$w
  
  ## Extract abundance matrix, average over last 50% of save time steps
  df_N <- model$N
  
  N_df <- melt(df_N)
  
  names(N_df) <- c("Time", "Predator_Species", "Predator_Size_Class", "Predator_Abundance")
  
  N_df$Predator_Species <- as.factor(N_df$Predator_Species)
  
  df_size_classes <- model$w # Create new object with list of size classes in log10 g
  df_size_classes <- as.data.frame(df_size_classes) # Convert to a dataframe
  names(df_size_classes) <- c("Predator_Size") # Create new name Size for the variable
  df_size_classes$Predator_Size_Class <- 1:length(model$w)
  
  df_species_predator <- model$param$groups$species # Create new object with species names from model
  df_species_predator <- as.data.frame((df_species_predator)) # Convert to a dataframe
  names(df_species_predator) <- c("Predator") # Rename variable
  df_species_predator$Predator_Species <- 1:length(model$param$groups$species) # Create new integer variable in ascending order with 'Species'
  df_species_predator$Predator_Species <- as.factor(df_species_predator$Predator_Species) # Convert to factor
  
  ## join df_N with df_size_class data in order to create log10 size classes rather then the size     class integer 1:184 that is in df_N
  df_combined <- N_df %>%
    left_join(df_species_predator,
              by="Predator_Species")
  
  df_combined <- left_join(df_combined, df_size_classes,
                           by="Predator_Size_Class")
  
  df_biomass  <- df_combined %>%
    mutate(Biomass = Predator_Abundance*Predator_Size) %>%
    filter(Biomass > 0)
  
  df_biomass$Predator <- as.factor(df_biomass$Predator)
  
  df_biomass_int  <- df_biomass %>%
    group_by(Time, Predator) %>%
    summarise(Biomass_sum = sum(Biomass*dx))
  
  
  return(df_biomass_int)
  
}

Compute_CV <- function(df, cutoff_time){
  
  tmax <- length(df$Annual_Biomass)/4
  
  names(df) <- c("Species", "Year", "Annual_Biomass","Species Order")
  
  df <- df %>%
    filter(Year >= cutoff_time, Year <= tmax)
  
  df_community <- df %>%
    group_by(Year) %>%
    summarise(Annual_Biomass = sum(Annual_Biomass))
  
  df_active <- df %>%
    filter(Species=="Active Cephs") %>%
    group_by(Year) %>%
    summarise(Annual_Biomass = sum(Annual_Biomass))
  
  df_inactive <- df %>%
    filter(Species=="Inactive Cephs") %>%
    group_by(Year) %>%
    summarise(Annual_Biomass = sum(Annual_Biomass))
  
  df_zoo <- df %>%
    filter(Species=="Zooplankton") %>%
    group_by(Year) %>%
    summarise(Annual_Biomass = sum(Annual_Biomass))
  
  df_fish <- df %>%
    filter(Species=="Fish") %>%
    group_by(Year) %>%
    summarise(Annual_Biomass = sum(Annual_Biomass))
  
  CV_community <- signif(cv(df_community$Annual_Biomass), 6)
  CV_Active    <- signif(cv(df_active$Annual_Biomass), 6)
  CV_Inactive  <- signif(cv(df_inactive$Annual_Biomass),6)
  CV_Zoo       <- signif(cv(df_zoo$Annual_Biomass),6)
  CV_Fish      <- signif(cv(df_fish$Annual_Biomass),6)
  
  # CV_list <- list(CV_community,CV_Active, CV_Inactive, CV_Zoo, CV_Fish)
  
  # return(CV_list)
  
  
  print(paste("Active Cephs CV =",CV_Active))
  print(paste("Inactive Cephs CV =",CV_Inactive))
  print(paste("Fish CV =",CV_Fish))
  print(paste("Zooplankton CV =",CV_Zoo))
  print(paste("Community CV =",CV_community))
  
}




Emergent_PPMR_Allometry <- function(model){
 
  # model <- Model
  # Average diet array for each predator size class for the final half of the projection
  dynamdiet = apply(model$dynam_diet_full[c(floor(0.5*dim(model$dynam_diet_full)[1]):dim(model$dynam_diet_full)[1]),,,,], c(2,3,4,5), mean)
  phytodiet = apply(model$phyto_diet_full[c(floor(0.5*dim(model$phyto_diet_full)[1]):dim(model$phyto_diet_full)[1]),,,], c(2,3,4), mean)
  
  tmax <- model$param$tmax
  tmax <- tmax-1
  time_cutoff <- tmax/2
  w <- model$w
  
  ## Extract abundance matrix, average over last 50% of save time steps
  df_N <- model$N
  N_ave <- apply(df_N[time_cutoff:tmax,,], c(2,3), FUN = mean, na.rm = TRUE)
  
  N_df <- melt(N_ave)
  
  names(N_df) <- c("Predator_Species", "Predator_Size_Class", "Predator_Abundance")
  
  N_df$Predator_Species <- as.factor(N_df$Predator_Species)
  
  
  curr_phyto_diet <- phytodiet # Current phyto diet
  curr_dynam_diet <- dynamdiet # Current heterotroph diet
  
  start_dynam_tl <- matrix(2, nrow = dim(dynamdiet)[1], ncol = dim(dynamdiet)[2]) # Start TL for dynamic size groups - start at 2 for all groups and sizes
  
  total_diet <- apply(curr_phyto_diet, c(1,2), sum) + apply(curr_dynam_diet, c(1,2), sum) # Total consumption, in grams wet weight, by pred group and pred size
  
  curr_phyto_frac <- apply(sweep(curr_phyto_diet, c(1,2), total_diet, '/'), c(1,2), sum) # Fraction of diet from phyto, by pred group and pred sizes
  curr_dynam_frac <- sweep(curr_dynam_diet, c(1,2), total_diet, '/') # Fraction of diet from each prey group and prey size, for each pred group and pred size
  
  curr_dynam_frac_melt <- melt(curr_dynam_frac)
  curr_phyto_frac_melt <- melt(curr_phyto_frac)
  
  # Assign names to all variables in df_diet
  names(curr_phyto_frac_melt) <- c("Predator_Species", "Predator_Size_Class", "Prey_Frac")
  
  # make predators and prey factor variables
  curr_phyto_frac_melt$Predator_Species <- as.factor(curr_phyto_frac_melt$Predator_Species)
  
  df_size_classes <- model$w # Create new object with list of size classes in log10 g
  df_size_classes <- as.data.frame(df_size_classes) # Convert to a dataframe
  names(df_size_classes) <- c("Predator_Size") # Create new name Size for the variable
  df_size_classes$Predator_Size_Class <- 1:length(model$w)
  
  df_species_predator <- model$param$groups$species # Create new object with species names from model
  df_species_predator <- as.data.frame((df_species_predator)) # Convert to a dataframe
  names(df_species_predator) <- c("Predator") # Rename variable
  df_species_predator$Predator_Species <- 1:length(model$param$groups$species) # Create new integer variable in ascending order with 'Species'
  df_species_predator$Predator_Species <- as.factor(df_species_predator$Predator_Species) # Convert to factor
  
  # df_phyto_size_classes <- model$w_phyto # Create new object with list of size classes in log10 g
  # df_phyto_size_classes <- as.data.frame(df_phyto_size_classes) # Convert to a dataframe
  # names(df_phyto_size_classes) <- c("Prey_Size") # Create new name Size for the variable
  # df_phyto_size_classes$Prey_Size_Class <- 1:length(model$w_phyto)
  # 
  ## join df_N with df_size_class data in order to create log10 size classes rather then the size     class integer 1:184 that is in df_N
  df_combined_frac <- left_join(curr_phyto_frac_melt, df_size_classes,
                                by="Predator_Size_Class") %>%
    
    left_join(df_species_predator,
              by="Predator_Species") # %>%
  # left_join(df_phyto_size_classes,
  #           by="Prey_Size_Class")
  
  df_combined_frac$Prey <- c("Phyto")
  
  ### Dynamic diet dataframe
  
  
  # names(dynamdiet)
  
  # Assign names to all variables in df_diet
  names(curr_dynam_frac_melt) <- c("Predator_Species", "Predator_Size_Class", "Prey_Species", "Prey_Size_Class", "Prey_Frac")
  
  # make predators and prey factor variables
  curr_dynam_frac_melt$Predator_Species <- as.factor(curr_dynam_frac_melt$Predator_Species)
  curr_dynam_frac_melt$Prey_Species <- as.factor(curr_dynam_frac_melt$Prey_Species)
  
  # df_size_classes <- model$w # Create new object with list of size classes in log10 g
  # df_size_classes <- as.data.frame(df_size_classes) # Convert to a dataframe
  # names(df_size_classes) <- c("Predator_Size") # Create new name Size for the variable
  # df_size_classes$Predator_Size_Class <- 1:length(model$w)
  
  df_species_prey <- model$param$groups$species # Create new object with species names from model
  df_species_prey <- as.data.frame((df_species_prey)) # Convert to a dataframe
  names(df_species_prey) <- c("Prey") # Rename variable
  df_species_prey$Prey_Species <- 1:length(model$param$groups$species) # Create new integer variable in ascending order with 'Species'
  df_species_prey$Prey_Species <- as.factor(df_species_prey$Prey_Species) # Convert to factor
  
  df_prey_size_classes <- model$w # Create new object with list of size classes in log10 g
  df_prey_size_classes <- as.data.frame(df_prey_size_classes) # Convert to a dataframe
  names(df_prey_size_classes) <- c("Prey_Size") # Create new name Size for the variable
  df_prey_size_classes$Prey_Size_Class <- 1:length(model$w)
  
  ## join df_N with df_size_class data in order to create log10 size classes rather then the size     class integer 1:184 that is in df_N
  df_combined_dynam_frac <- left_join(curr_dynam_frac_melt, df_size_classes,
                                      by="Predator_Size_Class") %>%
    
    left_join(df_species_predator,
              by="Predator_Species") %>%
    left_join(df_prey_size_classes,
              by="Prey_Size_Class") %>%
    left_join(df_species_prey,
              by="Prey_Species")
  
  df_comb_all <- df_combined_dynam_frac %>%
    group_by(Predator,Predator_Size, Prey) %>%
    summarise(Prey_Frac = sum(Prey_Frac)) 
  
  
  df_phyto_frac_v1 <- df_combined_frac %>%
    select(Predator, Prey, Predator_Size, Prey_Frac)
  
  df_test_v3 <- bind_rows(df_comb_all,df_phyto_frac_v1)
  
  df_test_v3$Prey <- as.factor(df_test_v3$Prey)
  
  df_test_v3 <- df_test_v3[Reduce(`&`, lapply(df_test_v3, function(x) !is.na(x)  & is.finite(x))),]
  
  df_cephs_prey_frac <- df_test_v3 %>%
    filter(Prey == "Active Cephs"|Prey == "Inactive Cephs") %>%
    group_by(Predator, Predator_Size) %>%
    summarise(Prey_Frac = sum(Prey_Frac))
  
  df_cephs_prey_frac$Prey <- "Cephs"
  
  df_rest_prey <- df_test_v3 %>%
    filter(Prey == "Fish"|Prey == "Zooplankton"|Prey == "Phyto")
  
  df_ceph_prey_comb <- rbind(df_cephs_prey_frac, df_rest_prey)
  
  df_ceph_prey_comb$Prey <- as.factor(df_ceph_prey_comb$Prey)
  
  
  #####
  
  df_d_frac <- df_comb_all %>%
    select(Predator, Prey, Predator_Size, Prey_Frac)
  df_pd_frac <- df_combined_frac %>%
    select(Predator, Prey, Predator_Size, Prey_Frac)
  
  df_d_frac <- df_d_frac[Reduce(`&`, lapply(df_d_frac, function(x) !is.na(x)  & is.finite(x))),]
  
  df_pd_frac <- df_pd_frac[complete.cases(df_pd_frac),]
  
  df_frac <- bind_rows(df_d_frac, df_pd_frac) # join data frames, one on top of the other
  
  # ggplot(df_frac, aes(x=log10(Predator_Size), y = Prey_Frac, fill = Prey)) +
  #   geom_area() +
  #   facet_wrap(~Predator, scales = "free_x") +
  #   scale_fill_manual(name = 'Prey',
  #                     values =c("#F8766D","#7CAE00", "#00BFC4","green2","#C77CFF"),
  #                     labels = c('Active', 'Fish','Inactive', "Phytoplankton",'Zooplankton')) +
  #   theme_classic() +
  #   xlab(expression(paste("Predator ",Size~(10^X~g)))) +
  #   ylab("Prey Proportion")
  

  
  phytodiet <- melt(phytodiet)
  # names(phytodiet)
  
  # Assign names to all variables in df_diet
  names(phytodiet) <- c("Predator_Species", "Predator_Size_Class", "Prey_Size_Class", "Prey_Abundance")
  
  # make predators and prey factor variables
  phytodiet$Predator_Species <- as.factor(phytodiet$Predator_Species)
  
  df_size_classes <- model$w # Create new object with list of size classes in log10 g
  df_size_classes <- as.data.frame(df_size_classes) # Convert to a dataframe
  names(df_size_classes) <- c("Predator_Size") # Create new name Size for the variable
  df_size_classes$Predator_Size_Class <- 1:length(model$w)
  
  df_species_predator <- model$param$groups$species # Create new object with species names from model
  df_species_predator <- as.data.frame((df_species_predator)) # Convert to a dataframe
  names(df_species_predator) <- c("Predator") # Rename variable
  df_species_predator$Predator_Species <- 1:length(model$param$groups$species) # Create new integer variable in ascending order with 'Species'
  df_species_predator$Predator_Species <- as.factor(df_species_predator$Predator_Species) # Convert to factor
  
  df_phyto_size_classes <- model$w_phyto # Create new object with list of size classes in log10 g
  df_phyto_size_classes <- as.data.frame(df_phyto_size_classes) # Convert to a dataframe
  names(df_phyto_size_classes) <- c("Prey_Size") # Create new name Size for the variable
  df_phyto_size_classes$Prey_Size_Class <- 1:length(model$w_phyto)
  
  ## join df_N with df_size_class data in order to create log10 size classes rather then the size     class integer 1:184 that is in df_N
  df_combined <- left_join(phytodiet, df_size_classes,
                           by="Predator_Size_Class") %>%
    
    left_join(df_species_predator,
              by="Predator_Species") %>%
    left_join(df_phyto_size_classes,
              by="Prey_Size_Class")
  
  df_combined <- df_combined %>%
    mutate(Prey_Biomass = Prey_Abundance*Prey_Size)
  
  df_combined$Prey <- c("Phyto")
  
  ### Dynamic diet dataframe
  
  dynamdiet <- melt(dynamdiet)
  
  # names(dynamdiet)
  
  # Assign names to all variables in df_diet
  names(dynamdiet) <- c("Predator_Species", "Predator_Size_Class", "Prey_Species", "Prey_Size_Class", "Prey_Abundance")
  
  # make predators and prey factor variables
  dynamdiet$Predator_Species <- as.factor(dynamdiet$Predator_Species)
  dynamdiet$Prey_Species <- as.factor(dynamdiet$Prey_Species)
  
  # df_size_classes <- model$w # Create new object with list of size classes in log10 g
  # df_size_classes <- as.data.frame(df_size_classes) # Convert to a dataframe
  # names(df_size_classes) <- c("Predator_Size") # Create new name Size for the variable
  # df_size_classes$Predator_Size_Class <- 1:length(model$w)
  
  df_species_prey <- model$param$groups$species # Create new object with species names from model
  df_species_prey <- as.data.frame((df_species_prey)) # Convert to a dataframe
  names(df_species_prey) <- c("Prey") # Rename variable
  df_species_prey$Prey_Species <- 1:length(model$param$groups$species) # Create new integer variable in ascending order with 'Species'
  df_species_prey$Prey_Species <- as.factor(df_species_prey$Prey_Species) # Convert to factor
  
  df_prey_size_classes <- model$w # Create new object with list of size classes in log10 g
  df_prey_size_classes <- as.data.frame(df_prey_size_classes) # Convert to a dataframe
  names(df_prey_size_classes) <- c("Prey_Size") # Create new name Size for the variable
  df_prey_size_classes$Prey_Size_Class <- 1:length(model$w)
  
  ## join df_N with df_size_class data in order to create log10 size classes rather then the size     class integer 1:184 that is in df_N
  df_combined_dynam <- left_join(dynamdiet, df_size_classes,
                           by="Predator_Size_Class") %>%
    
    left_join(df_species_predator,
              by="Predator_Species") %>%
    left_join(df_prey_size_classes,
              by="Prey_Size_Class") %>%
    left_join(df_species_prey,
              by="Prey_Species")
  
  
  ###
  
  df_d <- df_combined_dynam %>%
    select(Predator, Prey, Predator_Size, Prey_Size, Prey_Abundance)
  df_pd <- df_combined %>%
    select(Predator, Prey, Predator_Size, Prey_Size, Prey_Abundance)
  
  df <- rbind(df_d, df_pd) # join data frames, one on top of the other
  
  ## Combine abundance df with diet df

  df_ind_PPMR_prelim <- N_df %>%
    left_join(df_species_predator,
              by="Predator_Species")
  
  df_ind_PPMR <- left_join(df_ind_PPMR_prelim, df_size_classes,
                           by="Predator_Size_Class")
  
  df_ind_PPMR <- df_ind_PPMR %>%
    select(Predator, Predator_Size, Predator_Abundance)
  
  df_ind_PPMR_plot <- left_join(df, df_ind_PPMR)
  
  df_final <- df_ind_PPMR_plot %>%
    mutate(per_capita_diet_N = ((Prey_Abundance/Predator_Abundance)/Prey_Size))
  
  df_final_clean <- df_final[Reduce(`&`, lapply(df_final, function(x) !is.na(x)  & is.finite(x))),]
  
  # # reorder factors
  # df$Predator <- factor(df$Predator, 
  #                       levels = c("Zooplankton", "Fish", "Inactive Cephs", "Active Cephs"))
  # 
  # df$Prey <- factor(df$Prey, 
  #                   levels = c("Phyto","Zooplankton", "Fish", "Inactive Cephs", "Active Cephs"))
  # 
  # unique_sizes <- sort(unique(df$Prey_Size))
  # df_sizes <- tibble(
  #   Prey_Size = unique_sizes, 
  #   Prey_size_class = 1:length(unique_sizes)
  # )
  # 
  # df <- left_join(df, df_sizes, by = "Prey_Size") %>%
  #   arrange(Predator, Prey, Predator_Size, Prey_Size)
  # 
  # df_sizes <- tibble(
  #   Predator_Size = unique_sizes, 
  #   Predator_size_class = 1:length(unique_sizes)
  # )
  # 
  # df <- left_join(df, df_sizes, by = "Predator_Size") %>%
  #   arrange(Predator, Prey, Predator_Size, Prey_Size)
  # 
  # # rm(df_d, df_pd)
  # 
  # df_Pred_Diet_N_sum <- df %>%
  #   group_by(Predator, Predator_Size, Prey_Size) %>%
  #   summarise(Total_Prey_Abundance = sum(Prey_Abundance)) %>%
  #   filter(Total_Prey_Abundance > 0)
  # 
  # # reorder factors
  # df$Predator <- factor(df$Predator, 
  #                       levels = c("Zooplankton", "Fish", "Inactive Cephs", "Active Cephs"))
  # 
  # df$Prey <- factor(df$Prey, 
  #                   levels = c("Phyto", "Zooplankton", "Fish", "Inactive Cephs", "Active Cephs"))
  
  
  # create a data frame for plotting
  df_plot_comb <- df_final_clean %>%
    mutate(PySzPyAb = Prey_Size*per_capita_diet_N) %>%  # total size
    group_by(Predator, Predator_Size) %>% # groups of prey distributions
    summarise(
      sum_prey_abundance = sum(per_capita_diet_N), # total prey numbers
      sum_PySzPyAb       = sum(PySzPyAb), # total size
      mean_prey_size     = sum_PySzPyAb / sum_prey_abundance # mean size
    )
  
  # data for dashed line where expected prey size is 10% of predator size
  # df_ref <- tibble(
  #   Predator_Size = c(min(df_Pred_Diet_N_sum$Predator_Size), max(df_Pred_Diet_N_sum$Predator_Size))
  # )
  # df_ref$Expected_Prey_Size_PPMR100000 = 0.00001*df_ref$Predator_Size
  # df_ref$Expected_Prey_Size_PPMR10000 = 0.0001*df_ref$Predator_Size
  # df_ref$Expected_Prey_Size_PPMR1000 = 0.001*df_ref$Predator_Size
  # df_ref$Expected_Prey_Size_PPMR100 = 0.01*df_ref$Predator_Size
  # df_ref$Expected_Prey_Size_PPMR10 = 0.1*df_ref$Predator_Size


  # ggplot() +
  #   geom_line(data = df_ref, aes(x = Predator_Size, y = Expected_Prey_Size_PPMR10),
  #             linetype = "solid") +
  #   geom_line(data = df_ref, aes(x = Predator_Size, y = Expected_Prey_Size_PPMR100),
  #             linetype = "solid") +
  #   geom_line(data = df_ref, aes(x = Predator_Size, y = Expected_Prey_Size_PPMR1000),
  #             linetype = "longdash") +
  #   geom_line(data = df_ref, aes(x = Predator_Size, y = Expected_Prey_Size_PPMR10000),
  #             linetype = "dashed") +
  #   geom_line(data = df_ref, aes(x = Predator_Size, y = Expected_Prey_Size_PPMR100000),
  #             linetype = "dotted") +
  #   geom_line(data = df_plot_comb,
  #             aes(x = Predator_Size, y = mean_prey_size, color = Predator), size = 1) +
  #   scale_x_log10() +
  #   scale_y_log10() +
  #   # facet_wrap( ~ Predator) +
  #   labs(x = "Predator size", y = "Mean prey size") +
  #   theme_bw()
  
  
  return(list(df_plot_comb, df_frac, df_ceph_prey_comb))
  
}


Prey_Proportion_Plot_Cephs_Separate <- function(model){


# Average diet array for each predator size class for the final half of the projection
dynamdiet = apply(model$dynam_diet_full[c(floor(0.5*dim(model$dynam_diet_full)[1]):dim(model$dynam_diet_full)[1]),,,,], c(2,3,4,5), mean)
phytodiet = apply(model$phyto_diet_full[c(floor(0.5*dim(model$phyto_diet_full)[1]):dim(model$phyto_diet_full)[1]),,,], c(2,3,4), mean)

curr_phyto_diet <- phytodiet # Current phyto diet
curr_dynam_diet <- dynamdiet # Current heterotroph diet

start_dynam_tl <- matrix(2, nrow = dim(dynamdiet)[1], ncol = dim(dynamdiet)[2]) # Start TL for dynamic size groups - start at 2 for all groups and sizes

total_diet <- apply(curr_phyto_diet, c(1,2), sum) + apply(curr_dynam_diet, c(1,2), sum) # Total consumption, in grams wet weight, by pred group and pred size

curr_phyto_frac <- apply(sweep(curr_phyto_diet, c(1,2), total_diet, '/'), c(1,2), sum) # Fraction of diet from phyto, by pred group and pred sizes
curr_dynam_frac <- sweep(curr_dynam_diet, c(1,2), total_diet, '/') # Fraction of diet from each prey group and prey size, for each pred group and pred size

curr_dynam_frac_melt <- melt(curr_dynam_frac)
curr_phyto_frac_melt <- melt(curr_phyto_frac)

# Assign names to all variables in df_diet
names(curr_phyto_frac_melt) <- c("Predator_Species", "Predator_Size_Class", "Prey_Frac")

# make predators and prey factor variables
curr_phyto_frac_melt$Predator_Species <- as.factor(curr_phyto_frac_melt$Predator_Species)

df_size_classes <- model$w # Create new object with list of size classes in log10 g
df_size_classes <- as.data.frame(df_size_classes) # Convert to a dataframe
names(df_size_classes) <- c("Predator_Size") # Create new name Size for the variable
df_size_classes$Predator_Size_Class <- 1:length(model$w)

df_species_predator <- model$param$groups$species # Create new object with species names from model
df_species_predator <- as.data.frame((df_species_predator)) # Convert to a dataframe
names(df_species_predator) <- c("Predator") # Rename variable
df_species_predator$Predator_Species <- 1:length(model$param$groups$species) # Create new integer variable in ascending order with 'Species'
df_species_predator$Predator_Species <- as.factor(df_species_predator$Predator_Species) # Convert to factor

# df_phyto_size_classes <- model$w_phyto # Create new object with list of size classes in log10 g
# df_phyto_size_classes <- as.data.frame(df_phyto_size_classes) # Convert to a dataframe
# names(df_phyto_size_classes) <- c("Prey_Size") # Create new name Size for the variable
# df_phyto_size_classes$Prey_Size_Class <- 1:length(model$w_phyto)
# 
## join df_N with df_size_class data in order to create log10 size classes rather then the size     class integer 1:184 that is in df_N
df_combined_frac <- left_join(curr_phyto_frac_melt, df_size_classes,
                              by="Predator_Size_Class") %>%
  
  left_join(df_species_predator,
            by="Predator_Species") # %>%
# left_join(df_phyto_size_classes,
#           by="Prey_Size_Class")

df_combined_frac$Prey <- c("Phyto")

### Dynamic diet dataframe


# names(dynamdiet)

# Assign names to all variables in df_diet
names(curr_dynam_frac_melt) <- c("Predator_Species", "Predator_Size_Class", "Prey_Species", "Prey_Size_Class", "Prey_Frac")

# make predators and prey factor variables
curr_dynam_frac_melt$Predator_Species <- as.factor(curr_dynam_frac_melt$Predator_Species)
curr_dynam_frac_melt$Prey_Species <- as.factor(curr_dynam_frac_melt$Prey_Species)

# df_size_classes <- model$w # Create new object with list of size classes in log10 g
# df_size_classes <- as.data.frame(df_size_classes) # Convert to a dataframe
# names(df_size_classes) <- c("Predator_Size") # Create new name Size for the variable
# df_size_classes$Predator_Size_Class <- 1:length(model$w)

df_species_prey <- model$param$groups$species # Create new object with species names from model
df_species_prey <- as.data.frame((df_species_prey)) # Convert to a dataframe
names(df_species_prey) <- c("Prey") # Rename variable
df_species_prey$Prey_Species <- 1:length(model$param$groups$species) # Create new integer variable in ascending order with 'Species'
df_species_prey$Prey_Species <- as.factor(df_species_prey$Prey_Species) # Convert to factor

df_prey_size_classes <- model$w # Create new object with list of size classes in log10 g
df_prey_size_classes <- as.data.frame(df_prey_size_classes) # Convert to a dataframe
names(df_prey_size_classes) <- c("Prey_Size") # Create new name Size for the variable
df_prey_size_classes$Prey_Size_Class <- 1:length(model$w)

## join df_N with df_size_class data in order to create log10 size classes rather then the size     class integer 1:184 that is in df_N
df_combined_dynam_frac <- left_join(curr_dynam_frac_melt, df_size_classes,
                                    by="Predator_Size_Class") %>%
  
  left_join(df_species_predator,
            by="Predator_Species") %>%
  left_join(df_prey_size_classes,
            by="Prey_Size_Class") %>%
  left_join(df_species_prey,
            by="Prey_Species")

df_comb_all <- df_combined_dynam_frac %>%
  group_by(Predator,Predator_Size, Prey) %>%
  summarise(Prey_Frac = sum(Prey_Frac)) 
  
  #####

df_d_frac <- df_comb_all %>%
  select(Predator, Prey, Predator_Size, Prey_Frac)
df_pd_frac <- df_combined_frac %>%
  select(Predator, Prey, Predator_Size, Prey_Frac)

df_d_frac <- df_d_frac[Reduce(`&`, lapply(df_d_frac, function(x) !is.na(x)  & is.finite(x))),]

df_pd_frac <- df_pd_frac[complete.cases(df_pd_frac),]

df_frac <- bind_rows(df_d_frac, df_pd_frac) # join data frames, one on top of the other

# ggplot(df_frac, aes(x=log10(Predator_Size), y = Prey_Frac, fill = Prey)) +
#   geom_area() +
#   facet_wrap(~Predator, scales = "free_x") +
#   scale_fill_manual(name = 'Prey',
#                     values =c("#F8766D","#7CAE00", "#00BFC4","green2","#C77CFF"),
#                     labels = c('Active', 'Fish','Inactive', "Phytoplankton",'Zooplankton')) +
#   theme_classic() +
#   xlab(expression(paste("Predator ",Size~(10^X~g)))) +
#   ylab("Prey Proportion")

return(df_frac)

}

Prey_Proportion_Plot_Cephs_Combined <- function(model){
  
  # Average diet array for each predator size class for the final half of the projection
  dynamdiet = apply(model$dynam_diet_full[c(floor(0.5*dim(model$dynam_diet_full)[1]):dim(model$dynam_diet_full)[1]),,,,], c(2,3,4,5), mean)
  phytodiet = apply(model$phyto_diet_full[c(floor(0.5*dim(model$phyto_diet_full)[1]):dim(model$phyto_diet_full)[1]),,,], c(2,3,4), mean)
  
  curr_phyto_diet <- phytodiet # Current phyto diet
  curr_dynam_diet <- dynamdiet # Current heterotroph diet
  
  start_dynam_tl <- matrix(2, nrow = dim(dynamdiet)[1], ncol = dim(dynamdiet)[2]) # Start TL for dynamic size groups - start at 2 for all groups and sizes
  
  total_diet <- apply(curr_phyto_diet, c(1,2), sum) + apply(curr_dynam_diet, c(1,2), sum) # Total consumption, in grams wet weight, by pred group and pred size
  
  curr_phyto_frac <- apply(sweep(curr_phyto_diet, c(1,2), total_diet, '/'), c(1,2), sum) # Fraction of diet from phyto, by pred group and pred sizes
  curr_dynam_frac <- sweep(curr_dynam_diet, c(1,2), total_diet, '/') # Fraction of diet from each prey group and prey size, for each pred group and pred size
  
  curr_dynam_frac_melt <- melt(curr_dynam_frac)
  curr_phyto_frac_melt <- melt(curr_phyto_frac)
  
  # Assign names to all variables in df_diet
  names(curr_phyto_frac_melt) <- c("Predator_Species", "Predator_Size_Class", "Prey_Frac")
  
  # make predators and prey factor variables
  curr_phyto_frac_melt$Predator_Species <- as.factor(curr_phyto_frac_melt$Predator_Species)
  
  df_size_classes <- model$w # Create new object with list of size classes in log10 g
  df_size_classes <- as.data.frame(df_size_classes) # Convert to a dataframe
  names(df_size_classes) <- c("Predator_Size") # Create new name Size for the variable
  df_size_classes$Predator_Size_Class <- 1:length(model$w)
  
  df_species_predator <- model$param$groups$species # Create new object with species names from model
  df_species_predator <- as.data.frame((df_species_predator)) # Convert to a dataframe
  names(df_species_predator) <- c("Predator") # Rename variable
  df_species_predator$Predator_Species <- 1:length(model$param$groups$species) # Create new integer variable in ascending order with 'Species'
  df_species_predator$Predator_Species <- as.factor(df_species_predator$Predator_Species) # Convert to factor
  
  # df_phyto_size_classes <- model$w_phyto # Create new object with list of size classes in log10 g
  # df_phyto_size_classes <- as.data.frame(df_phyto_size_classes) # Convert to a dataframe
  # names(df_phyto_size_classes) <- c("Prey_Size") # Create new name Size for the variable
  # df_phyto_size_classes$Prey_Size_Class <- 1:length(model$w_phyto)
  # 
  ## join df_N with df_size_class data in order to create log10 size classes rather then the size     class integer 1:184 that is in df_N
  df_combined_frac <- left_join(curr_phyto_frac_melt, df_size_classes,
                                by="Predator_Size_Class") %>%
    
    left_join(df_species_predator,
              by="Predator_Species") # %>%
  # left_join(df_phyto_size_classes,
  #           by="Prey_Size_Class")
  
  df_combined_frac$Prey <- c("Phyto")
  
  ### Dynamic diet dataframe
  
  
  # names(dynamdiet)
  
  # Assign names to all variables in df_diet
  names(curr_dynam_frac_melt) <- c("Predator_Species", "Predator_Size_Class", "Prey_Species", "Prey_Size_Class", "Prey_Frac")
  
  # make predators and prey factor variables
  curr_dynam_frac_melt$Predator_Species <- as.factor(curr_dynam_frac_melt$Predator_Species)
  curr_dynam_frac_melt$Prey_Species <- as.factor(curr_dynam_frac_melt$Prey_Species)
  
  # df_size_classes <- model$w # Create new object with list of size classes in log10 g
  # df_size_classes <- as.data.frame(df_size_classes) # Convert to a dataframe
  # names(df_size_classes) <- c("Predator_Size") # Create new name Size for the variable
  # df_size_classes$Predator_Size_Class <- 1:length(model$w)
  
  df_species_prey <- model$param$groups$species # Create new object with species names from model
  df_species_prey <- as.data.frame((df_species_prey)) # Convert to a dataframe
  names(df_species_prey) <- c("Prey") # Rename variable
  df_species_prey$Prey_Species <- 1:length(model$param$groups$species) # Create new integer variable in ascending order with 'Species'
  df_species_prey$Prey_Species <- as.factor(df_species_prey$Prey_Species) # Convert to factor
  
  df_prey_size_classes <- model$w # Create new object with list of size classes in log10 g
  df_prey_size_classes <- as.data.frame(df_prey_size_classes) # Convert to a dataframe
  names(df_prey_size_classes) <- c("Prey_Size") # Create new name Size for the variable
  df_prey_size_classes$Prey_Size_Class <- 1:length(model$w)
  
  ## join df_N with df_size_class data in order to create log10 size classes rather then the size     class integer 1:184 that is in df_N
  df_combined_dynam_frac <- left_join(curr_dynam_frac_melt, df_size_classes,
                                      by="Predator_Size_Class") %>%
    
    left_join(df_species_predator,
              by="Predator_Species") %>%
    left_join(df_prey_size_classes,
              by="Prey_Size_Class") %>%
    left_join(df_species_prey,
              by="Prey_Species")
  
  df_comb_all <- df_combined_dynam_frac %>%
    group_by(Predator,Predator_Size, Prey) %>%
    summarise(Prey_Frac = sum(Prey_Frac)) 
  ###
  
  df_pd_frac <- df_combined_frac %>%
    select(Predator, Prey, Predator_Size, Prey_Frac)
  
  df_test_v3 <- bind_rows(df_comb_all,df_pd_frac)
  
  df_test_v3$Prey <- as.factor(df_test_v3$Prey)
  
  df_test_v3 <- df_test_v3[Reduce(`&`, lapply(df_test_v3, function(x) !is.na(x)  & is.finite(x))),]
  
  df_cephs_prey_frac <- df_test_v3 %>%
    filter(Prey == "Active Cephs"|Prey == "Inactive Cephs") %>%
    group_by(Predator, Predator_Size) %>%
    summarise(Prey_Frac = sum(Prey_Frac))
  
  df_cephs_prey_frac$Prey <- "Cephs"
  
  df_rest_prey <- df_test_v3 %>%
    filter(Prey == "Fish"|Prey == "Zooplankton"|Prey == "Phyto")
  
  df_ceph_prey_comb <- rbind(df_cephs_prey_frac, df_rest_prey)
  
  df_ceph_prey_comb$Prey <- as.factor(df_ceph_prey_comb$Prey)
  
  
  # ggplot(df_ceph_prey_comb, aes(x=log10(Predator_Size), y = Prey_Frac, fill = Prey)) +
  #   geom_area() +
  #   facet_wrap(~Predator, scales = "free_x") +
  #   scale_fill_manual(name = 'Prey',
  #                     values =c("#F8766D","#7CAE00", "green2","#C77CFF"),
  #                     labels = c('Cephs', 'Fish', "Phytoplankton",'Zooplankton')) +
  #   theme_classic() +
  #   xlab(expression(paste("Predator ",Size~(10^X~g)))) +
  #   ylab("Prey Proportion")
  
  return(df_ceph_prey_comb)
  
}


PPMR_Allometry <- function(projection){
  
  model <- Model
  
  # Average diet array for each predator size class for the final half of the projection
  dynamdiet = apply(model$dynam_diet_full[c(floor(0.5*dim(model$dynam_diet_full)[1]):dim(model$dynam_diet_full)[1]),,,,], c(2,3,4,5), mean)
  phytodiet = apply(model$phyto_diet_full[c(floor(0.5*dim(model$phyto_diet_full)[1]):dim(model$phyto_diet_full)[1]),,,], c(2,3,4), mean)
  
  phytodiet <- melt(phytodiet)
  # names(phytodiet)
  
  # Assign names to all variables in df_diet
  names(phytodiet) <- c("Predator_Species", "Predator_Size_Class", "Prey_Size_Class", "Prey_Abundance")
  
  # make predators and prey factor variables
  phytodiet$Predator_Species <- as.factor(phytodiet$Predator_Species)
  
  df_size_classes <- model$w # Create new object with list of size classes in log10 g
  df_size_classes <- as.data.frame(df_size_classes) # Convert to a dataframe
  names(df_size_classes) <- c("Predator_Size") # Create new name Size for the variable
  df_size_classes$Predator_Size_Class <- 1:length(model$w)
  
  df_species_predator <- model$param$groups$species # Create new object with species names from model
  df_species_predator <- as.data.frame((df_species_predator)) # Convert to a dataframe
  names(df_species_predator) <- c("Predator") # Rename variable
  df_species_predator$Predator_Species <- 1:length(model$param$groups$species) # Create new integer variable in ascending order with 'Species'
  df_species_predator$Predator_Species <- as.factor(df_species_predator$Predator_Species) # Convert to factor
  
  df_phyto_size_classes <- model$w_phyto # Create new object with list of size classes in log10 g
  df_phyto_size_classes <- as.data.frame(df_phyto_size_classes) # Convert to a dataframe
  names(df_phyto_size_classes) <- c("Prey_Size") # Create new name Size for the variable
  df_phyto_size_classes$Prey_Size_Class <- 1:length(model$w_phyto)
  
  ## join df_N with df_size_class data in order to create log10 size classes rather then the size     class integer 1:184 that is in df_N
  df_combined <- left_join(phytodiet, df_size_classes,
                           by="Predator_Size_Class") %>%
    
    left_join(df_species_predator,
              by="Predator_Species") %>%
    left_join(df_phyto_size_classes,
              by="Prey_Size_Class")
  
  df_combined <- df_combined %>%
    mutate(Prey_Biomass = Prey_Abundance*Prey_Size)
  
  df_combined$Prey <- c("Phyto")
  
  ### Dynamic diet dataframe
  
  dynamdiet <- melt(dynamdiet)
  
  # names(dynamdiet)
  
  # Assign names to all variables in df_diet
  names(dynamdiet) <- c("Predator_Species", "Predator_Size_Class", "Prey_Species", "Prey_Size_Class", "Prey_Abundance")
  
  # make predators and prey factor variables
  dynamdiet$Predator_Species <- as.factor(dynamdiet$Predator_Species)
  dynamdiet$Prey_Species <- as.factor(dynamdiet$Prey_Species)
  
  # df_size_classes <- model$w # Create new object with list of size classes in log10 g
  # df_size_classes <- as.data.frame(df_size_classes) # Convert to a dataframe
  # names(df_size_classes) <- c("Predator_Size") # Create new name Size for the variable
  # df_size_classes$Predator_Size_Class <- 1:length(model$w)
  
  df_species_prey <- model$param$groups$species # Create new object with species names from model
  df_species_prey <- as.data.frame((df_species_prey)) # Convert to a dataframe
  names(df_species_prey) <- c("Prey") # Rename variable
  df_species_prey$Prey_Species <- 1:length(model$param$groups$species) # Create new integer variable in ascending order with 'Species'
  df_species_prey$Prey_Species <- as.factor(df_species_prey$Prey_Species) # Convert to factor
  
  df_prey_size_classes <- model$w # Create new object with list of size classes in log10 g
  df_prey_size_classes <- as.data.frame(df_prey_size_classes) # Convert to a dataframe
  names(df_prey_size_classes) <- c("Prey_Size") # Create new name Size for the variable
  df_prey_size_classes$Prey_Size_Class <- 1:length(model$w)
  
  ## join df_N with df_size_class data in order to create log10 size classes rather then the size     class integer 1:184 that is in df_N
  df_combined_dynam <- left_join(dynamdiet, df_size_classes,
                                 by="Predator_Size_Class") %>%
    
    left_join(df_species_predator,
              by="Predator_Species") %>%
    left_join(df_prey_size_classes,
              by="Prey_Size_Class") %>%
    left_join(df_species_prey,
              by="Prey_Species")
  
  
  df_combined_dynam <- df_combined_dynam %>%
    mutate(Prey_Biomass = Prey_Abundance*Prey_Size)
  ###
  
  df_d <- df_combined_dynam %>%
    select(Predator, Prey, Predator_Size, Prey_Size, Prey_Abundance, Prey_Biomass)
  df_pd <- df_combined %>%
    select(Predator, Prey, Predator_Size, Prey_Size, Prey_Abundance, Prey_Biomass)
  
  df <- rbind(df_d, df_pd) # join data frames, one on top of the other
  
  
  # rm(df_d, df_pd)
  
  df_Pred_Diet_N_sum <- df %>%
    group_by(Predator, Predator_Size, Prey_Size) %>%
    summarise(Total_Prey_Abundance = sum(Prey_Abundance)) %>%
    filter(Total_Prey_Abundance > 0)
  # 
  # # reorder factors
  # df$Predator <- factor(df$Predator, 
  #                       levels = c("Zooplankton", "Fish", "Inactive Cephs", "Active Cephs"))
  # 
  # df$Prey <- factor(df$Prey, 
  #                   levels = c("Phyto", "Zooplankton", "Fish", "Inactive Cephs", "Active Cephs"))
  
  
  # create a data frame for plotting
  df_plot_comb <- df_Pred_Diet_N_sum %>%
    mutate(PySzPyAb = Prey_Size*Total_Prey_Abundance) %>%  # total size
    group_by(Predator, Predator_Size) %>% # groups of prey distributions
    summarise(
      sum_prey_abundance = sum(Total_Prey_Abundance), # total prey numbers
      sum_PySzPyAb       = sum(PySzPyAb), # total size
      mean_prey_size     = sum_PySzPyAb / sum_prey_abundance # mean size
    )
  
  # data for dashed line where expected prey size is 10% of predator size
  # df_ref <- tibble(
  #   Predator_Size = c(min(df_Pred_Diet_N_sum$Predator_Size), max(df_Pred_Diet_N_sum$Predator_Size))
  # )
  # df_ref$Expected_Prey_Size_PPMR100000 = 0.00001*df_ref$Predator_Size
  # df_ref$Expected_Prey_Size_PPMR10000 = 0.0001*df_ref$Predator_Size
  # df_ref$Expected_Prey_Size_PPMR1000 = 0.001*df_ref$Predator_Size
  # df_ref$Expected_Prey_Size_PPMR100 = 0.01*df_ref$Predator_Size
  # df_ref$Expected_Prey_Size_PPMR10 = 0.1*df_ref$Predator_Size
  
  
  # ggplot() +
  #   geom_line(data = df_ref, aes(x = Predator_Size, y = Expected_Prey_Size_PPMR10),
  #             linetype = "solid") +
  #   geom_line(data = df_ref, aes(x = Predator_Size, y = Expected_Prey_Size_PPMR100),
  #             linetype = "solid") +
  #   geom_line(data = df_ref, aes(x = Predator_Size, y = Expected_Prey_Size_PPMR1000),
  #             linetype = "longdash") +
  #   geom_line(data = df_ref, aes(x = Predator_Size, y = Expected_Prey_Size_PPMR10000),
  #             linetype = "dashed") +
  #   geom_line(data = df_ref, aes(x = Predator_Size, y = Expected_Prey_Size_PPMR100000),
  #             linetype = "dotted") +
  #   geom_line(data = df_plot_comb,
  #             aes(x = Predator_Size, y = mean_prey_size, color = Predator), size = 1) +
  #   scale_x_log10() +
  #   scale_y_log10() +
  #   # facet_wrap( ~ Predator) +
  #   labs(x = "Predator size", y = "Mean prey size") +
  #   theme_bw()
  
  
  return(df_plot_comb)
  
}

TL_Allom_Plot_Inde_Scaling_df<- function(projection){
  
  model <- projection
  
  # Average diet array for each predator size class for the final half of the projection
  dynamdiet = apply(model$dynam_diet_full[c(floor(0.5*dim(model$dynam_diet_full)[1]):dim(model$dynam_diet_full)[1]),,,,], c(2,3,4,5), mean)
  phytodiet = apply(model$phyto_diet_full[c(floor(0.5*dim(model$phyto_diet_full)[1]):dim(model$phyto_diet_full)[1]),,,], c(2,3,4), mean)
  
  phyto_tl <- 1 # TL for phytoplankton is 1
  start_dynam_tl <- matrix(2, nrow = dim(dynamdiet)[1], ncol = dim(dynamdiet)[2]) # Start TL for dynamic size groups - start at 2 for all groups and sizes
  
  curr_phyto_diet <- phytodiet # Current phyto diet
  curr_dynam_diet <- dynamdiet # Current heterotroph diet
  
  total_diet <- apply(curr_phyto_diet, c(1,2), sum) + apply(curr_dynam_diet, c(1,2), sum) # Total consumption, in grams wet weight, by pred group and pred size
  
  curr_phyto_frac <- apply(sweep(curr_phyto_diet, c(1,2), total_diet, '/'), c(1,2), sum) # Fraction of diet from phyto, by pred group and pred sizes
  curr_dynam_frac <- sweep(curr_dynam_diet, c(1,2), total_diet, '/') # Fraction of diet from each prey group and prey size, for each pred group and pred size
  
  pb = txtProgressBar(min = 0, max = 100, initial = 0, style = 3)
  
  for(j in 1:100){ # Gauss-Siedel iterative loop to calculate trophic levels
    setTxtProgressBar(pb, j)
    
    calc_dynam_tl = sweep(curr_dynam_frac, c(3,4), start_dynam_tl, '*')
    calc_dynam_tl[which(is.nan(calc_dynam_tl) == TRUE)] = 0 # Get rid of nans - these are entrys where there is no biomass for a given group
    #calc_dynam_tl[which(calc_dynam_tl == Inf)] = 0 # Get rid of infinite values, occurs with asymptotic size bins, because there is no biomass to have a diet in those bins
    start_dynam_tl = 1+phyto_tl*curr_phyto_frac + apply(calc_dynam_tl, c(1,2), sum) # Update trophic level matrix
  } # End Gauss-Siedel loop
  
  # Create long form diet dataframe
  start_dynam_tl <- melt(start_dynam_tl)
  
  # Assign names to all variables in df_diet
  names(start_dynam_tl) <- c("Predator_Species", "Predator_Size_Class", "RTP")
  
  # make predators and prey factor variables
  start_dynam_tl$Predator_Species <- as.factor(start_dynam_tl$Predator_Species)
  
  df_size_classes <- model$w # Create new object with list of size classes in log10 g
  df_size_classes <- as.data.frame(df_size_classes) # Convert to a dataframe
  names(df_size_classes) <- c("Predator_Size") # Create new name Size for the variable
  df_size_classes$Predator_Size_Class <- 1:length(model$w)
  
  df_species_predator <- model$param$groups$species # Create new object with species names from model
  df_species_predator <- as.data.frame((df_species_predator)) # Convert to a dataframe
  names(df_species_predator) <- c("Predator") # Rename variable
  df_species_predator$Predator_Species <- 1:length(model$param$groups$species) # Create new integer variable in ascending order with 'Species'
  df_species_predator$Predator_Species <- as.factor(df_species_predator$Predator_Species) # Convert to factor
  
  ## join df_N with df_size_class data in order to create log10 size classes rather then the size     class integer 1:184 that is in df_N
  df_combined <- left_join(start_dynam_tl, df_size_classes,
                           by="Predator_Size_Class") %>%
    
    left_join(df_species_predator,
              by="Predator_Species") 
  
  # Read in empirical trophic level data
  df_emp_SIA <- read_csv("df_emp_SIA.csv")
  # Make copy
  df_emp_SIA_rel <- df_emp_SIA
  
  #####
  # Rescale empirical RTP for active cephs
  df_emp_active <- df_emp_SIA_rel %>%
    filter(name == "Active Cephs")
  
  model_lm_emp_active <- lm(RTP ~ log10(Predator_Size), data = df_emp_active)
  
  df_active <- df_combined %>%
    filter(Predator == "Active Cephs")
  
  min_emp_RTP_active <- model_lm_emp_active$coefficients[1] + (model_lm_emp_active$coefficients[2] * log10(min(df_emp_active$Predator_Size)))
  
  model_lm_active <- lm(RTP ~ log10(Predator_Size), data = df_active)
  
  df_model_min_RTP <- df_active %>%
    filter(Predator_Size_Class == 111 | Predator_Size_Class == 112) %>%
    # filter(Predator_Size_Class == 220 | Predator_Size_Class == 221) %>%
    mutate(avg_min_RTP = mean(RTP))
  
  min_mod_RTP_active <- df_model_min_RTP$avg_min_RTP[1]
  
  rescale_RTP_active_val <- min_mod_RTP_active - min_emp_RTP_active
  
  pos_val_active <- sqrt(rescale_RTP_active_val * rescale_RTP_active_val)
  
  
  #####
  # Rescale empirical RTP for inactive group
  df_emp_inactive <- df_emp_SIA_rel %>%
    filter(name == "Inactive Cephs")
  
  df_inactive <- df_combined %>%
    filter(Predator == "Inactive Cephs")
  
  
  model_lm_emp_inactive <- lm(RTP ~ log10(Predator_Size), data = df_emp_inactive)
  
  min_emp_RTP_inactive <- model_lm_emp_inactive$coefficients[1] + (model_lm_emp_inactive$coefficients[2] * log10(min(df_emp_inactive$Predator_Size)))
  
  
  df_model_min_RTP_inactive <- df_inactive %>%
    filter(Predator_Size_Class == 111 | Predator_Size_Class == 112) %>%
    # filter(Predator_Size_Class == 220 | Predator_Size_Class == 221) %>%
    mutate(avg_min_RTP = mean(RTP))
  
  min_mod_RTP_inactive <- df_model_min_RTP_inactive$avg_min_RTP[1]
  
  rescale_RTP_inactive_val <- min_mod_RTP_inactive - min_emp_RTP_inactive
  
  pos_val_inactive <- sqrt(rescale_RTP_inactive_val * rescale_RTP_inactive_val)
  
  Act_rescale_val <-  pos_val_active
  
  Inact_rescale_val <- pos_val_inactive 
  
  df_emp_inactive$RTP <- df_emp_inactive$RTP + Inact_rescale_val
  
  df_emp_active$RTP <- df_emp_active$RTP + Act_rescale_val
  
  df_zoom <- df_combined %>%
    filter(Predator == "Inactive Cephs" | Predator ==  "Active Cephs")
  
  return(df_combined)
  
}

