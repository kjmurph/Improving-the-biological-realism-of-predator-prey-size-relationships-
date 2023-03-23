## An extension of the model published in Heneghan et al., (2016):
## Fixed spectrum: phytoplankton spectrum, driven by longhurst province environmental data
## Dynamic spectra: Models single zoo, single fish and however many ceph communities specified
## No temperature effect in this model
## Pre-maturation senesence factor 'ZSpre' increased from 0.02 to 0.1

## Set up Model Parameter List
params <- function(fileGroups, tmax, intercept, slope, phyto_min, phyto_max, f_mort){
  
  # Read in functional group specific parameters from file:
  groups = fileGroups
  #	groups = read.csv(fileGroups)
  nutrition = groups$carbon
  environ = list("a" = intercept, "b" = slope, "phyto_min" = phyto_min, "phyto_max" = phyto_max)
  environ = as.data.frame(environ)
  
  # Set up parameter list
  param = list(groups = groups, 
               environ = environ,
               nutrition = nutrition,
               
               # Model parameters
               ngrps = dim(groups)[1],		# no. of groups
               tmax = tmax,					# no. of years
               dt = 0.08333333, 		# monthly time step 
              # dt = 0.02, 				# time step (0.02 of a year)
               # dx = 0.05,         # reduced size difference in log10 weight step
               dx = 0.1,         # log10 weight step
               day = 12,          # day length (hours of each day in sun)
               gge_base = 0.25, # baseline gross growth efficiency
               w0 = 10^(min(groups$W0)),		# minimum dynamic size class
               wMax = 10^(max(groups$Wmax)),# maximum dynamic size class
               # Original value of 0.02, but Ryan suggested this was probably too low and increasing to 0.1 may reduce travelling waves
               ZSpre = 0.1, # senescence mortality prefactor
               ZSexp = 0.3,                 # senescence mortality exponent
               f_mort = f_mort, # fishing mortality
               w0_phyto = 10^(environ$phyto_min),		# minimum phytoplankton size class (1um)
               wMax_phyto = 10^environ$phyto_max		# maximum phytoplankton size class
               
  )
  # param$isave = 10	# how often to save results (every 'isave' time steps
  param$isave = 12	# how often to save results (every 'isave' time steps
  if(param$isave == 0){param$isave = 1}
  param$fish_grps = which(groups$community == "fish") # Which rows are fish
  param$ceph_grps = which(groups$community == "ceph") # Which rows are cephs
  param$zoo_grps = which(groups$community == "zoo") # Which rows are zooplankton
  param$num_zoo = length(which(groups$community == "zoo")) # How many zooplankton
  param$num_fish = length(which(groups$community == "fish"))  # How many fish
  param$num_ceph = length(which(groups$community == "ceph"))  # How many cephs
  return(param)
}

## Set up Model components
## Create a list of model components, and initial values for N and nPP. Call the list "model"
Setup <- function(param, restart, saved_N){
  # Pull out some useful parameters - just a shortcut
  grp = param$groups
  ngrps = param$ngrps
  dt = param$dt
  dx = param$dx
  environ = param$environ
  fish_grps = param$fish_grps
  zoo_grps = param$zoo_grps
  ceph_grps = param$ceph_grps
  num_zoo = param$num_zoo
  num_fish = param$num_fish
  num_ceph = param$num_ceph
  
  # Set up dynamic grid
  w <- 10^(seq(from = log10(param$w0), to =  log10(param$wMax), dx))
  ngrid <- length(w)
  
  # Set up phytoplankton size classes
  w_phyto <- 10^(seq(from = log10(param$w0_phyto), to = log10(param$wMax_phyto), dx))
  ngridPP <- length(w_phyto)
  
  # Number of time slots to save
  nsave   <- floor(param$tmax/(param$dt*param$isave)) 
  
  # Diel Vertical Migration - change availability of phyto to zoo
  # and zoo to fish based on slope of phytoplankton (calculated as
  # proportion of day searching for food and available for predation)
  
  dvm_max = param$day/24 # maximum proportion of day spent migrating
  ESD_sizes = 2*(3/(4*pi)*w)^(1/3) # convert g wet weight to ESD (cm)
  dvm <- dvm_max*(1.02*(ESD_sizes) - 0.02) # size-dependent amount of time spent away from surface
  dvm[which(w < 10^-5.4)] = 0 # Microzoo don't migrate (ESD < 0.02cm)
  dvm[which(w > 10^-0.3)] = dvm_max # Macrozoo don't migrate (ESD > 10cm)
  
  # Dynamic prey availability matrix: dim1 is predators, dim2 is predator size classes,
  # dim3 is prey groups, dim 4 is prey size classes.
  
  dvm_mat = matrix(dvm, nrow = ngrid, ncol = ngrid, byrow = TRUE)
  dvm_mat = 1 - dvm_mat
  dvm_mat[lower.tri(dvm_mat)] <- 0
  dvm_mat <- t(dvm_mat) + dvm_mat
  diag(dvm_mat) <- diag(dvm_mat)/2
  
  dynam_theta = array(1, dim = c(ngrps, ngrid, ngrps, ngrid))
  dynam_theta = sweep(dynam_theta, c(2,4), dvm_mat,"*")
  
  # Phyto availability matrix: rows are predators, columns are their size classes,
  # entries are time spent feeding on phytoplankton for the size class
  phyto_theta = matrix(1-dvm, nrow = ngrps, ncol = ngrid, byrow = TRUE)
  
  ### REMOVE DVM - I HAVE TURNED OFF DIEL VERTICAL MIGRATION AT THE MOMENT
  dynam_theta = array(1, dim = c(ngrps, ngrid, ngrps, ngrid))
  phyto_theta = matrix(1, nrow = ngrps, ncol = ngrid, byrow = TRUE)
  ###
  
  herb_grps = which(grp$type == 'H')
  carn_grps = which(grp$type == 'C')
  phyto_theta[carn_grps,] = 0 
  dynam_theta[herb_grps,,,] = 0 # Herbivorous groups only prey on phyto
  
  cc_phyto <- 0.15   # Carbon content of phytoplankton size classes
  
  ## Makes the model object, full of constant functions for model
  model <- list(
    param = param,
    environ = environ,
    ngrid = ngrid,
    ngridPP = ngridPP,
    
    # Phytoplankton abundance
    nPP = 10^(environ$a)*(w_phyto^(environ$b)), # phyto abundance spectrum
    
    # Grid parameters
    w = w,
    dx = dx,
    w_phyto = w_phyto,
    
    # Group parameters storage
    phyto_growthkernel = array(NA, dim = c(ngrps, ngrid, ngridPP)),
    # ngrps = Predators, ngrid = size classes of predators, ngridPP = Prey size classes
    phyto_diffkernel = array(NA, dim = c(ngrps, ngrid, ngridPP)), # predation on phytoplankton
    dynam_growthkernel =  array(NA, dim = c(ngrps, ngrid, ngrid)),
    dynam_diffkernel = array(NA, dim = c(ngrps, ngrid, ngrid)), # predation on dynamic component
    dynam_mortkernel = array(NA, dim = c(ngrps, ngrid, ngrid)), # mortality from predation on dynamic component
    M_sb = matrix(0, nrow = ngrps, ncol = ngrid), # senescence + background mortality
    psi = matrix(0, nrow = ngrps, ncol = ngrid), # maturation function
    fish_mort = matrix(0, nrow = ngrps, ncol = ngrid), # fishing mortality
    
    #### STORAGE FOR DIET KERNELS
    phyto_dietkernel = array(NA, dim = c(ngrps, ngrid, ngridPP)),
    dynam_dietkernel = array(NA, dim = c(ngrps, ngrid, ngrid)),
    
    # Output storage
    N = array(0, dim = c(nsave, ngrps, ngrid)),#, # dynamic abundance spectrum
   Z = array(0, dim = c(nsave, ngrps, ngrid)), # total mortality
   gg = array(0, dim = c(nsave, ngrps, ngrid)), # growth
   phyto_diet_full = array(0, dim = c(nsave, ngrps, ngrid)), # Save phyto diet, dim1 = save step, dim2 = pred group, dim3 = prey size
   dynam_diet_full = array(0, dim = c(nsave, ngrps, ngrid, ngrps, ngrid)) # Save dynam diet, dim1 = save step, dim2 = pred group, dim3 = pred size, dim4 = prey group, dim5 = prey size
   # diet = array(0, dim = c(nsave, ngrps, ngrid, ngrps, ngrid)), # diet (ngrps = predators, ngrid = predator size class, ngrps = prey, ngrid = prey size class)
   # Biomass = matrix(0, nrow = nsave, ncol = ngrps), # biomass of each group
   # Diff = array(0, dim = c(nsave, ngrps, ngrid)) # save diffusion
  )
  
  # GGE for different groups
  assim_phyto =  param$groups$alpha*cc_phyto/param$nutrition # Phytoplankton
  assim_dynam =  matrix(param$groups$alpha*param$nutrition, nrow = ngrps, ncol = ngrps, byrow = TRUE)/
    matrix(param$nutrition, nrow = ngrps, ncol = ngrps)
  
  #### INITIAL DYNAMIC POPULATION ABUNDANCES
  a_dynam = 10^(environ$a)*(w[1]^(environ$b+1)) # calculate coefficient for initial dynamic spectrum
  
  # Initial abundances form a continuation of the plankton spectrum
  tempN <- matrix(a_dynam*(w)^-1, nrow = ngrps, ncol = ngrid, byrow = TRUE) 
  props_z <- grp$prop[zoo_grps] # Zooplankton proportions 
  props_z <- as.numeric(props_z)
  tempN[zoo_grps,] <- props_z*tempN[zoo_grps,]

  props_ceph <- grp$prop[ceph_grps] # Ceph proportions
  props_ceph <- as.numeric(props_ceph) 
  tempN[ceph_grps,] <- props_ceph*tempN[ceph_grps,]
  
  # For each group, set densities at w > Winf and w < Wmin to 0
  tempN[unlist(tapply(round(log10(w), digits = 2), 1:length(w), function(wx,Winf) Winf < wx, Winf = (grp$Wmax)))] <- 0
  tempN[unlist(tapply(round(log10(w), digits = 2), 1:length(w), function(wx,Wmin) Wmin > wx, Wmin = (grp$W0)))] <- 0
  if (restart == TRUE) {
    model$N[1,,] <- tempN
  } else {
    model$N[1,,] <- saved_N
  }
  
  #model$N[1,,] <- tempN
  
  # Fishing mortality
  model$fish_mort[fish_grps, c(w >= 2.30103)] = param$f_mort
  
  # Untransformed body size intervals for dynamic spectrum
  dw <- c(diff(model$w),diff(model$w)[length(diff(model$w))])
  model$dw <- dw # Save dw for later
  
  ### MATRICES FOR LOG TRANSFORM OF EQUATION
  # Predators are rows, phyto prey weights are columns   
  gg_log_t_phyto = ((w^-1) %*% t(w_phyto))/log(10) # Growth
  diff_log_t_phyto = ((w^-1) %*% t(w_phyto^2))/log(10) # Diffusion
  diet_log_t_phyto = matrix(w_phyto, nrow = length(w), ncol = length(w_phyto), byrow = TRUE) # Diet
  
  # Predators are rows, dynam prey weights are columns
  gg_log_t_dynam = ((w^-1) %*% t(w))/log(10) # Growth
  diff_log_t_dynam = ((w^-1) %*% t(w^2))/log(10) # Diffusion
  diet_log_t_dynam = matrix(w, nrow = length(w), ncol = length(w), byrow = TRUE) # Diet
  
  ### PREDATION KERNELS FOR PHYTOPLANKTON SPECTRUM AND DYNAMIC SPECTRUM
  phyto_pred_weight_matrix = matrix(w, nrow = ngrid, ncol = ngridPP)
  dynam_pred_weight_matrix = matrix(w, nrow = ngrid, ncol = ngrid)
  phyto_prey_weight_matrix = matrix(w_phyto, nrow = ngrid, ncol = ngridPP, byrow = TRUE)
  dynam_prey_weight_matrix = matrix(w, nrow = ngrid, ncol = ngrid, byrow = TRUE)
  
  ## Search Volume storage
  SearchVol = matrix(NA, nrow = ngrps, ncol = ngrid) # Search volume
  
  # Simpson's Rule matrices for growth, diffusion and mortality integrals
  simp_phyto = array(1, dim = ngridPP)
  simp_phyto[c(seq(2,ngridPP-1,2))] = 4
  simp_phyto[c(seq(3,ngridPP-1,2))] = 2
  sm_phyto = matrix(simp_phyto, nrow = ngrid, ncol = ngridPP, byrow = TRUE)*(dx/3)
  
  simp_dynam = array(1, dim = ngrid)
  simp_dynam[c(seq(2,ngrid-1,2))] = 4
  simp_dynam[c(seq(3,ngrid-1,2))] = 2
  sm_dynam = matrix(simp_dynam, nrow = ngrid, ncol = ngrid, byrow = TRUE)*(dx/3)
  
  ## Temperature Effect Matrix
  # Effect of temperature on feeding rate and background mortality
 
  ### NO TEMPERATURE EFFECT - ALL SET TO 1
  temp_zoo <- rep(1, num_zoo) # exp(23.93 - 0.58/(8.62e-05*(273+environ$sst)))
  temp_ceph <- rep(1, num_ceph)
  temp_fish <- rep(1, num_fish) # exp(25.55 - 0.63/(8.62e-05*(273+environ$sst)))
  temp_effect <- matrix(c(temp_zoo, temp_ceph, temp_fish), nrow = ngrps, ncol = ngrid)
  
  #### CALCULATES CONSTANT BITS OF THE MODEL FUNCTIONS FOR EACH GROUP
  for(i in 1:ngrps){ 
    
    ## Maturation function
    curr_weights = c(log10(w) >= grp$W0[i] & log10(w) <= grp$Wmax[i])
    model$psi[i,curr_weights] <- (1+(w[curr_weights]/(0.25*max(w[curr_weights])))^(-10))^(-1)
    
    ## Senescence mortality
    model$M_sb[i,] <- param$ZSpre*(w/(10^(grp$Wmat[i])))^param$ZSexp
    model$M_sb[i, 10^(grp$Wmax[i]) < w] <- 0
    model$M_sb[i, 10^(grp$Wmat[i]) > w] <- 0	
    
    ### Search volume
    SearchVol[i,] <- (grp$carbon[i]*grp$gamma[i])*(w^(grp$q[i]))
    SearchVol[i, 10^(grp$Wmax[i]) < w] <- 0
    SearchVol[i, 10^(grp$W0[i]) > w] <- 0		
    
    ### Predation Kernels
    if(is.na(grp$m[i]) == FALSE){ # If group has an m-value (zooplankton)
      # Calculate PPMR for zooplankton, which changes according to body-size (Wirtz, 2012)	
      D.z <- 2*(3*w*1e12/(4*pi))^(1/3) # convert body mass g to ESD (um)
      betas =  (exp(0.02*log(D.z)^2 - grp$m[i] + 1.832))^3
      beta_mat_phyto = matrix(betas, nrow = ngrid, ncol = ngridPP)
      beta_mat_dynam = matrix(betas, nrow = ngrid, ncol = ngrid)
      
      sp_phyto_predkernel = exp(-0.5*(log((beta_mat_phyto*phyto_prey_weight_matrix)/
                                            phyto_pred_weight_matrix)/grp$sigma[i])^2)/
        sqrt(2*pi*grp$sigma[i]^2)
      sp_dynam_predkernel = exp(-0.5*(log((beta_mat_dynam*dynam_prey_weight_matrix)/
                                            dynam_pred_weight_matrix)/grp$sigma[i])^2)/
        sqrt(2*pi*grp$sigma[i]^2)
      
    } else { # If group does not have an m-value (fish)
      beta_mat_phyto = matrix(grp$beta[i], nrow = ngrid, ncol = ngridPP)
      beta_mat_dynam = matrix(grp$beta[i], nrow = ngrid, ncol = ngrid)
      
      sp_phyto_predkernel = exp(-0.5*(log((beta_mat_phyto*phyto_prey_weight_matrix)/
                                            phyto_pred_weight_matrix)/grp$sigma[i])^2)/
        sqrt(2*pi*grp$sigma[i]^2)
      sp_dynam_predkernel = exp(-0.5*(log((beta_mat_dynam*dynam_prey_weight_matrix)/
                                            dynam_pred_weight_matrix)/grp$sigma[i])^2)/
        sqrt(2*pi*grp$sigma[i]^2)
    }
    
    ### GROWTH INTEGRAL CONSTANTS
    # Predators are rows, prey are columns	
    model$phyto_growthkernel[i,,] = matrix(SearchVol[i,], nrow = ngrid, ncol = ngridPP)*
      sp_phyto_predkernel*gg_log_t_phyto*sm_phyto
    model$dynam_growthkernel[i,,] = matrix(SearchVol[i,], nrow = ngrid, ncol = ngrid)*
      sp_dynam_predkernel*gg_log_t_dynam*sm_dynam
    
    ### DIFFUSION INTEGRAL CONSTANTS
    # Predators are rows, prey are columns	
    model$phyto_diffkernel[i,,] = matrix(SearchVol[i,], nrow = ngrid, ncol = ngridPP)*
      sp_phyto_predkernel*diff_log_t_phyto*sm_phyto
    model$dynam_diffkernel[i,,] = matrix(SearchVol[i,], nrow = ngrid, ncol = ngrid)*
      sp_dynam_predkernel*diff_log_t_dynam*sm_dynam	
    
    ### DIET INTEGRAL CONSTANTS
    # Predators are rows, prey are columns
    model$phyto_dietkernel[i,,] = matrix(SearchVol[i,], nrow = ngrid, ncol = ngridPP)*
      sp_phyto_predkernel*diet_log_t_phyto*sm_phyto
    model$dynam_dietkernel[i,,] = matrix(SearchVol[i,], nrow = ngrid, ncol = ngrid)*
      sp_dynam_predkernel*diet_log_t_dynam*sm_dynam
    
    ### MORTALITY INTEGRAL CONSTANTS                                  
    # Prey are rows, predators are columns
    model$dynam_mortkernel[i,,] = matrix(SearchVol[i,], nrow = ngrid, ncol = ngrid, byrow = TRUE)*
      t(sp_dynam_predkernel)*sm_dynam 
  }
  
  #no_sen = which(grp$species == c("Flagellates", "Ciliates")) # no senescence mortality for flagellates and ciliates
  model$M_sb = temp_effect*model$M_sb 
  model$M_sb[c(ngrps)] = 0
  
  model$phyto_growthkernel = sweep(sweep(model$phyto_growthkernel, c(1,2), phyto_theta, "*"), 1, assim_phyto, "*")
  model$phyto_diffkernel = sweep(sweep(model$phyto_diffkernel, c(1,2), phyto_theta, "*"), 1, assim_phyto^2, "*")
  model$phyto_dietkernel =  sweep(sweep(model$phyto_dietkernel, c(1,2), phyto_theta, "*"), 1, 1, "*")
  
  model$dynam_growthkernel = sweep(sweep(sweep(dynam_theta, c(1,2,4), model$dynam_growthkernel, "*"), 
                                         c(1,3), assim_dynam, "*"), c(1,2), temp_effect, "*")
  model$dynam_diffkernel = sweep(sweep(sweep(dynam_theta, c(1,2,4), model$dynam_diffkernel, "*"), c(1,3), assim_dynam^2, "*"),
                                 c(1,2), temp_effect^2, "*")     
  model$dynam_mortkernel = sweep(aperm(sweep(aperm(dynam_theta, c(3,1,4,2)), 
                                             c(2,3,4), model$dynam_mortkernel, "*"), c(1,3,2,4)),
                                 c(3,4), temp_effect, "*")
  model$dynam_dietkernel = sweep(sweep(sweep(dynam_theta, c(1,2,4), model$dynam_dietkernel, "*"), 
                                       c(1,3), 1, "*"), c(1,2), temp_effect, "*")
  
  model$ingested_phyto = temp_effect*(rowSums(sweep(model$phyto_growthkernel, 3, model$nPP, "*"), dims = 2))
  model$diff_phyto = temp_effect^2*(rowSums(sweep(model$phyto_diffkernel, 3, model$nPP, "*"), dims = 2))
  
  #### ADDED BY R HENEGHAN - THU 19 DECEMBER 2019
  model$diet_phyto = temp_effect*(rowSums(sweep(model$phyto_dietkernel, 3, model$nPP, "*"), dims = 2)) # Diet of total phyto
  
  return(model)
} # End of Setup function

## Run model forward in time
Project <- function(model){
  
  # Pull out some useful parameters - just a shortcut
  param <- model$param
  grp <- param$group
  ngrid <- model$ngrid
  ngridPP <- model$ngridPP
  ngrps <- param$ngrps
  dt <- param$dt
  fish_grps <- param$fish_grps
  zoo_grps <- param$zoo_grps
  ceph_grps <- param$ceph_grps
  dx <- model$dx
  w <- model$w
  w_phyto <- model$w_phyto
  dw <- model$dw
  assims <- param$nutrition/0.1*(grp$alpha)
  # w0idx <- which(grp$W0 > min(grp$W0) & is.na(grp$prop) == FALSE) # JB: group index - which groups is not the smallest  minimim egg size group, taht also has a prop? This is becuase the smallest group was  hinged to phytoplankton, whereas others are not 
  # w0mins <- rep(0, length(w0idx))
  
  #props_ceph <- grp$prop[w0idx] # Ceph proportions - JB: changed to ceph and fish group proportions, coupld just be props if all groups are treated in same way
  props <- grp$prop # just a shortcut
  props <- as.numeric(props)
  w_mat = matrix(w, nrow = ngrps, ncol = ngrid, byrow = TRUE)
  
  # for(i in 1:length(w0idx)){
  #   w0mins[i] <- which(round(log10(w), digits = 2) == grp$W0[w0idx[i]])
  # }
  
  # JB: we still need this but change to number of dynamic groups rather than the subset of groups, also use a sweep
  # dont need a loop as which() works on a vector 
  w0mins <- which(round(log10(w), digits = 2) == grp$W0)
  
  # JB - this is already set to w0mins
  # Which size element should eggs be dumped into
  # w0idx <- as.vector(tapply(grp$W0,1:length(grp$W0),function(W0,wx) max(which(wx<=W0)),wx=log10(w)))
  
  
  
  # Handy stuff - 2:ngrid as it is skipping the 1st size class
  idx <- 2:ngrid
  itimemax  <- param$tmax / dt  #max index of time array
  
  # Matrices for solver
  A.iter <- C.iter <- S.iter <- A <- B <- C <- S <- matrix(0,nrow=ngrps,ncol=ngrid)
  
  # Temporary Matrices that get updated each time step
  # some of these saved for output
  N <- matrix(model$N[1,,], nrow = ngrps, ncol = ngrid)
  nPP <- model$nPP
  
  pb = txtProgressBar(min = 0, max = itimemax, initial = 1, style = 3)
  # BIG TIME LOOP
  for (itime in 1:itimemax)
  {
    setTxtProgressBar(pb, itime)
    ### Create an ngrps*ngrid*ngrps*ngrid array of abundances, to save time without sweeps
    # dim1 = pred groups, dim 2 = pred sizes, dim 3 = prey groups, dim 4 = prey sizes
    
    N_array <- aperm(replicate(ngrid, N), c(3,1,2))
    N_array <- aperm(replicate(ngrps, N_array), c(4,1,2,3))
    
    ### GROWTH
    gg <- (model$ingested_phyto + 
             rowSums(rowSums(model$dynam_growthkernel*N_array, dims = 3), dims = 2))     
    
    ### MORTALITY
    # Predation mortality
    M2 <- (rowSums(rowSums(model$dynam_mortkernel*N_array, dims = 3), dims = 2))
    
    # Total dynamic spectrum mortality
    Z = M2 + model$M_sb  + model$fish_mort
    
    
    ### DIFFUSION
    diff <- (model$diff_phyto + rowSums(rowSums(model$dynam_diffkernel*N_array, dims = 3), dims = 2))
    
    ### MvF WITH DIFFUSION ALGORITHM
    
    idx.iter <- 2:ngrid
    idx <- 2:(ngrid-1)
    
    # Numerical implementation matrices
    A.iter[,idx.iter] <- dt/dx*gg[,idx.iter-1] 
    C.iter[,idx.iter] <- 1 + dt*Z[,idx.iter] + dt/dx*gg[,idx.iter]
    S.iter[,idx.iter] <- N[,idx.iter]
    N.iter <- N
    
    A[,idx] <- dt/dx*(gg[,idx-1] + diff[,idx-1]/(2*dx))
    B[,idx] <- diff[,idx+1]*dt/(2*dx^2)
    C[,idx] <- 1 + dt*Z[,idx] + dt/dx*(gg[,idx] + diff[,idx]/dx)
    S[,idx] <- N[,idx]
    
    for(i in 1:ngrps){
      
      ## Set size range index for current group
      curr_min_size = which(round(log10(w), digits = 2) == param$groups$W0[i])
      curr_max_size = which(round(log10(w), digits = 2) == param$groups$Wmax[i])
      idx_curr = (curr_min_size+1):curr_max_size
      
      for(j in idx_curr){## Find the abundance at the next size class with standard MvF
        N.iter[i,j] <- (S.iter[i,j] + A.iter[i,j]*N[i,j-1])/(C.iter[i,j])
        if(j >= (idx_curr[1]+1)){ ## Find abundance with MvF with diffusion
          k = j - 1
          N[i,k] = (S[i,k] + A[i,k]*N[i,k-1] + B[i,k]*N.iter[i,k+1])/C[i,k]
        }
        # MvF without diffusion for last size class  
        if(j == idx_curr[length(idx_curr)]){ 
          N[i,j] = 0
          # N[i,curr_min_size] <- N.iter[i,curr_min_size] # Keep starting sizes constant
        }
      }
    }
    
    
    
    #### Keep smallest fish community size class as equal to equivalent zooplankton size class
    
    
    ### Keep smallest ceph size class abundnace 
    ### for each group locked to others in size spectrum
    ### JB: Kieran and I discussed and changed to reflect these as proportion of 
    ### what the fixed phytoplankton would be at each groups min size (assuming it would extend to the whole community) 
    ### this is a simpler assumption

    for(n in 1:ngrps){ 
    N[n, w0mins[n]] = props[n]*10^(param$environ$a)*w[w0mins[n]]^(param$environ$b)
    }
    
    # #for(i in 1:length(w0idx)){
    #   w_min_curr = w0mins[i]
    #   exclude_mins = w0idx[which(w0mins == w_min_curr)]
    #  #N[w0idx[i], w_min_curr] = props_ceph[i]*sum(N[-exclude_mins, w_min_curr])
    #   N[w0idx[i], w_min_curr] = props[i]*10^(param$environ$a)*w[w_min_curr]^(param$environ$b)
    #   # make a proportion of the fixed phyto abundance spectrum
    #   
    #    }
    # 
    # 
    # fish_mins = unlist(lapply(param$groups$W0[fish_grps], 
    #                           function(x){which(round(log10(model$w), digits = 2) == x)}))
    # 
    # if(length(fish_grps) > 1){
    #   N[fish_grps,fish_mins] = (1/3)*(colSums(N[-fish_grps,fish_mins]))
    # }else{
    #   N[fish_grps, fish_mins] = sum(N[-fish_grps, fish_mins])
    # }
    
    # if(fish_on == FALSE){
    #   N[fish_grps,] = 0 # switch off fish_groups
    # }
    
    # Save results:
    if((itime %% param$isave) == 0){
      isav<-itime/param$isave
      
      
      ##### NEW ADDITION - TUE 17 DEC 2019, BY RYAN HENEGHAN - SAVE PHYTO DIET AND DYNAM DIET,
      ##### PHYTO DIET DIM1 = SAVE STEP, DIM2 = PRED GROUP, DIM3 = PRED SIZE
      ##### DYNAM DIET DIM1 = SAVE STEP, DIM2 = PRED GROUP, DIM3 = PRED SIZE, DIM 4 = PREY GROUP, DIM 5 = PREY SIZE
      phyto_diet = apply(sweep(model$diet_phyto, c(1,2), N,  "*"),c(1,2), sum)
      model$phyto_diet_full[isav,,] = phyto_diet
      dynam_diet =  sweep(model$dynam_dietkernel*N_array, c(1,2), N, "*")
      model$dynam_diet_full[isav,,,,] = dynam_diet
      
      #phyto_diet = rowSums(model$diet_phyto*N)
      #dynam_diet =  sweep(model$dynam_dietkernel*N_array, c(1,2), N, "*")
      #model$diet[isav,,,,] = dynam_diet
      
      # model$gg[isav,,] <- gg # Save growth rates
      model$N[isav,,] <- N # Save abundance
      if(length(zoo_grps) > 1){
        model$N[isav,c(1,2),c(60,61)] <- 0
      }
      # model$Biomass[isav,] <- rowSums(model$N[isav,,] # Save biomass
                                      # *matrix(model$w, nrow = ngrps, ncol = ngrid, byrow = TRUE)) 
      model$Z[isav,,] <- Z # Save total mortality rates
      #	model$Diff[isav,,] <- diff # Save total diffusion rates
    }
    
  } # End of time loop
  
  return(model)
  
}

#ave_diet = colMeans(modelss$diet[(ceiling(0.5*dim(modelss$diet)[1])):(dim(modelss$diet)[1]),,,,], dim = 1)


#### PLOT AVERAGE SPECTRUM
Spectrum_Plot <- function(models, fish_on){
  
  model = models
  modelss = models
  param = modelss$param
  
  N_sav = modelss$N
  N_ave = colMeans(modelss$N[(ceiling(0.5*dim(N_sav)[1])):(dim(N_sav)[1]),,], dim = 1)
  zoo_groups =  param$zoo_grps
  ceph_groups = param$ceph_grps
  fish_groups = param$fish_grps
  num_zoo = param$num_zoo
  num_cephs = param$num_ceph
  num_fish = param$num_fish
  
  if(length(zoo_groups) > 1){
    tot_zoo = colSums(N_ave[zoo_groups,])
  }
  
  if(length(zoo_groups) == 1){
    tot_zoo = N_ave[zoo_groups,]  
  }
  
  if(length(ceph_groups) > 1){
    tot_cephs = colSums(N_ave[ceph_groups,])
  }
  
  if(length(ceph_groups) == 1){
    tot_cephs = N_ave[ceph_groups,]
  }
  
  if(num_fish > 0){
    if(num_fish > 1){
      tot_fish = colSums(N_ave[fish_groups,])
    }else{tot_fish = N_ave[fish_groups,]}
  }
  
  y_d_ref = which(round(log10(model$w), digits = 2) == 5)
  if(fish_on == TRUE){
    y_down = floor(log10(N_ave[dim(N_ave)[1],][y_d_ref]))
  }else{y_down = floor(log10(N_ave[num_zoo,][which(model$w == 10)]))}
  y_up = ceiling(log10(model$nPP[5]))
  
  ## PLOT PHYTO-ZOO-FISH TOTALS
  par(mar = c(4,4,4,2))
  plot(log10(model$w_phyto), log10(model$nPP), type= "l", col = "green",
       lwd = 2, xlim = c(log10(model$w_phyto[5]), log10(model$w[length(model$w)])), 
       ylim = c(y_down, y_up), 
       xlab = "", 
       ylab ="")
  lines(log10(model$w), log10(tot_zoo), col = "purple", lwd = 3)
  lines(log10(model$w), log10(tot_cephs), col = "yellow", lwd = 3.5)
  if(num_fish > 0){
    lines(log10(model$w), log10(tot_fish), col = "blue", lwd = 3)
  }
  mtext(expression(paste("log"[10], "(Abundance ", "# m "^-3, ")")), side = 2, line = 2.2, cex=1)
  mtext(expression(paste("log"[10], "(Body Weight, g)")), side = 1, line = 2.2, cex = 1)
  legend("bottomleft", legend = c("Phytoplankton", "Zooplankton Community", "Total Cephs", "Fish Community"),
         # col = c("green","red", "purple","blue"), lty = 1, lwd = 2, bty = "n", cex = 0.9)
         col = c("green","purple", "yellow","blue"), lty = 1, lwd = 2, bty = "n", cex = 0.9)
  title(main = paste("Phytoplankton Slope = " , round(modelss$environ$b, digits = 2)))#,
                    # "\nTemperature = ", round(modelss$environ$sst), "C"))
  ## PLOT CEPH GROUPS
  if(num_cephs > 1){
    coll = c("#F8766D", "#00BFC4")
    # coll = rainbow(num_cephs)
    
    for(i in 1:num_cephs){
      colll = coll[i]
      lines(log10(model$w), log10(N_ave[ceph_groups[i],]), lty = 2, col = colll, lwd = 2.5)
    }
    legend("topright", legend = as.character(param$groups$species[ceph_groups]),
           col = coll, lty = 2, lwd = 2.5,bty = "n", cex = 0.8)
  }
  ## PLOT FISH GROUPS
  if(num_fish > 0){
    for(i in 1:num_fish){
      lines(log10(model$w), log10(N_ave[fish_groups[i],]), lty = 2, col =  "blue", lwd = 1.5)
    }
  }
}

Spectrum_Plot_No_Cephs <- function(models, fish_on){
  
  model = models
  modelss = models
  param = modelss$param
  
  N_sav = modelss$N
  N_ave = colMeans(modelss$N[(ceiling(0.5*dim(N_sav)[1])):(dim(N_sav)[1]),,], dim = 1)
  zoo_groups =  param$zoo_grps
  # ceph_groups = param$ceph_grps
  fish_groups = param$fish_grps
  num_zoo = param$num_zoo
  # num_cephs = param$num_ceph
  num_fish = param$num_fish
  
  if(length(zoo_groups) > 1){
    tot_zoo = colSums(N_ave[zoo_groups,])
  }
  
  if(length(zoo_groups) == 1){
    tot_zoo = N_ave[zoo_groups,]  
  }
  
  # if(length(ceph_groups) > 1){
  #   tot_cephs = colSums(N_ave[ceph_groups,])
  # }
  
  # if(length(ceph_groups) == 1){
  #   tot_cephs = N_ave[ceph_groups,]
  # }
  
  if(num_fish > 0){
    if(num_fish > 1){
      tot_fish = colSums(N_ave[fish_groups,])
    }else{tot_fish = N_ave[fish_groups,]}
  }
  
  y_d_ref = which(round(log10(model$w), digits = 2) == 5)
  if(fish_on == TRUE){
    y_down = floor(log10(N_ave[dim(N_ave)[1],][y_d_ref]))
  }else{y_down = floor(log10(N_ave[num_zoo,][which(model$w == 10)]))}
  y_up = ceiling(log10(model$nPP[5]))
  
  ## PLOT PHYTO-ZOO-FISH TOTALS
  par(mar = c(4,4,4,2))
  plot(log10(model$w_phyto), log10(model$nPP), type= "l", col = "green",
       lwd = 2, xlim = c(log10(model$w_phyto[5]), log10(model$w[length(model$w)])), 
       ylim = c(y_down, y_up), 
       xlab = "", 
       ylab ="")
  lines(log10(model$w), log10(tot_zoo), col = "purple", lwd = 3)
  # lines(log10(model$w), log10(tot_cephs), col = "yellow", lwd = 3.5)
  if(num_fish > 0){
    lines(log10(model$w), log10(tot_fish), col = "blue", lwd = 3)
  }
  mtext(expression(paste("log"[10], "(Abundance ", "# m "^-3, ")")), side = 2, line = 2.2, cex=1)
  mtext(expression(paste("log"[10], "(Body Weight, g)")), side = 1, line = 2.2, cex = 1)
  legend("bottomleft", legend = c("Phytoplankton", "Zooplankton Community", "Fish Community"),
         # col = c("green","red", "purple","blue"), lty = 1, lwd = 2, bty = "n", cex = 0.9)
         col = c("green","purple","blue"), lty = 1, lwd = 2, bty = "n", cex = 0.9)
  title(main = paste("Phytoplankton Slope = " , round(modelss$environ$b, digits = 2)))#,
  # "\nTemperature = ", round(modelss$environ$sst), "C"))
  ## PLOT CEPH GROUPS
  # if(num_cephs > 1){
  #   coll = c("#F8766D", "#00BFC4")
  #   # coll = rainbow(num_cephs)
  #   
  #   for(i in 1:num_cephs){
  #     colll = coll[i]
  #     lines(log10(model$w), log10(N_ave[ceph_groups[i],]), lty = 2, col = colll, lwd = 2.5)
  #   }
  #   legend("topright", legend = as.character(param$groups$species[ceph_groups]),
  #          col = coll, lty = 2, lwd = 2.5,bty = "n", cex = 0.8)
  # }
  ## PLOT FISH GROUPS
  if(num_fish > 0){
    for(i in 1:num_fish){
      lines(log10(model$w), log10(N_ave[fish_groups[i],]), lty = 2, col =  "blue", lwd = 1.5)
    }
  }
}

#### PLOT Zoom AVERAGE SPECTRUM
Zoom_Spectrum_Plot <- function(models, fish_on){
  
  model = models
  modelss = models
  param = modelss$param
  
  N_sav = modelss$N
  N_ave = colMeans(modelss$N[(ceiling(0.5*dim(N_sav)[1])):(dim(N_sav)[1]),,], dim = 1)
  zoo_groups =  param$zoo_grps
  ceph_groups = param$ceph_grps
  fish_groups = param$fish_grps
  num_zoo = param$num_zoo
  num_cephs = param$num_ceph
  num_fish = param$num_fish
  
  if(length(zoo_groups) > 1){
    tot_zoo = colSums(N_ave[zoo_groups,])
  }
  
  if(length(zoo_groups) == 1){
    tot_zoo = N_ave[zoo_groups,]  
  }
  
  if(length(ceph_groups) > 1){
    tot_cephs = colSums(N_ave[ceph_groups,])
  }
  
  if(length(ceph_groups) == 1){
    tot_cephs = N_ave[ceph_groups,]
  }
  
  if(num_fish > 0){
    if(num_fish > 1){
      tot_fish = colSums(N_ave[fish_groups,])
    }else{tot_fish = N_ave[fish_groups,]}
  }
  
  y_d_ref = which(round(log10(model$w), digits = 2) == 5)
  if(fish_on == TRUE){
    y_down = floor(log10(N_ave[dim(N_ave)[1],][y_d_ref]))
  }else{y_down = floor(log10(N_ave[num_zoo,][which(model$w == 10)]))}
  y_up = ceiling(log10(model$nPP[5]))
  
  ## PLOT PHYTO-ZOO-FISH TOTALS
  par(mar = c(4,4,4,2))
  plot(log10(model$w_phyto), log10(model$nPP), type= "l", col = "green",
       lwd = 2, xlim = c(-4, 4), 
       ylim = c(y_down, 3), 
       xlab = "", 
       ylab ="")
  lines(log10(model$w), log10(tot_zoo), col = "purple", lwd = 3)
  lines(log10(model$w), log10(tot_cephs), col = "yellow", lwd = 3.5)
  if(num_fish > 0){
    lines(log10(model$w), log10(tot_fish), col = "blue", lwd = 3)
  }
  mtext(expression(paste("log"[10], "(Abundance ", "# m "^-3, ")")), side = 2, line = 2.2, cex=1)
  mtext(expression(paste("log"[10], "(Body Weight, g)")), side = 1, line = 2.2, cex = 1)
  legend("bottomleft", legend = c("Phytoplankton", "Zooplankton Community", "Total Cephs", "Fish Community"),
         col = c("green","purple", "yellow","blue"), lty = 1, lwd = 2, bty = "n", cex = 0.9)
  title(main = paste("Phytoplankton Slope = " , round(modelss$enviro$b, digits = 2)))#,
                     #"\nTemperature = ", round(enviro$sst), "C"))
  ## PLOT CEPH GROUPS
  if(num_cephs > 1){
    # coll = rainbow(num_cephs)
    coll = c("#F8766D", "#00BFC4")
    
    for(i in 1:num_cephs){
      colll = coll[i]
      lines(log10(model$w), log10(N_ave[ceph_groups[i],]), lty = 2, col = colll, lwd = 3)
    }
    legend("topright", legend = as.character(param$groups$species[ceph_groups]),
           col = coll, lty = 2, lwd = 2.5,bty = "n", cex = 1)
  }
  ## PLOT FISH GROUPS
  if(num_fish > 0){
    for(i in 1:num_fish){
      lines(log10(model$w), log10(N_ave[fish_groups[i],]), lty = 2, col = "blue", lwd = 1.5)
    }
  }
}

## PLOT BIOMASS OVER TIME
# Biomass_Plot <- function(models){
#   modelss <- models
#   param <- modelss$param
#   
#   zoo_groups =  param$zoo_grps
#   ceph_groups = param$ceph_grps
#   fish_groups = param$fish_grps
#   num_zoo = param$num_zoo
#   num_cephs = param$num_ceph
#   num_fish = param$num_fish
#   
#   if(num_fish > 0){
#     FB_total = modelss$Biomass[(ceiling(0.5*dim(modelss$Biomass)[1])):(dim(modelss$Biomass)[1]), fish_groups]
#   }
#   Z_groups = modelss$Biomass[(ceiling(0.5*dim(modelss$Biomass)[1])):(dim(modelss$Biomass)[1]), zoo_groups]
#   
#   C_groups = modelss$Biomass[(ceiling(0.5*dim(modelss$Biomass)[1])):(dim(modelss$Biomass)[1]), ceph_groups]
#   
#   if(num_cephs > 1){
#     C_total = rowSums(C_groups)
#   }else{
#     C_total = C_groups
#   }
#   
#   if(num_zoo > 1){
#     ZB_total = rowSums(Z_groups) 
#   }else{
#     ZB_total = Z_groups  
#   }
#   
#   if(num_fish > 0){
#     if(num_fish > 1){
#       FB_total = rowSums(FB_total)
#     }
#   }
#   
#   par(mar = c(4,4,4,10))
#   if(num_fish > 0){
#     plot(FB_total, type = "l", col = "blue", lwd = 2, ylim = c(0, (max(c(FB_total, ZB_total)))), xaxt = "n",
#          ylab = "", xlab = "")
#   }else{
#     plot(0, type = "n", xlim = c(0, length(ZB_total)), ylim = c(0, max(ZB_total)), xaxt = "n",
#          ylab = "", xlab = "")
#   }
#   lines(ZB_total, type = "l", col = "red", lwd = 2)
#   lines(C_total, type = "l", col = "purple", lwd = 2)
#   
#   mtext(text = expression(paste("Total Biomass ", "g m"^-3)), side = 2, line = 2.5)
#   mtext(text = "Time (years)", side = 1, line = 2.5)
#   axis(side = 1, at = seq(1,length(ZB_total), length.out = 10), labels = round(seq(0.1*length(ZB_total)/12.17,length(ZB_total)/12.17,length.out =10)), cex.axis = 1)
#   
#   if(length(zoo_groups) > 1){
#     coll = rainbow(num_zoo)
#     for(i in 1:length(zoo_groups)){
#       lines(Z_groups[,i], lty = 2, col = coll[i], lwd = 1)
#     }
#   }
#   
#   if(length(ceph_groups) > 1){
#     coll = rainbow(num_cephs)
#     for(i in 1:length(ceph_groups)){
#       lines(C_groups[,i], lty = 2, col = coll[i], lwd = 1)
#     }
#   }
#   
#   title(main = paste("Phytoplankton Slope = " , round(enviro$b, digits = 2),
#                      "\nTemperature = ", round(enviro$sst), "C"))
#   par(xpd = TRUE)
#   max_x = length(ZB_total)
#   
#   if(length(ceph_groups) > 1){
#       max_y = max(c(FB_total, ZB_total))
#       legend(x= (max_x + 2), y = max_y, legend = c("All Zooplankton", "All Cephs", "All Fish", as.character(param$groups$species[ceph_groups])),
#              col = c("Red", "Purple", "Blue", coll), lty = c(1,1,1, rep(2, length(ceph_groups))), lwd = c(2,2, rep(1.5, length(ceph_groups))),bty = "n", cex = 0.8)
#       par(xpd = FALSE)
#     }
# }


## PLOT BIOMASS CONTRIBUTIONS
# Bio_Cont_Plot <- function(models){
#   modelss = models
#   param = modelss$params
#   
#   w_min = which(round(log10(model$w), digits = 2) == -6.2)
#   w_max = which(round(log10(model$w), digits = 2) == 4)
#   N_sav = modelss$N
#   N_ave = colMeans(modelss$N[(ceiling(0.5*dim(N_sav)[1])):(dim(N_sav)[1]),,], dim = 1)
#   N_ave = N_ave[,(w_min:w_max)]
#   zoo_groups =  which(is.na(param$groups$prop) == FALSE & rowSums(N_ave) != 0)
#   N_ave = N_ave[zoo_groups,]
#   Biom_ave = sweep(N_ave, 2, model$w[w_min:w_max], "*")
#   Biom_props = sweep(Biom_ave,2, colSums(Biom_ave), "/")
#   
#   par(mar = c(4,4,4,8))
#   x = log10(model$w[w_min:w_max])
#   plot(x, rep(1,length(x)), type = "n", ylim = c(0.01,1), xlab = expression(paste("log"[10], "(Body Weight, g)")),
#        ylab = "% Biomass Contribution")
#   coll = rainbow(length(zoo_groups))
#   for(i in 1:length(zoo_groups)){
#     lines(x, Biom_props[i,], col = coll[i], lwd = 1.5, lty = 1)
#   }
#   title(main = paste("Phytoplankton Slope = " , round(enviro$b, digits = 2),
#                      "\nTemperature = ", round(enviro$sst), "C"))
#   par(xpd = TRUE)
#   max_x = x[length(x)]
#   max_y = 1
#   legend(x= (max_x + 0.1), y = max_y, legend = c(as.character(param$groups$species[zoo_groups])),
#          col = coll, lty = 1, lwd = 1.5,bty = "n", cex = 0.8)
#   par(xpd = FALSE)
# }


Ceph_Abund_Pie <- function(modelss, cut_point1, cut_point2){
  N_sav = modelss$N
  N_ave = colMeans(modelss$N[(ceiling(0.5*dim(N_sav)[1])):(dim(N_sav)[1]),,], dim = 1)
 
  ceph_groups = param$ceph_grps
  num_cephs = param$num_ceph

  weight_cut = which(modelss$w >= 10^cut_point1 & modelss$w <= 10^cut_point2)
  ceph_abunds = rowSums(N_ave[ceph_groups, weight_cut])
  
  par(mfrow = c(1,1))
  
  slices <- ceph_abunds
  lbls <- param$groups$species[ceph_groups]
  pct <- round(slices/sum(slices)*100)
  lbls <- paste(lbls, pct)
  lbls <- paste(lbls, "%", sep = "")
  pie(slices, labels = lbls, col = c("#F8766D", "#00BFC4"), 
  # pie(slices, labels = lbls, col = rainbow(length(lbls)), 
      main = "Ceph Abundances")
}

Ceph_Biom_Pie <- function(modelss, cut_point1, cut_point2){
  N_sav = modelss$N
  N_ave = colMeans(modelss$N[(ceiling(0.5*dim(N_sav)[1])):(dim(N_sav)[1]),,], dim = 1)
  ceph_groups = param$ceph_grps
  num_cephs = param$num_ceph
  
  weight_cut = which(modelss$w >= 10^cut_point1 & modelss$w <= 10^cut_point2)
  ceph_abunds = rowSums(N_ave[ceph_groups, weight_cut])
  weight_cut = which(modelss$w >= 10^cut_point1 & modelss$w <= 10^cut_point2)
  ceph_bioms = sweep(N_ave[ceph_groups,], 2, modelss$w, "*")
  ceph_biomss = rowSums(ceph_bioms[, weight_cut])
  
  par(mfrow = c(1,1))
  slices <- ceph_biomss
  lbls <- param$groups$species[ceph_groups]
  pct <- round(slices/sum(slices)*100)
  lbls <- paste(lbls, pct)
  lbls <- paste(lbls, "%", sep = "")
  pie(slices, labels = lbls, col = c("#F8766D", "#00BFC4"), 
  # pie(slices, labels = lbls, col = rainbow(length(lbls)), 
      main = "Ceph Biomass")
} 



###################################################################################
###################################################################################
##################### THESE FUNCTIONS ARE FOR EMERGENT DIETS OF THE FISH - THEY NEED TO BE 
##################### FIXED UP FOR CEPHS
###################################################################################

# Diet_Matrix <- function(models, params){
#   model = modelss
#   param = params
#   diet = model$diet
#   fish_diet = model$fish_diet
#   ngrps = dim(param$groups)[1]
#   
#   assim_dynam <- matrix(param$groups$alpha*param$nutrition, nrow = ngrps, ncol = ngrps, byrow = TRUE)/
#     matrix(param$nutrition, nrow = ngrps, ncol = ngrps)
#   
#   tot_phyto_prod = ((10^1.58)*(param$environ$chlo^1.29)*1e-2*365)
#   
#   diet_matt = matrix(0, nrow = ngrps, ncol = (ngrps+1))
#   
#   #  fish_s = seq(1,which(round(log10(model$w), digits = 2) == 2),1) # (10^-3g - 10^2g)
#   # fish_l = seq((length(fish_s)+1), length(model$w), 1) # (10^2g - 10^6g)
#   
#   # fish_diet_mat = colMeans(fish_diet[ceiling(0.5*dim(fish_diet)[1]):dim(fish_diet)[1],,], dim = 1)
#   diet_mat = colMeans(diet[ceiling(0.5*dim(diet)[1]):dim(diet)[1],,,], dim = 1)
#   diet_mat[1:ngrps, ,2:(ngrps+1)] <- sweep(diet_mat[1:ngrps,,2:(ngrps+1)], c(1,3), assim_dynam, "/")
#   diet_mat = aperm(diet_mat, c(2,1,3))
#   
#   diet_matt[(1:ngrps),1:(ngrps+1)] = colSums(diet_mat)
#   
#   ### DON'T NEED IF WE HAVE SMALL, MEDIUM, LARGE FISH
#   # diet_matt[10,1:11] = colSums(diet_mat[fish_s,10,])
#   #  diet_matt[11,1:11] = colSums(diet_mat[fish_l,10,])
#   #  diet_matt[10,11] = fish_diet_mat[1,1]
#   # diet_matt[11,11] = fish_diet_mat[2,1]
#   #  diet_matt[11,12] = fish_diet_mat[2,2]
#   
#   #  diet_mat[11,,] = diet_mat[10,,]*matrix(fish_large, nrow = length(fish_large), ncol = 11)  
#   
#   diet_matt = round(diet_matt, digits = 6)*1000
#   group_names = as.character(param$groups$species)
#   rownames(diet_matt) = as.character(param$groups$species)
#   colnames(diet_matt) = c("Phytoplankton", group_names)
#   return(list(diet_matt, tot_phyto_prod))
#   
# }
# 
# ## FISH DIET PLOTS
# Fish_Diet <- function(models, params){
#   modelss = models
#   param = params
#   
#   zoo_groups =  param$zoo_grps
#   fish_groups = param$fish_grps
#   num_zoo = param$num_zoo
#   num_fish = param$num_fish
#   
#   ngrps = dim(param$groups)[1]
#   w_min = which(log10(model$w) == -3)
#   w_max = which(log10(model$w) == 6)
#   assim_dynam <- matrix(param$groups$alpha*param$nutrition, nrow = ngrps, ncol = ngrps, byrow = TRUE)/
#     matrix(param$nutrition, nrow = ngrps, ncol = ngrps)
#   diet_sav = modelss$diet
#   diet_ave = colMeans(diet_sav[(ceiling(0.5*dim(diet_sav)[1])):(dim(diet_sav)[1]),,,], dim = 1)
#   diet_ave = diet_ave[,(w_min:w_max),]
#   fish_ident = which(is.na(param$groups$prop) == TRUE)
#   fish_diet = diet_ave[fish_ident,,]
#   fish_diet_ingest = sweep(fish_diet, 2, c(1,assim_dynam[ngrps,]), "/")
#   fish_diet_props = sweep(fish_diet_ingest, 1, rowSums(fish_diet_ingest), "/")
#   fish_diet_props = fish_diet_props[,-1]
#   zoo_groups =  which(apply(fish_diet_props, 2, max) > 0.03)
#   zoo_groups = zoo_groups[which(is.na(param$groups$prop[zoo_groups]) == FALSE )]
#   fish_diet_dom_zoo = fish_diet_props
#   zoo_names = as.character(param$groups$species[1:num_zoo])  
#   
#   
#   par(mar = c(4,4,4,2))
#   x = log10(model$w[w_min:w_max])
#   plot(x[-length(x)], rep(1, length(x[-length(x)])), type = "n", ylim = c(0, ceiling(max(fish_diet_dom_zoo[,1:9], na.rm = TRUE)*100)/100),
#        xlab = expression(paste("log"[10], "(Body Weight, g)")), ylab = "% of Diet")
#   
#   coll = rainbow(9)
#   for(i in 1:9){
#     lines(x[-length(x)], fish_diet_dom_zoo[-length(x),i], lty = 1, lwd = 2, col = coll[i])
#   }
#   legend("topright", legend = zoo_names, lty = 1, col = coll,
#          lwd = 2, bty = "n", cex = 0.8)
#   title(main = paste("Phytoplankton Slope = " , round(enviro$b, digits = 2),
#                      "\nTemperature = ", round(enviro$sst), "C"))
# }
# 
# 

#####################################################################################################

#####################################################################################################





# Summary_Results <- function(models){
#   
#    modelss = models
#    param = modelss$param
#   
#   N_sav = modelss$N
#   N_ave = colMeans(modelss$N[(ceiling(0.5*dim(N_sav)[1])):(dim(N_sav)[1]),,], dim = 1)
#   zoo_groups =  param$zoo_grps    #zoo_groups =  which(is.na(param$groups$prop) == FALSE)
#   num_zoo = length(zoo_groups)
#   tot_zoo = sum(N_ave[zoo_groups,])
#   ceph_groups = param$ceph_grps   # New line
#   num_cephs = length(ceph_groups) # New line
#   tot_ceph = colSums(N_ave[ceph_groups,]) # New line
#   fish_groups = param$fish_grps   # New line
#   num_fish = length(fish_groups)  # New line
#   tot_fish = sum(N_ave[fish_groups,]) # tot_fish = sum(N_ave[-zoo_groups])
#   
#   ## BIOMASS OF DIFFERENT GROUPS
#   zoo_bioms = sum(N_ave[zoo_groups,]*matrix(modelss$w, nrow = length(zoo_groups), 
#                                                 ncol = dim(N_ave)[2], byrow = TRUE))
#   ceph_bioms = rowSums(N_ave[ceph_groups,]*matrix(modelss$w, nrow = length(ceph_groups), 
#                                                 ncol = dim(N_ave)[2], byrow = TRUE))
#   fish_s = seq(1,which(round(log10(modelss$w), digits = 2) == 2),1) # (10^-3g - 10^2g)
#   fish_l = seq((length(fish_s)+1), length(modelss$w), 1) # (10^2g - 10^6g)
#   tot_small_fish = sum(N_ave[fish_groups, fish_s]*modelss$w[fish_s]) # tot_small_fish = sum(N_ave[-zoo_groups, fish_s]*modelss$w[fish_s])
#   tot_large_fish =  sum(N_ave[fish_groups, fish_l]*modelss$w[fish_l])
#   zoo_specs = as.character(param$groups$species)[zoo_groups]
#   ceph_specs = as.character(param$groups$species)[ceph_groups]
#   
#   
#   tot_bioms =  matrix(c(zoo_bioms, ceph_bioms, tot_small_fish, tot_large_fish), nrow = 8, ncol = 1)
#   rownames(tot_bioms) = c(zoo_specs, ceph_specs, "Small Fish", "Large Fish")
#   colnames(tot_bioms) = "Biomass g m^-3"
#   tot_bioms = round(tot_bioms, digits = 4)
#   
#   ## Calculate information tables for zooplankton and fish communities
#   zoo_groups =  param$zoo_grps
#   num_zoo = length(zoo_groups)
#   
#   FB_total = modelss$Biomass[(ceiling(0.5*dim(modelss$Biomass)[1])):(dim(modelss$Biomass)[1]), fish_groups]
#   ZB_total = modelss$Biomass[(ceiling(0.5*dim(modelss$Biomass)[1])):(dim(modelss$Biomass)[1]), zoo_groups]
#   #ZB_total = sum(Z_groups) 
#   C_groups = modelss$Biomass[(ceiling(0.5*dim(modelss$Biomass)[1])):(dim(modelss$Biomass)[1]), ceph_groups]
#   CB_total = rowSums(C_groups)
#   
#   # Average Biomass
#   ave_zoo = round(mean(ZB_total), digits = 2)
#   ave_ceph = round(mean(CB_total), digits = 2)
#   ave_fish = round(mean(FB_total), digits = 2)
#   ave_phyto = round(sum(model$nPP*model$w_phyto), digits = 2)
#   
#   # Average Slope
#   fish_start = which(model$w == 10^(param$groups$W0[dim(param$groups)[1]]))
#   fish_finish = which(model$w == 10^(param$groups$Wmat[dim(param$groups)[1]])) 
#   max_phyto = round(log10(param$wMax_phyto), digits = 2)
#   zoo_start =  which(round(log10(model$w), digits = 2) == -12)
#   zoo_finish = which(round(log10(model$w), digits = 2) == -2) 
#   zoo_slope2 = round(lm(log10(tot_zoo[zoo_start:zoo_finish])~log10(model$w[zoo_start:zoo_finish]))$
#                        coefficients[2], digits = 2)
#   ceph_start =  which(round(log10(model$w), digits = 2) == -2)
#   ceph_finish = which(round(log10(model$w), digits = 2) == 3) 
#   ceph_slope2 = round(lm(log10(tot_ceph[ceph_start:ceph_finish])~log10(model$w[ceph_start:ceph_finish]))$
#                        coefficients[2], digits = 2)
#   # zoo_slope =  round((log10(tot_zoo[zoo_finish]) - log10(tot_zoo[zoo_start]))/
#   #                      (log10(model$w[zoo_finish]) - log10(model$w[zoo_start])), digits = 2)
#   fish_slope2 = round(lm(log10(N_ave[dim(N_ave)[1], c(fish_start:fish_finish)])~
#                            log10(model$w[c(fish_start:fish_finish)]))$coefficients[2], digits = 2)
#   #fish_slope2 = round((log10(N_ave[dim(N_ave)[1],fish_finish]) - log10(N_ave[dim(N_ave)[1],fish_start]))/
#   #                      (log10(model$w[fish_finish]) - log10(model$w[fish_start])), digits = 2)
#   phyto_slope = param$environ$b
#   
#   # Proportions and PPMR
#   prop_zoo_biom = round(colMeans(Z_groups)/sum(colMeans(Z_groups)), digits = 2) # Biomass Proportions
#   zoo_ave_m = round(sum(colMeans(Z_groups)/sum(colMeans(Z_groups))* # Average m-value
#                           param$groups$m[zoo_groups]), digits = 2)
#   
#   prop_ceph_biom = round(colMeans(C_groups)/sum(colMeans(C_groups)), digits = 2) # Biomass Proportions
#   ceph_ave_m = round(sum(colMeans(C_groups)/sum(colMeans(C_groups))* # Average m-value
#                           param$groups$m[ceph_groups]), digits = 2)
#   
#   NN = rowSums(modelss$N[,,which(round(log10(modelss$w), digits = 2) == -4.5):dim(modelss$N)[3]], dim = 2) ## What is this?
#   
#   prop_zoo_abund =  round(colMeans(NN[(ceiling(0.5*dim(NN)[1])):(dim(NN)[1]), zoo_groups])/
#                             sum(colMeans(NN[(ceiling(0.5*dim(NN)[1])):(dim(NN)[1]), zoo_groups])), digits = 2)# Abundance Proportions
#   
#   prop_ceph_abund =  round(colMeans(NN[(ceiling(0.5*dim(NN)[1])):(dim(NN)[1]), ceph_groups])/
#                             sum(colMeans(NN[(ceiling(0.5*dim(NN)[1])):(dim(NN)[1]), ceph_groups])), digits = 2)# Abundance Proportions
#   
#   results <- matrix(c(ave_phyto, ave_zoo, ave_ceph, ave_fish, round(param$environ$b, digits = 2), zoo_slope2, ceph_slope2, fish_slope2, NA, zoo_ave_m, ceph_ave_m, NA),
#                     nrow = 4, ncol = 3)
#   colnames(results) <- c("Biomass", "Slope", "m")
#   rownames(results) <- c("Phytoplankton","Zooplankton", "Cephalopods", "Fish")
#   
#   zoo_results <- matrix(c(prop_zoo_abund, prop_zoo_biom), nrow = num_zoo, ncol = 2)
#   colnames(zoo_results) <- c("% Abundance Total", "% Biomass Total")
#   rownames(zoo_results) <- as.character(param$groups$species[zoo_groups])
#   
#   ceph_results <- matrix(c(prop_ceph_abund, prop_ceph_biom), nrow = num_cephs, ncol = 2)
#   colnames(ceph_results) <- c("% Abundance Total", "% Biomass Total")
#   rownames(ceph_results) <- as.character(param$groups$species[ceph_groups])
#   
#   return(list("General" = results, "Zoo_Groups" = zoo_results, "Ceph_Croups" = ceph_results, "Total_Bioms" = tot_bioms))
# }

TL_Allom_Plot_Fixed_Scaling_dx_reduced <- function(projection){
  
  
  model <- projection
  
  # Average diet array for each predator size class for the final quarter of the projection
  dynamdiet = apply(model$dynam_diet_full[c(floor(0.75*dim(model$dynam_diet_full)[1]):dim(model$dynam_diet_full)[1]),,,,], c(2,3,4,5), mean)
  phytodiet = apply(model$phyto_diet_full[c(floor(0.75*dim(model$phyto_diet_full)[1]):dim(model$phyto_diet_full)[1]),,], c(2,3), mean)
  
  phyto_tl <- 1 # TL for phytoplankton is 1
  start_dynam_tl <- matrix(2, nrow = dim(dynamdiet)[1], ncol = dim(dynamdiet)[2]) # Start TL for dynamic size groups - start at 2 for all groups and sizes
  
  curr_phyto_diet <- phytodiet # Current phyto diet
  curr_dynam_diet <- dynamdiet # Current heterotroph diet
  
  total_diet <- curr_phyto_diet + apply(curr_dynam_diet, c(1,2), sum) # Total consumption, in grams wet weight, by pred group and pred size
  
  curr_phyto_frac <- curr_phyto_diet/total_diet # Fraction of diet from phyto, by pred group and pred sizes
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
    # filter(Predator_Size_Class == 111 | Predator_Size_Class == 112) %>%
    filter(Predator_Size_Class == 220 | Predator_Size_Class == 221) %>%
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
    # filter(Predator_Size_Class == 111 | Predator_Size_Class == 112) %>%
    filter(Predator_Size_Class == 220 | Predator_Size_Class == 221) %>%
    mutate(avg_min_RTP = mean(RTP))
  
  min_mod_RTP_inactive <- df_model_min_RTP_inactive$avg_min_RTP[1]
  
  rescale_RTP_inactive_val <- min_mod_RTP_inactive - min_emp_RTP_inactive
  
  pos_val_inactive <- sqrt(rescale_RTP_inactive_val * rescale_RTP_inactive_val)
  
  rescale_val <- (pos_val_inactive + pos_val_active)/2
  
  df_emp_inactive$RTP <- df_emp_inactive$RTP + rescale_val
  
  df_emp_active$RTP <- df_emp_active$RTP + rescale_val
  
  df_zoom <- df_combined %>%
    filter(Predator == "Inactive Cephs" | Predator ==  "Active Cephs")
  
p1 <- df_combined %>%
    mutate(name = fct_relevel(Predator, 'Active Cephs', 'Inactive Cephs','Zooplankton','Fish')) %>%
    ggplot(aes(x=Predator_Size, y= RTP, colour=name)) +
    scale_colour_manual(values = c( "#F8766D", "#00BFC4","#C77CFF", "#7CAE00")) +
    geom_smooth(method = "lm", size=1.5, alpha = 0, linetype = "dashed") +
    #facet_wrap(~Species_Order, scales = "free_x") +
    geom_point(aes(x=Predator_Size, y= RTP, colour=name),data = df_emp_active, size = .5) +
    geom_smooth(aes(x=Predator_Size, y= RTP, colour=name), data = df_emp_active,  method = "lm", size = 1, alpha = 0.5) +
    
    geom_point(aes(x=Predator_Size, y= RTP, colour=name),data = df_emp_inactive, size = .5) +
    geom_smooth(aes(x=Predator_Size, y= RTP, colour=name), data = df_emp_inactive,  method = "lm", size = 1, alpha = 0.5) +
    scale_x_log10() +
    ggtitle(label = "Community-level RTP Fixed-rescaling: Ceph FGs empirical fit (dots and solid lines) vs model output (dashed lines)") +
    theme_bw()

p2 <- df_zoom %>%
  mutate(name = fct_relevel(Predator, 'Active Cephs', 'Inactive Cephs')) %>%
  ggplot(aes(x=Predator_Size, y= RTP, colour=name)) +
  scale_colour_manual(values = c("#F8766D","#00BFC4")) +
  geom_smooth(method = "lm", size=1.5, alpha = 0, linetype = "dashed") +
  #facet_wrap(~Species_Order, scales = "free_x") +
  geom_point(aes(x=Predator_Size, y= RTP, colour=name),data = df_emp_active, size = .5) +
  geom_smooth(aes(x=Predator_Size, y= RTP, colour=name), data = df_emp_active,  method = "lm", size = 1, alpha = 0.25) +
  geom_point(aes(x=Predator_Size, y= RTP, colour=name),data = df_emp_inactive, size = .5) +
  geom_smooth(aes(x=Predator_Size, y= RTP, colour=name), data = df_emp_inactive,  method = "lm", size = 1, alpha = 0.25) +
  scale_x_log10() +
  ggtitle(label = "Ceph FGs RTP Fixed-rescaling: Ceph FGs empirical fit (dots and solid lines) vs model output (dashed lines)") +
  theme_bw()

list(p1, p2)
  
}

TL_Allom_Plot_Fixed_Scaling<- function(projection){
  
  
  model <- projection
  
  # Average diet array for each predator size class for the final quarter of the projection
  dynamdiet = apply(model$dynam_diet_full[c(floor(0.75*dim(model$dynam_diet_full)[1]):dim(model$dynam_diet_full)[1]),,,,], c(2,3,4,5), mean)
  phytodiet = apply(model$phyto_diet_full[c(floor(0.75*dim(model$phyto_diet_full)[1]):dim(model$phyto_diet_full)[1]),,], c(2,3), mean)
  
  phyto_tl <- 1 # TL for phytoplankton is 1
  start_dynam_tl <- matrix(2, nrow = dim(dynamdiet)[1], ncol = dim(dynamdiet)[2]) # Start TL for dynamic size groups - start at 2 for all groups and sizes
  
  curr_phyto_diet <- phytodiet # Current phyto diet
  curr_dynam_diet <- dynamdiet # Current heterotroph diet
  
  total_diet <- curr_phyto_diet + apply(curr_dynam_diet, c(1,2), sum) # Total consumption, in grams wet weight, by pred group and pred size
  
  curr_phyto_frac <- curr_phyto_diet/total_diet # Fraction of diet from phyto, by pred group and pred sizes
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
  
  rescale_val <- (pos_val_inactive + pos_val_active)/2
  
  df_emp_inactive$RTP <- df_emp_inactive$RTP + rescale_val
  
  df_emp_active$RTP <- df_emp_active$RTP + rescale_val
  
  df_zoom <- df_combined %>%
    filter(Predator == "Inactive Cephs" | Predator ==  "Active Cephs")
  
  p1 <- df_combined %>%
    mutate(name = fct_relevel(Predator, 'Active Cephs', 'Inactive Cephs','Zooplankton','Fish')) %>%
    ggplot(aes(x=Predator_Size, y= RTP, colour=name)) +
    scale_colour_manual(values = c( "#F8766D", "#00BFC4","#C77CFF", "#7CAE00")) +
    geom_smooth(method = "lm", size=1.5, alpha = 0, linetype = "dashed") +
    #facet_wrap(~Species_Order, scales = "free_x") +
    geom_point(aes(x=Predator_Size, y= RTP, colour=name),data = df_emp_active, size = .5) +
    geom_smooth(aes(x=Predator_Size, y= RTP, colour=name), data = df_emp_active,  method = "lm", size = 1, alpha = 0.5) +
    
    geom_point(aes(x=Predator_Size, y= RTP, colour=name),data = df_emp_inactive, size = .5) +
    geom_smooth(aes(x=Predator_Size, y= RTP, colour=name), data = df_emp_inactive,  method = "lm", size = 1, alpha = 0.5) +
    scale_x_log10() +
    ggtitle(label = "Community-level RTP Fixed-rescaling: Ceph FGs empirical fit (dots and solid lines) vs model output (dashed lines)") +
    theme_bw()
  
  p2 <- df_zoom %>%
    mutate(name = fct_relevel(Predator, 'Active Cephs', 'Inactive Cephs')) %>%
    ggplot(aes(x=Predator_Size, y= RTP, colour=name)) +
    scale_colour_manual(values = c("#F8766D","#00BFC4")) +
    geom_smooth(method = "lm", size=1.5, alpha = 0, linetype = "dashed") +
    #facet_wrap(~Species_Order, scales = "free_x") +
    geom_point(aes(x=Predator_Size, y= RTP, colour=name),data = df_emp_active, size = .5) +
    geom_smooth(aes(x=Predator_Size, y= RTP, colour=name), data = df_emp_active,  method = "lm", size = 1, alpha = 0.25) +
    geom_point(aes(x=Predator_Size, y= RTP, colour=name),data = df_emp_inactive, size = .5) +
    geom_smooth(aes(x=Predator_Size, y= RTP, colour=name), data = df_emp_inactive,  method = "lm", size = 1, alpha = 0.25) +
    scale_x_log10() +
    ggtitle(label = "Ceph FGs RTP Fixed-rescaling: Ceph FGs empirical fit (dots and solid lines) vs model output (dashed lines)") +
    theme_bw()
  
  list(p1, p2)
  
}

####
TL_Allom_Plot_Fixed_Scaling_Simple<- function(projection){
  
  
  model <- projection
  
  # Average diet array for each predator size class for the final quarter of the projection
  dynamdiet = apply(model$dynam_diet_full[c(floor(0.75*dim(model$dynam_diet_full)[1]):dim(model$dynam_diet_full)[1]),,,,], c(2,3,4,5), mean)
  phytodiet = apply(model$phyto_diet_full[c(floor(0.75*dim(model$phyto_diet_full)[1]):dim(model$phyto_diet_full)[1]),,], c(2,3), mean)
  
  phyto_tl <- 1 # TL for phytoplankton is 1
  start_dynam_tl <- matrix(2, nrow = dim(dynamdiet)[1], ncol = dim(dynamdiet)[2]) # Start TL for dynamic size groups - start at 2 for all groups and sizes
  
  curr_phyto_diet <- phytodiet # Current phyto diet
  curr_dynam_diet <- dynamdiet # Current heterotroph diet
  
  total_diet <- curr_phyto_diet + apply(curr_dynam_diet, c(1,2), sum) # Total consumption, in grams wet weight, by pred group and pred size
  
  curr_phyto_frac <- curr_phyto_diet/total_diet # Fraction of diet from phyto, by pred group and pred sizes
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
  
  
  
  df_combined %>%
    mutate(name = fct_relevel(Predator, 'Zooplankton','Fish')) %>%
    ggplot(aes(x=Predator_Size, y= RTP, colour=name)) +
    scale_colour_manual(values = c( "#C77CFF", "#7CAE00")) +
    geom_smooth(method = "lm", size=1.5, alpha = 0, linetype = "dashed") +
    #facet_wrap(~Species_Order, scales = "free_x") +
    # geom_point(aes(x=Predator_Size, y= RTP, colour=name),data = df_emp_active, size = .5) +
    # geom_smooth(aes(x=Predator_Size, y= RTP, colour=name), data = df_emp_active,  method = "lm", size = 1, alpha = 0.5) +
    # 
    # geom_point(aes(x=Predator_Size, y= RTP, colour=name),data = df_emp_inactive, size = .5) +
    # geom_smooth(aes(x=Predator_Size, y= RTP, colour=name), data = df_emp_inactive,  method = "lm", size = 1, alpha = 0.5) +
    scale_x_log10() +
    ggtitle(label = "Community-level Emergent RTP: Simple Model output (dashed lines)") +
    theme_bw()
  
  # p2 <- df_zoom %>%
  #   mutate(name = fct_relevel(Predator, 'Active Cephs', 'Inactive Cephs')) %>%
  #   ggplot(aes(x=Predator_Size, y= RTP, colour=name)) +
  #   scale_colour_manual(values = c("#F8766D","#00BFC4")) +
  #   geom_smooth(method = "lm", size=1.5, alpha = 0, linetype = "dashed") +
  #   #facet_wrap(~Species_Order, scales = "free_x") +
  #   geom_point(aes(x=Predator_Size, y= RTP, colour=name),data = df_emp_active, size = .5) +
  #   geom_smooth(aes(x=Predator_Size, y= RTP, colour=name), data = df_emp_active,  method = "lm", size = 1, alpha = 0.25) +
  #   geom_point(aes(x=Predator_Size, y= RTP, colour=name),data = df_emp_inactive, size = .5) +
  #   geom_smooth(aes(x=Predator_Size, y= RTP, colour=name), data = df_emp_inactive,  method = "lm", size = 1, alpha = 0.25) +
  #   scale_x_log10() +
  #   ggtitle(label = "Ceph FGs RTP Fixed-rescaling: Ceph FGs empirical fit (dots and solid lines) vs model output (dashed lines)") +
  #   theme_bw()
  # 
  # list(p1, p2)
  
}
####

TL_Allom_Plot_Fixed_Scaling_from_df<- function(df){
  
  
  df_combined <- df
  
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
  
  rescale_val <- (pos_val_inactive + pos_val_active)/2
  
  df_emp_inactive$RTP <- df_emp_inactive$RTP + rescale_val
  
  df_emp_active$RTP <- df_emp_active$RTP + rescale_val
  
  df_zoom <- df_combined %>%
    filter(Predator == "Inactive Cephs" | Predator ==  "Active Cephs")
  
  p1 <- df_combined %>%
    mutate(name = fct_relevel(Predator, 'Active Cephs', 'Inactive Cephs','Zooplankton','Fish')) %>%
    ggplot(aes(x=Predator_Size, y= RTP, colour=name)) +
    scale_colour_manual(values = c( "#F8766D", "#00BFC4","#C77CFF", "#7CAE00")) +
    geom_smooth(method = "lm", size=1.5, alpha = 0, linetype = "dashed") +
    #facet_wrap(~Species_Order, scales = "free_x") +
    geom_point(aes(x=Predator_Size, y= RTP, colour=name),data = df_emp_active, size = .5) +
    geom_smooth(aes(x=Predator_Size, y= RTP, colour=name), data = df_emp_active,  method = "lm", size = 1, alpha = 0.5) +
    
    geom_point(aes(x=Predator_Size, y= RTP, colour=name),data = df_emp_inactive, size = .5) +
    geom_smooth(aes(x=Predator_Size, y= RTP, colour=name), data = df_emp_inactive,  method = "lm", size = 1, alpha = 0.5) +
    scale_x_log10() +
    ggtitle(label = "Community-level RTP Fixed-rescaling: Ceph FGs empirical fit (dots and solid lines) vs model output (dashed lines)") +
    theme_bw()
  
  p2 <- df_zoom %>%
    mutate(name = fct_relevel(Predator, 'Active Cephs', 'Inactive Cephs')) %>%
    ggplot(aes(x=Predator_Size, y= RTP, colour=name)) +
    scale_colour_manual(values = c("#F8766D","#00BFC4")) +
    geom_smooth(method = "lm", size=1.5, alpha = 0, linetype = "dashed") +
    #facet_wrap(~Species_Order, scales = "free_x") +
    geom_point(aes(x=Predator_Size, y= RTP, colour=name),data = df_emp_active, size = .5) +
    geom_smooth(aes(x=Predator_Size, y= RTP, colour=name), data = df_emp_active,  method = "lm", size = 1, alpha = 0.25) +
    geom_point(aes(x=Predator_Size, y= RTP, colour=name),data = df_emp_inactive, size = .5) +
    geom_smooth(aes(x=Predator_Size, y= RTP, colour=name), data = df_emp_inactive,  method = "lm", size = 1, alpha = 0.25) +
    scale_x_log10() +
    ggtitle(label = "Ceph FGs RTP Fixed-rescaling: Ceph FGs empirical fit (dots and solid lines) vs model output (dashed lines)") +
    theme_bw()
  
  list(p1, p2)
  
}

####
TL_Allom_Fixed_Scaling_df<- function(projection){
  
  
  model <- projection
  
  # Average diet array for each predator size class for the final quarter of the projection
  dynamdiet = apply(model$dynam_diet_full[c(floor(0.5*dim(model$dynam_diet_full)[1]):dim(model$dynam_diet_full)[1]),,,,], c(2,3,4,5), mean)
  phytodiet = apply(model$phyto_diet_full[c(floor(0.5*dim(model$phyto_diet_full)[1]):dim(model$phyto_diet_full)[1]),,], c(2,3), mean)
  
  phyto_tl <- 1 # TL for phytoplankton is 1
  start_dynam_tl <- matrix(2, nrow = dim(dynamdiet)[1], ncol = dim(dynamdiet)[2]) # Start TL for dynamic size groups - start at 2 for all groups and sizes
  
  curr_phyto_diet <- phytodiet # Current phyto diet
  curr_dynam_diet <- dynamdiet # Current heterotroph diet
  
  total_diet <- curr_phyto_diet + apply(curr_dynam_diet, c(1,2), sum) # Total consumption, in grams wet weight, by pred group and pred size
  
  curr_phyto_frac <- curr_phyto_diet/total_diet # Fraction of diet from phyto, by pred group and pred sizes
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
  
  # # Read in empirical trophic level data
  # df_emp_SIA <- read_csv("df_emp_SIA.csv")
  # # Make copy
  # df_emp_SIA_rel <- df_emp_SIA
  # 
  # #####
  # # Rescale empirical RTP for active cephs
  # df_emp_active <- df_emp_SIA_rel %>%
  #   filter(name == "Active Cephs")
  # 
  # model_lm_emp_active <- lm(RTP ~ log10(Predator_Size), data = df_emp_active)
  # 
  # df_active <- df_combined %>%
  #   filter(Predator == "Active Cephs")
  # 
  # min_emp_RTP_active <- model_lm_emp_active$coefficients[1] + (model_lm_emp_active$coefficients[2] * log10(min(df_emp_active$Predator_Size)))
  # 
  # model_lm_active <- lm(RTP ~ log10(Predator_Size), data = df_active)
  # 
  # df_model_min_RTP <- df_active %>%
  #   filter(Predator_Size_Class == 111 | Predator_Size_Class == 112) %>%
  #   # filter(Predator_Size_Class == 220 | Predator_Size_Class == 221) %>%
  #   mutate(avg_min_RTP = mean(RTP))
  # 
  # min_mod_RTP_active <- df_model_min_RTP$avg_min_RTP[1]
  # 
  # rescale_RTP_active_val <- min_mod_RTP_active - min_emp_RTP_active
  # 
  # pos_val_active <- sqrt(rescale_RTP_active_val * rescale_RTP_active_val)
  # 
  # 
  # #####
  # # Rescale empirical RTP for inactive group
  # df_emp_inactive <- df_emp_SIA_rel %>%
  #   filter(name == "Inactive Cephs")
  # 
  # df_inactive <- df_combined %>%
  #   filter(Predator == "Inactive Cephs")
  # 
  # 
  # model_lm_emp_inactive <- lm(RTP ~ log10(Predator_Size), data = df_emp_inactive)
  # 
  # min_emp_RTP_inactive <- model_lm_emp_inactive$coefficients[1] + (model_lm_emp_inactive$coefficients[2] * log10(min(df_emp_inactive$Predator_Size)))
  # 
  # 
  # df_model_min_RTP_inactive <- df_inactive %>%
  #   filter(Predator_Size_Class == 111 | Predator_Size_Class == 112) %>%
  #   # filter(Predator_Size_Class == 220 | Predator_Size_Class == 221) %>%
  #   mutate(avg_min_RTP = mean(RTP))
  # 
  # min_mod_RTP_inactive <- df_model_min_RTP_inactive$avg_min_RTP[1]
  # 
  # rescale_RTP_inactive_val <- min_mod_RTP_inactive - min_emp_RTP_inactive
  # 
  # pos_val_inactive <- sqrt(rescale_RTP_inactive_val * rescale_RTP_inactive_val)
  # 
  # rescale_val <- (pos_val_inactive + pos_val_active)/2
  # 
  # df_emp_inactive$RTP <- df_emp_inactive$RTP + rescale_val
  # 
  # df_emp_active$RTP <- df_emp_active$RTP + rescale_val
  
  return(df_combined)
  
  
}
####
Total_Consumption <- function(projection){
  
  model <- projection
  
  # Average diet array for each predator size class for the final quarter of the projection
  dynamdiet = apply(model$dynam_diet_full[c(floor(0.75*dim(model$dynam_diet_full)[1]):dim(model$dynam_diet_full)[1]),,,,], c(2,3,4,5), mean)
  phytodiet = apply(model$phyto_diet_full[c(floor(0.75*dim(model$phyto_diet_full)[1]):dim(model$phyto_diet_full)[1]),,], c(2,3), mean)
  
  phyto_tl <- 1 # TL for phytoplankton is 1
  start_dynam_tl <- matrix(2, nrow = dim(dynamdiet)[1], ncol = dim(dynamdiet)[2]) # Start TL for dynamic size groups - start at 2 for all groups and sizes
  
  curr_phyto_diet <- phytodiet # Current phyto diet
  curr_dynam_diet <- dynamdiet # Current heterotroph diet
  
  total_diet <- curr_phyto_diet + apply(curr_dynam_diet, c(1,2), sum) # Total consumption, in grams wet weight, by pred group and pred size
  
  # Create long form diet dataframe
  total_diet_long <- melt(total_diet)

  # Assign names to all variables in df_diet
  names(total_diet_long) <- c("Predator_Species", "Predator_Size_Class", "Consumption")

  # make predators and prey factor variables
  total_diet_long$Predator_Species <- as.factor(total_diet_long$Predator_Species)

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
  df_combined <- left_join(total_diet_long, df_size_classes,
                           by="Predator_Size_Class") %>%

    left_join(df_species_predator,
              by="Predator_Species")
  
  return(df_combined)
  
  
}

####

N_to_Biomass <- function(projection){
  
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
  
  p1 <- ggplot(df_Species_Biomass, aes(x=Time, y= Annual_Biomass, colour=Species_Order)) +
    geom_point() +
    scale_colour_manual(values = c("#F8766D", "#00BFC4", "#7CAE00", "#C77CFF")) +
    xlab("Year") +
    facet_wrap(~Species_Order) +
    ylab(expression(Biomass~(g~m^3))) +
    geom_smooth(method = "loess", alpha=F) +
    #geom_smooth(method = "lm", alpha=F, linetype = "dashed") +
    #scale_y_log10() +
    # ggtitle(paste("Species Biomass Time Series: Model", scenario)) +
    theme_bw() +
    theme(legend.position = "none") +
    theme(plot.title = element_text(size = 9))
  
 p2 <-  ggplot(df_Species_Biomass, aes(x=Time, y= Annual_Biomass, colour=Species_Order)) +
    geom_point() +
    scale_colour_manual(values = c("#F8766D", "#00BFC4", "#7CAE00", "#C77CFF")) +
    xlab("Year") +
    facet_wrap(~Species_Order, scales = "free_y") +
    ylab(expression(Biomass~(g~m^3))) +
    geom_smooth(method = "loess", alpha=F) +
    #geom_smooth(method = "lm", alpha=F, linetype = "dashed") +
    #scale_y_log10() +
    # ggtitle(paste("Species Biomass Time Series: Model", scenario)) +
    theme_bw() +
    theme(legend.position = "none") +
    theme(plot.title = element_text(size = 9))
 
 list(p1, p2)

}
###
N_to_Biomass_Simple <- function(projection){
  
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
                                             levels = c(
                                                      'Zooplankton',
                                                        'Fish'),
                                             ordered = TRUE)
  
  p1 <- ggplot(df_Species_Biomass, aes(x=Time, y= Annual_Biomass, colour=Species_Order)) +
    geom_point() +
    scale_colour_manual(values = c( "#7CAE00", "#C77CFF")) +
    xlab("Year") +
    facet_wrap(~Species_Order) +
    ylab(expression(Biomass~(g~m^3))) +
    geom_smooth(method = "loess", alpha=F) +
    #geom_smooth(method = "lm", alpha=F, linetype = "dashed") +
    #scale_y_log10() +
    # ggtitle(paste("Species Biomass Time Series: Model", scenario)) +
    theme_bw() +
    theme(legend.position = "none") +
    theme(plot.title = element_text(size = 9))
  
  p2 <-  ggplot(df_Species_Biomass, aes(x=Time, y= Annual_Biomass, colour=Species_Order)) +
    geom_point() +
    scale_colour_manual(values = c("#7CAE00", "#C77CFF")) +
    xlab("Year") +
    facet_wrap(~Species_Order, scales = "free_y") +
    ylab(expression(Biomass~(g~m^3))) +
    geom_smooth(method = "loess", alpha=F) +
    #geom_smooth(method = "lm", alpha=F, linetype = "dashed") +
    #scale_y_log10() +
    # ggtitle(paste("Species Biomass Time Series: Model", scenario)) +
    theme_bw() +
    theme(legend.position = "none") +
    theme(plot.title = element_text(size = 9))
  
  list(p1, p2)
  
}

###

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

#####
Mortality_df<- function(projection){
  
  param <- projection$param
  tmax <- projection$param$tmax
  time_cutoff <- tmax/2
  
  df_Z <- projection$Z
  df_Z <- as.data.frame(df_Z)
  
  df_Z <- melt(projection$Z) # Extract adundance from model projection and use melt to convert 3D                                  array into long form data
  #glimpse(df_N) # Show a glimpse of the new df created
  names(df_Z) <- c("Time", "Species_Names", "Size_Class", "Z") # Rename variable names
  #glimpse(df_N) # Show a glimpse of the renamed variables in df_N
  df_Z$Species_Names <- factor(df_Z$Species_Names) # Change Species_Names to a factor
  #glimpse(df_N) # Check it's done the job you think it ha
  
  #eval(parse(text=paste0("projection"))) <<- df_N
  #list2env(df_N, envir = paste(projection),.GlobalEnv)
  
  df_size_classes <- projection$w # Create new object with list of size classes in log10 g
  df_size_classes <- as.data.frame(df_size_classes) # Convert to a dataframe
  names(df_size_classes) <- c("Size") # Create new name Size for the variable 
  df_size_classes$Size_Class <- 1:length(projection$w) # Create new variable 1:184 in ascending order with log10 g                                             size classes, which will be used to left join later
  
  df_species <- projection$param$groups$species # Create new object with species names from model
  df_species <- as.data.frame((df_species)) # Convert to a dataframe
  names(df_species) <- c("Species") # Rename variable
  df_species$Species_Names <- 1:length(param$groups$species) # Create new integer variable in ascending order with 'Species'
  df_species$Species_Names <- as.factor(df_species$Species_Names) # Convert to factor
  
  
  
  ## join df_N with df_size_class data in order to create log10 size classes rather then the size     class integer 1:184 that is in df_N
  df_combined <- left_join(
    select(df_Z, Time, Species_Names, Size_Class, Z),
    select(df_size_classes, Size, Size_Class),
    by="Size_Class") %>%
    # then join df_species to df_combined so that there is a variable with actual species names rather than a numeric variable 1:n (in this case 1:7)
    left_join(df_species,
              by="Species_Names")
  
  
  # # Create a filtered subset of df_combined that only takes the last year
  # df_plot <- df_combined %>%
  #   filter(Time == tmax)   %>% 
  #   group_by(Species, Size) %>% # Then group this subset by species and size
  #   summarise(Mean_N = mean(N)) %>% # For each subset of species and size, caluclate mean N
  #   filter(Mean_N > 0) # Get rid of any zero values
  
  
  df_Z_TimeAvg <- df_combined %>%
    filter(Time >= time_cutoff, Time <= tmax)   %>% 
    group_by(Species, Size) %>% # Then group this subset by species and size
    summarise(Mean_Z = mean(Z)) #%>% # For each subset of species and size, caluclate mean N
  #filter(Mean_N > 0)
  
  #write_csv(df_plot, paste("Model Simulation Output/",file_name, "Abundance_df.csv"))
  return(df_Z_TimeAvg)
  
}


PB_df <- function(projection){
  
  param <- projection$param
  tmax <- projection$param$tmax
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

