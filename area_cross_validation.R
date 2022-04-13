# PACKAGES ####
library(sp)
library(raster)
library(ctmm)

# ANALYSIS ####
# Read job number from command line
# args = commandArgs(trailingOnly=TRUE)
# sim_no <- args[1]
# print(sim_no)

# Define the OUF model parameters ###
# Length of a day in seconds
ds <- 86400

# Spatial variance in m^2
sig <- 200000

# True 95% range area
trueRngArea <- -2*log(0.05)*pi*sig

# Specify an OUF model for simulation
mod <- ctmm(tau=c(ds,ds/6), isotropic=TRUE, sigma=sig, mu=c(0,0))

# SIMULATION ####

# Sampling frequencies to quantify
samp <- 24

# Create an empty data.frame for saving results
# name_df <- c("sim_no","samp_freq", "wrsf_coef", "wrsf_lcl", "wrsf_ucl", "runtime")
# df_sims <- array(rep(NaN), dim = c(0, length(name_df)))
# colnames(df_sims) <- name_df

# Create raster
r1 <- raster(nrows = 100, ncols = 100, xmn = -5000, xmx = 5000, ymn = -5000, ymx = 5000, 
             vals = c(rep(1,9999),1.0000001)) # c(rep(1,800),rep(0,800))))
projection(r1) <- "+proj=aeqd +lon_0=0 +lat_0=0 +datum=WGS84"

# Record start time to monitor how long replicates take to compute
sTime <- Sys.time()
print(sTime)

# Loop over sampling frequencies (samp)
  
  # Specify variables to manipulate sampling frequency while holding duration constant
  nd <- 400 # number of days
  pd <- samp # [i] # number of sampled points per day
  
  # Sampling schedule
  st <- 1:(nd*pd)*(ds/pd) 
  
  # Simulate from the movement model ###
  sim <- simulate(mod, t=st, complete = TRUE)
  df <- as.data.frame(sim)
  
  # Create training and test sets ###
  tl <- length(st)
  fh <- tl/2
  sh <- tl/2+1
  train <- sim[1:fh,]
  test <- sim[sh:tl,]
  
  # Fit the movement model to the simulated data
  fit <- ctmm.fit(train, CTMM=mod, control=list(method="pNewton")) #
  print("Fitted movement model")
  
  # Calculate the UDs ###
  ud <- akde(train, fit, weights=TRUE)
  print("UD created")
  
  # Fit the iRSF ###
  rsf <- ctmm:::rsf.fit(train, UD=ud, R=list(test=r1), debias=TRUE, error=0.1)
  summary(rsf)
  print("Fitted iRSF") 
  
  # Fit the conventional RSF ###
  rsf2 <- ctmm:::rsf.fit(train, UD=ud, R=list(test=r1), integrated = FALSE, debias=TRUE, error=0.1)
  summary(rsf2)
  print("Fitted cRSF")
  
  # Check out AGDE of iRSF ###
  AGDE <- agde(rsf, grid=r1)
  # plot(AGDE)
  # hist(AGDE@.Data[[1]][AGDE@.Data[[1]]!=0])
  
  # Create raster of test set density ###
  sp <- SpatialPoints(test[2:3], proj4string=CRS("+proj=aeqd +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))
  countr <- rasterize(sp, r1, fun='count')
  countr2 <- overlay(countr, fun = function(x){ x / length(test@.Data[[1]]) })
  raster::plot(countr2)
  
  
  AGDE2 <- agde(rsf2, R=list(test=r1), grid=r1)
  plot(AGDE2)
  
  check <- data.frame(agde = as.vector(AGDE@.Data[[1]]/sum(AGDE@.Data[[1]])), agde2 = as.vector(AGDE2@.Data[[1]]/sum(AGDE2@.Data[[1]])), count = countr2@data@values)
  check$count = ifelse(is.na(check$count), 0, check$count)
  check$log_agde <- log(check$agde)
  check$log_agde2 <- log(check$agde2)
  check$log_count <- log(check$count)
  quantile(check$count, probs = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))
  plot(check$count ~ check$agde, xlim=c(0,.012), ylim=c(0,.012), main = "Integrated",
       xlab="Probability Density (AGDE)", ylab="Probability Density (Movement Track)")
  abline(a=0, b=1)
  plot(check$count ~ check$agde2, xlim=c(0,.012), ylim=c(0,.012), main = "Conventional",
       xlab="Probability Density (Uniform from AKDE)", ylab="Probability Density (Movement Track)")
  abline(a=0, b=1)
  cor(check$agde, check$count, method = "pearson")
  cor(check$agde2, check$count, method = "pearson")
  LaplacesDemon::KLD(check$agde, check$count)$sum.KLD.px.py
  LaplacesDemon::KLD(check$agde2, check$count)$sum.KLD.px.py
  exp(-LaplacesDemon::KLD(check$agde, check$count)$sum.KLD.px.py/LaplacesDemon::KLD(check$agde2, check$count)$sum.KLD.px.py)
  
  suit <- suitability(rsf, R=list(), grid=ud)
  table(suit@data@values)
  
  eTime <- Sys.time()
  
  # Extract variables of interest ###
  sim_no <- sim_no
  samp_freq <- pd
  wrsf_coef <- summary(rsf)$CI[1,2]
  wrsf_lcl <- summary(rsf)$CI[1,1]
  wrsf_ucl <- summary(rsf)$CI[1,3]
  runtime <- difftime(eTime, sTime, units="mins")
  
  #################################
  # Vector of results to return
  x <- data.frame(sim_no, samp_freq, wrsf_coef, wrsf_lcl, wrsf_ucl, runtime)
  
  # Store results in data.frame
  write.table(x, 'results/final/wrsf_sim_results_lo_wrsf.csv', append=TRUE, row.names=FALSE, col.names=FALSE, sep=',') 
  
  # Print indicators of progress
  print(pd)
  print(eTime)
  print(eTime - sTime)
