# PACKAGES ####
library(sp)
library(raster)
library(ctmm)

# ANALYSIS ####
# Read job number from command line
args = commandArgs(trailingOnly=TRUE)
sim_no <- args[1]
print(sim_no)

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
# Create an empty data.frame for saving results
# name_df <- c("sim_no","samp_freq", "wrsf_coef", "wrsf_lcl", "wrsf_ucl", "runtime")
# df_sims <- array(rep(NaN), dim = c(0, length(name_df)))
# colnames(df_sims) <- name_df

# Create raster
r1 <- raster(nrows = 100, ncols = 100, xmn = -5000, xmx = 5000, ymn = -5000, ymx = 5000, 
             vals = rep(1,10000))
projection(r1) <- "+proj=aeqd +lon_0=0 +lat_0=0 +datum=WGS84"

# Record start time to monitor how long replicates take to compute
sTime <- Sys.time()
print(sTime)

# Loop over sampling frequencies (samp)
  
# Specify variables to manipulate sampling frequency while holding duration constant
nd <- 1000 # number of days
pd <- 24   # number of sampled points per day
 
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
irsf <- ctmm:::rsf.fit(train, UD=ud, R=list(), debias=TRUE, error=0.1)
summary(irsf)
print("Fitted iRSF") 
  
# Fit the conventional RSF ###
crsf <- ctmm:::rsf.fit(train, UD=ud, R=list(), integrated = FALSE, debias=TRUE, error=0.1)
summary(crsf)
print("Fitted cRSF")
  
# Check out AGDE of iRSF ###
agde_irsf <- agde(irsf, grid=r1)
agde_irsf_r <- raster(agde_irsf, DF="PMF")
# plot(agde_irsf)
# hist(agde_irsf@.Data[[1]][agde_irsf@.Data[[1]]!=0])
  
# Create raster of test set density ###
sp <- SpatialPoints(test[2:3], proj4string=CRS("+proj=aeqd +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))
countr <- rasterize(sp, r1, fun='count')
countr2 <- overlay(countr, fun = function(x){ x / length(test@.Data[[1]]) })
# raster::plot(countr2)
  
# Check out AGDE of cRSF
agde_crsf <- agde(crsf, R=list(test=r1), grid=r1)
agde_crsf_r <- raster(agde_crsf, DF="PMF")
# plot(agde_crsf)

check <- data.frame(sim_no = rep(sim_no, length(agde_irsf@data@values)),
                    prob_irsf = agde_irsf_r@data@values), 
                    prob_crsf = agde_crsf_r@data@values), 
                    prob_emp = countr2@data@values)
check$prob_emp = ifelse(is.na(check$prob_emp), 0, check$prob_emp)
# quantile(check$count, probs = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))
# plot(check$count ~ check$agde, xlim=c(0,.012), ylim=c(0,.012), main = "Integrated",
#      xlab="Probability Density (AGDE)", ylab="Probability Density (Movement Track)")
# abline(a=0, b=1)
# plot(check$count ~ check$agde2, xlim=c(0,.012), ylim=c(0,.012), main = "Conventional",
#     xlab="Probability Density (Uniform from AKDE)", ylab="Probability Density (Movement Track)")
# abline(a=0, b=1)
eTime <- Sys.time()
  
# Extract variables of interest ###
# sim_no <- sim_no
cor_irsf <- cor(check$prob_irsf, check$prob_emp, method = "pearson")
cor_crsf <-cor(check$prob_crsf, check$prob_emp, method = "pearson")
kld_irsf <- LaplacesDemon::KLD(check$prob_irsf, check$prob_emp)$sum.KLD.px.py
kld_crsf <- LaplacesDemon::KLD(check$prob_crsf, check$prob_emp)$sum.KLD.px.py
kld_r2 <- exp(-LaplacesDemon::KLD(check$prob_irsf, check$prob_emp)$sum.KLD.px.py/LaplacesDemon::KLD(check$prob_crsf, check$prob_emp)$sum.KLD.px.py)
  
#################################
# Vector of results to return
results <- data.frame(sim_no, cor_irsf, cor_crsf, kld_irsf, kld_crsf, kld_r2)
  
# Store results in data.frame
write.table(results, 'results/area_cv_summary_sims.csv', append=TRUE, row.names=FALSE, col.names=FALSE, sep=',') 
write.table(check, 'results/area_cv_data_sims.csv', append=TRUE, row.names=FALSE, col.names=FALSE, sep=',')
  
# Print indicators of progress
print(sim_no)
print(eTime)
print(eTime - sTime)
