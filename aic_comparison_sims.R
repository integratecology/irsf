lCKAGES ####
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

# Create an empty data.frame for saving results
name_df <- c("sim_no", "irsf_coef", "irsf_lcl", "irsf_ucl", "crsf_coef", "crsf_lcl", "crsf_ucl", "runtime")
df_sims <- array(rep(NaN), dim = c(0, length(name_df)))
colnames(df_sims) <- name_df

# Create raster
r1 <- raster(nrows = 1000, ncols = 1000, 
             xmn = -0.05, xmx = 0.05, ymn = -0.05, ymx = 0.05,
             vals = as.logical(rep(0:1, 500000)))
projection(r1) <- "+proj=longlat +datum=WGS84"

# Record start time to monitor how long replicates take to compute
sTime <- Sys.time()
print(sTime)

# Loop over sampling frequencies (samp)
 
# Specify variables to manipulate sampling frequency while holding duration constant
set.seed(sim_no) # unique seed for each simulation run
nd <- 800 # number of days
pd <- 96 # number of sampled points per day
  
# Sampling schedule
st <- 1:(nd*pd)*(ds/pd) 
    
# Simulate from the movement model ###
sim <- simulate(mod, t=st, complete = TRUE)

# Create a new track with habitat selection (shift half of points in movement track right by 0.0001 degrees)
sim_sub <- sim
sim_sub$longitude <- ifelse(sim_sub$longitude > 0 & floor(sim_sub$longitude*10000) %% 2 == 0, 
                            sim_sub$longitude + 0.0001, sim_sub$longitude)

df2 <- data.frame(sim_sub)
pts2 <- df2[,6:7]
sp <- SpatialPoints(pts2, proj4string = CRS(projection(r1)))

# Check number of points in each habitat
df2$habitat <- raster::extract(r1, sp)
sum(df2$habitat)

# Create test and training sets 
tl <- length(sim_sub$timestamp)
fh <- round(tl/2)
sh <- round(fh)+1
train <- sim_sub[1:fh,]
test <- sim_sub[sh:tl,]

# Export the test data as an sp object
test_sp <- SpatialPoints.telemetry(test)
   
# Fit the movement model to the simulated data
fit <- ctmm.fit(train, CTMM=mod, control=list(method="pNewton")) #
summary(fit)
  
# Calculate the UD ###
ud <- akde(train, fit, weights=TRUE)
summary(ud)
  
# Fit and compare the RSFs ###
irsf <- ctmm:::rsf.fit(train, UD=ud, R=list(test=r1), debias=TRUE, error=0.01)
summary(irsf)

agde_irsf <- agde(irsf, grid=r1)  
summary(agde_irsf)

crsf - ctmm:::rsf.fit(train, UD=ud, R=list(test=1), debias=TRUE, error=0.01, integrated=FALSE)
summary(crsf)

agde_crsf <- agde(crsf, grid=r1)
summary(agde_irsf)

aic_comp <- summary(list(irsf=irsf, crsf=crsf))

# Cross-validated log-likelihoods
ud_irsf_r <- raster(agde_irsf, DF = "PDF")
ud_crsf_r <- raster(agde_crsf, DF = "PDF")

cv_ll_i <- raster::extract(ud_irsf_r, test_sp)
cv_ll_irsf <- sum(cv_ll_i)

cv_ll_c <- raster::extract(ud_crsf_r, test_sp)
cv_ll_crsf <- sum(cv_ll_c)

eTime <- Sys.time()
  
# Extract variables of interest ###
sim_no <- sim_no
irsf_est <- summary(irsf)$CI[1,2]
irsf_lcl <- summary(irsf)$CI[1,1]
irsf_ucl <- summary(irsf)$CI[1,3]
crsf_est <- summary(crsf)$CI[1,2]
crsf_lcl <- summary(crsf)$CI[1,1]
crsf_ucl <- summary(crsf)$CI[1,3]
best_mod <- row.names(aic_comp)[1]
delta_aic <- aic_comp[2,1]
runtime <- difftime(eTime, sTime, units="mins")
habitat1 <- sum(df2$habitat)
habitat2 <- length(df2$habitat)-habitat1
  
#################################
# Vector of results to return
x <- data.frame(sim_no, irsf_est, irsf_lcl, irsf_ucl, crsf_est, crsf_lcl, crsf_ucl, best_mod, delta_aic, habitat1, habitat2, cv_ll_irsf, cv_ll_crsf, runtime)
  
# Store results in data.frame
write.table(x, 'results/aic_comparison_sims.csv', append=TRUE, row.names=FALSE, col.names=FALSE, sep=',') 
  
# Print indicators of progress
print(pd)
print(eTime)
print(eTime - sTime)

