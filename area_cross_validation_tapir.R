# PACKAGES ####
library(sp)
library(raster)
library(ctmm)

# ANALYSIS ####
# Read job number from command line
args = commandArgs(trailingOnly=TRUE)
ind_file <- args[1]
print(ind_file)

# Load data and create telemetry object
df <- fread(ind_file)
aid <- df$individual.local.identifier[1]

l <- as.telemetry(df)

# Record start time to monitor how long replicates take to compute
sTime <- Sys.time()
print(sTime)

# Create training and test sets ###
tl <- length(l$timestamp)
fh <- tl/2
sh <- tl/2+1
train <- l[1:fh,]
test <- l[sh:tl,]

# ANALYSIS ####
# Fit the movement model to the test data
svf <- variogram(train)
guess <- ctmm.guess(train, variogram=svf, interactive=FALSE)
fit <- ctmm.select(train, guess, trace=2)
summary(fit) 
print("Fitted movement model")
  
# Calculate the UD ###
ud <- akde(train, fit, weights=TRUE)
print("UD created")

# Fit the RSFs ###
irsf <- ctmm:::rsf.fit(train, UD=ud, R=list(), debias=TRUE, error=0.01)
summary(irsf)
print("Fitted iRSF") 

crsf <- ctmm:::rsf.fit(train, UD=ud, R=list(), integrated = FALSE, debias=TRUE, error=0.01)
summary(crsf)
print("Fitted cRSF")
  
# Check out AGDEs of fitted RSFs ###
# Create a dummy raster for the AGDE grid #
r1 <- raster(nrows = 100, ncols = 100, xmn = -5000, xmx = 5000, ymn = -5000, ymx = 5000, 
                 vals = rep(1,10000))
projection(r1) <- ctmm:::projection(irsf)

# Integrated RSF #
agde_irsf <- agde(irsf, grid=r1)
# plot(agde_irsf)
# hist(agde_irsf@.Data[[1]][agde_irsf@.Data[[1]]!=0])

# Conventional RSF #
agde_crsf <- agde(crsf, grid=r1)
# plot(agde_crsf)
  
# Create raster of test set density ###
sp <- SpatialPoints(test[6:7], proj4string=CRS(ctmm:::projection(irsf))
countr <- rasterize(sp, r1, fun='count')
countr2 <- overlay(countr, fun = function(x){ x / length(test@.Data[[1]]) })
# raster::plot(countr2)
  
check <- data.frame(aid = rep(aid, length(as.vector(agde_irsf@.Data[[1]]))),
                    prob_irsf = as.vector(agde_irsf@.Data[[1]]/sum(agde_irsf@.Data[[1]])), 
                    prob_crsf = as.vector(agde_crsf@.Data[[1]]/sum(agde_crsf@.Data[[1]])), 
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
cor_irsf <- cor(check$prob_irsf, check$prob_emp, method = "pearson")
cor_crsf <-cor(check$prob_crsf, check$prob_emp, method = "pearson")
kld_irsf <- LaplacesDemon::KLD(check$prob_irsf, check$prob_emp)$sum.KLD.px.py
kld_crsf <- LaplacesDemon::KLD(check$prob_crsf, check$prob_emp)$sum.KLD.px.py
kld_r2 <- exp(-LaplacesDemon::KLD(check$prob_irsf, check$prob_emp)$sum.KLD.px.py/LaplacesDemon::KLD(check$prob_crsf, check$prob_emp)$sum.KLD.px.py)
  
# Create vector of results
results <- data.frame(aid, ind_file, cor_irsf, cor_crsf, kld_irsf, kld_crsf, kld_r2)
  
# Store results in data.frame
write.table(results, 'results/area_cv_summary_tapir.csv', append=TRUE, row.names=FALSE, col.names=FALSE, sep=',') 
write.table(check, 'results/area_cv_data_tapir.csv', append=TRUE, row.names=FALSE, col.names=FALSE, sep=',')
  
# Print indicators of progress
print(aid)
print(eTime)
print(eTime - sTime)
