# PACKAGES ####
library(ctmm)
library(data.table)
library(raster)

# DATA ####
setwd("/home/jalston/proj/wrsf")

mongoose <- read.csv("mongoose.csv")
serval <- read.csv("streicher_servals.csv")
caracal <- read.csv("streicher_caracals.csv")

names(mongoose) <- c("species", "sex", "individual.local.identifier", "timestamp", "location.lat", "location.long")
mongoose$timestamp <- as.POSIXct(mongoose$timestamp, format="%m/%d/%Y %H:%M")

serval$timestamp <- as.POSIXct(serval$timestamp, format="%m/%d/%Y %H:%M")
caracal$timestamp <- as.POSIXct(caracal$timestamp, format="%m/%d/%Y %H:%M")

caracal4 <- caracal[caracal$individual.local.identifier==4,]
caracal4 <- caracal4[1:504,]
caracal4 <- caracal4[,-2]

caracal_t <- as.telemetry(caracal)

fwrite(caracal4, "data/caracal_example.csv")

mongoose_t <- as.telemetry(mongoose)
mongoose_l <- mongoose_t[[10]]

fwrite(mongoose_l, "data/mongoose_example.csv")

serval2 <- serval[serval$individual.local.identifier=="Nine",]
tdiff <- serval_t[[13]]$t - lag(serval_t[[13]]$t,1)
serval2 <- serval2[tdiff!=0,]

serval_t <- as.telemetry(serval)
serval_1 <- serval_t[[1]]
serval_l <- serval_l[!duplicated(serval_l$timestamp),]

fwrite(serval_l, "data/serval_example.csv")


DATA <- mongoose_l

SVF <- variogram(DATA)
GUESS <- ctmm.guess(DATA, variogram=SVF, ctmm(error=TRUE, isotropic=TRUE), interactive=FALSE)
FIT <- ctmm.select(DATA, GUESS, trace=2)
summary(FIT)
zoom(SVF, FIT)

UD <- akde(DATA, FIT, weights=TRUE)
plot(DATA, UD)
summary(UD)

# Create raster and check where data falls on it
r1 <- raster("KZN_2017_Landuse_latlong.tif")
# projection(r1) <- "+proj=aeqd +lon_0=0 +lat_0=0 +datum=WGS84 +units=m"
raster::plot(r1)
points(DATA$latitude~DATA$longitude)

coords <- DATA[,2:3]
sp <- SpatialPoints(coords)

check <- raster::extract(r1, sp)
table(check)


e <- extent(min(DATA$longitude) - 0.2, max(DATA$longitude) + 0.2,
            min(DATA$latitude) - 0.2, max(DATA$latitude) + 0.2)
lc2 <- crop(r1, e)
raster::plot(lc2)
points(DATA$latitude ~ DATA$longitude)
table(lc2@data@values)

# Manipulating rasters
m1 <- c(1, 1, 0, 2, 2, 1, 3, 43, 0, 44, 45, 1, 46, 48,0)
rclmat1 <- matrix(m1, ncol = 3, byrow = TRUE)
plantation <- reclassify(lc2, rclmat1, include.lowest=TRUE, right = NA)
plantation@data@values <- as.logical(plantation@data@values)
raster::plot(plantation)
points(DATA$latitude ~ DATA$longitude)

m2 <- c(1, 11, 0, 12, 12, 1,13, 13, 0, 14, 14, 1, 15, 33, 0, 34, 35, 1, 36, 41, 0, 42, 43, 1, 44, 48, 0)
rclmat2 <- matrix(m2, ncol = 3, byrow = TRUE)
urban <- reclassify(lc2, rclmat2, include.lowest=TRUE, right = NA)
urban@data@values <- as.logical(urban@data@values)
raster::plot(urban)
points(DATA$latitude ~ DATA$longitude)

m3 <- c(1, 17, 0, 18, 20, 1, 21, 24, 0, 25, 25, 1, 26, 44, 0)
rclmat3 <- matrix(m3, ncol = 3, byrow = TRUE)
forest <- reclassify(lc2, rclmat3, include.lowest=TRUE, right = NA)
forest@data@values <- as.logical(forest@data@values)
raster::plot(forest)
points(DATA$latitude ~ DATA$longitude)

m4 <- c(1, 22, 0, 23, 23, 1, 24, 26, 0, 27, 27, 1, 28, 44, 0)
rclmat4 <- matrix(m4, ncol = 3, byrow = TRUE)
grassland <- reclassify(lc2, rclmat4, include.lowest=TRUE, right = NA)
grassland@data@values <- as.logical(grassland@data@values)
raster::plot(grassland)
points(DATA$latitude ~ DATA$longitude)

# Resource selection function with two habitat covariates
RSF.W <- ctmm:::rsf.fit(DATA, UD=UD, R=list(plantation=plantation,urban=urban), 
                              debias=TRUE, error=0.1)
summary(RSF.W)

# Look at a plot of the distribution of availability
# plot(agde(RSF.HABITAT))

# now consider the IID case (much slower because N<<n)
IID <- ctmm.fit(DATA,CTMM=ctmm(isotropic=TRUE))
UD.IID <- akde(DATA,IID)
plot(DATA, UD.IID)

RSF.IID <- ctmm:::rsf.fit(DATA, UD=UD.IID, R=list(urban=urban, plantation=plantation), 
                          debias=TRUE, error=0.25)
summary(RSF.IID)
# compare to CIs of
summary(RSF.W)
# confidence intervals are much wider

