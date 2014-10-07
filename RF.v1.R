library(randomForest)
library(ncdf)
library(RNetCDF)
library(ncdf4)
library(fields)
library(raster)

#########################################################################
# Section for individual adjustment
# Please specify all variables in this section individually for every single parameter
#########################################################################
# define the number of trees for the random forest
NTree <- 500

#BE SURE THAT ALL THE FOLDERS CONTAIN THE SAME YEARS!!!
#How many years are contained in each folder?
Nyears <- 35

#define the paths of the predictors:
path1 <- "//TEC-JULIA-FRA/Users/julia_francesca/Documents/Data/2_Anomalies/ERA_Interim/monthly/1_Geopotential_500/3deg"
path2 <- "//TEC-JULIA-FRA/Users/julia_francesca/Documents/Data/2_Anomalies/ERA_Interim/monthly/2_SLP/3deg"
path3 <- "//TEC-JULIA-FRA/Users/julia_francesca/Documents/Data/2_Anomalies/ERA_Interim/monthly/3_T2m/3deg"

#define the path of the variable to be predicted:
path4 <- "//TEC-JULIA-FRA/Users/julia_francesca/Documents/Data/2_Anomalies/ERA_Interim/monthly/4_TP/3deg"

# what are the dimension of the grid?
geoLat <- 60
geoLon <- 120

OutputVar <- "RF_TP"
########################################################################
fileList1 <- list.files(path1, pattern=".nc")
fileList2 <- list.files(path2, pattern=".nc")
fileList3 <- list.files(path3, pattern=".nc")
fileList4 <- list.files(path4, pattern=".nc")

nrOfTimeSteps <- Nyears*12
nrOfPixels <- geoLat*geoLon

vector1 <- vector(mode = "numeric", length = nrOfTimeSteps*nrOfPixels/4)
vector2 <- vector(mode = "numeric", length = nrOfTimeSteps*nrOfPixels/4)
vector3 <- vector(mode = "numeric", length = nrOfTimeSteps*nrOfPixels/4)
vector4 <- vector(mode = "numeric", length = nrOfTimeSteps*nrOfPixels/4)

#predictor 1
setwd(path1)
dataMatrix <- matrix(0, nrOfTimeSteps, nrOfPixels/4)
t <- 1  #initialize a count variable
for(m in fileList1){

        GP_anomaly <- open.ncdf(m, write=FALSE)
        GP_var <- get.var.ncdf(GP_anomaly, "GP500_anomaly")
        for(z in 1:12){
                tempRaster <- GP_var[,,z]
                tempRaster[tempRaster=="NaN"] <- NA
                tempRaster <- raster(tempRaster)
                tempRaster <- aggregate(tempRaster, fact=2, FUN=mean, na.rm=TRUE)
                tempRaster <- as.matrix(tempRaster)
                dataMatrix[t,] <- tempRaster
                t <- t+1 
        }
}
vector1 <- as.vector(mode="numeric", dataMatrix)

#predictor 2
setwd(path2)
dataMatrix <- matrix(0, nrOfTimeSteps, nrOfPixels/4)
t <- 1  #initialize a count variable
for(m in fileList2){
        
        GP_anomaly <- open.ncdf(m, write=FALSE)
        GP_var <- get.var.ncdf(GP_anomaly, "SLP_anomaly")
        for(z in 1:12){
                tempRaster <- GP_var[,,z]
                tempRaster[tempRaster=="NaN"] <- NA
                tempRaster <- raster(tempRaster)
                tempRaster <- aggregate(tempRaster, fact=2, FUN=mean, na.rm=TRUE)
                tempRaster <- as.matrix(tempRaster)
                dataMatrix[t,] <- tempRaster
                t <- t+1 
        }
}
vector2 <- as.vector(mode="numeric", dataMatrix)

#predictor 3
setwd(path3)
dataMatrix <- matrix(0, nrOfTimeSteps, nrOfPixels/4)
t <- 1  #initialize a count variable
for(m in fileList3){
        
        GP_anomaly <- open.ncdf(m, write=FALSE)
        GP_var <- get.var.ncdf(GP_anomaly, "T2m_anomaly")
        for(z in 1:12){
                tempRaster <- GP_var[,,z]
                tempRaster[tempRaster=="NaN"] <- NA
                tempRaster <- raster(tempRaster)
                tempRaster <- aggregate(tempRaster, fact=2, FUN=mean, na.rm=TRUE)
                tempRaster <- as.matrix(tempRaster)
                dataMatrix[t,] <- tempRaster
                t <- t+1 
        }
}
vector3 <- as.vector(mode="numeric", dataMatrix)

#output
setwd(path4)
dataMatrix <- matrix(0, nrOfTimeSteps, nrOfPixels/4)
t <- 1  #initialize a count variable
for(m in fileList4){
        
        GP_anomaly <- open.ncdf(m, write=FALSE)
        GP_var <- get.var.ncdf(GP_anomaly, "TP_anomaly")
        for(z in 1:12){
                tempRaster <- GP_var[,,z]
                tempRaster[tempRaster=="NaN"] <- NA
                tempRaster <- raster(tempRaster)
                tempRaster <- aggregate(tempRaster, fact=2, FUN=mean, na.rm=TRUE)
                tempRaster <- as.matrix(tempRaster)
                dataMatrix[t,] <- tempRaster
                t <- t+1 
        }
}
vector4 <- as.vector(mode="numeric", dataMatrix)

set.seed(20)

DF <- data.frame(vector1, vector2, vector3)

RF <- randomForest(x=DF, y=vector4, ntree=NTree)
