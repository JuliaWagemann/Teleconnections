# Necessary libraries
library(raster)
library(ncdf)
library(fields)
library(ncdf4)

# Load functions needed for anomaly calculation --> if necessary adjust path to source code
source("F:/Julia/Data/anomalies_functions_CCI.R")

# Set working directory
<<<<<<< HEAD
path <- "F:/Julia/Data/CCI/Aerosol/"
=======
<<<<<<< HEAD
path <- "F:/Julia/Data/ERA_interim/1979_2013/GP/GP_500/"
=======
<<<<<<< HEAD
path <- "F:/Julia/Data/ERA_interim/1979_2013/GP/GP_500/"
=======
path <- "F:/Julia/Data/CCI/SM/"
>>>>>>> b76b7e7e5fb454b4e9d6d19421fdbdd137cf5c78
>>>>>>> 0f8c816104918abd5389b14768c06abe324633c8
>>>>>>> origin/master
setwd(path)


###############################################################################
# Section for individual adjustment
# Please specify the all the variables in this section individually for anomaly
# calculation and netcdf writing
###############################################################################
<<<<<<< HEAD
nrOfYears <- 8 # Please adjust here the amount of years over which the anomalies are calculated
varName <- "AOD_anomaly" # Name of the variable name as indicated in the netcdf file
varUnit <- "" # Unit of the variable stored in the netcdf file
varDescription <- "Atmospheric optical depth anomaly" # Description of the variable in the netcdf file
fileName <- "AOD_anomaly_" # Beginning of the individual netcdf file name --> for each year, a running number will be concatenated
=======
<<<<<<< HEAD
=======
<<<<<<< HEAD
>>>>>>> 0f8c816104918abd5389b14768c06abe324633c8
nrOfYears <- 35 # Please adjust here the amount of years over which the anomalies are calculated
varName <- "GP500_anomaly" # Name of the variable name as indicated in the netcdf file
varUnit <- "m**2 s**-2" # Unit of the variable stored in the netcdf file
varDescription <- "Geopotential at 500 mbar anomaly" # Description of the variable in the netcdf file
fileName <- "GP500_anomaly_" # Beginning of the individual netcdf file name --> for each year, a running number will be concatenated
<<<<<<< HEAD
=======
=======
nrOfYears <- 32 # Please adjust here the amount of years over which the anomalies are calculated
varName <- "SM_anomaly" # Name of the variable name as indicated in the netcdf file
varUnit <- "Volumetric soil moisture [m3/m3]" # Unit of the variable stored in the netcdf file
varDescription <- "Soil moisture anomaly" # Description of the variable in the netcdf file
fileName <- "SM_anomaly_" # Beginning of the individual netcdf file name --> for each year, a running number will be concatenated
>>>>>>> b76b7e7e5fb454b4e9d6d19421fdbdd137cf5c78
>>>>>>> 0f8c816104918abd5389b14768c06abe324633c8
>>>>>>> origin/master
###############################################################################

# All files in the working directory are stored in a list
# Required file structure: 1 netcdf file per year with monthly values of the parameter
fileList <- list.files(path,pattern="*.nc")

# NetCDF Output preparation
lat <- getLat(fileList[1])
lon <- getLon(fileList[1])

# Define the dimensions of the resulting netCDF file
x <- ncdim_def("lon", "degrees",lon)
y <- ncdim_def("lat", "degrees", lat)
time <- ncdim_def("Time", "Month", c(1:12))

# Define Variable which is written to the netCDF file
netCDF_varDef <- ncvar_def(varName, units= varUnit,dim=list(y,x,time), missval=NA,
                           longname=varDescription, prec= 'double', verbose=FALSE)

# Routine which restructures the stack based on month to stack based on years
# containing the monthly anomalies per year
stackVector <- createStackVector(fileList)
for(j in 1:nrOfYears){
  tempStack <- stack()
  monthVector <- c(01,02,03,04,05,06,07,08,09,10,11,12)
  for(i in stackVector){
    tempStack <- addLayer(tempStack,i[[j]])
  }  
  # Create the actual netCDF output file
<<<<<<< HEAD
  netCDFFile <- nc_create(paste(fileName,j,".nc",sep=""), netCDF_varDef, force_v4=FALSE,verbose=FALSE)
=======
<<<<<<< HEAD
  netCDFFile <- nc_create(paste(fileName,j,"_mm.nc",sep=""), netCDF_varDef, force_v4=FALSE,verbose=FALSE)
=======
<<<<<<< HEAD
  netCDFFile <- nc_create(paste(fileName,j,"_mm.nc",sep=""), netCDF_varDef, force_v4=FALSE,verbose=FALSE)
=======
  netCDFFile <- nc_create(paste(fileName,j,"_test.nc",sep=""), netCDF_varDef, force_v4=FALSE,verbose=FALSE)
>>>>>>> b76b7e7e5fb454b4e9d6d19421fdbdd137cf5c78
>>>>>>> 0f8c816104918abd5389b14768c06abe324633c8
>>>>>>> origin/master
  
  # fill the created netCDF file with the actual content --> calculated anomalies
  ncvar_put(netCDFFile, varName, as.matrix(brick(tempStack)),verbose=FALSE)

  # close the created netCDF file
  nc_close(netCDFFile)
  rm(tempStack,i)
}

#######################################################################################
# Testing environment for the created files
# Opening the netcdf, getting the variable stored and plot individual month of it
<<<<<<< HEAD
test <- open.ncdf("SST_anomaly_3.nc",write=FALSE)
testVariable <- get.var.ncdf(test,"SST_anomaly")
image(testVariable[,,5])
=======
<<<<<<< HEAD
test <- open.ncdf("TP_anomaly_1.nc",write=FALSE)
testVariable <- get.var.ncdf(test,"TP_anomaly")
raster(testVariable)
image(testVariable[,,1])
=======
<<<<<<< HEAD
test <- open.ncdf("TP_anomaly_1.nc",write=FALSE)
testVariable <- get.var.ncdf(test,"TP_anomaly")
raster(testVariable)
image(testVariable[,,1])
=======
test <- open.ncdf(fileList[2],write=FALSE)
testVariable <- get.var.ncdf(test,"01")
raster(testVariable)
image(testVariable)
>>>>>>> b76b7e7e5fb454b4e9d6d19421fdbdd137cf5c78
>>>>>>> 0f8c816104918abd5389b14768c06abe324633c8

breakpoints <- c(-0.6, -0.2, 0, 0.2, 0.6)
colors <- c("red", "orange","white","blue","darkblue")
>>>>>>> origin/master

close.ncdf(test)
