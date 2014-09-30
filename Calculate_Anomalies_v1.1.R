# Necessary libraries
library(raster)
library(ncdf)
library(fields)
library(ncdf4)

# Load functions needed for anomaly calculation
source("C:/Users/Julia Wagemann/Desktop/Teleconnections/anomalies_functions.R")

# Set working directory
path <- "//TEC-JULIA-FRA/Users/julia_francesca/Documents/Data/1_Raw_Data_NCDF/ERA_Interim/Tair_2m"
setwd(path)


###############################################################################
# Section for individual adjustment
# Please specify the all the variables in this section individually for anomaly
# calculation and netcdf writing
###############################################################################
nrOfYears <- 24 # Please adjust here the amount of years over which the anomalies are calculated
varName <- "t2m_anomaly" # Name of the variable name as indicated in the netcdf file
varUnit <- "K" # Unit of the variable stored in the netcdf file
varDescription <- "Air temperature 2m anomalies" # Description of the variable in the netcdf file
fileName <- "Tair_2m_anomaly_" # Beginning of the individual netcdf file name --> for each year, a running number will be concatenated
NA_val <- -32767 # NA value of the input file --> will be set to NAs
parName <- "t2m"
nrOfYears <- 35 # Please adjust here the amount of years over which the anomalies are calculated
varName <- "GP500_anomaly" # Name of the variable name as indicated in the netcdf file
varUnit <- "m**2 s**-2" # Unit of the variable stored in the netcdf file
varDescription <- "Geopotential at 500 mbar anomalies" # Description of the variable in the netcdf file
fileName <- "GP500_anomaly_" # Beginning of the individual netcdf file name --> for each year, a running number will be concatenated
NA_val <- -32767 # NA value of the input file --> will be set to NAs
parName <- "z"
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

stackVector <- createStackVector(fileList)
# Routine which restructures the stack based on month to stack based on years
# containing the monthly anomalies per year
for(j in 1: nrOfYears){
  tempStack <- stack()
  stackVector <- createStackVector(fileList)
  for(i in stackVector){
    tempStack <- addLayer(tempStack,i[[j]])
  }
  # make a copy of each yearly stack to the current work space
  assign(paste("t2m_",j,sep=""),tempStack)  
  assign(paste("GP500_",j,sep=""),tempStack)  
  
  # Create the actual netCDF output file
  netCDFFile <- nc_create(paste(fileName,j,".nc",sep=""), netCDF_varDef, force_v4=FALSE,verbose=FALSE)
  
  # fill the created netCDF file with the actual content --> calculated anomalies
  ncvar_put(netCDFFile, varName, as.matrix(brick(tempStack)),verbose=FALSE)
  # close the created netCDF file
  nc_close(netCDFFile)
  rm(tempStack,i)
}

#######################################################################################
# Testing environment for the created files
# Opening the netcdf, getting the variable stored and plot individual month of it
test <- open.ncdf("SLP_anomaly_2.nc",write=FALSE)
testVariable <- get.var.ncdf(test,"SLP_anomaly")
test <- open.ncdf("GP500_anomaly_1.nc",write=FALSE)
testVariable <- get.var.ncdf(test,"GP500_anomaly")
image(testVariable[,,1])

close.ncdf(test)
