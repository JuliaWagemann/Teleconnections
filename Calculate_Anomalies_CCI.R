# Necessary libraries
library(raster)
library(ncdf)
library(fields)
library(ncdf4)

# Load functions needed for anomaly calculation --> if necessary adjust path to source code
source("F:/Julia/Data/anomalies_functions_CCI.R")

# Set working directory
path <- "F:/Julia/Data/CCI/Aerosol/"
setwd(path)


###############################################################################
# Section for individual adjustment
# Please specify the all the variables in this section individually for anomaly
# calculation and netcdf writing
###############################################################################
nrOfYears <- 8 # Please adjust here the amount of years over which the anomalies are calculated
varName <- "AOD_anomaly" # Name of the variable name as indicated in the netcdf file
varUnit <- "" # Unit of the variable stored in the netcdf file
varDescription <- "Atmospheric optical depth anomaly" # Description of the variable in the netcdf file
fileName <- "AOD_anomaly_" # Beginning of the individual netcdf file name --> for each year, a running number will be concatenated
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
test <- open.ncdf("SST_anomaly_3.nc",write=FALSE)
testVariable <- get.var.ncdf(test,"SST_anomaly")
image(testVariable[,,5])

close.ncdf(test)
