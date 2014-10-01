


# Necessary libraries
library(raster)
library(ncdf)
library(fields)
library(ncdf4)

# Load functions needed for anomaly calculation
source("C:/Users/Julia Wagemann/Desktop/Teleconnections/anomalies_functions_v1.1.R")

# Set working directory
path <- "//TEC-JULIA-FRA/Users/julia_francesca/Documents/Data/1_Raw_Data_NCDF/ERA_Interim/Total_precipitation"
setwd(path)

fileList <- list.files(path)

lat <- getLat(fileList[1])
lon <- getLon(fileList[1])

######################################################################################################
# Call of function "anomalizeToNCDF"
# Following parameters have to be specified:
# fileList      -       list of netCDF files to be processed
# fun           -       specify aggregation method for calculation of seasonal anomalies (mean, sum)
# yearsVec      -       specify over what time period the anomalies should base on
# seasonMonth   -       specify if you want to calculate monthly or seasonal anomalies
# lat/lon       -       provide vector of latitudes / longitudes
# timeUnit      -       specify the unit of the time dimension in the resulting netCDF file
# varName       -       specify variable name in the resulting netCDF file
# varUnit       -       specify variable unit in the resulting netCDF file
# varDescription-       provide detailed description of netCDF variable
# timeStepVec_ncdf      provide number of layers in resulting netCDF file, e.g. seasonal anomalies of ERA-interim (1:35) and
#                       monthly anomalies of ERA-interim (1:12) 
# fileName      -       provide file name of resulting netCDF file
#####################################################################################################
anomalizeToNCDF(fileList=fileList,fun=sum,yearsVec=1:35,seasonMonth="month",lat=lat,lon=lon,
                timeUnit="month",varName="TP_anomaly",varUnit="m",
                varDescription="Total precipitation anomalies",
                timeStepVec_ncdf=c(1:12),fileName="TP_anomaly_monthly")

