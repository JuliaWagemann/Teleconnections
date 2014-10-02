################################################################################
# Name:         Calculate_Anomalies_v1.2.R
# version:      v1.2
# Author:       Julia Wagemann
# Date:         2014/10/01
################################################################################

# Necessary libraries
library(raster)
library(ncdf)
library(fields)
library(ncdf4)

# Load functions needed for anomaly calculation
source("C:/Users/Julia Wagemann/Desktop/Teleconnections/anomalies_functions_v1.1.R")

# Set working directory
path <- "//TEC-JULIA-FRA/Users/julia_francesca/Documents/Data/1_Raw_Data_NCDF/CCI/SLA/3deg/"
setwd(path)

fileList <- list.files(path)

lat <- sort(getLat(fileList[1]),decreasing=TRUE)
lon <- getLon(fileList[1])

######################################################################################################
# Call of function "anomalizeToNCDF"
# Function calculates monthly or seasonal anomalies based on an input netCDF file list.
# Following parameters have to be specified:
# fileList      -       list of netCDF files to be processed
# fun           -       specify aggregation method for calculation of seasonal anomalies (mean, sum)
# timeStepVec   -       specify the month of the season (winter - c(12,1,2), spring - c(3:5), summer - c(6:8),
#                       autumn - c(9:11))
# yearsVec      -       specify over what time period the anomalies should base on
# seasonMonth   -       specify if you want to calculate monthly or seasonal anomalies, with either "month" or "season"
# lat/lon       -       provide vector of latitudes / longitudes
# timeUnit      -       specify the unit of the time dimension in the resulting netCDF file
# varName       -       specify variable name in the resulting netCDF file
# varUnit       -       specify variable unit in the resulting netCDF file
# varDescription-       provide detailed description of netCDF variable
# timeStepVec_ncdf      provide number of layers in resulting netCDF file, e.g. seasonal anomalies of ERA-interim (1:35) and
#                       monthly anomalies of ERA-interim (1:12) 
# fileName      -       provide file name of resulting netCDF file
#####################################################################################################
anomalizeToNCDF(fileList=fileList,timeStepVec=c(9:11),fun=mean,yearsVec=1:18,seasonMonth="season",lat=lat,lon=lon,
                timeUnit="year",varName="SLA_anomaly_SON",varUnit="m",
                varDescription="Autumn mean sea level height anomalies in m",
                timeStepVec_ncdf=c(1:18),fileName="SLA_anomaly_3deg_SON")

