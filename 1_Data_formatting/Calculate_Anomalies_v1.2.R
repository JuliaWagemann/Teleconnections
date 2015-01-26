%%R # to integrate this script into a IPython notebook 

################################################################################
# Name:         Calculate_Anomalies_v1.2.R
# Description:	This script calculates anomalies and writes them to a netCDF
#				output.
# version:      v1.2
# Author:       Julia Wagemann
# Date:         2014/01/23
################################################################################

# Necessary libraries
library(raster)
library(ncdf)
library(fields)
library(ncdf4)

# Load functions needed for anomaly calculation
source("H:/1_Data_formatting/anomalies_functions_v1.1.R")

# Set working directory
path <- "P:/1_Raw_Data_NCDF/ERA_Interim/Tair_2m/3deg/1990-2011"
setwd(path)

fileList <- list.files(path)

lat <- sort(getLat(fileList[1]),decreasing=TRUE)
lon <- getLon(fileList[1])

######################################################################################################
# Call of function "anomalizeToNCDF"
# Function calculates monthly or seasonal anomalies based on an input netCDF file list.
# Args:
# 	fileList:	list of netCDF files to be processed
#	fun:	specify aggregation method for calculation of seasonal anomalies (mean, sum)
# 	timeStepVec:	specify the month of the season (winter - c(12,1,2), spring - c(3:5), summer - c(6:8),
#                   autumn - c(9:11))
# 	yearsVec:	specify over what time period the anomalies should base on
# 	seasonMonth:	specify if you want to calculate monthly or seasonal anomalies, with either "month" or "season"
# 	lat/lon:	provide vector of latitudes / longitudes
# 	timeUnit:	specify the unit of the time dimension in the resulting netCDF file
# 	varName:	specify variable name in the resulting netCDF file
# 	varUnit:	specify variable unit in the resulting netCDF file
# 	varDescription:	provide detailed description of netCDF variable
# 	timeStepVec_ncdf:	provide number of layers in resulting netCDF file, e.g. seasonal anomalies of ERA-interim (1:35) and
#                       monthly anomalies of ERA-interim (1:12) 
# 	fileName:	provide file name of resulting netCDF file
#####################################################################################################
anomalizeToNCDF(fileList=fileList,timeStepVec=c(1:12),fun=mean,yearsVec=1:22,seasonMonth="month",lat=lat,lon=lon,
                timeUnit="year",varName="Tair_2m_anomaly",varUnit="degC",
                varDescription="Air temperature at 2m anomalies in degC",
                timeStepVec_ncdf=c(1:12),fileName="Tair2m_1990-2011_anomaly")

