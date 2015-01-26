###############################################################################
# Copyright 2015
# Author: Julia Wagemann
# Date: 2015-01-23
#
# This script uses EOT analysis to identify the main spatial patterns of one
# given spatial dataset. 
#
# Input: spatial data field in netCDF format (e.g. SST)
# Output: netCDF file of EOT output
#         1. stack of response EOTs as netCDF
#         2. .txt containing EOT time series for each mode
#         3. .txt containing for each mode EV and lat/lon coordinates for base
#            point
#
# Requirements: anomalies_functions_v1.1.R for helper functions
#               libraries: remote, ncdf, MASS
# TODO:
###############################################################################

# load R-Script for helper functions
source("H:anomalies_functions_v1.1.R")

# load libraries
library(remote)
library(ncdf)
library(MASS)

# Path to spatial data set
path_data <- "P:/2_Anomalies/CCI/monthly/SST/3deg/"
# Path to folder where output shall be stored
path_EOT <- "P:/4_ML_Results/2_EOT/"

nameVector <- c("Jan_", "Feb_", "Mar_", "Apr_", "May_", "Jun_", "Jul_","Aug_","Sep_","Oct_","Nov_","Dec_" )

# get lat/lon information from response field
setwd(path_data)
fileList <- list.files(pattern=".nc")

#get lat/lon information
lat_resp <- sort(getLat(fileList[1]),decreasing=TRUE)
lon_resp <- getLon(fileList[1])

# create big stack out of all input data
st_data <- createBigStack(path_data,nameVector,1)
layerNr <- nlayers(st_data)

# denoise the data
st_data.dns <- denoise(st_data,expl.var=0.9,weighted=FALSE)

# EOT analysis of the input data set --> amount of modes can be specified with n = ...
modes <- eot(st_data.dns,resp=NULL,n=6,standardized=FALSE,print.console=TRUE)

# write calculated EOT modes to netCDF file
writeEOTOutput(path_data,modes,driver=FALSE,layerNr,pred="SST",resp="SST",lat_resp=lat_resp,lon_resp=lon_resp,add="_noDriver")