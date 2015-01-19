###############################################################################
# Copyright 2015
# Author: Julia Wagemann
#
# This script uses EOT analysis to identify the main spatial drivers between 
# two given spatial data sets. Teleconnections of two given spatial datasets
# can therefore be calculated. The two variables are opposed to each other up
# to a lag time of 12 month (1 year).
#
# Input: spatial data fields in netCDF format (e.g. SST and Tair2m)
# Output: for each lag time combination in total 4 output files are generated
#         1. stack of response EOTs as netCDF
#         2. stack of predictor EOTs as netCDF
#         3. .txt containing EOT time series for each mode
#         4. .txt containing for each mode EV and lat/lon coordinates for base
#            point
#
# Requirements: anomalies_functions_v1.1.R for helper functions
#               libraries: remote, ncdf
# TODO:
###############################################################################

# set WD to load R-Script for helper functions
setwd("//TEC-JULIA-FRA/Users/julia_francesca/Desktop/Teleconnections")
source("anomalies_functions_v1.1.R")

# load libraries
library(remote)
library(ncdf)

# Path to response data set (equivalent to y-variable in a regression)
path_resp <- "P:/2_Anomalies/ERA_Interim/monthly/4_TP/3deg"
# Path to predictor data set (equivalent to x-variable in a regression)
path_pred <- "P:/2_Anomalies/CCI/monthly/SST/3deg"

# Path, where output files shall be stored
path_EOT <- "P:/4_ML_Results/2_EOT/predictor_response/"

# Name vector to name raster layers in a coherent and understandable way
nameVector <- c("Jan_", "Feb_", "Mar_", "Apr_", "May_", "Jun_", "Jul_","Aug_","Sep_","Oct_","Nov_","Dec_" )

# get lat/lon information from response field
setwd(path_resp)
fileList <- list.files(pattern=".nc")
lat_resp <- getLat(fileList[1])
lon_resp <- getLon(fileList[1])

# get lat/lon information from predictor field
setwd(path_pred)
fileList <- list.files(pattern=".nc")
lat_pred <- getLat(fileList[1])
lon_pred <- getLon(fileList[1])

# put all netCDF files of a fileList into a big stack and apply a conversion
# factor if applicable
st_resp <- createBigStack(path_resp,nameVector,1000)
st_pred <- createBigStack(path_pred,nameVector,1)

# define lag vector
lag_vec <- c(0:12)
# Loop through lag vector
for(i in lag_vec){
  lagged <- lagalize(st_pred,st_resp,lag=i,freq=12)
  
  # denoise the data, means apply a principal component filter to the data
  st_resp.dns <- denoise(lagged[[2]],expl.var=0.9)
  st_pred.dns <- denoise(lagged[[1]],expl.var=0.9)
  
  # amount of layers in the stack
  nrLayers <- nlayers(st_resp.dns)
  
  # calculate eot modes --> indicate amount of modes with n = ...
  modes <- eot(x = st_pred, y = st_resp.dns, 
               n = 3, standardised = FALSE, 
               reduce.both = FALSE, print.console = TRUE)
  
  # file name addition
  add <- paste("_lag_",i,sep="")
  
  # write calculated results as output 
  writeEOTOutput(path_EOT,modes,driver=TRUE,nrLayers,"SST","TP",lat_pred,lon_pred,lat_resp,lon_resp,add)
}
