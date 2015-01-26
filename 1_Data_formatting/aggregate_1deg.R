################################################################################
# Name:         aggregate_1deg.R
# Description:	Raster datasets in netCDF format will be aggregated to a coarser
#				resolution.
# version:      v1.0
# Author:       Julia Wagemann
# Date:         2014/10/01
################################################################################
library(ncdf)
library(RNetCDF)
library(ncdf4)
library(fields)
library(raster)

source("H:/anomalies_functions_v1.1.R")

path <-"P:/1_Raw_Data_NCDF/CCI/SLA/1deg/"
setwd(path)

fileList <- list.files(path,pattern=".nc")

######################################################################################################
# Call of function "aggregate1"
# Function aggregates raster datasets (given by a netCDF file list) to a coarser resolution
# Following parameters have to be specified:
# fileList      -       list of netCDF files to be processed
# fun           -       specify aggregation method (in general: "mean", e.g. "sum" for precipitation fields)
# aggVal        -       specify the aggregation value, factor of aggregation
# timeUnit      -       specify the unit of the time dimension in the resulting netCDF file
# varName       -       specify variable name in the resulting netCDF file
# varUnit       -       specify variable unit in the resulting netCDF file
# varDescription-       provide detailed description of netCDF variable
# timeStepVec   -       provide number of layers in resulting netCDF file (in general: c(1:12))
# degrees       -       aggregation value for lat/lon vectors
# ncdfMode      -       how is the structure of the input netCDF files - e.g "multiple", if the netcdf file has
#                       one variable, with various layers and "single", if the netcdf file has various layers
# addfileName   -       specify what is added to the resulting file name
#####################################################################################################
aggregate1(fileList=fileList,fun=mean,aggVal=3,timeUnit="month",varName="SLA",varUnit="m",
           varDescription="Monthly sea level anomalies in m",timeStepVec=c(1:12),degrees=3,ncdfMode="multiple",
           addFileName="_3deg")
