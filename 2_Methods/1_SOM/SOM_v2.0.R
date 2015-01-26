################################################################################
# Name:         SOM_v2.0.R
# Description:	This script calls the function SOMtoNCDF that calculates SOM for
#				for a given netCDF raster data file list. The output SOM object
#				is written to a netCDF file.
# version:      v2.0
# Author:       Julia Wagemann
# Date:         2014/12/02
################################################################################
# load necessary functions
source("H:/2_Methods/2_SOM/SOM_functions_1.1.R")

# load necessary packages
library(kohonen)
library(ncdf)
library(raster)

# path for input netcdf files
path <- "P:/2_Anomalies/ERA_Interim/seasonal/3_summer/1_GP500/"
# path for output files
path1 <- "P:/4_Results/SOM/ERA_interim/GP500/"

setwd(path)
fileList <- list.files(path, pattern="_1deg.nc")

###############################################################################
# Call of function "SOMToNCD"
# Function applies the machine learning technique "Self-organizing-maps" to netcdf files of an input
# file list and writes the resulting SOM nodes into a netCDF file.
#
# Args:
#	fileList:	list of raster time series data for SOM calculation stored in netCDF format
#	gridSize1:	Nr. of rows of resulting SOM object
#	gridSize2:	Nr. of cols of resulting SOM object
#   neighbourhoodRadius:	radius of the SOM nodes (e.g. 5)
#   trainingRate:	training rate -> how often is the SOM presented (e.g. 500)
#   path:	path of the output files
#   freqFileName:	file name of the .txt file containing SOM node frequencies
#	timeUnit:	time unit of the resulting ncdf file
#   varName:	name of the variable of the ncdf file
#   varUnit:	unit of the ncdf file variable
#   varDescription:	description of the ncdf variable more in detail
#   ncdfFileName:	name of the ouptut ncdf file
# Returns:
#	no object is returned --> results will be written to a netCDF file 
###############################################################################
SOMToNCDF(fileList,gridSize1=3, gridSize2=4,neighbourhoodRadius=5,trainingRate=500,
          path=path1,freqFileName="SOM_MAM_freqs.txt",timeUnit="SOM_nodes_MAM",
          varName="SOMs",varUnit="",varDescription="SOM nodes of ERA-interim spring TP anomalies",
          ncdfFileName="TP_SOM_MAM_1deg")