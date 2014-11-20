################################################################################
# Name:         SOM_v2.0R
# version:      v2.0
# Author:       Julia Wagemann
# Date:         2014/10/03
################################################################################
library(kohonen)
library(ncdf)
library(RNetCDF)
library(ncdf4)
library(fields)
library(raster)
# load necessary functions
source("C:/Users/julia_francesca/Desktop/Teleconnections/SOM_functions_1.1.R")

# path for input netcdf files
path <- "C:/Users/julia_francesca/Documents/Data/2_Anomalies/ERA_Interim/seasonal/3_summer/1_GP500/"
# path for output files
path1 <- "C:/Users/julia_francesca/Documents/Data/4_Results/SOM/ERA_interim/GP500/"

setwd(path)
fileList <- list.files(path, pattern="_1deg.nc")

######################################################################################################
# Call of function "SOMToNCDF"
# Function applies the machine learning technique "Self-organizing-maps" to netcdf files of an input
# file list and writes the resulting SOM nodes into a netCDF file.
#
# Following parameters have to be specified:
# fileList      -       list of netCDF files to be processed
# gridSize 1/2  -       specify the size of the resulting SOM grid
# neighbourhoodRadius   specify the radius of the neighbourhood for every SOM node (e.g. 5)
# trainingRate  -       specify the amount of training repetition of the SOM procedure (e.g. 500)
# path          -       specify the output path of the resulting files
# freqFileName  -       provide file name of .txt file for calculated SOM node frequencies
# timeUnit      -       specify the unit of the time dimension in the resulting netCDF file
# varName       -       specify variable name in the resulting netCDF file
# varUnit       -       specify variable unit in the resulting netCDF file
# varDescription-       provide detailed description of netCDF variable
# ncdfFileName  -       provide file name of resulting netCDF file
#####################################################################################################
SOMToNCDF(fileList,gridSize1=3, gridSize2=4,neighbourhoodRadius=5,trainingRate=500,
          path=path1,freqFileName="SOM_MAM_freqs.txt",timeUnit="SOM_nodes_MAM",
          varName="SOMs",varUnit="",varDescription="SOM nodes of ERA-interim spring TP anomalies",
          ncdfFileName="TP_SOM_MAM_1deg")
