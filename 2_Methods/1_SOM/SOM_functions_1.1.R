################################################################################
# Name:         SOM_functions_v1.1.R
# Description:	This script contains various helper functions for formatting and
#				manipulating netCDF raster data sets in order to make a SOM
#				analysis
# version:      v1.1
# Author:       Julia Wagemann
# Date:         2014/10/03
################################################################################


getLat <- function(file){
###############################################################################
# Get-function to retrieve the values for the latitude dimension of the input 
# netcdf-file
#
# Args:
#	file:	netCDF file
#
# Returns:
#	latitude values of a specified netCDF file
##############################################################################
        temp <- open.ncdf(file,write=FALSE)
        lat <- temp$dim$lat$vals
        close.ncdf(temp)
        return (lat)
}

getLon <- function(file){
###############################################################################
# Get-function to retrieve the values for the longitude dimension of the input 
# netcdf-file
#
# Args:
#	file:	netCDF file
#
# Returns:
#	longitude values of a specified netCDF file
##############################################################################
        temp <- open.ncdf(file,write=FALSE)
        lon <- temp$dim$lon$vals
        close.ncdf(temp)
        return (lon)
}

getLayers <- function(file){
###############################################################################
# Get-function to retrieve the values for the third
# (time) dimension of the input netcdf-files
#
# Args:
#	file:	input netCDF file
#
# Returns:
#	number of layers (third dimension) of an netCDF file
##############################################################################
        temp <- open.ncdf(file,write=FALSE)
        nlayers <- temp$dim$Time$vals
        close.ncdf(temp)
        return(nlayers)
}

getSomDataMatrix <- function(fileList,nrOfPixels,nrOfTimeSteps,nlayers){
###############################################################################
# This function prepares the input raster files of a given fileList for the
# calculation of SOM nodes. Each raster file will be reshaped to a 1D-row-vector.
# All row vectors will be put together and a matrix is built, with number of pixels
# as the amount of columns and number of raster files (time steps) of the entire
# fileList as number of rows.
#
# Args:
#	fileList:	
#   nrOfPixels:	nrOfPixels of raster files
#   nrOfTimeSteps:	length of file list * number of layers
# 					in each ncdf file
#   nlayers:	amount of layers in ncdf file
#
# Returns: 
#	matrix with all 1D-row-vectors
###############################################################################
        t <- 1
        dataMatrix <- matrix(ncol=nrOfPixels, nrow=nrOfTimeSteps)
        for(i in fileList){              
                tempNCDF <- open.ncdf(paste(path,i,sep=""), write=FALSE)
                parName <- names(tempNCDF$var)
                NCDFVar <- get.var.ncdf(tempNCDF, parName)
                for(z in 1:nlayers){
                        tempRaster <- NCDFVar[,,z]
                        dataMatrix[t,] <- tempRaster
                        t <- t+1 
                }
        }
        return(dataMatrix)
}


scaleDataMatrix <- function(dataMatrix,nrOfPixels){
###############################################################################
# This function normalizes the data matrix of the anomalies with the help of the
# standard deviation
#
# Args:
#	dataMatrix:
#	nrOfPixels:
# Parameters required:  - dataMatrix -> as resulting from function getSomDataMatrix
#                       - nrOfPixels of raster files
#
# Returns: 
#	matrix with all 1D-row-vectors
###############################################################################
        for (k in 1:nrOfPixels){
                std=sd(dataMatrix[,k], na.rm=TRUE)
                dataMatrix[,k]<-dataMatrix[,k]/std
        }
        return(dataMatrix)
}

writeFreqs <- function(SOM,path,fileName){
###############################################################################
# This function calculates the frequency occurences of each SOM node pattern
# and writes the results into a .txt file
#
# Args:
#	SOM:	SOM object
#	path:	path for the output .txt file
#	fileName:	file name of the .txt file
#
# Returns:
# no object is returned --> results will be written to a .txt output file
###############################################################################
        outputFile <- file(paste(path,fileName),open="w")
        freq <- table(SOM$unit.classif)/length(SOM$unit.classif)
        writeLines(as.character(freq),outputFile,sep="\n")
        close(outputFile)
}

reshapeToStack <- function(somCodes,lat,lon){
###############################################################################
# This function reshapes the 1D-row vectors based on the code variable of the SOM
# object to better interpretable matrix raster files.
#
# Args:
#	somCodes:	code vectors of SOM object
#	lat:	nrow of resulting raster data set
#	lon:	ncol of resulting raster data set
#
# Returns: 
#	stack of all reshaped code vectors of the SOM object
###############################################################################
        st <- stack()
        for(i in 1:nrow(somCodes)){
                mat <- matrix(somCodes[i,],nrow=lat,ncol=lon)
                r <- raster(mat)
                st <- addLayer(st,r)
        }
        return(st)
}


writeStackToNCDF <-function(stack,lon,lat,timeUnit,varName,varUnit,varDescription,timeStepVec,fileName){
###############################################################################
# This function writes a given raster stack to a netCDF file
#
# Args:
#	stack:	stack to write as .nc file
#	lat/lon:	lat/lon values for the output raster files
#	timeUnit:	character identifying the unit of the time dimension
#   varUnit:	character identifying the unit of the variable
#   varDescription:	is a sentence describing the variables more in specific
#   timeStepVec:	is a vector indicating the number of layers to be put to the netCDF file
#   fileName:	specifies the name of the resulting netCDF file
# Returns:
# 	no object is returned --> results will be written to a netCDF file 
###############################################################################
        # Define the dimensions of the resulting netCDF file
        x <- ncdim_def("lon", "degrees",lon)
        y <- ncdim_def("lat", "degrees", lat)
        time <- ncdim_def("Time", timeUnit, timeStepVec)
        
        # Define Variable which is written to the netCDF file
        netCDF_varDef <- ncvar_def(varName, units= varUnit,dim=list(y,x,time), missval=NA,
                                   longname=varDescription, prec= 'double', verbose=FALSE)
        
        netCDFFile <- nc_create(paste(fileName,".nc",sep=""), netCDF_varDef, force_v4=FALSE,verbose=FALSE)
        
        # fill the created netCDF file with the actual content --> calculated anomalies
        ncvar_put(netCDFFile, varName, as.matrix(brick(stack)),verbose=FALSE)
        # close the created netCDF file
        nc_close(netCDFFile)
}

SOMToNCDF <- function(fileList,gridSize1,gridSize2,neighbourhoodRadius,trainingRate,path,freqFileName,
                      timeUnit,varName,varUnit,varDescription,ncdfFileName){
###############################################################################
# This function is the main function to calculate Self-Organizing-Map nodes based
# on an input fileList. All the results are written to a ncdf file.
#
# Args:
#	fileList:	list of raster time series data for SOM calculation
#	gridSize1:	Nr. of rows of resulting SOM object
#	gridSize2:	Nr. of cols of resulting SOM object
#   neighbourhoodRadius:	radius of the SOM nodes
#   trainingRate:	training rate -> how often is the SOM presented
#   path:	path of the output files
#   freqFileName:	file name of the .txt file containing SOM frequencies
#	timeUnit:	time unit of the resulting ncdf file
#   varName:	name of the variable of the ncdf file
#   varUnit:	unit of the ncdf file variable
#   varDescription:	description of the ncdf variable more in detail
#   ncdfFileName:	name of the ouptut ncdf file
# Returns:
#	no object is returned --> results will be written to a netCDF file 
###############################################################################
        lat <- getLat(fileList[1])
        lon <- getLon(fileList[1])
        nlayers <- nrow(getLayers(fileList[1]))
        
        nrOfPixels <- nrow(lat)*nrow(lon)
        nrOfTimeSteps <- length(fileList)*nlayers
        
        dataMatrix <- getSomDataMatrix(fileList,nrOfPixels=nrOfPixels,nrOfTimeSteps=nrOfTimeSteps,nlayers=nlayers)
        dataMatrix_scale <- scaleDataMatrix(dataMatrix,nrOfPixels)
        set.seed(20)
        SOM <- som(data=dataMatrix_scale, grid=somgrid(gridSize1,
                                                       gridSize2,"rectangular"), 
                   rlen=trainingRate, alpha=c(0.05, 0.01), radius=neighbourhoodRadius, toroidal=FALSE, 
                   n.hood="circular", keep.data=TRUE)
        
        writeFreqs(SOM,path,freqFileName)
        codes <- SOM$codes
        
        codeStack <- reshapeToStack(somCodes=codes,lat=nrow(lon),lon=nrow(lat))
        writeStackToNCDF(stack=codeStack,lat=lat,lon=lon,timeUnit=timeUnit,varName=varName,varUnit=varUnit,
                         varDescription=varDescription,
                         timeStepVec=c(1:nrow(codes)),fileName=ncdfFileName)        
}