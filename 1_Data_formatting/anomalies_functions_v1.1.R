################################################################################
# Name:         anomalies_functions_v1.1.R
# Description:	This script contains various helper functions for formatting netCDF
#				data, calculating anomalies, aggregating netCDF or writing raster data
#				to a netCDF output.
# version:      v1.1
# Author:       Julia Wagemann
# Date:         2014/01/23
################################################################################

# load required libraries
library(ncdf4)
library(raster)
library(MASS)

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
# netcdf-files
#
# Args:
#	file: 	netCDF file 
#
# Returns:
#	longitude values of a specified netCDF file
##############################################################################
        temp <- open.ncdf(file,write=FALSE)
        lon <- temp$dim$lon$vals
        close.ncdf(temp)
        return (lon)
}

getMonthStack <- function(fileList, month){
###############################################################################
# This function creates a stack of the specified input data for one specific
# month (e.g. 1 for January, ..., 12 for December) 
#
# Args:
# 	fileList: a list of ncdf files
#   month:	integer indicating the month of a year that you wish to create a stack
#			from
# Returns: 
#	stack of monthly raster data sets
###############################################################################
  st <- stack()
  for(i in fileList){
    temp <- open.ncdf(i, write=FALSE)
    parName <- names(temp$var)
    NA_val <- temp$var[[1]]$missval
    
    tempRaster <- get.var.ncdf(temp,parName) 
    z <- raster(tempRaster[,,month])
    z[z==NA_val] <- NA
    st <- addLayer(st,z)
    close.ncdf(temp)
  }
  return(st)
}

getSeasonalStack <- function(fileList,season,fun){
###############################################################################
# This function creates a stack of seasonal mean values of raster fields.
# It can be distinguished between winter (DJF), spring (MAM), summer (JJA) or
# autumn (SON)
# If winter season is specified, the resulting stack has one layer less, as the
# season mean is calculated on a multi-year basis, eg. Dec 1993, Jan-Feb 1994
#
# Args:
# 	fileList:	a list of netCDF files
#   season:	vector for season (either winter - c(12,1,2),
#           spring - c(3:5), summer - c(6:8), autumn - c(9:11))
#	fun:	function that should be applied to calculate the seasonal values, e.g
#			mean or sum
#
# Returns: 
#	stack of seasonal raster data sets
###############################################################################
        st <- stack()
        if(isTRUE(all.equal(season,c(12,1,2)))){
                for(i in 2:length(fileList)){
                        temp1 <- open.ncdf(fileList[i-1],write=FALSE)
                        NA_val <- temp1$var[[1]]$missval
                        tempRaster1 <- stack(fileList[i-1],bands=12)
                        tempStack <- stack(fileList[i],bands=c(1:2),varname=names(temp1$var))
                        tempStack <- addLayer(tempStack,tempRaster1)
                        seasonalRas <- stackApply(tempStack,indices=c(1),fun=fun)
                        seasonalRas[seasonalRas==NA_val] <- NA
                        st <- addLayer(st,seasonalRas)
                }
        } else {
                for(i in fileList){
                        temp <- open.ncdf(i, write=FALSE)
                        NA_val <- temp$var[[1]]$missval                       
                        seasonalRas <- stackApply(stack(i,bands=season,varname=names(temp$var)),indices=c(1),fun=fun)
                        seasonalRas[seasonalRas==NA_val] <- NA
                        st <- addLayer(st,seasonalRas)        
                }
        }
        return(st)       
}


calculateAnomaly <- function(meanStack,yearsVec){
###############################################################################
# This function calculates anomalies based on an input stack and a vector
# indicating the years of the climatology
#
# Args:
#	meanStack:	either monthly or seasonal raster stack 
#	yearsVec:	vector with years, e.g. 1:35 for year 1 to year 35
#
# Returns: 
#	stack of monthly or seasonal anomalies
###############################################################################
        climatology <- mean(meanStack[[yearsVec]])
        st <- stack()
        j <- 1
        for (j in 1:nlayers(meanStack)){
                actualAn <- raster(meanStack,layer=j)-climatology
                st <- addLayer(st,actualAn)
                }
        return(st)
}


writeStackToNCDF <-function(stack,lat,lon,timeUnit,varName,varUnit,varDescription,timeStepVec,fileName){
###############################################################################
# Writes a given raster stack to a netCDF file
#
# Args:  
#	stack:	stack that should be written to a .nc file
#   lat/lon:	lat/lon values for the output raster files
#   timeUnit:	timeUnit is a character identifying the unit of the time dimension
#   varUnit:	varUnit is a character identifying the unit of the variable
#   varDescription:	varDescription is a sentence describing the variables more specifically
#   timeStepVec:	timeStepVec is a vector indicating the number of layers to be put to the netCDF file
#   fileName:	fileName specifies the name of the resulting netCDF file
#
# Returns:
#	nothing is returned --> output is written to a netCDF file
###############################################################################
        # Define the dimensions of the resulting netCDF file
        x <- ncdim_def("lon", "degrees",lon)
        y <- ncdim_def("lat", "degrees", lat)
        time <- ncdim_def("Time", timeUnit, timeStepVec)
        
        # Define Variable which is written to the netCDF file
        netCDF_varDef <- ncvar_def(varName, units= varUnit,dim=list(x,y,time), missval=NA,
                                   longname=varDescription, prec= 'double', verbose=FALSE)
        
        netCDFFile <- nc_create(paste(fileName,".nc",sep=""), netCDF_varDef, force_v4=FALSE,verbose=FALSE)
        
        # fill the created netCDF file with the actual content --> calculated anomalies
        ncvar_put(netCDFFile, varName, as.matrix(brick(stack)),verbose=FALSE)
        # close the created netCDF file
        nc_close(netCDFFile)
}


anomalizeToNCDF <- function(fileList,timeStepVec,fun,yearsVec,seasonMonth,lat,lon,timeUnit,
                            varName,varUnit,varDescription,timeStepVec_ncdf,fileName){
###############################################################################
# This function combines several functions to calculate anomalies of raster fields, given
# by an input file list as .nc files and writes the calculated anomalies into new netCDF files
# It can be specified if monthly or seasonal anomalies shall be calculated
#
# Args:  
#	fileList:	list of files of which the anomalies shall be calculated
#   timeStepVec:	optional, only necessary for seasonal anomaly calculations
#   fun:	specify function for anomaly calculation, e.g mean or sum(e.g for precipitation)
#   yearsVec:	vector of years on which anomalies base on (time series of climatology)
#   seasonMonth:	character specifying if seasonal or monthly anomalies shall be calculated
#   lat/lon:	spatial dimension of the resulting ncdf file
#   timeUnit:	character identifying the unit of the time dimension
#   varUnit:	character identifying the unit of the variable
#   varDescription:	sentence describing the variables more specifically
#   timeStepVec_ncdf:	vector indicating the number of layers to be put to the netCDF file
#   fileName:	specifies the name of the resulting netCDF file
#
# Returns: 
#	nothing is returned --> output is written to a netCDF file
###############################################################################
        st <- stack()
        if(seasonMonth=="season"){
                st <- getSeasonalStack(fileList,timeStepVec,fun)
                st_an <- calculateAnomaly(st,yearsVec)
                writeStackToNCDF(st_an,lat,lon,timeUnit,varName,varUnit,varDescription,timeStepVec_ncdf,fileName)
                
        } else {
                stackVector <- c()
                for(i in 1:12){
                        st <- getMonthStack(fileList,i)
                        temp <- calculateAnomaly(st,yearsVec)
                        stackVector[i] <- c(temp)
                }
                for(j in 1:nlayers(temp)){
                        tempStack <- stack()
                        for(k in stackVector){
                                tempStack <- addLayer(tempStack,t(k[[j]]))
                        }
                        fileNameTemp <- paste(fileName,"_",j,sep="")
                        writeStackToNCDF(tempStack,lat,lon,timeUnit,varName,varUnit,varDescription,timeStepVec_ncdf,fileNameTemp)                        
                        rm(tempStack,k)
                }  
        }        
}

aggregate1 <- function(fileList,fun,aggVal,timeUnit,varName,varUnit,varDescription,timeStepVec,degrees,ncdfMode,addFileName){
######################################################################################################
# Aggregation of raster datasets (given by a netCDF file list) to a coarser resolution
#
# Args:
#	fileList:	list of netCDF files to be processed
#   fun:	specify aggregation method (in general: "mean",e.g. for temperature, or "sum" for precipitation 
#			fields)
#   aggVal:	specify the aggregation value, factor of aggregation, e.g from 1 deg to 3 deg is factor 3
#   timeUnit:	specify the unit of the time dimension in the resulting netCDF file
#   varName:	specify variable name in the resulting netCDF file
#   varUnit:	specify variable unit in the resulting netCDF file
#   varDescription:	provide detailed description of netCDF variable
#   timeStepVec:	provide number of layers in resulting netCDF file (in general: c(1:12))
#   degrees:	aggregation value for lat/lon vectors
#   ncdfMode:	how is the structure of the input netCDF files - e.g "multiple", if the netcdf file has
#   			one variable with various layers (e.g. ERA-interim) and "single", if the netcdf file 
#				has various layers (e.g. ESA CCI data)
#   addfileName:	specify the addition to the resulting file name
#
# Returns:
#	nothing is returned --> output is written to a netCDF file 
#####################################################################################################
        for (i in fileList){
                ncdf <- open.ncdf(paste(path,i,sep=""),write=FALSE)
                
                lat <- ncdf$dim$lat$vals
                lon <- ncdf$dim$lon$vals
                lat_first <- mean(lat[1:aggVal])
                lat_last <- mean(lat[(length(lat)-(aggVal-1)):length(lat)])
                lon_first <- mean(lon[1:aggVal])
                lon_last <- mean(lon[(length(lon)-(aggVal-1)):length(lon)])
                lon_new <- seq(lon_first,lon_last,degrees)
                lat_new <- seq(lat_first,lat_last,degrees)

                vec <- names(ncdf$var)
                st <- stack()
                if(ncdfMode=="single"){
                        for(z in vec){
                                ncdf_var <- raster(get.var.ncdf(ncdf, z))
                                aggRaster <- aggregate(ncdf_var,fact=aggVal,fun=fun,expand=FALSE,na.rm=FALSE)
                                st <- addLayer(st,t(aggRaster))
                        }                       
                } else {
                        for(z in 1:12){
                                ncdf_var <- get.var.ncdf(ncdf,vec[1])
                                aggRaster <- aggregate(raster(ncdf_var[,,z]),fact=aggVal,fun=fun,expand=FALSE,na.rm=FALSE)
                                st <- addLayer(st,t(aggRaster))
                        }
                }
                fileName <- paste(substr(i,1,nchar(i)-3),addFileName,sep="")
                writeStackToNCDF(st,lat_new,lon_new,timeUnit,varName,varUnit,varDescription,timeStepVec, fileName)
        }
}

createBigStack <- function(path,nameVector,convFactor){
######################################################################################################
# Function puts all netCDF files of a given fileList into one rasterStack and applies, if needed, a
# conversion factor to the stack
#
# Args:
# 	path: path to netCDF files to be processed
#   nameVector: specify vector of layer names
#   convFactor: specify, if needed, a converstion factor the stack is multiplied with
#
# Returns: 
#	one big stack with all every layer of all ncdf files
#####################################################################################################
        setwd(path)
        fileList <- list.files(path,pattern=".nc")
        st <- stack()
        j <- 1
       
        for(i in fileList){
                temp <- stack(i)
                names(temp) <- paste(nameVector,j,sep="")
                st <- addLayer(st,temp)               
                j <- j+1
        }
        st <- st*convFactor
        return (st)
}

writeEOTOutput <- function(path,EOTStack,driver=FALSE,nrLayers,pred,resp,lat_pred,lon_pred,lat_resp,lon_resp,add){
######################################################################################################
# Function to write results of an EOT analysis (regression-based or not) as output. For each
# analysed lag time step two netCDF files containing predictor and response EOT modes, one .txt file with
# EOT time series of each mode and one .txt file containing EV (explained variance) and lat/lon coordinates
# of the identified base point will be written as output.
# 
# Args:
#	path: path where the output should be written
#   EOTStack: resulting object of an EOT analysis
#   driver:	EOT caclulation mode - if driver=FALSE, the regression-based approach
#           is used
#   nrLayers: number of layers of the created big stack
#   pred: string of predictor variable, e.g. "SST"
#   resp: string of response variable, e.g. "Tair2m"
#   lat_pred / lon_pred: lat / lon data from predictor field
#   lat_resp / lon_resp: lat / lon data from response field
#   add: supplement for file naming extension
#
# Returns:
#	nothing is returned --> output is written to a netCDF / .txt files
#####################################################################################################
  
  # Initialize response stack
  responseStack <- stack() 
  nrOfModes <- nmodes(EOTStack)
  
  # Initialize two matrix for .txt file outputs
  matrix1 <- matrix(nrow=nrLayers,ncol=nmodes(EOTStack))
  matrix2 <- matrix(nrow=nrOfModes,ncol=3)

  # Loop through amount of modes in EOTStack
  for(i in 1:nrOfModes){
    # extract first mode of EOTStack
    tempMode <- EOTStack@modes[[i]]
    if(driver==TRUE){
      predictStack <- stack()
      # Extract predictor EOT
      predictStack <- addLayer(predictStack,tempMode@rsq_predictor)
    } else
    # Extract response EOT
    responseStack <- addLayer(responseStack,tempMode@r_response)
    # add EOT ts to first matrix
    matrix1[,i] <- tempMode@eot
    # add EOT bp and EV to 2nd matrix
    matrix2[i,] <- c(tempMode@cum_exp_var,tempMode@coords_bp)  
  }
  if(exists("predictStack")==TRUE){
    # write predictor stack to netCDF file
    writeStackToNCDF(stack=predictStack,lat=lat_pred,lon=lon_pred,
                     timeUnit="EOTmodes",varName=paste("EOT_predictor_",pred,sep=""),
                     varUnit="",
                     varDescription='coefficient of determination between base point and each pixel of the predictor domain',
                     timeStepVec=c(1:nrOfModes),
                     fileName=paste(path,"predictor_",pred,"_",resp,add,sep=""))
  } else{}
    # write response stack to netCDF file
  
  writeStackToNCDF(stack=responseStack,lat=lat_resp,lon=lon_resp,
                     timeUnit="EOTmodes",varName="EOT_response",
                     varUnit="",
                     varDescription='correlation coefficients between base point and each pixel of the predictor domain',
                     timeStepVec=c(1:nrOfModes),
                     fileName=paste(path,"response_",pred,"_",resp,add,sep=""))    
 
  # write both matrices to file
  write.matrix(matrix1,file=paste(path,"EOT_ts_",pred,"_",resp,add,".txt",sep=""))
  write.matrix(matrix2,file=paste(path,"EOT_EV_BP_",pred,"_",resp,add,".txt",sep=""))  
}