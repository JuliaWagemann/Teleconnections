# Function "getLat" is a help function to retrieve the values for the latitude
# dimension of the input netcdf-files
getLat <- function(file){
        temp <- open.ncdf(file,write=FALSE)
        lat <- temp$dim$latitude$vals
        close.ncdf(temp)
        return (lat)
}

# Function "getLon" is a help function to retrieve the values for the longitude
# dimension of the input netcdf-files
getLon <- function(file){
        temp <- open.ncdf(file,write=FALSE)
        lon <- temp$dim$lon$vals
        close.ncdf(temp)
        return (lon)
}

# Function "getMonthStack"creates a stack of input data for one specific month,
# based on a fileList
# This function is called in the function calculateAnomaly
getMonthStack <- function(fileList, month){
  st <- stack()
  for(i in fileList){
    temp <- open.ncdf(i, write=FALSE)
    parName <- names(temp$var)
    NA_val <- temp$var[[1]]$missval
    
    tempRaster <- raster(get.var.ncdf(temp,names(temp$var)[month]))  
#    z <- raster(tempRaster[,,month])
#    z[z==NA_val] <- NA
    tempRaster[tempRaster==NA_val] <- NA
#    st <- addLayer(st,z)
    st <- addLayer(st,tempRaster)
    close.ncdf(temp)
  }
  return(st)
}


getSeasonalStack <- function(fileList,season,fun){
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

# Function "calculateAnomaly" calculates monthly anomalies based on the 
# climatology given by an input fileList and returns a stack based on the
# individual months
calculateAnomaly <- function(meanStack,yearsVec){
        climatology <- mean(meanStack[[yearsVec]])
        st <- stack()
        j <- 1
        for (j in 1:nlayers(meanStack)){
                actualAn <- raster(meanStack,layer=j)-climatology
                st <- addLayer(st,actualAn)
                }
        return(st)
}

anomalizeToNCDF <- function(fileList,timeStepVec,fun,yearsVec,seasonMonth,lat,lon,timeUnit,
                            varName,varUnit,varDescription,timeStepVec_ncdf,fileName){
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

writeStackToNCDF <-function(stack,lat,lon,timeUnit,varName,varUnit,varDescription,timeStepVec,fileName){
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

aggregate <- function(fileList,fun,aggVal,timeUnit,varName,varUnit,varDescription,timeStepVec,degrees){
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
                for(z in vec){
                        ncdf_var <- raster(get.var.ncdf(ncdf, z))
                        aggRaster <- aggregate(ncdf_var,fact=4,fun=mean,expand=FALSE,na.rm=FALSE)
                        st <- addLayer(st,aggRaster)
                }
                fileName <- paste(i,"_1deg.nc",sep="")
                writeStackToNCDF(st,lat_new,lon_new,timeUnit,varName,varUnit,varDescription,timeStepVec, fileName)
        }
}

aggregate(fileList=fileList,fun=fun,aggVal=4,timeUnit="month",varName="TP",
          varDescription="Monthly total precipitation in m",timeStepVec=c(1:12),degrees=1)
