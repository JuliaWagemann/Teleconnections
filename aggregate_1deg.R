library(ncdf)
library(RNetCDF)
library(ncdf4)
library(fields)
library(raster)


path <-"//TEC-JULIA-FRA/Users/julia_francesca/Documents/Data/1_Raw_Data_NCDF/ERA_Interim/Total_precipitation/"
setwd(path)

fileList <- list.files(path,pattern=".nc")

varName <- "SLP_anomaly_3deg"
varUnit <- "Pa"
varDesc <- "Sea level pressure anomalies in Pa - 3 deg spatial resolution"

lat <- rev(seq(-89,90,3))
lon <- seq(-179,180,3)

aggregate <- function(fileList,aggregateFactor,fun,lat,lon,)
for(i in fileList_precip){
        ncdf <- open.ncdf(paste(precip_path,i,sep=""), write=FALSE)

        ncdf_var <- get.var.ncdf(ncdf, "SLP_anomaly")
        st <- stack()
        for(z in 1:12){
                tempRaster <- raster(ncdf_var[,,z])
                aggRaster <- aggregate(tempRaster,fact=3,fun=mean,expand=FALSE,na.rm=FALSE)
                st <- addLayer(st,aggRaster)        
        }

        ncdfName <- paste(precip_path,substr(i,1,nchar(i)-3),"_3deg.nc",sep="")
             
        # Define the dimensions of the resulting netCDF file
        x <- ncdim_def("lon", "degrees",lon)
        y <- ncdim_def("lat", "degrees", lat)
        anomalies <- ncdim_def("Month", " ", c(1:12))
        
        # Define Variable which is written to the netCDF file
        netCDF_varDef <- ncvar_def(varName, units=varUnit,dim=list(x,y,anomalies), missval=NA,
                                   longname=varDesc, prec= 'double', verbose=FALSE)
        
        # Create the actual netCDF output file
        netCDFFile <- nc_create(ncdfName, netCDF_varDef, force_v4=FALSE,verbose=FALSE)
        
        # fill the created netCDF file with the actual content --> calculated anomalies
        ncvar_put(netCDFFile, netCDF_varDef, as.matrix(brick(st)),verbose=FALSE)
        
        # close the created netCDF file
        nc_close(netCDFFile)           
}
