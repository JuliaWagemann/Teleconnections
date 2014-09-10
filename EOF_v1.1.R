library(wq)
library(remote)
library(ncdf)
library(RNetCDF)
library(ncdf4)
library(fields)
library(raster)

#########################################################################
# Section for individual adjustment
# Please specify all variables in this section individually for every single parameter
#########################################################################

#how many eof you want to obtain?
NofEOFtoRetain <- 12

#what's the size of the input data? 
geoLat <- 181
geoLon <- 360

#If the dataset is too big, the RAM might not be enough to perform the calculation.
#You may want to interpolate the data to reduce the dataset size. What factor you want to use?
Intepolation_factor <- 3

#Where are the input netCDF files?
path <- "C:/Users/francesca cecinati/Documents/Teleconnection_Analysis/Data/2_Anomalies/ERA_Interim/Tair_2m"

#what is the variable name in the input netCDF files?
varNcdfName <- "t2m_anomaly"

#what name you want to give to the variable in the output netCDF file? What unit does it have? what description?
varEOFName <- "T2M.eof"
varUnit <- "deg"
varDesc <- "EOF of Temperature at 2m anomalies"

#where do you want to save the output? what name do you want to give to the netCDF and the associated csv?
outputPath <- "C:/Users/francesca cecinati/Documents/Teleconnection_Analysis/Data/5_EOF"
outputNCDFname <- "t2m_eof.nc"
csvName <- "T2M_eof.csv"
########################################################################

setwd(path)
   
fileList <- list.files(path, pattern=".nc")
nrOfTimeSteps <- length(fileList)*12

# for all files in fileList, the monthly entries are taken and written to a stack
st <- stack()
j <- 1  #initialize a count variable
for(i in 1:length(fileList)){
        anomaly <- open.ncdf(fileList[i], write=FALSE)
        var <- get.var.ncdf(anomaly, varid=varNcdfName)
        
        for(k in 1:12){
                temp <- raster(var[,,k])
                wtemp <- geoWeight(temp)
                atemp <- aggregate(wtemp, fact=Intepolation_factor, fun=mean, na.rm=TRUE)
                st <- addLayer(st, atemp)
               
        }
}

#if interpolation has been performed, the longitude and latitude are reduced
newx <- st@layers[[1]]@ncols
newy <- st@layers[[1]]@nrows


# the stack is transformed into a matrix for computational purposes
nrOfPixels <- length(st[[1]]@data@values)
unstacked <- unstack(st)
datamatrix <- matrix(nrow=nrOfTimeSteps, ncol=nrOfPixels, dimnames = list(c(1:nrOfTimeSteps), c(1:nrOfPixels)))

for (k in 1:nrOfTimeSteps){
        datamatrix[k,] <- unstacked[[k]]@data@values
}
                
# dataMatrix is scaled (normalized) 
datMat.sc <- scale(datamatrix)
                
# EOF calculation
result <- eof(datMat.sc, NofEOFtoRetain)
      
#the calculated rotated eof modes are stored into a stack
reof_st <- stack()

for (l in 1:NofEOFtoRetain){
        temp_eof_vect <- result$REOF[[l+1]]
        temp_eof_mat <- matrix(temp_eof_vect, newx, newy)
        temp_eof_mat <- temp_eof_mat[,ncol(temp_eof_mat):1]
        reof_st <- addLayer(reof_st, raster(temp_eof_mat))
}

#the calculated amplitude time series are stored in a data frame
eofDF <- data.frame(time_series=result$amplitude)
eofDF <- rbind(c(0, result$eigen.pct[1:NofEOFtoRetain]), eofDF)

#preparing to write the output files
setwd(outputPath)

#the amplitude time series are written in a csv
write.csv(eofDF, csvName)

# The EOF are written in a NetCDF file

#preparation of longitude and latitude dimensions
lat <- 1:Interpolation_factor:newy*Interpolation_factor
lat<- lat-9-Interpolation_factor
lon <- 1:Interpolation_factor:newx*Interpolation_factor
lon <- lon-Interpolation_factor
                
# Define the dimensions of the resulting netCDF file
x <- ncdim_def("lon", "degrees",lon)
y <- ncdim_def("lat", "degrees", lat)
eofs <- ncdim_def("eof", " ", c(1:NofEOFtoRetain))
                
# Define Variable which is written to the netCDF file
netCDF_varDef <- ncvar_def(varEOFName, units=varUnit,dim=list(y,x,eofs), missval=NA,
                 longname=varDesc, prec= 'double', verbose=FALSE)
                
# Create the actual netCDF output file
netCDFFile <- nc_create(outputNCDFname, netCDF_varDef, force_v4=FALSE,verbose=FALSE)
                
# fill the created netCDF file with the actual content
ncvar_put(netCDFFile, netCDF_varDef, as.matrix(brick(reof_st)),verbose=FALSE)
                
# close the created netCDF file
nc_close(netCDFFile)




