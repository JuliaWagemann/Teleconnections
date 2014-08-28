#load necessary libaries

library(kohonen)
library(Reot)
library(ncdf)
library(RNetCDF)
library(ncdf4)
library(fields)
library(raster)
library(rworldmap)
library(RColorBrewer)

#set path
path <- "C:/Users/francesca cecinati/Documents/Teleconnection_Analysis/Data/2_Anomalies/ERA_Interim/Tair_2m"
setwd(path)

fileList <- list.files(path, pattern=".nc")
filenumber<-length(fileList)
set.seed(30)

#import from the ncdf

dataMatrix <- matrix(nrow=12*filenumber,ncol=65160)
count<-1

for (files in 1:filenumber){
        t2m_anomaly <- open.ncdf(fileList[files], write=FALSE)
        t2m_var <- get.var.ncdf(t2m_anomaly, "t2m_anomaly")
        
        for(i in 1:12){
                  
                temp <- t2m_var[,,i]
                dataMatrix[count,]<-temp
                count<-count+1
        }
}

#normalize anomalies
for (k in 1:65160){
        std=sd(dataMatrix[,k], na.rm=TRUE)
        dataMatrix[,k]<-dataMatrix[,k]/std
}


#define the dimensions of the som grid
somrows<-3
somcols<-3
somelements<-somrows*somcols

#compute the som
t2m.som <- som(data=dataMatrix, grid=somgrid(somrows,somcols,"rectangular"), rlen=500, alpha=c(0.05, 0.01),
               radius=5,toroidal=FALSE, n.hood="circular", keep.data=FALSE)

codes=t2m.som$codes

#reshape the results in a raster stack
st<-stack()

for(j in 1:somelements){
        mat<-matrix(codes[j,], 181, 360)
        mat<-mat[nrow(mat):1,] 
        r<-raster(mat)
        st<-addLayer(st, r)
}

#change path
path<- "C:/Users/francesca cecinati/Documents/Teleconnection_Analysis/Data/3_SOM/ERA_Interim/t2m_som"
setwd(path)

# NetCDF Output preparation
lat <- 1:181
lat<- lat-90
lon <- 1:360

# Define the dimensions of the resulting netCDF file
x <- ncdim_def("lon", "degrees",lon)
y <- ncdim_def("lat", "degrees", lat)
soms <- ncdim_def("soms", " ", c(1:somelements))

# Define Variable which is written to the netCDF file
netCDF_varDef <- ncvar_def("t2m.som", units= "deg",dim=list(x,y,soms), missval=NA,
                           longname="t2m.som", prec= 'double', verbose=FALSE)

# Create the actual netCDF output file
netCDFFile <- nc_create("som_norm.nc", netCDF_varDef, force_v4=FALSE,verbose=FALSE)

# fill the created netCDF file with the actual content --> calculated anomalies
ncvar_put(netCDFFile, netCDF_varDef, as.matrix(brick(st)),verbose=FALSE)

# close the created netCDF file
nc_close(netCDFFile)



# raster_test <- raster(GP_var[,,1])
# lat <- GP$dim$latitude$vals
# lon <- GP$dim$longitude$vals
# 
# offset = c(-max(GP$dim$longitude$vals)/2,min(GP$dim$latitude$vals))
# cellsize = c(abs(diff(GP$dim$longitude$vals[1:2])),abs(diff(GP$dim$latitude$vals[1:2])))
# offset = offset + cellsize/2
# 
# cells.dim = c(GP$dim$longitude$len, GP$dim$latitude$len)
# gt <- GridTopology(cellcentre.offset = offset, cellsize= cellsize, cells.dim=cells.dim)
# catMethod=seq(from=48000, to=58000,by=500)
# colourPalette=c(brewer.pal(20,"RdBu"))
# gridVals <- data.frame(as.vector(GP_jan))
# sGDF <-SpatialGridDataFrame(gt,data=gridVals)
# 
# mapParams <- mapGriddedData(sGDF, nameColumnToPlot='att', catMethod=catMethod,
#                             colourPalette=colourPalette,addLegend=FALSE)





# data('wines')
# 
# data('nir')
# set.seed(7)
# 
# nir.som <- som(data=spectra, grid=somgrid(3,3,"rectangular"), rlen=500, alpha=c(0.05, 0.01),
#                radius=5,toroidal=TRUE, n.hood="square", keep.data=TRUE)
# 
# image(nir.som$distances)

