library(kohonen)
library(ncdf)
library(RNetCDF)
library(ncdf4)
library(fields)
library(raster)
library(rworldmap)
library(RColorBrewer)

path <- "F:/Julia/Data/ERA_interim/Anomalies/GP500/"
setwd(path)

fileList <- list.files(path, pattern=".nc")

set.seed(300)

GP_anomaly <- open.ncdf(fileList[1], write=FALSE)
GP_var <- get.var.ncdf(GP_anomaly, "GP500_anomaly")
dataMatrix <- matrix(nrow=12,ncol=65160)

matrixList <- list()
for(i in 1:12){
#       
        tempRaster <- matrix(GP_var[,,i])
        matrixList <- list(matrixList,tempRaster)

#dataMatrix[i,] <- temp
}


nir.som <- supersom(data=matrixList, grid=somgrid(3,3,"rectangular"), rlen=500, alpha=c(0.05, 0.01),
               radius=5,toroidal=FALSE, n.hood="circular", keep.data=FALSE)











raster_test <- raster(GP_var[,,1])
lat <- GP$dim$latitude$vals
lon <- GP$dim$longitude$vals

offset = c(-max(GP$dim$longitude$vals)/2,min(GP$dim$latitude$vals))
cellsize = c(abs(diff(GP$dim$longitude$vals[1:2])),abs(diff(GP$dim$latitude$vals[1:2])))
offset = offset + cellsize/2

cells.dim = c(GP$dim$longitude$len, GP$dim$latitude$len)
gt <- GridTopology(cellcentre.offset = offset, cellsize= cellsize, cells.dim=cells.dim)
catMethod=seq(from=48000, to=58000,by=500)
colourPalette=c(brewer.pal(20,"RdBu"))
gridVals <- data.frame(as.vector(GP_jan))
sGDF <-SpatialGridDataFrame(gt,data=gridVals)

mapParams <- mapGriddedData(sGDF, nameColumnToPlot='att', catMethod=catMethod,
                            colourPalette=colourPalette,addLegend=FALSE)





data('wines')

data('nir')
set.seed(7)

nir.som <- som(data=spectra, grid=somgrid(3,3,"rectangular"), rlen=500, alpha=c(0.05, 0.01),
               radius=5,toroidal=TRUE, n.hood="square", keep.data=TRUE)

image(nir.som$distances)

