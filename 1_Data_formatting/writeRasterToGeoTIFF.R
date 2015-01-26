library(raster)
library(ncdf)
library(ggplot2)
library(kohonen)
library(fields)
library(rgdal)

path <- "//TEC-JULIA-FRA/Users/julia_francesca/Documents/Data/2_Anomalies/CCI/monthly/SLA/1deg"
setwd(path)
fileList <- list.files(path,pattern="_1deg.nc")
test <- stack(fileList[1],bands=c(1:12),varname='SLA')

i <- 1
j ,- 1
for(i in 1:12){
        j <- writeRaster(test[[i]],fileName=paste("test_",i,".tif",sep=""),format="Gtiff",overwrite=TRUE)
}