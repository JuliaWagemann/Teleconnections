library(R.matlab)
library(raster)

# load R-Script with additional funtions
source("C:/Users/julia_francesca/Desktop/Teleconnections/anomalies_functions_v1.1.R")

# Path of general working directory
path <- "C:/Users/julia_francesca/Documents/3_Data/1_Raw_Data_NCDF/ET_5km/"

# Set WD to .mat files with lat/lon information
setwd(paste(path,"Info/",sep=""))
fileList_extent <- list.files(pattern=".mat")
# read .mat files
lat <- readMat(fileList_extent[1])
lon <- readMat(fileList_extent[2])

# Set WD to .mat files containing data
setwd(paste(path,"matFiles/",sep=""))
fileList_data <- list.files()

j = 1 # running variable
st <- stack()

# Loop through data list --> .mat file data are accumulated in a raster stack and every 12 month (each year)
#                            the raster stack is written as output
for(i in fileList_data){
  # read .mat file
  temp <- readMat(i)
  # rasterize loaded .mat file data
  tempRast <- raster(temp$ETmon,xmn=min(lon$Longitude),xmx=max(lon$Longitude),
                     ymn=min(lat$Latitude),ymx=max(lat$Latitude))
  if(maxValue(tempRast)>400 & minValue(tempRast)<0){
    tempRast[tempRast > 400] <- NA
    tempRast[tempRast < 0] <- NA
  }
  # rename data in raster file
  names(tempRast) <- substring(i,1,6)
  
  # if one year is complete, write stack to output
  if(j>12){
    res <- 0.05
    latVals <- sort(seq(min(lat$Latitude),max(lat$Latitude)-res,by=res),decreasing=TRUE)
    lonVals <- seq(min(lon$Longitude),max(lon$Longitude)-res,by=res)
    writeStackToNCDF(st,lat=latVals,lon=lonVals,timeUnit="Month",varName="ET",varUnit="mm",
                     varDescription="Global Evapotranspiration (ET) with 5km resolution",timeStepVec=c(1:12),
                     fileName=paste(path,"results/",substring(t,1,4),sep=""))
    rm(st) # remove stack due to capacity issues
    st <- stack()
    j <- 1 # set running variable back to one
    st <- addLayer(st,tempRast)
    j <- j+1
  } else {
    st <- addLayer(st,tempRast)
    j <- j+1
    t <- i # helper to store the right year information in order to name the ncdf file appropriate
  }
}




