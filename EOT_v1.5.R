source("H:anomalies_functions_v1.1.R")

library(remote)
library(ncdf)
library(MASS)

# Path to response data set (equivalent to y-variable in a regression)
path_data <- "P:/2_Anomalies/CCI/monthly/SST/3deg/"
path_EOT <- "P:/4_ML_Results/2_EOT/"

nameVector <- c("Jan_", "Feb_", "Mar_", "Apr_", "May_", "Jun_", "Jul_","Aug_","Sep_","Oct_","Nov_","Dec_" )

# get lat/lon information from response field
setwd(path_data)
fileList <- list.files(pattern=".nc")
lat_resp <- sort(getLat(fileList[1]),decreasing=TRUE)
lon_resp <- getLon(fileList[1])

st_data <- createBigStack(path_data,nameVector,1)
layerNr <- nlayers(st_data)
st_data.dns <- denoise(st_data,expl.var=0.9,weighted=FALSE)

modes <- eot(st_data.dns,resp=NULL,n=6,standardized=FALSE,print.console=TRUE)

writeEOTOutput(path_data,modes,driver=FALSE,layerNr,pred="SST",resp="SST",lat_resp=lat_resp,lon_resp=lon_resp,add="_noDriver")
