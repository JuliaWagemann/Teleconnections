%%R

library(remote)
library(latticeExtra)
library(gridExtra)
library(ncdf)
library(MASS)


setwd("//TEC-JULIA-FRA/julia_francesca/Desktop/Teleconnections")
source("anomalies_functions_v1.1.R")

path_resp <- "//TEC-JULIA-FRA/Users/julia_francesca/Documents/3_Data/2_Anomalies/ERA_Interim/monthly/4_TP/3deg"
path_pred <- "//TEC-JULIA-FRA/Users/julia_francesca/Documents/3_Data/2_Anomalies/ERA_Interim/monthly/2_SLP/3deg"

path_EOT <- "//TEC-JULIA-FRA/Users/julia_francesca/Documents/3_Data/4_Results/2_EOT/ERA_Interim/predictor_response/"
nameVector <- c("Jan_", "Feb_", "Mar_", "Apr_", "May_", "Jun_", "Jul_","Aug_","Sep_","Oct_","Nov_","Dec_" )


setwd(path_resp)
fileList <- list.files(pattern=".nc")
lat_resp <- getLat(fileList[1])
lon_resp <- getLon(fileList[1])

setwd(path_pred)
fileList <- list.files(pattern=".nc")
lat_pred <- getLat(fileList[1])
lon_pred <- getLon(fileList[1])

st_resp <- createBigStack(path_resp,nameVector,1000)
st_pred <- createBigStack(path_pred,nameVector,0.01)

lag_vec <- c(1:12)
for(i in lag_vec){
  lagged <- lagalize(st_pred,st_resp,lag=6,freq=12)
  
  st_resp.dns <- denoise(lagged[[2]],expl.var=0.9)
  st_pred.dns <- denoise(lagged[[1]],expl.var=0.9)

  modes <- eot(x = st_pred.dns, y = st_resp.dns, 
               n = 3, standardised = FALSE, 
               reduce.both = FALSE, print.console = TRUE)
  
  add <- paste("lag_",i,sep="")
  writeEOTOutput(path_EOT,modes,"SLP","TP",lat_pred,lon_pred,lat_resp,lon_resp,add)
}



#plot(modes, 1, 
#        show.eot.loc = TRUE, 
#        arrange = "lon")
