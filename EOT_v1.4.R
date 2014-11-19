library(remote)
library(latticeExtra)
library(gridExtra)
library(ncdf4)
library(raster)

source("C:/Users/Julia Wagemann/Desktop/Teleconnections/anomalies_functions_v1.1.R")

data("australiaGPCP")
data("pacificSST")

path <- "//TEC-JULIA-FRA/Users/julia_francesca/Documents/3_Data/2_Anomalies/ERA_Interim/monthly/4_TP/3deg"
path <- "//TEC-JULIA-FRA/Users/julia_francesca/Documents/3_Data/2_Anomalies/ERA_Interim/monthly/2_SLP/3deg"




fileList_SLP <- list.files(path,pattern=".nc")
fileList_TP <- list.files(path,pattern=".nc")
anomaly <- open.ncdf(fileList_SLP[1], write=FALSE)
var <- raster(get.var.ncdf(anomaly, "SLP_anomaly"))
stack(var)

st_SLP<- stack()
j <- 1
nameVector <- c("Jan_", "Feb_", "Mar_", "Apr_", "May_", "Jun_", "Jul_","Aug_","Sep_","Oct_","Nov_","Dec_" )



st.dns <- denoise(st,expl.var=0.9)
st_TP.dns <- denoise(st_TP,expl.var=0.9)
st_SLP.dns <- denoise(st_SLP,expl.var=0.9)



sst.pred <- deseason(pacificSST, cycle.window = 12)
gpcp.resp <- deseason(australiaGPCP, cycle.window = 12)

sst.pred.dns <- denoise(sst.pred, expl.var = 0.9)
gpcp.resp.dns <- denoise(gpcp.resp,expl.var = 0.9)




modes <- eot(x = st_SLP.dns, y = st_TP.dns, 
             n = 3, standardised = FALSE, 
             reduce.both = FALSE, print.console = TRUE)

plot(modes, 1, 
        show.eot.loc = TRUE, 
        arrange = "lon")




opar <- par(mfrow=c(1,2))
image(gpcp.resp[[1]],main="original")
image(gpcp.resp.dns[[1]],main="deseasoned")
par(opar)