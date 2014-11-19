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
# how many modes do you want to retain?
nOfEOT <- 12

#what are the dimensions of the input netCDF file?
geoLat <- 60
geoLon <- 120
nrOfPixels <- geoLat*geoLon

#what are the input and output paths?
inputPath <- "//TEC-JULIA-FRA/Users/julia_francesca/Documents/3_Data/2_Anomalies/ERA_Interim/monthly/2_SLP/3deg"
outputPath <- "//TEC-JULIA-FRA/Users/julia_francesca/Documents/EOT"

#what is the input variable name and units?
varNcdfName <- "Tair2m_anomaly_3deg"
varNcdfUnits <- "degrees"

#name and description of the output netCDF file:
ncdfName <- "T2M_eot.nc"
csvName <- "T2M_eot.csv"
varDesc <- "EOT of Temperature at 2m anomalies"
########################################################################

# move to the input folder and list the ncdf files
setwd(inputPath)
fileList <- list.files(inputPath, pattern=".nc")

nrOfTimeSteps <- length(fileList)*12

#import from the ncdf

dataStack<- stack()

for (files in 1:length(fileList)){
        anomaly <- open.ncdf(fileList[files], write=FALSE)
        var <- get.var.ncdf(anomaly, varNcdfName)
        
        for(i in 1:12){
                temp<-var[,,i]
                temp_r<-raster(temp)
                temp_w<-geoWeight(temp_r)
                dataStack <- addLayer(dataStack, temp_w)
        }
}

#compute the eot
var.eot <- eot(dataStack, resp = NULL, n=nOfEOT, standardized=FALSE, write.out=TRUE, path.out=outputPath, prefix="remote_test")

modes <- var.eot@modes

# writes the nOfEOT value layers, EOT vectors and variance values in a stack, a matrix and a vector

modes_st<-stack()
eot_matrix <- matrix (0, nrOfTimeSteps, nOfEOT)
variance_vector <- vector("numeric", nOfEOT)
names <- vector("character",nOfEOT)
for (k in 1:nOfEOT){
        
        string1 <- paste("vals <- var.eot@modes$mode_", formatC(k, width=2, flag=0), "@r_predictor@data@values", sep="")
        string2 <- paste("eots <- var.eot@modes$mode_", formatC(k, width=2, flag=0), "@eot", sep="")
        string3 <- paste("cum_var <- var.eot@modes$mode_", formatC(k, width=2, flag=0), "@cum_exp_var", sep="")
        eval(parse(text=string1))
        eval(parse(text=string2))
        eval(parse(text=string3))
        m <- matrix(vals, 120, 60)
        m <- m[,ncol(m):1]
        m <- raster(m)
        modes_st <- addLayer(modes_st, m)
        eot_matrix[,k] <- eots
        variance_vector[k] <- cum_var
        names[k]=paste("eot_", toString(k),sep="")
                
}

#writing output
setwd(outputPath)

#preparation of a csv
eotDF <-data.frame(rownames=names, Variance=variance_vector, time_step=t(eot_matrix))
write.csv(eotDF, csvName)

# NetCDF Output preparation
lat <- anomaly$dim$lat$vals
lon <- anomaly$dim$lon$vals



# Define the dimensions of the resulting netCDF file
x <- ncdim_def("lon", "degrees",lon)
y <- ncdim_def("lat", "degrees", lat)
nrEOT <- ncdim_def("n of EOT", " ", c(1:12))

# Define Variable which is written to the netCDF file
netCDF_modes <- ncvar_def("modes_st", units="", dim=list(x,y,nrEOT), missval=NA,
                           longname="EOT modes", prec= 'double', verbose=FALSE)

# Create the actual netCDF output file
netCDFFile <- nc_create(ncdfName, netCDF_modes, force_v4=FALSE,verbose=FALSE)

# fill the created netCDF file with the actual content
ncvar_put(netCDFFile, netCDF_modes, as.matrix(brick(modes_st)),verbose=FALSE)

# close the created netCDF file
nc_close(netCDFFile)



