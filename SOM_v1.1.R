library(kohonen)
library(ncdf)
library(RNetCDF)
library(ncdf4)
library(fields)
library(raster)

path <- "F:/Julia/Data/ERA_interim/Anomalies/GP500/"
setwd(path)

fileList <- list.files(path, pattern=".nc")

#########################################################################
# Section for individual adjustment
# Please specify all variables in this section individually for every single parameter
#########################################################################
nrOfPixels <- 65160
nrOfTimeSteps <- length(fileList)*12
somRows <- 4
somCols <- 3

trainingRate <- 500
neighbourhoodRadius <- 5

geoLat <- 181
geoLon <- 360

varName <- "GPA.som"
varUnit <- "hPa"
varDesc <- "SOM nodes of Geopotential 500 anomalies"

ncdfName <- "som_GPA500.nc"
########################################################################

set.seed(300)
dataMatrix <- matrix(ncol=nrOfPixels, nrow=nrOfTimeSteps)


# for all files in fileList, the monthly entries are taken and written to a dataMatrix
# the resulting dataMatrix TxP contains the time values as rows (T) and the number of pixels
# as columns (P)
j <- 1  #initialize a count variable
for(i in 1:length(fileList)){
        GP_anomaly <- open.ncdf(fileList[i], write=FALSE)
        GP_var <- get.var.ncdf(GP_anomaly, "GP500_anomaly")
        for(k in 1:12){
                tempRaster <- GP_var[,,k]
                dataMatrix[j,] <- tempRaster
                j <- j+1  
        }
}

# dataMatrix is scaled based on the number of columns of the matrix
for (k in 1:nrOfPixels){
        std=sd(dataMatrix[,k], na.rm=TRUE)
        dataMatrix[,k]<-dataMatrix[,k]/std
}

# Training of the SOM for the entire amount of input data
# Parameters to be carefully chosen: 
# - size of the somgrid:        objective identification of the ideal size of the somgrid with q- and t-error
#                               (refer to Rousi et al. 2014)
# - rlen:       number of times the complete dataset will be presented to the network
# - alpha:      learning rate, a vector of two numbers indicating the amount of change
#               small value leads to a slow and smooth learning process, high value leads to a fast and unstable learning
# - radius:     radius of the neighbourhood
# - toroidal:   if TRUE --> edges of the neighbourhood map are joined
# - n.hood:     shape of the neighbourhood - either "circular" or "square"
GPA.som <- som(data=datMat.sc, grid=somgrid(somCols,somRows,"rectangular"), rlen=trainingRate, alpha=c(0.05, 0.01),
               radius=neighbourhoodRadius,toroidal=FALSE, n.hood="circular", keep.data=TRUE)

# most important element of the som-object: codebook vectors that are stored in rows
# --> for each somgrid node one codebood vector will be produced with the amount of pixels
# these codebook vectors can then be plotted as a matrix with the amount of rows / columns
codes <- GPA.som$codes

#reshape the results in a raster stack
st<-stack()
somElements <- somRows*somCols
for(j in 1:somElements){
        mat<-matrix(codes[j,], geoLat, geoLon)
        mat<-mat[nrow(mat):1,] 
        r<-raster(mat)
        st<-addLayer(st, r)
}

# Various plot opportunities offered by the kohonen package 
# shows the mean distance of the closest codebook vector during training
# --> this is a good opportunity to evaluate the training process
plot(GPA.som,type="changes")
# shows the different codebook vectors
plot(GPA.som,type="codes")
# shows the number of objects mapped to the individual units
plot(GPA.som,type="counts")
# shows the sum of the distances to all immediate neighbours --> also known as a U-matrix plot
plot(GPA.som,type="dist.neighbours")
# shows where objects are mapped. Uses the classif argument and needs labels or pchs
plot(GPA.som$unit.classif,type="mapping")


# Each codebook vector is extracted from the som-object, rearranged in a NxM matrix and then plotted as
# raster field
par(mfrow=c(4,3))
som.node1 <- codes[1,]
image(matrix(som.node1,ncol=360,nrow=181))
som.node2 <- codes[2,]
image(matrix(som.node2,ncol=360,nrow=181))
som.node3 <- codes[3,]
image(matrix(som.node3,ncol=360,nrow=181))
som.node4 <- codes[4,]
image(matrix(som.node4,ncol=360,nrow=181))
som.node5 <- codes[5,]
image(matrix(som.node5,ncol=360,nrow=181))
som.node6 <- codes[6,]
image(matrix(som.node6,ncol=360,nrow=181))
som.node7 <- codes[7,]
image(matrix(som.node7,ncol=360,nrow=181))
som.node8 <- codes[8,]
image(matrix(som.node8,ncol=360,nrow=181))
som.node9 <- codes[9,]
image(matrix(som.node9,ncol=360,nrow=181))
som.node10 <- codes[10,]
image(matrix(som.node10,ncol=360,nrow=181))
som.node11 <- codes[11,]
image(matrix(som.node11,ncol=360,nrow=181))
som.node12 <- codes[12,]
image(matrix(som.node12,ncol=360,nrow=181))


# NetCDF Output preparation
lat <- 1:geoLat
lat<- lat-90
lon <- 1:geoLon

# Define the dimensions of the resulting netCDF file
x <- ncdim_def("lon", "degrees",lon)
y <- ncdim_def("lat", "degrees", lat)
soms <- ncdim_def("soms", " ", c(1:nrow(codes)))

# Define Variable which is written to the netCDF file
netCDF_varDef <- ncvar_def(varName, units=varUnit,dim=list(x,y,soms), missval=NA,
                           longname=varDesc, prec= 'double', verbose=FALSE)

# Create the actual netCDF output file
netCDFFile <- nc_create(ncdfName, netCDF_varDef, force_v4=FALSE,verbose=FALSE)

# fill the created netCDF file with the actual content --> calculated anomalies
ncvar_put(netCDFFile, netCDF_varDef, as.matrix(brick(st)),verbose=FALSE)

# close the created netCDF file
nc_close(netCDFFile)



