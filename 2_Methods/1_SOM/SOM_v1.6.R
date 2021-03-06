################################################################################
# Name:         SOM_v1.6.R
# Description:	This script gives some examples how the R package "kohonen" needs
#				needs to be set up in order to calculate "Self-organizing-maps"
#				for a given time series of raster data sets 
# version:      v1.6
# Author:       Julia Wagemann
# Date:         2014/10/03
################################################################################
library(kohonen)
library(ncdf)
library(RNetCDF)
library(ncdf4)
library(fields)
library(raster)

path <- "C:/Users/julia_francesca/Documents/Data/2_Anomalies/CCI/SST_anomaly/1deg_aggregated/"
setwd(path)

path1 <- "C:/Users/julia_francesca/Documents/Data/4_Results/SOM/CCI/SST/"

fileList <- list.files(path, pattern=".nc")

#########################################################################
# Section for individual adjustment
# Please specify all variables in this section individually for every single parameter
#########################################################################
nrOfTimeSteps <- length(fileList)*12

trainingRate <- 500
neighbourhoodRadius <- 5

gridSize_1 <- 3
gridSize_2 <- 4

geoLat <- 180
geoLon <- 360
nrOfPixels <- geoLat*geoLon

varName <- "SST.som"
varUnit <- "degC"
varDesc <- "SOM nodes of sea surface temperature anomalies"

#ncdfName <- "som_SLPA_rlen_1000_radius_3.nc"
########################################################################



matrixList <- list()
dataMatrix <- matrix(ncol=nrOfPixels, nrow=nrOfTimeSteps)
# for all files in fileList, the monthly entries are taken and written to a dataMatrix
# the resulting dataMatrix TxP contains the time values as rows (T) and the number of pixels
# as columns (P)
t <- 1  #initialize a count variable
for(m in fileList){

        GP_anomaly <- open.ncdf(m, write=FALSE)
        GP_var <- get.var.ncdf(GP_anomaly, "SST_anomaly_1deg")
        for(z in 1:12){
                tempRaster <- GP_var[,,z]
                tempRaster[tempRaster=="NaN"] <- NA
                dataMatrix[t,] <- tempRaster
                t <- t+1 
        }
}

matrixList[[1]] <- dataMatrix

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
outputFile <- file(paste(path1,"SST_SOM_evaluation.txt"),open="w")
outputFile_2 <- file(paste(path1,"SST_freq.txt"),open="w")
writeLines(c("SOM_size,radius,q-error,topo-error"),outputFile,sep="\n")

set.seed(20)
ncol <- gridSize_1
nrow <- gridSize_2
GPA.som <- supersom(data=matrixList, grid=somgrid(gridSize_1,
                                             gridSize_2,"rectangular"), 
               rlen=trainingRate, alpha=c(0.05, 0.01), radius=neighbourhoodRadius, toroidal=FALSE, 
               n.hood="circular", keep.data=TRUE)

# Calculate q- and topo error for calculated som object
q.error <- mean(GPA.som$distances)
t.error <- topo.error(GPA.som,"bmu")

fileLine <- c(paste(gridSize_1,"x",gridSize_2,
                    ",",neighbourhoodRadius,",",q.error,",",t.error,sep=""))
writeLines(fileLine,outputFile,sep="\n")




# most important element of the som-object: codebook vectors that are stored in rows
# --> for each somgrid node one codebood vector will be produced with the amount of pixels
# these codebook vectors can then be plotted as a matrix with the amount of rows / columns
codes <- GPA.som$codes[[1]]
newCodes <- matrix(ncol=nrOfPixels,nrow=12)
for(i in 1:12){
  temp <- replace(dataMatrix[22,],!is.na(dataMatrix[22,]),codes[i,])
  newCodes[i,] <- temp
}

freq <- table(GPA.som$unit.classif)/length(GPA.som$unit.classif)
writeLines(as.character(freq),outputFile_2,sep="\n")

#reshape the results in a raster stack
st<-stack()
somElements <- gridSize_1*gridSize_2
for(k in 1:somElements){
  mat<-matrix(newCodes[k,], geoLat, geoLon)
  #      mat<-mat[nrow(mat):1,] 
  r<-raster(mat)
  st<-addLayer(st, r)
}

# Various plot opportunities offered by the kohonen package 
# shows the mean distance of the closest codebook vector during training
# --> this is a good opportunity to evaluate the training process
tiff(filename=paste(path1,"training_process_somgrid_",gridSize_1,"x",
                    gridSize_2,"_radius_",neighbourhoodRadius,".tiff",sep=""), 
     width=800,height=700,units="px",pointsize=12,res=100)
plot(GPA.som,type="changes")
dev.off()

ncdfName <- paste(path1,"som_SST_somgrid_",gridSize_1,"x",
                  gridSize_2,"_radius_",neighbourhoodRadius,".nc",sep="")

# NetCDF Output preparation
lat <- GP_anomaly$dim$lat$vals
lon <- GP_anomaly$dim$lon$vals

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

close(outputFile)
close(outputFile_2)





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





