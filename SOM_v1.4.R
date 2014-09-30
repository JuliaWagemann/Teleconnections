library(kohonen)
library(ncdf)
library(RNetCDF)
library(ncdf4)
library(fields)
library(raster)



path <-"//TEC-JULIA-FRA/Users/julia_francesca/Documents/Data/2_Anomalies/ERA_Interim/Tair_2m/"
setwd(path)

path1 <- "//TEC-JULIA-FRA/Users/julia_francesca/Documents/SOM_exploration/Tair_2m/"

fileList <- list.files(path, pattern=".nc")

#########################################################################
# Section for individual adjustment
# Please specify all variables in this section individually for every single parameter
#########################################################################
nrOfPixels <- 65160
nrOfTimeSteps <- length(fileList)*12

trainingRate <- 500
neighbourhoodRadius <- c(1,3,5,10)
neighbourhoodRadius <- 5

gridSize_1 <- c(3,3,4,4,5)
gridSize_2 <- c(3,4,4,5,5)
gridSize <- list(firstCol=gridSize_1,secondCol=gridSize_2)


geoLat <- 181
geoLon <- 360

varName <- "Tair.som"
varUnit <- "degrees Celsius"
varDesc <- "SOM nodes of air temperature anomalies at 2m"

#ncdfName <- "som_SLPA_rlen_1000_radius_3.nc"
########################################################################


dataMatrix <- matrix(ncol=nrOfPixels, nrow=nrOfTimeSteps)


# for all files in fileList, the monthly entries are taken and written to a dataMatrix
# the resulting dataMatrix TxP contains the time values as rows (T) and the number of pixels
# as columns (P)
t <- 1  #initialize a count variable
for(m in fileList){
        GP_anomaly <- open.ncdf(m, write=FALSE)
        GP_var <- get.var.ncdf(GP_anomaly, "t2m_anomaly")
        for(z in 1:12){
                tempRaster <- GP_var[,,z]
                dataMatrix[t,] <- tempRaster
                t <- t+1 
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
outputFile <- file(paste(path1,"gridSizeAssessment.txt"),open="w")
writeLines(c("SOM_size,radius,q-error,topo-error"),outputFile,sep="\n")
for(i in 1:length(gridSize_1)){
  for (j in neighbourhoodRadius){
    set.seed(20)
    ncol <- gridSize$firstCol[i]
    nrow <- gridSize$secondCol[i]
    GPA.som <- som(data=dataMatrix, grid=somgrid(gridSize$firstCol[i],
                                                 gridSize$secondCol[i],"rectangular"), 
                   rlen=trainingRate, alpha=c(0.05, 0.01), radius=j, toroidal=FALSE, 
                   n.hood="circular", keep.data=TRUE)
    
    # Calculate q- and topo error for calculated som object
    q.error <- mean(GPA.som$distances)
    t.error <- topo.error(GPA.som,"bmu")
    
    qList <- c()
    tList <- c()
    qList[i] <- q.error
    tList[i] <- t.error
    
    fileLine <- c(paste(gridSize$firstCol[i],"x",gridSize$secondCol[i],
                        ",",j,",",qList[i],",",tList[i],sep=""))
    writeLines(fileLine,outputFile,sep="\n")
    
    
    
    # most important element of the som-object: codebook vectors that are stored in rows
    # --> for each somgrid node one codebood vector will be produced with the amount of pixels
    # these codebook vectors can then be plotted as a matrix with the amount of rows / columns
    codes <- GPA.som$codes
    
    #reshape the results in a raster stack
    st<-stack()
    somElements <- gridSize$firstCol[i]*gridSize$secondCol[i]
    for(k in 1:somElements){
      mat<-matrix(codes[k,], geoLat, geoLon)
#      mat<-mat[nrow(mat):1,] 
      r<-raster(mat)
      st<-addLayer(st, r)
    }
    
    # Various plot opportunities offered by the kohonen package 
    # shows the mean distance of the closest codebook vector during training
    # --> this is a good opportunity to evaluate the training process
    tiff(filename=paste(path1,"training_process_somgrid_",gridSize$firstCol[i],"x",
                        gridSize$secondCol[i],"_radius_",j,".tiff",sep=""), 
         width=800,height=700,units="px",pointsize=12,res=100)
    plot(GPA.som,type="changes")
    dev.off()
    
    ncdfName <- paste(path1,"som_Tair_somgrid_,",gridSize$firstCol[i],"x",
                      gridSize$secondCol[i],"_radius_",j,"_2.nc",sep="")
    
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
  }
}
close(outputFile)






# shows the different codebook vectors
plot(GPA.som,type="codes")
# shows the number of objects mapped to the individual units
plot(GPA.som,type="counts")
# shows the sum of the distances to all immediate neighbours --> also known as a U-matrix plot
plot(GPA.som,type="dist.neighbours")
# shows where objects are mapped. Uses the classif argument and needs labels or pchs
plot(GPA.som$unit.classif)

freq <- c(GPA.som$unit.classif)

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





