# Function "getMonthStack"creates a stack of input data for one specific month,
# based on a fileList
# This function is called in the function calculateAnomaly
getMonthStack <- function(fileList, month,parName,NA_val){
  st <- stack()
  for(i in fileList){
    temp <- open.ncdf(i, write=FALSE)
    tempRaster <- get.var.ncdf(temp,parName)
    z <- raster(tempRaster[,,month])
    z[z==NA_val] <- NA
    st <- addLayer(st,z)
    close.ncdf(temp)
  }
  return(st)
}

# Function "calculateAnomaly" calculates monthly anomalies based on the 
# climatology given by an input fileList and returns a stack based on the
# individual months
calculateAnomaly <- function(month,fileList){
  monthGP <- getMonthStack(fileList,month,parName,NA_val)
  climatology <- mean(monthGP)
  st <- stack()
  j <- 1
  for (j in 1:nlayers(monthGP)){
    actualGP <- raster(monthGP,layer=j)-climatology
    st <- addLayer(st,actualGP)
  }
  return(st)
}

# Function "createStackVector" takes the fileList as input, calculates the
# anomalies and creates a stack of all the monthly stacks
createStackVector <- function(fileList){
  stackVector <- c()
  for (i in 1:12){
    temp <- calculateAnomaly(i, fileList)
    stackVector[i] <- c(temp)
  }
  return(stackVector)
}

# Function "getLat" is a help function to retrieve the values for the latitude
# dimension of the input netcdf-files
getLat <- function(file){
  temp <- open.ncdf(file,write=FALSE)
  lat <- temp$dim$latitude$vals
  close.ncdf(temp)
  return (lat)
}

# Function "getLat" is a help function to retrieve the values for the longitude
# dimension of the input netcdf-files
getLon <- function(file){
  temp <- open.ncdf(file,write=FALSE)
  lon <- temp$dim$longitude$vals
  close.ncdf(temp)
  return (lon)
}


