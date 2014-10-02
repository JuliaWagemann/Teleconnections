library(raster)
library(ggplot2)
library(sp)
library(ncdf)
library(RColorBrewer)
library(rasterVis)
library(rgdal)
library(fields)

########################################################################################
# Section for individual adjustment
# Please specify all variables in this section individually for every single parameter
########################################################################################

#path where the shapefiles are:
shp_path <- "//TEC-JULIA-FRA/Users/julia_francesca/Documents/Data/3_shapefiles"

#eot netCDF files folder:
eot_path <- "//TEC-JULIA-FRA/Users/julia_francesca/Documents/Data/4_Results/EOT/ERA_Interim/3deg/GP500_monthly_3deg_eot.nc"

#name of the variable in the netCDF file:
vname <- "modes_st"

#Path and name you want to give to the plot tiff
tiff_path <- "//TEC-JULIA-FRA/Users/julia_francesca/Documents/Data/4_Results/EOT/ERA_Interim/3deg/GP500_monthly_3deg_EOT_modes.tiff"

#######################################################################################


# Shapefiles from Natural Earth (http://www.naturalearthdata.com/features/) for World
# boundary layers and a raster grid of 30 degrees
wmap <- readOGR(dsn=shp_path, layer="ne_110m_land")
grid <- readOGR(dsn=shp_path,,layer="ne_110m_graticules_30")

my.colors <- colorRampPalette(c("darkblue","blue","lightblue","beige","tomato","red","darkred"))

# Store netCDF data directly into a raster stack
st <- stack(eot_path,bands=c(1:12),varname=vname)

# ncdf data are stored from longitude 0 to 360 --> therefore the raster objects stored 
# in the rasterstack have to be rearranged to the extent -180 to 180 longitude to match with
# world boundary layer
# Within the loop, each raster layer within the raster stack is divided into two parts and
# rearranged - a new raster stack is built
# ext1 <- extent(c(-91.5,91.5, -1.5,178.5))
# ext2 <- extent(c(-91.5,91.5, 178.5,358.5))
# list <- c(st$X1,st$X2, st$X3)
st_new <- stack()
for(i in 1:12){
#         tempRast_1 <- crop(st[[i]],ext1)
#         tempRast_2 <- crop(st[[i]],ext2)
#         tempRast_1 <- flip(tempRast_1, direction='x')
#         tempRast_2 <- flip(tempRast_2, direction='x')
#         tempRast_1 <- flip(tempRast_1, direction='y')
#         tempRast_2 <- flip(tempRast_2, direction='y')
#         tempRast_11 <- shift(tempRast_1, y=-178.5)
#         tempRast_22 <- shift(tempRast_2, y=-178.5)
#         newRast <- merge(tempRast_22,tempRast_11)
        newRast <- st[[i]]                            
#         newRast <- flip(newRast, direction='x')       
#         newRast <- flip(newRast, direction='y')       
        st_new <- addLayer(st_new,newRast)        

}


# Naming of the raster layers within the raster stack is changed
names(st_new) <- c("SOM1", "SOM2", "SOM3", "SOM4", 
                   "SOM5", "SOM6", "SOM7", "SOM8", "SOM9", "SOM10", "SOM11","SOM12")


###############################################################################
# Plotting
###############################################################################

#define individual colours for plot
colourPalette=c(rev(brewer.pal(10,"RdBu")))
# define the breaks for the raster colour bar
# brks <- c(-4,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,4)
brks <- c(-1,-0.5,-0.2,-0.1,-0.05,0,0.05,0.1,0.2,0.5, 1)

tiff(filename=tiff_path,
     width=1700,height=1650, units="px",
     pointsize=10, res=300)

par(mfrow=c(4,3),mar=c(0.5,0.5,0.5,0.5), oma=c(4,2,0,2),bty="n")
#1
plot(st_new$SOM1,col=colourPalette,ext=extent(c(-180,180,-90,90)),legend=FALSE,
     axes=FALSE,bigplot=c(0.05,0.9,0.2,1),breaks=brks)
plot(wmap,add=TRUE)
axis(side=1,at=c(-179, -90, 0, 90, 180),labels=c(-180,-90,0,90,180),line=-1.1)
axis(side=2,at=c(-90,-45, 0,45,90),labels=c(-90,-45,0,45,90),line=0)
#2
plot(st_new$SOM2,col=colourPalette,ext=extent(c(-180,180,-90,90)),legend=FALSE,
     axes=FALSE,bigplot=c(0.05,0.9,0.2,1), breaks=brks)
plot(wmap,add=TRUE)
axis(side=1,at=c(-179, -90, 0, 90, 180),labels=c(-180,-90,0,90,180),line=-1.1)
axis(side=2,at=c(-90,-45, 0,45,90),labels=c(-90,-45,0,45,90),line=0)
#3
plot(st_new$SOM3,col=colourPalette,ext=extent(c(-180,180,-90,90)),legend=FALSE,
     axes=FALSE,bigplot=c(0.05,0.9,0.2,1),breaks=brks)
plot(wmap,add=TRUE)
axis(side=1,at=c(-179, -90, 0, 90, 180),labels=c(-180,-90,0,90,180),line=-1.1)
axis(side=2,at=c(-90,-45, 0,45,90),labels=c(-90,-45,0,45,90),line=0)
#4
plot(st_new$SOM4,col=colourPalette,ext=extent(c(-180,180,-90,90)),legend=FALSE,
     axes=FALSE,bigplot=c(0.05,0.9,0.2,1),breaks=brks)
plot(wmap,add=TRUE)
axis(side=1,at=c(-179, -90, 0, 90, 180),labels=c(-180,-90,0,90,180),line=-1.1)
axis(side=2,at=c(-90,-45, 0,45,90),labels=c(-90,-45,0,45,90),line=0)
#5
plot(st_new$SOM5,col=colourPalette,ext=extent(c(-180,180,-90,90)),legend=FALSE,
     axes=FALSE,bigplot=c(0.05,0.9,0.2,1),breaks=brks)
plot(wmap,add=TRUE)
axis(side=1,at=c(-179, -90, 0, 90, 180),labels=c(-180,-90,0,90,180),line=-1.1)
axis(side=2,at=c(-90,-45, 0,45,90),labels=c(-90,-45,0,45,90),line=0)
#6
plot(st_new$SOM6,col=colourPalette,ext=extent(c(-180,180,-90,90)),legend=FALSE,
     axes=FALSE,bigplot=c(0.05,0.9,0.2,1),breaks=brks)
plot(wmap,add=TRUE)
axis(side=1,at=c(-179, -90, 0, 90, 180),labels=c(-180,-90,0,90,180),line=-1.1)
axis(side=2,at=c(-90,-45, 0,45,90),labels=c(-90,-45,0,45,90),line=0)
#7
plot(st_new$SOM7,col=colourPalette,ext=extent(c(-180,180,-90,90)),legend=FALSE,
     axes=FALSE,bigplot=c(0.05,0.9,0.2,1),breaks=brks)
plot(wmap,add=TRUE)
axis(side=1,at=c(-179, -90, 0, 90, 180),labels=c(-180,-90,0,90,180),line=-1.1)
axis(side=2,at=c(-90,-45, 0,45,90),labels=c(-90,-45,0,45,90),line=0)
#8
plot(st_new$SOM8,col=colourPalette,ext=extent(c(-180,180,-90,90)),legend=FALSE,
     axes=FALSE,bigplot=c(0.05,0.9,0.2,1),breaks=brks)
plot(wmap,add=TRUE)
axis(side=1,at=c(-179, -90, 0, 90, 180),labels=c(-180,-90,0,90,180),line=-1.1)
axis(side=2,at=c(-90,-45, 0,45,90),labels=c(-90,-45,0,45,90),line=0)
#9
plot(st_new$SOM9,col=colourPalette,ext=extent(c(-180,180,-90,90)),legend=FALSE,
     axes=FALSE,bigplot=c(0.05,0.9,0.2,1),breaks=brks)
plot(wmap,add=TRUE)
axis(side=1,at=c(-179, -90, 0, 90, 180),labels=c(-180,-90,0,90,180),line=-1.1)
axis(side=2,at=c(-90,-45, 0,45,90),labels=c(-90,-45,0,45,90),line=0)
#10
plot(st_new$SOM10,col=colourPalette,ext=extent(c(-180,180,-90,90)),legend=FALSE,
     axes=FALSE,bigplot=c(0.05,0.9,0.2,1),breaks=brks)
plot(wmap,add=TRUE)
axis(side=1,at=c(-179, -90, 0, 90, 180),labels=c(-180,-90,0,90,180),line=-1.1)
axis(side=2,at=c(-90,-45, 0,45,90),labels=c(-90,-45,0,45,90),line=0)
#11
plot(st_new$SOM11,col=colourPalette,ext=extent(c(-180,180,-90,90)), 
     axes=FALSE,legend=FALSE,bigplot=c(0.05,0.9,0.2,1),breaks=brks)
plot(wmap,add=TRUE)
axis(side=1,at=c(-179, -90, 0, 90, 180),labels=c(-180,-90,0,90,180),line=-1.1)
axis(side=2,at=c(-90,-45, 0,45,90),labels=c(-90,-45,0,45,90),line=0)

plot(st_new$SOM11,legend.only=TRUE,zlim=c(-7,-4,-1,0,1,4,7),col=colourPalette, horizontal=TRUE,
     legend.mar=0,legend.width=1,nlevels=12, smallplot=c(0, 1, 0, 0.06))
#12
plot(st_new$SOM12,col=colourPalette,ext=extent(c(-180,180,-90,90)),legend=FALSE,
     axes=FALSE,bigplot=c(0.05,0.9,0.2,1),breaks=brks)
plot(wmap,add=TRUE)
axis(side=1,at=c(-179, -90, 0, 90, 180),labels=c(-180,-90,0,90,180),line=-1.1)
axis(side=2,at=c(-90,-45, 0,45,90),labels=c(-90,-45,0,45,90),line=0)

dev.off()
