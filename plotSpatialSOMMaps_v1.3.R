library(raster)
library(ggplot2)
library(sp)
library(ncdf)
library(RColorBrewer)
library(rasterVis)
library(rgdal)
library(fields)

path <- "F:/Julia/Data/"
setwd(path)

path_TECPC <- "//TEC-JULIA-FRA/Users/julia_francesca/Documents/SOM_exploration/Tair_2m/"
fileList <- list.files(path_TECPC,pattern=".nc")
# Shapefiles from Natural Earth (http://www.naturalearthdata.com/features/) for World
# boundary layers and a raster grid of 30 degrees
wmap <- readOGR(dsn=paste(path,"shapefiles",sep=""), layer="ne_110m_land")
grid <- readOGR(dsn=paste(path,"shapefiles",sep=""),layer="ne_110m_graticules_30")

my.colors <- colorRampPalette(c("darkblue","blue","lightblue","beige","tomato","red","darkred"))


for(i in fileList){
        fileName <- substr(i,10,nchar(i)-3)
        # Store netCDF data directly into a raster stack
        st <- stack(paste(path_TECPC,
                          i,sep=""),
                    bands=c(1:20),varname="Tair.som")
        
        # ncdf data are stored from longitude 0 to 360 --> therefore the raster objects stored 
        # in the rasterstack have to be rearranged to the extent -180 to 180 longitude to match with
        # world boundary layer
        # Within the loop, each raster layer within the raster stack is divided into two parts and
        # rearranged - a new raster stack is built
        ext1 <- extent(c(-0.5,179.5,-90.5,90.5))
        ext2 <- extent(c(179.5,359.5,-90.5,90.5))
        list <- c(st$X1,st$X2, st$X3)
        st_new <- stack()
        for(i in 1:20){
                tempRast_1 <- crop(st[[i]],ext1)
                tempRast_2 <- crop(st[[i]],ext2)
                tempRast_11 <- shift(tempRast_1, x=0.5)
                tempRast_22 <- shift(tempRast_2, x=-359.5)
                newRast <- merge(tempRast_22,tempRast_11)
                st_new <- addLayer(st_new,newRast)        
        }
        
        # Naming of the raster layers within the raster stack is changed
        names(st_new) <- c("SOM1", "SOM2", "SOM3", "SOM4", 
                           "SOM5", "SOM6", "SOM7", "SOM8", "SOM9" , "SOM10", "SOM11","SOM12",
                           "SOM13","SOM14","SOM15","SOM16","SOM17", "SOM18", "SOM19", "SOM20")
                           
        ,
                           "SOM21", "SOM22", "SOM23", "SOM24", "SOM25")
        
        ###############################################################################
        # Plotting
        ###############################################################################
        
        #define individual colours for plot
        colourPalette=c(rev(brewer.pal(10,"RdBu")))
        colourPalette=my.colors(100)
        # define the breaks for the raster colour bar
        brks <- c(-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5)
        brks <- seq(-3,3,length.out=100)
        
        tiff(filename=paste(path_TECPC,"SOM_nodes_",fileName,".tiff",sep=""),
             width=1700,height=1650, units="px",
             pointsize=10, res=300)
        
        par(mfrow=c(5,4),mar=c(0,2,0,0.5), oma=c(3,2,0,1),bty="n")
        #1
        plot(st_new$SOM1,col=colourPalette,ext=extent(c(-180,180,-90,90)),legend=FALSE,
             axes=FALSE,bigplot=c(0.05,0.9,0,1),breaks=brks,asp=1.2)
        plot(wmap,add=TRUE)
        axis(side=1,at=c(-179, -90, 0, 90, 180),labels=c(-180,-90,0,90,180),line=-1.3)
        axis(side=2,at=c(-90,-45, 0,45,90),labels=c(-90,-45,0,45,90),line=0)
        #2
        plot(st_new$SOM2,col=colourPalette,ext=extent(c(-180,180,-90,90)),legend=FALSE,
             axes=FALSE,bigplot=c(0.05,0.9,0,1), breaks=brks,asp=1.2)
        plot(wmap,add=TRUE)
        axis(side=1,at=c(-179, -90, 0, 90, 180),labels=c(-180,-90,0,90,180),line=-1.3)
#        axis(side=2,at=c(-90,-45, 0,45,90),labels=c(-90,-45,0,45,90),line=0)
        #3
        plot(st_new$SOM3,col=colourPalette,ext=extent(c(-180,180,-90,90)),legend=FALSE,
             axes=FALSE,bigplot=c(0.05,0.9,0,1),breaks=brks,asp=1.2)
        plot(wmap,add=TRUE)
        axis(side=1,at=c(-179, -90, 0, 90, 180),labels=c(-180,-90,0,90,180),line=-1.3)
#        axis(side=2,at=c(-90,-45, 0,45,90),labels=c(-90,-45,0,45,90),line=0)
        #4
        plot(st_new$SOM4,col=colourPalette,ext=extent(c(-180,180,-90,90)),legend=FALSE,
             axes=FALSE,bigplot=c(0.05,0.9,0,1),breaks=brks,asp=1.2)
        plot(wmap,add=TRUE)
        axis(side=1,at=c(-179, -90, 0, 90, 180),labels=c(-180,-90,0,90,180),line=-1.3)
#        axis(side=2,at=c(-90,-45, 0,45,90),labels=c(-90,-45,0,45,90),line=0)
        #5
#         plot(st_new$SOM5,col=colourPalette,ext=extent(c(-180,180,-90,90)),legend=FALSE,
#              axes=FALSE,bigplot=c(0.05,0.9,0,1),breaks=brks,asp=1.5)
#         plot(wmap,add=TRUE)
#         axis(side=1,at=c(-179, -90, 0, 90, 180),labels=c(-180,-90,0,90,180),line=-1.3)
# #        axis(side=2,at=c(-90,-45, 0,45,90),labels=c(-90,-45,0,45,90),line=0)
        #6
        plot(st_new$SOM5,col=colourPalette,ext=extent(c(-180,180,-90,90)),legend=FALSE,
             axes=FALSE,bigplot=c(0.05,0.9,0,1),breaks=brks,asp=1.2)
        plot(wmap,add=TRUE)
        axis(side=1,at=c(-179, -90, 0, 90, 180),labels=c(-180,-90,0,90,180),line=-1.3)
        axis(side=2,at=c(-90,-45, 0,45,90),labels=c(-90,-45,0,45,90),line=0)
        #7
        plot(st_new$SOM6,col=colourPalette,ext=extent(c(-180,180,-90,90)),legend=FALSE,
             axes=FALSE,bigplot=c(0.05,0.9,0,1),breaks=brks,asp=1.2)
        plot(wmap,add=TRUE)
        axis(side=1,at=c(-179, -90, 0, 90, 180),labels=c(-180,-90,0,90,180),line=-1.3)
 #       axis(side=2,at=c(-90,-45, 0,45,90),labels=c(-90,-45,0,45,90),line=0)
        #8
        plot(st_new$SOM7,col=colourPalette,ext=extent(c(-180,180,-90,90)),legend=FALSE,
             axes=FALSE,bigplot=c(0.05,0.9,0,1),breaks=brks,asp=1.2)
        plot(wmap,add=TRUE)
        axis(side=1,at=c(-179, -90, 0, 90, 180),labels=c(-180,-90,0,90,180),line=-1.3)
 #       axis(side=2,at=c(-90,-45, 0,45,90),labels=c(-90,-45,0,45,90),line=0)
        #9
        plot(st_new$SOM8,col=colourPalette,ext=extent(c(-180,180,-90,90)),legend=FALSE,
             axes=FALSE,bigplot=c(0.05,0.9,0,1),breaks=brks,asp=1.2)
        plot(wmap,add=TRUE)
        axis(side=1,at=c(-179, -90, 0, 90, 180),labels=c(-180,-90,0,90,180),line=-1.3)
 #       axis(side=2,at=c(-90,-45, 0,45,90),labels=c(-90,-45,0,45,90),line=0)
        #10
#         plot(st_new$SOM10,col=colourPalette,ext=extent(c(-180,180,-90,90)),legend=FALSE,
#              axes=FALSE,bigplot=c(0.05,0.9,0,1),breaks=brks,asp=1.5)
#         plot(wmap,add=TRUE)
#         axis(side=1,at=c(-179, -90, 0, 90, 180),labels=c(-180,-90,0,90,180),line=-1.3)
# #        axis(side=2,at=c(-90,-45, 0,45,90),labels=c(-90,-45,0,45,90),line=0)
        #11
        plot(st_new$SOM9,col=colourPalette,ext=extent(c(-180,180,-90,90)),legend=FALSE,
             axes=FALSE,bigplot=c(0.05,0.9,0,1),breaks=brks,asp=1.1)
        plot(wmap,add=TRUE)
        axis(side=1,at=c(-179, -90, 0, 90, 180),labels=c(-180,-90,0,90,180),line=-1.5)
        axis(side=2,at=c(-90,-45, 0,45,90),labels=c(-90,-45,0,45,90),line=0)
        #12
        plot(st_new$SOM10,col=colourPalette,ext=extent(c(-180,180,-90,90)),legend=FALSE,
             axes=FALSE,bigplot=c(0.05,0.9,0,1),breaks=brks,asp=1.1)
        plot(wmap,add=TRUE)
        axis(side=1,at=c(-179, -90, 0, 90, 180),labels=c(-180,-90,0,90,180),line=-1.5)
#       axis(side=2,at=c(-90,-45, 0,45,90),labels=c(-90,-45,0,45,90),line=0)
        #13
        plot(st_new$SOM11,col=colourPalette,ext=extent(c(-180,180,-90,90)),legend=FALSE,
             axes=FALSE,bigplot=c(0.05,0.9,0,1),breaks=brks,asp=1.1)
        plot(wmap,add=TRUE)
        axis(side=1,at=c(-179, -90, 0, 90, 180),labels=c(-180,-90,0,90,180),line=-1.5)
#       axis(side=2,at=c(-90,-45, 0,45,90),labels=c(-90,-45,0,45,90),line=0)
        #14
        plot(st_new$SOM12,col=colourPalette,ext=extent(c(-180,180,-90,90)),legend=FALSE,
             axes=FALSE,bigplot=c(0.05,0.9,0,1),breaks=brks,asp=1.1)
        plot(wmap,add=TRUE)
        axis(side=1,at=c(-179, -90, 0, 90, 180),labels=c(-180,-90,0,90,180),line=-1.5)
#       axis(side=2,at=c(-90,-45, 0,45,90),labels=c(-90,-45,0,45,90),line=0)
        #15
#         plot(st_new$SOM15,col=colourPalette,ext=extent(c(-180,180,-90,90)),legend=FALSE,
#              axes=FALSE,bigplot=c(0.05,0.9,0,1),breaks=brks,asp=1.5)
#         plot(wmap,add=TRUE)
#         axis(side=1,at=c(-179, -90, 0, 90, 180),labels=c(-180,-90,0,90,180),line=-1.3)
# #        axis(side=2,at=c(-90,-45, 0,45,90),labels=c(-90,-45,0,45,90),line=0)
        #16
        plot(st_new$SOM13,col=colourPalette,ext=extent(c(-180,180,-90,90)),legend=FALSE,
             axes=FALSE,bigplot=c(0.05,0.9,0.1,1),breaks=brks,asp=1.1)
        plot(wmap,add=TRUE)
        axis(side=1,at=c(-179, -90, 0, 90, 180),labels=c(-180,-90,0,90,180),line=-1.1)
        axis(side=2,at=c(-90,-45, 0,45,90),labels=c(-90,-45,0,45,90),line=0)
        #17
        plot(st_new$SOM14,col=colourPalette,ext=extent(c(-180,180,-90,90)),legend=FALSE,
             axes=FALSE,bigplot=c(0.05,0.9,0.1,1),breaks=brks,asp=1.1)
        plot(wmap,add=TRUE)
        axis(side=1,at=c(-179, -90, 0, 90, 180),labels=c(-180,-90,0,90,180),line=-1.1)
#        axis(side=2,at=c(-90,-45, 0,45,90),labels=c(-90,-45,0,45,90),line=0)
        #18
        plot(st_new$SOM15,col=colourPalette,ext=extent(c(-180,180,-90,90)),legend=FALSE,
             axes=FALSE,bigplot=c(0.05,0.9,0.1,1),breaks=brks,asp=1.1)
        plot(wmap,add=TRUE)
        axis(side=1,at=c(-179, -90, 0, 90, 180),labels=c(-180,-90,0,90,180),line=-1.1)
#        axis(side=2,at=c(-90,-45, 0,45,90),labels=c(-90,-45,0,45,90),line=0)
        #19
        plot(st_new$SOM16,col=colourPalette,ext=extent(c(-180,180,-90,90)), 
             axes=FALSE,legend=FALSE,bigplot=c(0.05,0.9,0.1,1),breaks=brks,asp=1.1)
        plot(wmap,add=TRUE)
        axis(side=1,at=c(-179, -90, 0, 90, 180),labels=c(-180,-90,0,90,180),line=-1.1)
#        axis(side=2,at=c(-90,-45, 0,45,90),labels=c(-90,-45,0,45,90),line=0)
        #20
#         plot(st_new$SOM20,col=colourPalette,ext=extent(c(-180,180,-90,90)),legend=FALSE,
#              axes=FALSE,bigplot=c(0.05,0.9,0.1,1),breaks=brks,asp=1.4)
#         plot(wmap,add=TRUE)
#         axis(side=1,at=c(-179, -90, 0, 90, 180),labels=c(-180,-90,0,90,180),line=-1)
# #        axis(side=2,at=c(-90,-45, 0,45,90),labels=c(-90,-45,0,45,90),line=0)
        #21
        plot(st_new$SOM17,col=colourPalette,ext=extent(c(-180,180,-90,90)),legend=FALSE,
             axes=FALSE,bigplot=c(0.05,0.9,0.25,1),breaks=brks,asp=1.1)
        plot(wmap,add=TRUE)
        axis(side=1,at=c(-179, -90, 0, 90, 180),labels=c(-180,-90,0,90,180),line=-0.2)
        axis(side=2,at=c(-90,-45, 0,45,90),labels=c(-90,-45,0,45,90),line=0) 
        #22
        plot(st_new$SOM18,col=colourPalette,ext=extent(c(-180,180,-90,90)),legend=FALSE,
             axes=FALSE,bigplot=c(0.05,0.9,0.25,1),breaks=brks,asp=1.1)
        plot(wmap,add=TRUE)
        axis(side=1,at=c(-179, -90, 0, 90, 180),labels=c(-180,-90,0,90,180),line=-0.2)
#        axis(side=2,at=c(-90,-45, 0,45,90),labels=c(-90,-45,0,45,90),line=0)
        #23
        plot(st_new$SOM19,col=colourPalette,ext=extent(c(-180,180,-90,90)),legend=FALSE,
             axes=FALSE,bigplot=c(0.05,0.9,0.25,1), breaks=brks,asp=1.1)
        plot(wmap,add=TRUE)
        axis(side=1,at=c(-179, -90, 0, 90, 180),labels=c(-180,-90,0,90,180),line=-0.2)
#        axis(side=2,at=c(-90,-45, 0,45,90),labels=c(-90,-45,0,45,90),line=0)
        plot(st_new$SOM19,legend.only=TRUE,zlim=c(-3, -2, -1.5,-1,-0.5,0,0.5,1,1.5,2,3),col=colourPalette, horizontal=TRUE,
             legend.mar=0,legend.width=3,nlevels=10,smallplot=c(0.05,0.9,0,0.1))
        #24
        plot(st_new$SOM20,col=colourPalette,ext=extent(c(-180,180,-90,90)),legend=FALSE,
             axes=FALSE,bigplot=c(0.05,0.9,0.25,1),breaks=brks,asp=1.1)
        plot(wmap,add=TRUE)
        axis(side=1,at=c(-179, -90, 0, 90, 180),labels=c(-180,-90,0,90,180),line=-0.2)
#        axis(side=2,at=c(-90,-45, 0,45,90),labels=c(-90,-45,0,45,90),line=0)
        #25
#         plot(st_new$SOM25,col=colourPalette,ext=extent(c(-180,180,-90,90)),legend=FALSE,
#              axes=FALSE,bigplot=c(0.05,0.9,0.25,1),breaks=brks,asp=1.3)
#         plot(wmap,add=TRUE)
#         axis(side=1,at=c(-179, -90, 0, 90, 180),labels=c(-180,-90,0,90,180),line=-0.5)
#         #        axis(side=2,at=c(-90,-45, 0,45,90),labels=c(-90,-45,0,45,90),line=0)

        dev.off()
}

