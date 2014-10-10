library(raster)
library(ggplot2)
library(sp)
library(ncdf)
library(RColorBrewer)
library(rasterVis)
library(rgdal)
library(fields)
library(base)

path <- "//TEC-JULIA-FRA/Users/julia_francesca/Documents/Data/4_Results/EOF/"
setwd(path)

csv_file <- read.csv("TP_eof.csv", header = TRUE, sep = ",", quote = "\"",
         dec = ".", fill = TRUE)

###############################################################################
# Plotting
###############################################################################

tiff(filename=paste(path,"TP_EOF_ts.tiff",sep=""),
     width=3400,height=1650, units="px",
     pointsize=10, res=300)

par(mfrow=c(6,2),mar=c(1,1,1,1), oma=c(4,2,0,2),bty="n")
#1

t <- seq(from <- as.Date("1979/1/1"), by = "month", to = as.Date("2013/12/31"))
label <- seq(from <-1979, to = 2013)

for (i in 1:12){
        y <- as.matrix(csv_file[i+2])
        variance <- y[1]
        y=y[2:length(y)] 
        title=paste("Mode", toString(i), "- Variance:", toString(variance))
        plot(x=t, y=y, type="l")
        axis(side=1,at=label,  labels=label, line=-1.1)
        title(main=title, line=-0.7)
}

dev.off()
