################################################################################
# Name:         visualize_EOT_TS.R
# Description:	This script retrieves time series values of calculated
#				EOT analysis visualizes them in time.
#				The result shows the time series of the identifies base point for
#				the specific EOT mode.
# version:      v1.0
# Author:       Julia Wagemann
# Date:         2014/01/23
################################################################################
path <- "P:/4_ML_results/2_EOT/predictor_response/ts/"
setwd(path)

# File with time series
dataFile <- "EOT_ts_SST_T2m_lag_0.txt"
ts <- read.csv(dataFile,header=FALSE,sep=";")

lagTime <- 0 #specify lag time to adjust vector for x-Axis plotting

dataVec <- seq(as.Date("1990/01/01"),as.Date("2011/12/01"),by="month")
minor <- seq(as.Date("01/01/1990", format = "%d/%m/%Y"), tail(dataVec, 1),
             by = "12 months")


plot(dataVec[1:(length(dataVec)-lagTime)],ts[,1],type="l",lwd=2,xaxt="n",xlab="",
     ylab="",col="blue")
abline(h=0,col="black")
abline(h=-2:2,v=as.numeric(minor),col="grey",lty="dotted")
axis.Date(side=1,dataVec[1:(length(dataVec)-0)],format="%m/%Y",at=minor,las=2)
