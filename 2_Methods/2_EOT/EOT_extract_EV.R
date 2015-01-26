################################################################################
# Name:         EOT_extract_EV.R
# Description:	This script retrieves explained variance (EV) values of calculated
#				EOT analysis with various lag times and visualizes them in time.
#				The result shows the change in EV with various lagged times of
#				two different parameter combinations.
# version:      v1.0
# Author:       Julia Wagemann
# Date:         2014/01/23
################################################################################

# Path to EOT outputs
path1 <- "P:/4_ML_results/2_EOT/predictor_response/EV_BP_SST_TP/"
path2 <- "P:/4_ML_results/2_EOT/predictor_response/EV_BP_SST_T2M/"

setwd(path1)
fileList1 <- list.files()

matr1 <- matrix(ncol=1,nrow=length(fileList1)) # Initialize matrix 1
j <- 1
# store all EV values to the matrix 1
for(i in fileList1){
        temp <- read.table(i,header=FALSE,sep=";")
        matr1[j,] <- temp[1,1]*100
        j <- j+1
}

setwd(path2)
fileList2 <- list.files()

matr2 <- matrix(ncol=1,nrow=length(fileList2)) # Initialize matrix 2
j <- 1
# store all EV values to the matrix 2
for(i in fileList2){       
        temp <- read.table(i,header=FALSE,sep=";")
        matr2[j,] <- temp[1,1]*100
        j <- j+1
}

# help vector for x-axis
vec <- seq(0,12)

# Visualize both matrices in one plot
plot(vec,matr1,type="b",xlab='lag time',ylab='EV in %',col="blue",lwd=2)
lines(vec,matr2,type="b", col="red",lwd=2)
legend(6.6,9,c("SST-Precipitation","SST-Temperature"),lty=c(1,1),lwd=c(2,2),col=c('blue','red'))