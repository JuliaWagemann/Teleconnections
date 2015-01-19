path1 <- "P:/4_ML_results/2_EOT/predictor_response/EV_BP_SST_TP/"
path2 <- "P:/4_ML_results/2_EOT/predictor_response/EV_BP_SST_T2M/"
setwd(path1)
fileList1 <- list.files()
matr1 <- matrix(ncol=1,nrow=length(fileList1))
j <- 1
for(i in fileList1){
        temp <- read.table(i,header=FALSE,sep=";")
        matr1[j,] <- temp[1,1]*100
        j <- j+1
}

setwd(path2)
fileList2 <- list.files()
matr2 <- matrix(ncol=1,nrow=length(fileList2))
j <- 1
for(i in fileList2){       
        temp <- read.table(i,header=FALSE,sep=";")
        matr2[j,] <- temp[1,1]*100
        j <- j+1
}

vec <- seq(0,12)
plot(vec,matr1,type="b",xlab='lag time',ylab='EV in %',col="blue",lwd=2)
lines(vec,matr2,type="b", col="red",lwd=2)
legend(6.6,9,c("SST-Precipitation","SST-Temperature"),lty=c(1,1),lwd=c(2,2),col=c('blue','red'))