rm(list=ls(all=TRUE))  # clean up left overs from previous run, though shouldn't need to.
# This code is meant to be run on a ROCKS 5.1 cluster
if (.Platform$OS.type=="windows") home="/users/radivot" else home="/home/radivot"
source(paste(home,"/start.r",sep=""))  # define host="machineName" in this file
setwd(home)  # place where models and results will go
if (1) library(ccems) else { # if 0 source in the package 
  pkgNms=dir(paste(home,"/case/active/ccems/ccems/R",sep=""),full.names=TRUE)
  for (i in pkgNms) source(i)  # source in all of the R files of the package
  RNR=read.table(file=paste(home,"/case/active/rnr/datasets/RNR.txt",sep=""),header=TRUE)
}

if (!is.null(dev.list())) for (i in 2:max(dev.list())) dev.off(i);
windows(width = 8, height = 4,restoreConsole = TRUE) 
par(mfrow=c(1,2),mar=c(4,4,2,1)+.1)

d=subset(RNR,(year==2002)&(fg==1)&(X>0),select=c(R,X,m,year))
names(d)[1:2]=c("RT","XT")
d=d[-(1:1),]  # remove bogus first data point since 90 is built into model anyway, or it should be estimated freely
d
with(d,plot(XT,m,type="p",pch=1, xlab="[ATP] (uM)", ylab="Mass (kDa)",log="x"))

hill4<-nls(m~M1+(M6-M1)*(XT/S50)^h/(1+(XT/S50)^h),d,start=list(M6=540,M1=90,S50=400,h=3))
hill4
hill3<-nls(m~M1+(5*M1)*(XT/S50)^h/(1+(XT/S50)^h),d,start=list(M1=90,S50=400,h=3))
hill3
mtext(paste("M1 = ",format(hill4$m$getPars()["M1"],digits=3),sep=""),line=-4.2,side=1,font=1,cex=0.7,adj=0.9)
#mtext(paste("M1 = ",format(hill3$m$getPars()["M1"],digits=3),sep=""),line=-4.2,side=1,font=1,cex=0.7,adj=0.1)
mtext(paste("M6 = ",format(hill4$m$getPars()["M6"],digits=3),sep=""),line=-3.2,side=1,font=1,cex=0.7,adj=0.9)
mtext(paste("S50 = ",format(hill4$m$getPars()["S50"],digits=3),sep=""),line=-2.2,side=1,font=1,cex=0.7,adj=0.9)
mtext(paste("h = ",format(hill4$m$getPars()["h"],digits=3),sep=""),line=-1.2,side=1,font=1,cex=0.7,adj=0.9)
mtext(paste("N = ",length(d$m),sep=""),line=-5.2,side=1,font=1,cex=0.7,adj=0.9)
lgx=log(d$XT)
upr=range(lgx)[2]
lwr=range(lgx)[1]
del=(upr-lwr)/50
fineX=exp(seq(lwr,upr,by=del))
lines(fineX,predict(hill4,list(XT=fineX)),col="black",lwd=1)
lines(fineX,predict(hill3,list(XT=fineX)),col="black",lwd=1)

plot(hill4$m$fitted(),hill4$m$resid(),xlab="Fitted Value",
    ylab="Residual",mar=c(2,2,0,1)+.1)
#points(hill3$m$fitted(),hill3$m$resid(),col="red",lwd=1)




