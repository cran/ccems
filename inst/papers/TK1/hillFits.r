rm(list=ls(all=TRUE))  # clean up left overs from previous run, though shouldn't need to.
if (.Platform$OS.type=="windows") home="/users/radivot" else home="/home/radivot"
setwd(home)  # place where models and results will go
source(paste(home,"/start.r",sep=""))  # define host="machineName" in this file
if (0) library(ccems) else { # if 0 source in the package to save install time 
  pkgNms=dir(paste(home,"/case/active/ccems/ccems/R",sep=""),full.names=TRUE)
  for (i in pkgNms) source(i)  # source in all of the R files of the package
  load(paste(home,"/case/active/papers/TK1/data/TK1.rda",sep=""))
}
kill3=TRUE
kill3=FALSE
kill1=TRUE
kill1=FALSE

#TK1=subset(TK1,subset=(S>.09)&(S<1.3))
#TK1

#kill1=FALSE
## Warning: the next line clears all existing figures!!
if (!is.null(dev.list())) for (i in 2:max(dev.list())) dev.off(i);
for (j in 1:2) {
  hparams=NULL
  if (.Platform$OS.type=="windows") 
    windows(width = 10, height = 4,restoreConsole = TRUE,ypos=ifelse(j==2,-50,0)) else X11(width=10,height=4)
  par(mfcol=c(2,5),mar=c(4,4,2,1)+.1)
  for (i in 1:5) {
    d=subset(TK1,index==i,select=c(E,S,k,frstAut,year))
    if ((i==1)&(kill3)) d=d[-(14:16),]
    if ((i==3)&(kill1)) d=d[-(1),]
    if (j==1) hillda<-nls(k~kmax*(S/S50)^h/(1+(S/S50)^h),d,start=list(kmax=5,S50=1,h=1))
    if (j==2) hillda<-nls(k~kmax*(S/S50)^h/(1+(S/S50)^h),d,start=list(kmax=5,S50=1,h=1),weights=1/k^2)
    print(hillda)
    plot(d$S,d$k,xlab="Total [dT]",log="xy", ylab="k (1/sec)", 
        main=paste(d[1,"frstAut"],d[1,"year"]))
    mtext(paste("kmax = ",format(hillda$m$getPars()["kmax"],digits=3),sep=""),line=-2.2,side=1,font=1,cex=0.7)
    mtext(paste("S50 = ",format(hillda$m$getPars()["S50"],digits=3),sep=""),line=-1.2,side=1,font=1,cex=0.7)
    mtext(paste("h = ",format(hillda$m$getPars()["h"],digits=3),sep=""),line=-3.2,side=1,font=1,cex=0.7)
    mtext(paste("N = ",length(d$k),sep=""),line=-4.2,side=1,font=1,cex=0.7)
    ## Note that the specific activity is ~16 fold higher in non-Birringer data
    lgx=log(d$S)
    upr=range(lgx)[2]
    lwr=range(lgx)[1]
    del=(upr-lwr)/50
    fineX=exp(seq(lwr,upr,by=del))
    lines(fineX,predict(hillda,list(S=fineX)),col="black",lwd=1)
    plot(hillda$m$fitted(),hillda$m$resid(),xlab="Fitted Value",
        ylab=paste(ifelse(j==2,"Relative",""),"Residual"),mar=c(2,2,0,1)+.1)
    ## Note that variance increases with the mean in non-Birringer data
    ## and that the 2000 and 1993 Hill fits are poor at low k (and [S])
    hparams=rbind(hparams,hillda$m$getPars())
  }
  print("***************************************")
  print(hparams)
  cat("The literature median is",apply(hparams,2,median),"\n\n\n")
}

