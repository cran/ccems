rm(list=ls(all=TRUE))  # clean up left overs from previous run, though shouldn't need to.
if (.Platform$OS.type=="windows") home="/users/radivot" else home="/home/radivot"
setwd(home)  # place where models and results will go
source(paste(home,"/start.r",sep=""))  # define host="machineName" in this file
if (1) library(ccems) else { # if 0 source in the package to save install time 
  pkgNms=dir(paste(home,"/case/active/ccems/ccems/R",sep=""),full.names=TRUE)
  for (i in pkgNms) source(i)  # source in all of the R files of the package
  load(paste(home,"/case/active/papers/TK1/data/TK1.rda",sep=""))
}

topology <- list(  
    heads=c("E1S0"), # E1S0 = substrate free E
    sites=list(                    
        c=list(    # c for catalytic site  
            t=c("E1S1","E1S2","E1S3","E1S4")   
        ) # t for tetramer 
    )
)   # in transform below, TK1 is 25kDa => 25mg/umole
g <-mkg(topology, activity=TRUE,TCC=FALSE)

KS=c("E1S0_S","E1S1_S","E1S2_S","E1S3_S")
kS=c("kE1S1","kE1S2","kE1S3","kE1S4")
## Warning: the next line clears all existing figures!!
if (!is.null(dev.list())) for (i in 2:max(dev.list())) dev.off(i);
# this is the K trumps k section
if (.Platform$OS.type=="windows")  # now create a window for ccems fits 
#  windows(width = 9.9, height = 6.6,restoreConsole = TRUE,ypos=100) else X11(width=10,height=2)
  windows(width = 8.4, height = 5.6,restoreConsole = TRUE,ypos=100) else X11(width=10,height=2)
#par(mfrow=c(2,3),mar=c(4,4,2,1)+.1,cex=1.1)
par(mfrow=c(2,3),mar=c(3,3,1,0)+.1,cex=1.1,oma=c(1,1,0,0))

getSlope<-function(K) {
  x=1:4
  yK=as.numeric(K)
  yK=yK/mean(yK)
  lmK=lm(yK~x)
  coef(lmK)["x"]
}


mymed<-function(x) sapply(x,median)

REP=FALSE
REP=TRUE
FULL=TRUE # make full fit full space
FULL=FALSE # make Lit Ave (LA5) fit full space
LA5=TRUE
LA5=FALSE
LA4=FALSE
LA4=TRUE
if (FULL) {
  load(ifelse(REP,"case/results/repFullouts","case/results/fullouts"))
  Kic=c( 2,1,0.5,0.1)
  kic=c(1,2,3,4)
  cp=1.5
} 
if (LA5){
  load(ifelse(REP,"case/results/repLA5outs","case/results/LA5outs"))
  cp=.55
  Kic=c( 0.9,0.7,0.7,0.55)
  kic=c(3.9,4.2,4.2,4.2) 
}
if (LA4){
  load(ifelse(REP,"case/results/repLA4outs","case/results/LA4outs"))
  cp=.55
  Kic=c( 0.9,0.7,0.7,0.55)
  kic=c(4.1,4.1,4.1,4.1) 
}
getMA<-function(df) {
  M=as.matrix(df[,c(KS,kS)])
  wgts=df$wgts/sum(df$wgts)
  exp(wgts%*%log(M)) # average in space of gibbs free energy changes
}

outs[[1]]$ma
getMA(outs[[1]]$df)


truepnt=c(getSlope(Kic),getSlope(kic))
globMApnt=c(getSlope(outs[[1]]$ma[1:4]),getSlope(outs[[1]]$ma[5:8]))
gains=c(1,0.7)
mypch=c(22,15)
mypch=c(1,22)
mypch=c(1,19)
for (ip in 3:8) {
  aK=NULL
  ak=NULL
  run=NULL
  bdf=NULL
  for (i in 1:length(outs)) {
    df=subset(outs[[i]]$df,subset=(numP==ip))
#  df=subset(outs[[i]]$df,subset=wgts>1e-6)
    rownames(df)<-paste(rownames(df),i,sep="")
    n=dim(df)[1]
#    print(n)
#    print(df)
    run=c(run,rep(i,n))
    for (j in 1:n) {
      aK=c(aK,getSlope(df[j,5:8]))
      ak=c(ak,getSlope(df[j,9:12]))
    }
    bdf=rbind(bdf,df)
  }
  df=transform(bdf,aK=aK,ak=ak,run=run)
  if (ip==3) df=df[-n,]
  print(df)
  #  df=data.frame(df,aK=aK,ak=ak)
  ur=dim(subset(df,subset=(ak>0)&(aK>0)))[1]
  lr=dim(subset(df,subset=(ak<0)&(aK>0)))[1]
  ll=dim(subset(df,subset=(ak<0)&(aK<0)))[1]
  ul=dim(subset(df,subset=(ak>0)&(aK<0)))[1]
  nruns=length(unique(df$run))
  with(df,plot(aK,ak,ylim=c(-cp,cp),pch=mypch[run],cex=gains[run],xlim=c(-cp,cp),col=1,cex.lab=1.2,cex.main=1,
#          ylab="k slope",xlab="K slope",main=paste(ip,"parameters  ",dim(df)[1],"models")))
#        ylab=NA,xlab=NA,cex.main=0.7,main=paste(ip,"parameters  ",dim(df)[1]/nruns,ifelse(ip<8,"models  ","model  "),nruns,"runs")))
          ylab=NA,xlab=NA,cex.main=0.7,main=paste(ip,"parameters  ",dim(df)[1],ifelse(ip<8,"models  ","model  "))))
#with(df,points(aK,ak,pch=1,cex=5*wgts/sum(wgts)))
#ylab=ifelse((ip==3)|(ip==6),"k slope",NA),xlab=NA,cex.main=0.9,main=paste(ip,"parameters   ",dim(df)[1],"models")))
  #  if (ip==3) 
  title(ylab="k slope",xlab="          K slope",outer=TRUE,line=-0.5,cex.lab=1.5)
  points(truepnt[1],truepnt[2],col=1,pch=4,cex=2,lwd=2)
  points(globMApnt[1],globMApnt[2],col=1,pch=4,cex=2,lwd=1)
  curMA=getMA(df)
  print(curMA)
  curMApnt=c(getSlope(curMA[1:4]),getSlope(curMA[5:8]))
  points(curMApnt[1],curMApnt[2],col=1,pch=1,cex=2,lwd=2)
  #  s<-by(df,run,mymed)
#  print(s)
#  for (i in 1:length(outs))  # i.e. i over runs
#      points(s[[i]]["aK"],s[[i]]["ak"],col=1,pch=mypch[i],cex=1.2)
  abline(h=0,v=0)
  text(.9*cp,.9*cp,ur)
  text(.9*cp,-.9*cp,lr)
  text(-.9*cp,-.9*cp,ll)
  text(-.9*cp,.9*cp,ul)
#  print(df[1:4,])
}  # loop on parameter number ip




