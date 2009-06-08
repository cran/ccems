rm(list=ls(all=TRUE))  # clean up left overs from previous run, though shouldn't need to.
if (.Platform$OS.type=="windows") home="/users/radivot" else home="/home/radivot"
setwd(home)  # place where models and results will go
source(paste(home,"/start.r",sep=""))  # define host="machineName" in this file
if (0) library(ccems) else { # if 0 source in the package to save install time 
  pkgNms=dir(paste(home,"/case/active/ccems/ccems/R",sep=""),full.names=TRUE)
  for (i in pkgNms) source(i)  # source in all of the R files of the package
  load(paste(home,"/case/active/papers/TK1/data/TK1.rda",sep=""))
}
# in S-phase TK1 is 4 ug/ml 
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
fineX=seq(.1,1.2,by=.1)
lgx=log(fineX)
upr=range(lgx)[2]
lwr=range(lgx)[1]
del=(upr-lwr)/50
fineX=exp(seq(lwr,upr,by=del))
predict <- data.frame(ET = rep(.0001,length(fineX)), ST = fineX)
testCases=list(
    list(Kic=0.6*c( 1,1,1,1),  kic=c(4,3,2,1)),
    list(Kic=0.6*c( 1,1,1,1),  kic=c(2,1.6,1.4,1)),
    list(Kic=0.6*c( 1,1,1,1),   kic=c(1.3,1.2,1.1,1)),
    list(Kic=0.6*c( 1,1,1,1), kic=c(1,1,1,1)),
    list(Kic=0.6*c( 1,1,1,1),  kic=c(.4,.6,.8,1)),
    list(Kic=0.6*c( 0.6,0.8,1.2,1.4),  kic=c( 1,1,1,1)),
    list(Kic=0.6*c( 1.4,1.2,.8,.6),    kic=c( 1,1,1,1)),
    list(Kic=0.6*c( .1,.5,1.5,10),  kic=c(1,1,1,1)),
    list(Kic=0.6*c( 1.4,1.2,.8,.6),    kic=c(1,1,1,1)),
    list(Kic=0.6*c( 1,1,1,1), kic=c(1,1,1,1))
)

Kmapping = mkKd2Kj(g)
## Warning: the next line clears all existing figures!!
if (!is.null(dev.list())) for (i in 2:max(dev.list())) dev.off(i);

if (.Platform$OS.type=="windows")  # now create a window for ccems fits 
  windows(width = 12.5, height = 5,restoreConsole = TRUE,ypos=-25) else X11(width=12.5,height=5)
par(mfcol=c(2,5),mar=c(4,4,2,1)+.1)
# set up MM reference line
Kic=0.6*c( 1,1,1,1)  
kic=c(1,1,1,1)  
names(Kic)<-KS
names(kic)<-kS
mdl = mkModel(g,"MM",Kdparams=Kic, Kd2KjLst=Kmapping,kparams=kic)
dfmm <- simulateData(mdl,predict=predict,typeYP="k")$predict  
names(dfmm)[3]<-"k"

for (i in 1:5) {
  Kic=testCases[[i]]$Kic  
  kic=testCases[[i]]$kic  
  names(Kic)<-KS
  names(kic)<-kS
  Kstrn="         K = "
  for (j in KS) 
    Kstrn=paste(Kstrn,sprintf("%3.2f",Kic[j]),sep=ifelse(j!="E1S0_S",", ",""))
  Kstrn 
  kstrn="         k = "
  for (j in kS) 
    kstrn=paste(kstrn,sprintf("%3.2f",kic[j]),sep=ifelse(j!="kE1S1",", ",""))
  kstrn
  
  mdl = mkModel(g,"full",Kdparams=Kic, Kd2KjLst=Kmapping,kparams=kic)
  df <- simulateData(mdl,predict=predict,typeYP="k")$predict  
  names(df)[3]<-"k"
  plot(dfmm$S,dfmm$k,xlab="[ST]",log="x", ylab="k",lty=1,type="l",lwd=2)
  lines(df$S,df$k,lty=1,type="l")
  mtext(Kstrn,line=-1.2,side=1,font=1,cex=0.7,adj=1)
  mtext(kstrn,line=-2.2,side=1,font=1,cex=0.7,adj=1)
  plot(log(dfmm$S),log(dfmm$k/(1-dfmm$k)),xlab="log [ST]",log="", ylab="log(k/(kmax-k))",lty=1,type="l",lwd=2)
  lines(log(df$S),log(df$k/(1-df$k)),lty=1,type="l")
}


