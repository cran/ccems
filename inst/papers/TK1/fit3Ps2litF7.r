rm(list=ls(all=TRUE))  # clean up left overs from previous run, though shouldn't need to.
if (.Platform$OS.type=="windows") home="/users/radivot" else home="/home/radivot"
setwd(home)  # place where models and results will go
if (0) library(ccems) else { # if 0 source in the package to save install time 
  pkgNms=dir(paste(home,"/case/active/ccems/ccems/R",sep=""),full.names=TRUE)
  for (i in pkgNms) source(i)  # source in all of the R files of the package
  load(paste(home,"/case/active/papers/TK1/data/TK1.rda",sep=""))
}
TK1

# This file does literature model average fits to the data. 
if (.Platform$OS.type=="windows")  # now create a window for ccems fits 
  windows(width = 12.5, height = 5,restoreConsole = TRUE) else X11(width=10,height=4)
topology <- list(  
    heads=c("E1S0"), # E1S0 = substrate free E
    sites=list(                    
        c=list(    # c for catalytic site  
            t=c("E1S1","E1S2","E1S3","E1S4")   
        ) # t for tetramer 
    )
)   # in transform below, TK1 is 25kDa => 25mg/umole
g <-mkg(topology, activity=TRUE,TCC=FALSE)

#KS=c("E1S0_S","E1S1_S","E1S2_S","E1S3_S")
#kS=c("kE1S1","kE1S2","kE1S3","kE1S4")
if (!is.null(dev.list())) for (i in 2:max(dev.list())) dev.off(i);
if (.Platform$OS.type=="windows")  # now create a window for ccems fits 
#  windows(width = 12.5, height = 7.5,restoreConsole = TRUE,ypos=100) else X11(width=10,height=6)
  windows(width = 12.5, height = 9.6,restoreConsole = TRUE,ypos=10) else X11(width=10,height=6)

ng=mkGrids(g,minTotalPs=3,maxTotalPs=3,KIC=0.6,kIC=4,m1=-100)
chunk=ng$chunk
#chunk=chunk[c(1,2,5,6),]
chunk
Kmapping=mkKd2Kj(g)
Keqs=ng$Keqs
keqs=ng$keqs
mdlNames=rownames(chunk)
lmdlNames=strsplit(mdlNames,split=".",fixed=TRUE)
names(lmdlNames)<-mdlNames
labs=mdlNames

labs=c("DFFF.DDDD","DDDM.DDDD","DDDD.DDLL","DDDD.DDDM")
options(digits=3)
par(mfrow=c(6,5),mar=c(2,2,1,0)+.1,cex=1,oma=c(2,2,0,0))
par(mfrow=c(1,1),mar=c(2,2,1,0)+.1,cex=1,oma=c(2,2,0,0))
aves=NULL
for (j in "DDDD.DDDM") 
#  for (j in mdlNames) 
{
  gk=rep(0,4) # global average
  gK=rep(0,4)
  ngk=rep(0,4) # global average
  ngK=rep(0,4)
  mngk=0
  mngK=0
  for (i in 1:1) {
      d=subset(TK1,index==i,select=c(E,S,k,frstAut,year))
      if (i==3) d=d[-1,]
      if (i==1) d=d[-(14:16),] 
      plot(d$S,d$k,xlab="Total [dT]",log="xy", ylab="k (1/sec)",cex.main=.7, 
          main=paste(d[1,"frstAut"],d[1,"year"]))
      names(d)[1:2]= c("ET","ST")
      kic=chunk[j,paste("k",g$Z,sep="")]
      if (i==5) {kic=rep(.25,4); names(kic)<-g$kS}
      mdl=mkModel(g,j,d,Kdparams=chunk[j,2:(g$nZ+1),drop=FALSE], 
          Keq=Keqs[[lmdlNames[[j]][1]]],Kd2KjLst=Kmapping,
          pparams=chunk[j,c("p","m1"),drop=FALSE],indx=chunk[j,"indx"],nParams=chunk[j,"nParams"],
          kparams=kic,keq=keqs[[lmdlNames[[j]][2]]])
    fmdl=fitModel(mdl)
    print(j)
    lines(fmdl$d$ST,fmdl$d$EY,lty=1,lwd=1)
    Kstrn="         K = "
    for (jj in 1:4) Kstrn=paste(Kstrn,sprintf("%3.2f",as.numeric(as.character(fmdl$report[g$KS[jj],"final"]))),sep=ifelse(jj>1,", ",""))
    Kstrn 
    kstrn="         k = "
    for (jj in 5:8) kstrn=paste(kstrn,sprintf("%3.2f",as.numeric(as.character(fmdl$report[g$kS[jj-4],"final"]))),sep=ifelse(jj>5,", ",""))
    kstrn
    mtext(Kstrn,line=-1.1,side=1,font=1,cex=0.55,adj=0.9)
    mtext(kstrn,line=-1.7,side=1,font=1,cex=0.55,adj=0.9)
    if (i==1) mtext(j,line=3.2,side=2,font=1,cex=0.7,adj=0.5)
    Kic=as.numeric(as.character(fmdl$report[g$KS,"final"]))
    kic=as.numeric(as.character(fmdl$report[g$kS,"final"]))
    gk=gk+kic
    gK=gK+Kic
    mngk=mngk+mean(kic)
    mngK=mngK+mean(Kic)
    ngk=ngk+kic/mean(kic)
    ngK=ngK+Kic/mean(Kic)
  }
  title(ylab="k (1/sec)",xlab="Total [dT]",outer=TRUE,line=0,cex.lab=1.3)
  gk=gk/5
  gK=gK/5
  gk
  gK
  ngk=(mngk/5)*ngk/5  # average amplitudes and shapes separately
  ngK=(mngK/5)*ngK/5
  aves[[j]]=c(ngK,ngk)
  Kstrn="         K = ("
  for (jj in 1:4) Kstrn=paste(Kstrn,sprintf("%3.2f",ngK[jj]),sep=ifelse(jj>1,", ",""))
  Kstrn=paste(Kstrn,")")
  kstrn="         k = ("
  for (jj in 1:4) kstrn=paste(kstrn,sprintf("%3.2f",ngk[jj]),sep=ifelse(jj>1,", ",""))
  kstrn=paste(kstrn,")")
  print(Kstrn)
  print(kstrn)
}
aves
