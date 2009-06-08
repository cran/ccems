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
  windows(width = 10, height = 4,restoreConsole = TRUE,ypos=100) else X11(width=10,height=6)
labs=c("L3","L4","L5")
labs=c("DDLL.DDDD","DDDD.DFFF","DFFM.DDDD","DFFM.DFFF")
labs=labs[1:2]
options(digits=3)
par(mfrow=c(2,5),mar=c(2,2,1,0)+.1,cex=1,oma=c(2,2,0,0))
for (imod in 1:length(labs) ) {
  gk=rep(0,4) # global average
  gK=rep(0,4)
  ngk=rep(0,4) # global average
  ngK=rep(0,4)
  mngk=0
  mngK=0
  for (i in 1:5) {
    if (imod==1) {Keq=c(E1S3_S="E1S2_S",E1S1_S="E1S0_S"); keq=c(kE1S2="kE1S1",kE1S3="kE1S1",kE1S4="kE1S1")}
    if (imod==2) {Keq=c(E1S3_S="E1S0_S",E1S2_S="E1S0_S",E1S1_S="E1S0_S"); keq=c(kE1S3="kE1S2",kE1S4="kE1S2")}
    if (imod==3) {Keq=c(E1S2_S="E1S1_S");                 keq=c(kE1S2="kE1S1",kE1S3="kE1S1",kE1S4="kE1S1")}
    if (imod==4) {Keq=c(E1S2_S="E1S1_S");                 keq=c(kE1S3="kE1S2",kE1S4="kE1S2")}
    d=subset(TK1,index==i,select=c(E,S,k,frstAut,year))
    if (i==3) d=d[-1,]
    if (i==1) d=d[-(14:16),]
#    d=transform(d,E=E/4) # doesn't matter in this limit
    plot(d$S,d$k,xlab="Total [dT]",log="xy", ylab="k (1/sec)",cex.main=.7, 
        main=paste(d[1,"frstAut"],d[1,"year"]))
    names(d)[1:2]= c("ET","ST")

    Kmapping=mkKd2Kj(g)
    Kic=rep(1,4)
    names(Kic)<-g$KS
#    names(Kic)<-c("E1S0_S","E1S1_S","E1S2_S","E1S3_S")
    kic=rep(5,4)
#  if (i==3) kic=rep(2,4)
    if (i==5) kic=rep(.25,4)
#    names(kic)<-c("kE1S1","kE1S2","kE1S3","kE1S4")
    names(kic)<-g$kS
    mdl=mkModel(g,"LA",d,Kdparams=Kic,Keq=Keq, Kd2KjLst=Kmapping,kparams=kic,keq=keq)
    fmdl=fitModel(mdl)
    lgx=log(d$ST)
    upr=range(lgx)[2]
    lwr=range(lgx)[1]
    del=(upr-lwr)/50
    fineX=exp(seq(lwr,upr,by=del))
    predict <- data.frame(ET = rep(d$ET[1],length(fineX)), ST = fineX)
    df <- simulateData(fmdl,predict=predict,typeYP="k")$predict  
    lines(df$ST,df$EY) 
    options(digits=2)
    Kstrn="         K = "
    for (jj in 1:4) Kstrn=paste(Kstrn,sprintf("%3.2f",as.numeric(as.character(fmdl$report[g$KS[jj],"final"]))),sep=ifelse(jj>1,", ",""))
    Kstrn 
    kstrn="         k = "
    for (jj in 5:8) kstrn=paste(kstrn,sprintf("%3.2f",as.numeric(as.character(fmdl$report[g$kS[jj-4],"final"]))),sep=ifelse(jj>5,", ",""))
    kstrn
    mtext(Kstrn,line=-1.1,side=1,font=1,cex=0.55,adj=0.9)
    mtext(kstrn,line=-1.7,side=1,font=1,cex=0.55,adj=0.9)
    if (i==1) mtext(labs[imod],line=3.2,side=2,font=1,cex=0.7,adj=0.5)
#      plot(fmdl$d$EY,fmdl$tres,xlab="Fitted Value",ylab="Residual")
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
  Kstrn="         K = ("
  for (jj in 1:4) Kstrn=paste(Kstrn,sprintf("%3.2f",ngK[jj]),sep=ifelse(jj>1,", ",""))
  Kstrn=paste(Kstrn,")")
  kstrn="         k = ("
  for (jj in 1:4) kstrn=paste(kstrn,sprintf("%3.2f",ngk[jj]),sep=ifelse(jj>1,", ",""))
  kstrn=paste(kstrn,")")
  print(Kstrn)
  print(kstrn)
}

