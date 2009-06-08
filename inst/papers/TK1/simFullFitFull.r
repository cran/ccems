rm(list=ls(all=TRUE))  # clean up left overs from previous run, though shouldn't need to.
if (.Platform$OS.type=="windows") home="/users/radivot" else home="/home/radivot"
setwd(home)  # place where models and results will go
source(paste(home,"/start.r",sep=""))  # define host="machineName" in this file
if (0) library(ccems) else { # if 0 source in the package to save install time 
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
if (!is.null(dev.list())) for (i in 2:max(dev.list())) dev.off(i);
options(digits=3)
SEQ=seq(.01,12,by=.01)
SEQ=c(0.001,0.002,0.005,0.01,0.02,0.05, .1, .2, .5, 1, 2, 5,10,20 )
lgx=log(c(.05,1.2))
upr=range(lgx)[2]
lwr=range(lgx)[1]
del=(upr-lwr)/100
SEQ=c(0.01,0.02,0.05, .1, .2, .5, 1, 2, 5,10 )
SEQ=exp(seq(lwr,upr,by=del))
fineX=rep(SEQ,each=1) 
predict <- data.frame(ET = rep(.0001,length(fineX)), ST = fineX)
Kmapping = mkKd2Kj(g)
 windows(width = 5, height = 5,restoreConsole = TRUE,ypos=100) 
 par(mfrow=c(1,1),mar=c(5,4,1,1)+.1)
 modS="DFLM.DFLM"
 for (j in c("A")) {
#   for (j in c("A","B")) {
    lgx=log(c(1e-8,1e7))
    upr=range(lgx)[2]
    lwr=range(lgx)[1]
    del=(upr-lwr)/100
    SEQ=exp(seq(lwr,upr,by=del))
    fineX=rep(SEQ,each=1) 
    predict <- data.frame(ET = rep(.0001,length(fineX)), ST = fineX)
    Keq=NULL
    keq=NULL
    Kt=c( 2,1,0.5,0.1)
    kt=c(1,2,3,4)
    Kt=c(1e-6,1e-3,1,1e3)
    kt=c(4,2,4/3,1)
    kt=c(100,10,40,10)
   if (j=="B") kt=sort(kt)
  names(Kt)<-KS
  names(kt)<-kS
  mdl = mkModel(g,"full",Kdparams=Kt, Kd2KjLst=Kmapping,Keq=Keq,kparams=kt,keq=keq)
  df <- simulateData(mdl,predict=predict,typeYP="k")$predict  
  names(df)[3]<-"k"
  plot(df$ST,df$k,xlab="Total [S]",log="x", ylab="k (1/sec)",cex.main=1,  #ylim=c(0,ifelse(j=="A",31,101)),
      cex.lab=1.2,
      pch=19,main=modS,cex=0.7)
#  if (j=="B") text(3e-8,95,"B",cex=1.2) else text(3e-8,28,"A",cex=1.2)
  Kic= c(0.6, .6,.6,.6);          kic=c(4,4,4,4)
  Kic= Kt/2 
  kic=kt/2
  names(Kic)<-KS
  names(kic)<-kS
  abline(h=kt*(1:4)/4,v=Kt*c(.25,.67,1.5,4))#,col="red")
  mdl2 = mkModel(g,"full",df,Kdparams=Kic, Kd2KjLst=Kmapping,Keq=Keq,kparams=kic,keq=keq,transform="none")
  fmdl=fitModel(mdl2)
  print(modS)
  x<-rbind(true=c(Kt,kt))
  print(rbind(x,t(fmdl$report[c(KS,kS),c("initial","final"),drop=FALSE])))
  lines(fmdl$d$ST,fmdl$d$EY)
}







