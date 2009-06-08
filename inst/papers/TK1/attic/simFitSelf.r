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
labs=c("DDLL.DDDD","DDDD.DFFF","DFFM.DDDD","DFFM.DFFF","DFLM.DFLM")
labs=labs[5]
if (length(labs)==1) windows(width = 6, height = 6,restoreConsole = TRUE,ypos=100) else
   windows(width =2.5*length(labs), height = 2,restoreConsole = TRUE,ypos=100) 
par(mfrow=c(1,length(labs)),mar=c(5,4,2,1)+.1)
# par(mfrow=c(1,1),mar=c(5,4,1,0)+.1)
 for (modS in labs ) {
  if (modS=="DDLL.DDDD") {
    Keq=c(E1S3_S="E1S2_S",E1S1_S="E1S0_S"); 
    keq=c(kE1S2="kE1S1",kE1S3="kE1S1",kE1S4="kE1S1")
    Kt = c(0.9, 0.9, 0.5, 0.5) 
    kt = c(4,4,4,4)
  }
  if (modS=="DDDD.DFFF") {
    Kt = c(0.77, 0.77, 0.77, 0.77) 
    kt = c(3, 4.6,   4.6, 4.6) 
    Keq=c(E1S3_S="E1S0_S",E1S2_S="E1S0_S",E1S1_S="E1S0_S"); 
    keq=c(kE1S3="kE1S2",kE1S4="kE1S2")
  }
  if (modS=="DFFM.DDDD") {
    Keq=c(E1S2_S="E1S1_S");                 
    keq=c(kE1S2="kE1S1",kE1S3="kE1S1",kE1S4="kE1S1")
    Kt=c( 0.9,0.7,0.7,0.55)
    kt=c(4.1,4.1,4.1,4.1) 
  }
  if (modS=="DFFM.DFFF") {
    lgx=log(c(1e-2,1e1))
    upr=range(lgx)[2]
    lwr=range(lgx)[1]
    del=(upr-lwr)/100
    SEQ=exp(seq(lwr,upr,by=del))
#    SEQ=c(1e-6,1e-5,1e-4,SEQ,1e2,1e3,1e4)
    fineX=rep(SEQ,each=1) 
    predict <- data.frame(ET = rep(.0001,length(fineX)), ST = fineX)
    Keq=c(E1S2_S="E1S1_S"); 
    keq=c(kE1S3="kE1S2",kE1S4="kE1S2")
    Kt=c( 0.9,0.7,0.7,0.55)
    kt=c(3.7,4.2,4.2,4.2) 
  }
  if (modS=="DFLM.DFLM") {
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
#    kt=sort(kt)
  }
  names(Kt)<-KS
  names(kt)<-kS
  mdl = mkModel(g,modS,Kdparams=Kt, Kd2KjLst=Kmapping,Keq=Keq,kparams=kt,keq=keq)
  df <- simulateData(mdl,predict=predict,typeYP="k")$predict  
  names(df)[3]<-"k"
#  df=transform(df,k=rnorm(k,k,0.05*k))
#df=transform(df,weights=1000)
  plot(df$ST,df$k,xlab="Total [S]",log="x", ylab="k (1/sec)",cex.main=ifelse(length(labs)==1,1,.7),
      cex.lab=ifelse(length(labs)==1,1.2,1),
      pch=19,main=modS,cex=0.7)
#  lines(df$ST,df$k)
  Kic= c(0.6, .6,.6,.6);          kic=c(4,4,4,4)
  if (modS=="DFLM.DFLM") { Kic= Kt/2;          kic=kt/2}
  names(Kic)<-KS
  names(kic)<-kS
  abline(h=kt*(1:4)/4,v=Kt*c(.25,.67,1.5,4))#,col="red")
#  abline(h=c(25,5,100,25),v=Kt)
  mdl2 = mkModel(g,modS,df,Kdparams=Kic, Kd2KjLst=Kmapping,Keq=Keq,kparams=kic,keq=keq,transform="none")
  fmdl=fitModel(mdl2)
  print(modS)
  x<-rbind(true=c(Kt,kt))
  print(rbind(x,t(fmdl$report[c(KS,kS),c("initial","final"),drop=FALSE])))
  lines(fmdl$d$ST,fmdl$d$EY)
}







