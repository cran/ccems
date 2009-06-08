rm(list=ls(all=TRUE))  # clean up left overs from previous run, though shouldn't need to.
if (.Platform$OS.type=="windows") home="/users/radivot" else home="/home/radivot"
setwd(home)  # place where models and results will go
if (0) library(ccems) else { # if 0 source in the package to save install time 
  pkgNms=dir(paste(home,"/case/active/ccems/ccems/R",sep=""),full.names=TRUE)
  for (i in pkgNms) source(i)  # source in all of the R files of the package
  load(paste(home,"/case/active/papers/TK1/data/TK1.rda",sep=""))
}
TK1

# This file does full model fits to the data. This happens anyway when running makeOuts.r 
# Here the idea was to just do the full models fits to the data. The results are uninteresting noisy parameters. 

if (.Platform$OS.type=="windows")  # now create a window for ccems fits 
  windows(width = 12.5, height = 7.5,restoreConsole = TRUE) else X11(width=10,height=4)
topology <- list(  
    heads=c("E1S0"), # E1S0 = substrate free E
    sites=list(                    
        c=list(    # c for catalytic site  
            t=c("E1S1","E1S2","E1S3","E1S4")   
        ) # t for tetramer 
    )
)   # in transform below, TK1 is 25kDa => 25mg/umole
g <-mkg(topology, activity=TRUE,TCC=FALSE)

par(mfcol=c(3,5),mar=c(4,4,2,1)+.1)
CV=NULL
for (i in 1:5) {
  d=subset(TK1,index==i,select=c(E,S,k,frstAut,year))
  d=transform(d,E=E/4) # doesn't matter in this limit
  plot(d$S,d$k,xlab="Total [dT]",log="xy", ylab="k (1/sec)", 
      main=paste(d[1,"frstAut"],d[1,"year"]))
  names(d)[1:2]= c("ET","ST")
  Kmapping=mkKd2Kj(g)
  Kic=rep(1,4)
  names(Kic)<-c("E1S0_S","E1S1_S","E1S2_S","E1S3_S")
  kic=rep(5,4)
#  if (i==3) kic=rep(2,4)
  if (i==5) kic=rep(.25,4)
  names(kic)<-c("kE1S1","kE1S2","kE1S3","kE1S4")
  mdl=mkModel(g,"DFLM.DFLM",d,Kdparams=Kic, Kd2KjLst=Kmapping,kparams=kic)
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
  for (i in 1:4) Kstrn=paste(Kstrn,sprintf("%3.2f",as.numeric(as.character(fmdl$report[i,"pointEstimate"]))),sep=ifelse(i>1,", ",""))
  Kstrn 
  kstrn="         k = "
  for (i in 5:8) kstrn=paste(kstrn,sprintf("%3.2f",as.numeric(as.character(fmdl$report[i,"pointEstimate"]))),sep=ifelse(i>5,", ",""))
  kstrn
  mtext(Kstrn,line=-2,side=1,font=1,cex=0.7)
  mtext(kstrn,line=-3,side=1,font=1,cex=0.7)
  plot(fmdl$d$EY,fmdl$tres,xlab="Fitted Values",ylab="Transformed Residuals")
  plot(fmdl$d$EY,fmdl$res,xlab="Fitted Values",ylab="Residuals")
#  mtext(paste("CV = ",format(cv<-sd(fmdl$res)/fmdl$d$EY,digits=2)),line=-1,side=3,font=1,cex=0.7,adj=0.05)
  mtext(paste("CV = ",format(cv<-mean(abs(fmdl$res-mean(fmdl$res))/fmdl$d$EY),digits=2)),line=-1,side=3,font=1,cex=0.7,adj=0.05)
#  mtext(paste("<e> = ",format(mean(fmdl$res),digits=2)),line=-2,side=3,font=1,cex=0.7,adj=0.05)
  CV=c(CV,cv)
  }
  print(mean(CV))
  print(mean(CV[2:5]))
