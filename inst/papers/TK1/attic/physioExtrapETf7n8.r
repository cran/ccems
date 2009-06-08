# this file shows all the data in the .1 to 1.2 physiological range
# model averages are plotted and literature averages taken

# errors created in extrapolations of literature average to high [ET] are also covered here

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
if (!is.null(dev.list())) for (i in 2:max(dev.list())) dev.off(i);
if (.Platform$OS.type=="windows")  # now create a window for ccems fits 
  windows(width = 8, height = 4,restoreConsole = TRUE) else X11(width=4,height=4)
par(mfcol=c(1,2),mar=c(4,4,2,1)+.1)
load("case/results/outs31")
KS=c("E1S0_S","E1S1_S","E1S2_S","E1S3_S")
kS=c("kE1S1","kE1S2","kE1S3","kE1S4")
hill<-function(S,kmax,S50,h) kmax*(S/S50)^h/(1+(S/S50)^h) 
nms=NULL
gk=rep(0,4) # global average
gK=rep(0,4)

ngk=rep(0,4) # global average
ngK=rep(0,4)
mngk=0
mngK=0
fineX=seq(.1,1.2,by=.1)
predict <- data.frame(ET = rep(.0001,length(fineX)), ST = fineX)

Kmapping=mkKd2Kj(g)

for (i in 1:5) {
  d=subset(TK1,(index==i)&(S>0.09)&(S<1.3),select=c(E,S,k,frstAut,year))
  if (i==1) plot(d$S,d$k,pch=i,xlab="Total [dT]",log="xy", ylab="k (1/sec)",ylim=c(0.02,7),xlim=c(0.095,1.2)) else
      points(d$S,d$k,pch=i)   
  Kic=outs[[i]]$ma[1:4]
  names(Kic)<-KS
  nms=c(nms,paste(d[1,"frstAut"],d[1,"year"]))
  kic=outs[[i]]$ma[,5:8]
  names(kic)<-kS
  mdl=mkModel(g,"ma",Kdparams=Kic, Kd2KjLst=Kmapping,kparams=kic)
  df <- simulateData(mdl,predict=predict,typeYP="k")$predict  
  lines(df$ST,df$EY,lty=i)
  gk=gk+kic
  gK=gK+Kic
  mngk=mngk+mean(kic)
  mngK=mngK+mean(Kic)
  ngk=ngk+kic/mean(kic)
  ngK=ngK+Kic/mean(Kic)
}
legend(.1,.3,bty="n",legend=nms,lty=1:5,pch=1:5,cex=0.6)
text(.1,6,"A")
gk=gk/5
gK=gK/5
options(digits=2)
gk
gK
ngk=(mngk/5)*ngk/5  # average amplitudes and shapes separately
ngK=(mngK/5)*ngK/5
ngk
ngK
mdl=mkModel(g,"LA",Kdparams=ngK, Kd2KjLst=Kmapping,kparams=ngk) # literature average
df <- simulateData(mdl,predict=predict,typeYP="k")$predict  
lines(df$ST,df$EY,lty=1,lwd=2)
names(df)[3]="k"
#  *******  END Figure 7A above

## ***** Start Figure 7B below *****
#if (.Platform$OS.type=="windows")  # now create a window for ccems fits 
#  windows(width = 4, height = 4,restoreConsole = TRUE,ypos=-25) else X11(width=4,height=4)
#par(mfcol=c(1,1),mar=c(4,4,0,1)+.1)
plot(df$ST,df$k,pch=1,xlab="Total [dT]",log="", ylab="k (1/sec)",ylim=c(0.2,2.7),xlim=c(.1,1.2))
lines(df$ST,df$k,lty=1)
text(.13,2.63,"B")


ng=mkGrids(g,minTotalPs=3,maxTotalPs=3,KIC=0.6,kIC=4)
chunk=ng$chunk
Keqs=ng$Keqs
keqs=ng$keqs
mdlNames=rownames(chunk)
lmdlNames=strsplit(mdlNames,split=".",fixed=TRUE)
names(lmdlNames)<-mdlNames
#  print(mdlNames)
#  print(lmdlNames)
  models=NULL
  for (j in mdlNames) 
     models[[j]]=mkModel(g,j,df,Kdparams=chunk[j,2:(g$nZ+1),drop=FALSE], 
          Keq=Keqs[[lmdlNames[[j]][1]]],Kd2KjLst=Kmapping,
          pparams=chunk[j,"p",drop=FALSE],indx=chunk[j,"indx"],nParams=chunk[j,"nParams"],
          kparams=chunk[j,paste("k",g$Z,sep="")],keq=keqs[[lmdlNames[[j]][2]]])
models
fmodels=lapply(models,fitModel)

for (j in mdlNames) 
  lines(fmodels[[j]]$d$ST,fmodels[[j]]$d$EY,lty=2,lwd=0.5)

hillda<-nls(k~kmax*(ST/S50)^h/(1+(ST/S50)^h),df,start=list(kmax=4,S50=1,h=1))
#lines(fineX,predict(hillda,list(S=fineX)),col="red",lwd=2) # to show it is perfect
hillda
#  end top plot at low ET

g <-mkg(topology, activity=TRUE,TCC=TRUE)
mdl=mkModel(g,"LA",Kdparams=ngK, Kd2KjLst=Kmapping,kparams=ngk)
for (j in mdlNames) 
{
  fmodels[[j]]$TCC=TRUE
  fmodels[[j]]$code=g$code
}

ET1=0.1
predict <- data.frame(ET = rep(ET1,length(fineX)), ST = fineX)
df <- simulateData(mdl,predict=predict,typeYP="k")$predict  
points(df$ST,df$EY,pch=2,col=1)
lines(df$ST,df$EY,lty=1)
#k=1
#for (j in mdlNames) {
#  df <- simulateData(fmodels[[j]],predict=predict,typeYP="k")$predict  
#  lines(df$ST,df$EY,lty=1)
##  points(df$ST,df$EY,pch=k,col=k)
#  k=k+1
#  locator(n=1)
#}
for (j in mdlNames) {
  df <- simulateData(fmodels[[j]],predict=predict,typeYP="k")$predict  
  lines(df$ST,df$EY,lty=2,lwd=0.5)
}
legend(.15,2.5,legend=c("[ET]=.0001 uM","[ET]=0.1 uM","[ET]=0.6 uM"),pch=1:3,cex=0.7,bty="n")
#  end middle plot at ET = 0.1 uM

ET1=0.6
predict <- data.frame(ET = rep(ET1,length(fineX)), ST = fineX)
df <- simulateData(mdl,predict=predict,typeYP="k")$predict  
points(df$ST,df$EY,pch=3,col=1)
lines(df$ST,df$EY,lty=1)
for (j in mdlNames) {
  df <- simulateData(fmodels[[j]],predict=predict,typeYP="k")$predict  
  lines(df$ST,df$EY,lty=2,lwd=0.5)
}

t(fmodels[[3]]$report[c(KS,kS),"final",drop=FALSE])
t(fmodels[[4]]$report[c(KS,kS),"final",drop=FALSE])



