# this file shows all the data in the .1 to 1.2 physiological range
# model averages are plotted and literature averages taken

# errors created in extrapolations of literature average to high [ET] are also covered here

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
#ng=mkGrids(g,minTotalPs=3,maxTotalPs=3,KIC=0.6,kIC=4)
#ng$Keq
#
#chunk=ng$chunk
#mkEq(g,chunk)




if (!is.null(dev.list())) for (i in 2:max(dev.list())) dev.off(i);
if (.Platform$OS.type=="windows")  # now create a window for ccems fits 
  windows(width = 8, height = 4,restoreConsole = TRUE) else X11(width=4,height=4)
par(mfcol=c(1,2),mar=c(4,4,2,1)+.1)

threeOnly<-function(outs,gibbs=TRUE) {
  for (i in 1:length(outs) ) {
    outs[[i]]$df=outs[[i]]$df[outs[[i]]$df$numP==3,]
    outs[[i]]$df$wgts=outs[[i]]$df$wgts/sum(outs[[i]]$df$wgts)
    M=as.matrix(outs[[i]]$df[,c(g$KS,g$kS)])
    if (gibbs)     ma=exp(outs[[i]]$df$wgts%*%log(M)) # average in space of gibbs free energy changes
    if (!gibbs)     ma=outs[[i]]$df$wgts%*%M # average in straight params
    outs[[i]]$ma=ma
  } 
  outs
}

#load("case/results/outs31")
#for (i in 1:5) print(outs[[i]]$ma)
#load("case/results/outs31G")
#for (i in 1:5) print(outs[[i]]$ma)
#
#load("case/results/outs31")
#outs=threeOnly(outs,gibbs=FALSE)
#for (i in 1:5) print(outs[[i]]$ma)
#load("case/results/outs31G")
#outs=threeOnly(outs)
#for (i in 1:5) print(outs[[i]]$ma)


avesWithin3PMods<-function(outs,normalize=TRUE) {
  BG=NULL
  nms=row.names(outs[[1]]$df)
#  print(nms)
  for (j in nms )
  {W=0
    print(j)
    gK=rep(0,4)
    gk=rep(0,4)
    ngk=rep(0,4) # global average
    mngk=0
  for (i in 1:length(outs) ) {
    W=W+outs[[i]]$df[j,"wgts"]
    k=as.numeric(outs[[i]]$df[j,g$kS])
    K=as.numeric(outs[[i]]$df[j,g$KS])
    names(K)<-g$KS
    names(k)<-g$kS
    gK=gK+K
    gk=gk+k
    mngk=mngk+mean(k)
    ngk=ngk+k/mean(k)
  } 
  if (normalize) BG[[j]]<-c("W"=W,gK/5,(mngk/5)*(ngk/5)) else  # average amplitudes and shapes separately
  BG[[j]]<-c("W"=W,gK/5,gk/5)  # average amplitudes and shapes separately
}
t(as.data.frame(BG))[c("DFFF.DDDD","DDLL.DDDD","DDDM.DDDD","DDDD.DFFF","DDDD.DDLL","DDDD.DDDM"),]
}


mkContribTab<-function(outs) {
  BG=NULL
  nms=names(outs)
#  print(nms)
  for (j in 1:length(nms) )
  {
    df=outs[[j]]$df[outs[[j]]$df[,"wgts"]>0.05,c("wgts",g$KS,g$kS)]
    names(df)<-c("Weight",paste("K",1:4,sep=""),paste("k",1:4,sep=""))
    df=data.frame(Dataset=j,Model=row.names(df),df)
#    df2=as.list(c("Dataset"=j,"Model"="average",Weights="",outs[[j]]$ma))
#    print(df2)
#    print(j)
    BG=rbind(BG,df)
#    BG=rbind(BG,df2)
  }
  row.names(BG)<-NULL
#  colnames(BG)[g$KS]=
#  colnames(BG)[g$kS]=paste("k",1:4,sep="")
  BG
}

mkMaTabs<-function(outs) {
  BG=NULL
  nms=names(outs)
  for (j in 1:length(nms) )
  {
    df=outs[[j]]$ma
    BG=rbind(BG,df)
  }
#  row.names(BG)<-NULL
  BG=data.frame(Dataset=1:length(nms),BG)
  names(BG)[2:9]<-c(paste("K",1:4,sep=""),paste("k",1:4,sep=""))
  BG
}



#load("case/results/outs01G")
load("case/results/outs31G")  # uncomment for 3Ps only where last 3 are out
outs=threeOnly(outs)   # comment to get lit ave of Fig. 4, leave in to get it for 3Ps only

library(R2HTML)
tab3P=avesWithin3PMods(outs)
colnames(tab3P)<-c("Weight",paste(" K",1:4,sep=""),paste(" k",1:4,sep=""))
tab3P

W=tab3P[,1]
M=tab3P[,-1]
LA3=t(t(M)%*%W/sum(W))
LA3
tabMAs=mkMaTabs(outs)
tabMAs

printMA<-function(x) {
 ss= sprintf("K = (%3.2f,%3.2f,%3.2f,%3.2f) k = (%3.2f,%3.2f,%3.2f,%3.2f)",x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9])
 print(ss)
}

printMA(tabMAs[3,])


tabContrib<-mkContribTab(outs)
tabContrib


target <- HTMLInitFile(home,filename="results/tab3P", BackGroundColor="#BBBBEE")
HTML(tab3P,file=target,Border=0,row.names=FALSE)
HTMLEndFile()
target <- HTMLInitFile(home,filename="results/tabMAs", BackGroundColor="#BBBBEE")
HTML(tabMAs,file=target,row.names=FALSE)
HTMLEndFile()
target <- HTMLInitFile(home,filename="results/tabContrib", BackGroundColor="#BBBBEE")
HTML(tabContrib,file=target,row.names=FALSE)
HTMLEndFile()




KS=c("E1S0_S","E1S1_S","E1S2_S","E1S3_S")
kS=c("kE1S1","kE1S2","kE1S3","kE1S4")
hill<-function(S,kmax,S50,h) kmax*(S/S50)^h/(1+(S/S50)^h) 
nms=NULL
gK=rep(0,4)
gk=rep(0,4)
ngk=rep(0,4) # global average
mngk=0
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
#  smdl <- simulateData(mdl,predict=predict,typeYP="k") 
#  print(smdl)
#  df=smdl$predict
  df <- simulateData(mdl,predict=predict,typeYP="k")$predict  
  lines(df$ST,df$EY,lty=i)
  gK=gK+Kic
  gk=gk+kic
  mngk=mngk+mean(kic)
  ngk=ngk+kic/mean(kic)
}
legend(.1,.3,bty="n",legend=nms,lty=1:5,pch=1:5,cex=0.6)
text(.1,6,"A")
gK=gK/5
gk=gk/5
options(digits=2)
gK
gk
ngk=(mngk/5)*ngk/5  # average amplitudes and shapes separately
ngk
mdl=mkModel(g,"LA",Kdparams=gK, Kd2KjLst=Kmapping,kparams=ngk) # literature average
df <- simulateData(mdl,predict=predict,typeYP="k")$predict  
lines(df$ST,df$EY,lty=1,lwd=2)
names(df)[3]="k"





#  *******  END Figure 7A above

## ***** Start Figure 7B below *****
#if (.Platform$OS.type=="windows")  # now create a window for ccems fits 
#  windows(width = 4, height = 4,restoreConsole = TRUE,ypos=-25) else X11(width=4,height=4)
#par(mfcol=c(1,1),mar=c(4,4,0,1)+.1)
plot(df$ST,df$k,pch=1,xlab="Total [dT]",log="", ylab="k (1/sec)",ylim=c(0,2.7),xlim=c(.1,1.2))
lines(df$ST,df$k,lty=1)
text(.13,2.63,"B")

df

hillda<-nls(k~kmax*(ST/S50)^h/(1+(ST/S50)^h),df,start=list(kmax=4,S50=1,h=1))
#lines(fineX,predict(hillda,list(S=fineX)),col="red",lwd=2) # to show it is perfect
hillda


ng=mkGrids(g,minTotalPs=3,maxTotalPs=3,KIC=0.6,kIC=4)
chunk=ng$chunk
chunk

Keqs=ng$Keqs
keqs=ng$keqs
Keqs
keqs


mdlNames=rownames(chunk)
lmdlNames=strsplit(mdlNames,split=".",fixed=TRUE)
names(lmdlNames)<-mdlNames
print(mdlNames)
print(lmdlNames)

models=NULL
for (j in mdlNames) 
  models[[j]]=mkModel(g,j,df,Kdparams=chunk[j,2:(g$nZ+1),drop=FALSE], 
#      Keq=Keqs[[mdlNames[j]]],Kd2KjLst=Kmapping,
      Keq=Keqs[[lmdlNames[[j]][1]]],Kd2KjLst=Kmapping,
      pparams=chunk[j,"p",drop=FALSE],indx=chunk[j,"indx"],nParams=chunk[j,"nParams"],
#      kparams=chunk[j,paste("k",g$Z,sep="")],keq=keqs[[mdlNames[j]]])
kparams=chunk[j,paste("k",g$Z,sep="")],keq=keqs[[lmdlNames[[j]][2]]])
#models
fmodels=lapply(models,fitModel)

for (j in mdlNames) 
{
  if (j %in% c("DDDD.DDLL","DDDD.DDDM","DDDD.DFFF")) lty=3 else lty=2   
  lines(fmodels[[j]]$d$ST,fmodels[[j]]$d$EY,lty=2,lwd=0.5)
}  

g <-mkg(topology, activity=TRUE,TCC=TRUE)
mdl=mkModel(g,"LA",Kdparams=gK, Kd2KjLst=Kmapping,kparams=ngk)
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
for (j in mdlNames) {
  df <- simulateData(fmodels[[j]],predict=predict,typeYP="k")$predict  
  if (j %in% c("DDDD.DDLL","DDDD.DDDM","DDDD.DFFF")) lty=3 else lty=2   
  lines(df$ST,df$EY,lty=lty,lwd=0.5)
}
legend(.15,2.5,legend=c("[ET]=.0001 uM","[ET]=0.1 uM","[ET]=0.6 uM"),pch=1:3,cex=0.7,bty="n")
#  end middle plot at ET = 0.1 uM
legend(.6,.5,legend=c("k mechanism models","K mechanism models"),lty=c(3,2),cex=0.7,bty="n")

ET1=0.6
predict <- data.frame(ET = rep(ET1,length(fineX)), ST = fineX)
df <- simulateData(mdl,predict=predict,typeYP="k")$predict  
points(df$ST,df$EY,pch=3,col=1)
lines(df$ST,df$EY,lty=1)
kk=1
for (j in mdlNames) {
  df <- simulateData(fmodels[[j]],predict=predict,typeYP="k")$predict  
  if (j %in% c("DDDD.DDLL","DDDD.DDDM","DDDD.DFFF")) lty=3 else lty=2   
  lines(df$ST,df$EY,lty=lty,lwd=0.5)# ,col=kk); kk=kk+1
}

t(fmodels[[3]]$report[c(KS,kS),"final",drop=FALSE])
t(fmodels[[4]]$report[c(KS,kS),"final",drop=FALSE])

print(ngk)
print(gK)

options(digits=5)
for (i in 1:5) print(sum(outs[[i]]$df$wgts[outs[[i]]$df$numP==3]))
for (i in 1:5) print(sum(outs[[i]]$df$wgts[outs[[i]]$df$numP==4]))
for (i in 1:5) {d=TK1[TK1$index==i,]
  if (i==3) d=d[-1,]
  print(r1<-range(d$S)); print(r1[2]/r1[1])}
options(digits=3)

for (i in 1:5) print(sum(outs[[i]]$df$wgts))

windows(width = 4, height = 4,restoreConsole = TRUE,ypos=-50)
par(mfcol=c(1,1),mar=c(4,4,0,1)+.1)
plot(df$ST,df$k,type="n",xlab="Total [dT]",log="", ylab="k (1/sec)",ylim=c(0,1.8),xlim=c(.1,1.2))
kk=1
rats=NULL
for (j in mdlNames) {
  df <- simulateData(fmodels[[j]],predict=predict,typeYP="k")$predict  
  lines(df$ST,df$EY,lty=kk,col=kk); 
  points(df$ST,df$EY,pch=kk,col=kk); 
  rats[kk]=fmodels[[kk]]$report[kS[1],"final"]#/fmodels[[kk]]$report[KS[1],"final"]
  lines(df$ST,rats[kk]*df$ST/4,lty=kk,col=kk); 
  kk=kk+1
}
legend(.15,1.8,legend=paste(mdlNames,format(rats,digits=2)),pch=1:6,lty=1:6,col=1:6,cex=0.7,bty="n")

for (kk in 1:6)
   print(t(fmodels[[kk]]$report[c(KS,kS),"final",drop=FALSE]))

 