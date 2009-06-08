rm(list=ls(all=TRUE))  # clean up left overs from previous run, though shouldn't need to.
if (.Platform$OS.type=="windows") home="/users/radivot" else home="/home/radivot"
setwd(home)  # place where models and results will go
source(paste(home,"/start.r",sep=""))  # define host="machineName" in this file
if (0) library(ccems) else { # if 0 source in the package to save install time 
  pkgNms=dir(paste(home,"/case/active/ccems/ccems/R",sep=""),full.names=TRUE)
  for (i in pkgNms) source(i)  # source in all of the R files of the package
  load(paste(home,"/case/active/papers/TK1/data/TK1.rda",sep=""))
}

cpusPerHost=c("localhost" = 4,"compute-0-0"=4,"compute-0-1"=4,"compute-0-2"=4,"compute-0-3"=4,"compute-0-4"=4) # for tk2
if ((host=="tk1")|(host=="stdn")) cpusPerHost=cpusPerHost[1]
if ((host=="dck")|(host=="rnrClust")) cpusPerHost=cpusPerHost[1:2]
if (host=="ATP") cpusPerHost=c("localhost" = 1)
if (host=="tk2") cpusPerHost=cpusPerHost[1:2]
if (host=="dNTP") cpusPerHost=c("localhost" = 3)

topology <- list(  
    heads=c("E1S0"), # E1S0 = substrate free E
    sites=list(                    
        c=list(    # c for catalytic site  
            t=c("E1S1","E1S2","E1S3","E1S4")   
        ) # t for tetramer 
    )
)   # in transform below, TK1 is 25kDa => 25mg/umole
g <-mkg(topology, activity=TRUE,TCC=TRUE)

getKk <- function(x) {t(x$report[c(paste("E1S",0:3,"_S",sep=""),
              paste("kE1S",1:4,sep="")),"final",drop=FALSE])}
getAIC <- function(x) { x$report["AIC","final"]}
getSSE <- function(x) { x$report["SSE","final"]}
getNumP <- function(x) { x$nOptParams}
getNames <- function(x) { x$mid}
options(digits=3)

fineX=seq(.1,1.2,.1)
SEQE=c(0.001,0.1, 0.6)
modS="ma3"
Kt = c(0.73,0.51,0.49,0.49) 
kt = c(1.97,3.16,3.90,3.89)
names(Kt)<-g$KS
names(kt)<-g$kS
Kmapping = mkKd2Kj(g)
runs=1  # to get distribution of many such experiments
mdl = mkModel(g,"TK1",Kdparams=Kt, Kd2KjLst=Kmapping,kparams=kt)
initPredict <- data.frame(ET = rep(.001,each=length(fineX)), ST = fineX)
grid <- data.frame(ET = rep(SEQE,each=length(fineX)), ST = fineX)
d <- simulateData(mdl,predict=initPredict,typeYP="k")$predict  
names(d)[3]<-"k"
names(d)[1:2]= c("ET","ST")
d=transform(d,k=rnorm(k,k,0.06*k))
d
curd=d
alld=d
ntop=6
adaptPnts=NULL
outs=NULL
for (i in 1:runs) {
  tops=ems(alld,g,maxTotalPs=3,minTotalPs=3,cpusPerHost=cpusPerHost,ptype="SOCK",
      doSpurs=FALSE,topN=6,kIC=rep(4,4),KIC=rep(.6,4))  # prime the pump with one full batch
  for (iBatch in 1:1) {  # run through 3 batches of 12 points
    expCond=grid  # to set up "row names"
    for (itop in 1:ntop)   # now tack on 6 predictions to each grid point/row
      expCond=cbind(expCond,simulateData(tops[[itop]],predict=grid,typeYP="k")$predict$EY)
    names(expCond)[3:(3-1+ntop)]=paste("M",1:ntop,sep="")
    mat=as.matrix(expCond[,3:(3-1+ntop)])
    sig=apply(mat,1,sd)
    denom=apply(mat,1,mean)
    expCond=cbind(expCond,cv=sig/denom)
    I=order(expCond$cv,decreasing=TRUE)
    sdf=expCond[I,]
    best=sdf[1:12,]   # best next 12 data points
    par(mfrow=c(1,1),mar=c(5.1,4.1,1.5,0.5))
    plot(best$ST,best$ET,type="p",
        ylab="[ET]",xlab="[ST]",log="xy",ylim=c(.001,1),xlim=c(.1,2))
    adaptPnts=rbind(adaptPnts,best[,1:2])
    curd <- simulateData(mdl,predict=best[,1:2],typeYP="k")$predict  
    names(curd)[3]<-"k"
    names(curd)[1:2]= c("ET","ST")
    curd=transform(curd,k=rnorm(k,k,0.06*k))
    alld=rbind(alld,curd)
    tops=ems(alld,g,maxTotalPs=3,minTotalPs=3,cpusPerHost=cpusPerHost,ptype="SOCK",
        doSpurs=FALSE,topN=64,kIC=rep(4,4),KIC=rep(.6,4))  # prime the pump with one full batch
  } 
  Kk=lapply(tops,getKk)
  nms=sapply(tops,getNames)
  rowList=data.frame(NULL)
  for (j in 1:length(nms))  rowList=rbind(rowList,Kk[[j]])
  rownames(rowList)<-nms
  aic=sapply(tops,getAIC)
  sse=sapply(tops,getSSE)
  numP=sapply(tops,getNumP)
  eDelAIC=exp(-(aic-min(aic)))
  wgts=eDelAIC/sum(eDelAIC)
  print(sum(wgts))
  df=data.frame(numP,sse,aic,wgts,rowList)
  M=as.matrix(rowList)
  ma=exp(wgts%*%log(M)) # average in space of gibbs free energy changes
  dataID=paste("run",i,sep="")
  print(dataID)
  outs[[dataID]]$df=df
  outs[[dataID]]$ma=ma
}
outs$trueValues=c(Kt,kt)
outs$trueName=modS
save(outs,file="case/results/adaptOuts") 



