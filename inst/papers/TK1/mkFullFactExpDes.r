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

#SEQS=c(0.01,0.02,0.05, .1, .2, .5, 1, 2, 5,10 )
#SEQE=c(0.001,0.01,0.1, 1)
SEQS=seq(.1,1.2,.1)
SEQE=c(0.001,0.1, 0.6)
modS="ma3"
Kt = c(0.73,0.51,0.49,0.49) 
kt = c(1.97,3.16,3.90,3.89)
Kmapping = mkKd2Kj(g)
names(Kt)<-g$KS
names(kt)<-g$kS
reps=2  # within an experiment
runs=100  # to get distribution of many such experiments
mdl = mkModel(g,"TK1",Kdparams=Kt, Kd2KjLst=Kmapping,kparams=kt)
fineX=rep(SEQS,each=reps) 
predict <- data.frame(ET = rep(SEQE,each=length(fineX)), ST = fineX)
d <- simulateData(mdl,predict=predict,typeYP="k")$predict  
names(d)[3]<-"k"
names(d)[1:2]= c("ET","ST")
outs=NULL
for (i in 1:runs) {
  df=transform(d,k=rnorm(k,k,0.06*k))
  tops=ems(df,g,maxTotalPs=3,minTotalPs=3,cpusPerHost=cpusPerHost,ptype="SOCK",
      doSpurs=FALSE,topN=64,
      kIC=rep(4,4),KIC=rep(.6,4))
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
save(outs,file="case/results/FFOutsR2") 



