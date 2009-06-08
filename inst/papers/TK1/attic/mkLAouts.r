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

getKk <- function(x) {t(x$report[c(paste("E1S",0:3,"_S",sep=""),
              paste("kE1S",1:4,sep="")),"final",drop=FALSE])}
getAIC <- function(x) { x$report["AIC","final"]}
getSSE <- function(x) { x$report["SSE","final"]}
getNumP <- function(x) { x$nOptParams}
getNames <- function(x) { x$mid}

KS=c("E1S0_S","E1S1_S","E1S2_S","E1S3_S")
kS=c("kE1S1","kE1S2","kE1S3","kE1S4")
options(digits=3)
REP=TRUE # ten times as many piled up our spread between
REP=FALSE 
if (REP) fineX=rep(seq(.1,1.2,by=.1),each=10) else fineX=seq(.1,1.2,by=.0092)
predict <- data.frame(ET = rep(.0001,length(fineX)), ST = fineX)
#predict

FULL=TRUE # make full and fit full space to it
FULL=FALSE  
LA5=TRUE# make 5 param Lit Ave (LA5)
LA5=FALSE
LA4=TRUE # make 4 param Lit Ave (LA5)
if (LA5) {Kic=c( 0.9,0.7,0.7,0.55);kic=c(3.9,4.2,4.2,4.2)} # 5 paramters
if (FULL) {Kic=c( 2,1,0.5,0.1);kic=c(1,2,3,4)} 
if (LA4){ Kic=c( 0.9,0.7,0.7,0.55);kic=c(4.1,4.1,4.1,4.1) } # 4 parameters
Kmapping = mkKd2Kj(g)
names(Kic)<-KS
names(kic)<-kS
mdl = mkModel(g,"TK1",Kdparams=Kic, Kd2KjLst=Kmapping,kparams=kic)
d <- simulateData(mdl,predict=predict,typeYP="k")$predict  
names(d)[3]<-"k"
names(d)[1:2]= c("ET","ST")
d=transform(d,k=rnorm(k,k,0.05*k))
#d

#outs=list(NULL)  # leaves a leading null element
outs=NULL
cpusPerHost=c("localhost" = 4,"compute-0-0"=4,"compute-0-1"=4,"compute-0-2"=4,"compute-0-3"=4,"compute-0-4"=4) # for tk2
if ((host=="tk1")|(host=="stdn")) cpusPerHost=cpusPerHost[1]
if ((host=="dck")|(host=="rnrClust")) cpusPerHost=cpusPerHost[1:4]
if (host=="ATP") cpusPerHost=c("localhost" = 1)
for (i in 1:1) {
  tops=ems(d,g,maxTotalPs=8,minTotalPs=3,cpusPerHost=cpusPerHost,ptype="SOCK",
      doSpurs=FALSE,topN=64,fullGrid=TRUE,
#      kIC=(kic<-exp(rnorm(1,log(4),.1*log(4)))),KIC=(Kic<-exp(rnorm(1,log(.6),.1*log(.6)))))
      kIC=(kic<-rep(4,4)),KIC=(Kic<-rep(.6,4)))
#  kIC=(kic<-rnorm(1,4,.2)),KIC=(Kic<-rnorm(1,.6,.03)))
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
  outs[[dataID]]$ics=c(KIC=Kic,kic=kic)
}
#outs=outs[-1] # remove leading NULL if initialiszed by list(NULL) rather than NULL
#print(outs)   # compare model averages across datasets
if (FULL)  if (REP) save(outs,file="case/results/repFullouts") else save(outs,file="case/results/fullouts")
if (LA5)  if (REP) save(outs,file="case/results/repLA5outs") else save(outs,file="case/results/LA5outs")
if (LA4)  if (REP) save(outs,file="case/results/repLA4outs") else save(outs,file="case/results/LA4outs")




