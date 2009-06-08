rm(list=ls(all=TRUE))  # clean up left overs from previous run, though shouldn't need to.
if (.Platform$OS.type=="windows") home="/users/radivot" else home="/home/radivot"
setwd(home)  # place where models and results will go
source(paste(home,"/start.r",sep=""))  # define host="machineName" in this file
if (0) library(ccems) else { # if 0 source in the package to save install time 
  pkgNms=dir(paste(home,"/case/active/ccems/ccems/R",sep=""),full.names=TRUE)
  for (i in pkgNms) source(i)  # source in all of the R files of the package
  load(paste(home,"/case/active/papers/TK1/data/TK1.rda",sep=""))
}
TK1

# This block is straight from the paper
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

getKk <- function(x) {t(x$report[c(paste("E1S",0:3,"_S",sep=""),
              paste("kE1S",1:4,sep="")),"final",drop=FALSE])}
getAIC <- function(x) { x$report["AIC","final"]}
getSSE <- function(x) { x$report["SSE","final"]}
getNumP <- function(x) { x$nOptParams}
getNames <- function(x) { x$mid}
outs=NULL

cpusPerHost=c("localhost" = 4,"compute-0-0"=4,"compute-0-1"=4,"compute-0-2"=4,"compute-0-3"=4,"compute-0-4"=4) # for tk2
if ((host=="tk1")|(host=="stdn")) cpusPerHost=cpusPerHost[1]
if ((host=="dck")|(host=="rnrClust")) cpusPerHost=cpusPerHost[1:4]
if (host=="ATP") cpusPerHost=c("localhost" = 1)
if (host=="dck") cpusPerHost=c("localhost" = 4)

kill3=FALSE
kill3=TRUE

kill1=FALSE
kill1=TRUE
killL1=TRUE
killL1=FALSE

gibbs=FALSE

par(mfcol=c(2,5),mar=c(4,4,2,1)+.1)
for (i in 1:5) {
  d=subset(TK1,index==i,select=c(E,S,k,frstAut,year))
#  d=subset(TK1,select=c(E,S,k,frstAut,year),subset=(index==i)&(S>.09)&(S<2.1))
  if ((i==1)&(kill3)) d=d[-(14:16),]
  if ((i==3)&(kill1)) d=d[-1,]
  if ((i==4)&(killL1)) d=d[-14,]
  if ((i==2)&(killL1)) d=d[-15,]
  #  d=transform(d,E=E/4) # now built into the frp function
   plot(d$S,d$k,xlab="Total [dT]",log="xy", ylab="k (1/sec)", 
   main=paste(d[1,"frstAut"],d[1,"year"]))
   names(d)[1:2]= c("ET","ST")
   tops=ems(d,g,maxTotalPs=8,minTotalPs=3,cpusPerHost=cpusPerHost,ptype="SOCK",
    doSpurs=FALSE,topN=64,fullGrid=TRUE,KIC=0.6,kIC=ifelse(i==5,.5,4))# takes ~15 sec for each dataset
    lgx=log(d$ST)
    upr=range(lgx)[2]
    lwr=range(lgx)[1]
    del=(upr-lwr)/50
    fineX=exp(seq(lwr,upr,by=del))
    predict <- data.frame(ET = rep(d$ET[1],length(fineX)), ST = fineX)
    df <- simulateData(tops[[1]],predict=predict,typeYP="k")$predict  
    lines(df$ST,df$EY) 
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
    if (!gibbs)    ma=wgts%*%M # average in straight params
if (gibbs)    ma=exp(wgts%*%log(M)) # average in space of gibbs free energy changes
    dataID=paste(d[1,"frstAut"],d[1,"year"],sep="")
    outs[[dataID]]$df=df
    outs[[dataID]]$ma=ma
    plot(tops[[1]]$d$EY,tops[[1]]$res,xlab="Fitted Value",
        ylab="Residual",main=tops[[1]]$mid)
  }
#  outs=outs[-1] # remove leading NULL
  print(outs)   # compare model averages across datasets
  par(mfrow=c(1,1))
if ((!kill1)&(!kill3)&(!killL1)) save(outs,file="case/results/outs00")
if ((kill1)&(kill3)&(!killL1)) save(outs,file=paste("case/results/outs31",ifelse(gibbs,"G",""),sep=""))
if ((kill1)&(!kill3)&(!killL1)) save(outs,file=paste("case/results/outs01",ifelse(gibbs,"G",""),sep=""))
if (killL1) save(outs,file="case/results/outsL1")
#load("case/results/outs")
  
