rm(list=ls(all=TRUE))  # clean up left overs from previous run, though shouldn't need to.
# This code is meant to be run on a ROCKS 5.1 cluster
if (.Platform$OS.type=="windows") home="/users/radivot" else home="/home/radivot"
source(paste(home,"/start.r",sep=""))  # define host="machineName" in this file
setwd(home)  # place where models and results will go
if (0) library(ccems) else { # if 0 source in the package 
  pkgNms=dir(paste(home,"/case/active/ccems/ccems/R",sep=""),full.names=TRUE)
  for (i in pkgNms) source(i)  # source in all of the R files of the package
  RNR=read.table(file=paste(home,"/case/active/rnr/datasets/RNR.txt",sep=""),header=TRUE)
}
d1=subset(RNR,(year==2001)&(fg==1)&(G==0)&(t>0),select=c(R,t,m,year))
#d2=subset(RNR,year==2006,select=c(R,t,m,year)) 
#d=rbind(d1,d2) 
d=d1 
names(d)[1:2]=c("RT","tT")

topology=list(  
    heads=c("R1t0","R2t0"),
    sites=list(
        s=list(
            m=c("R1t1"),
            d=c("R2t1","R2t2")
        )
    )
) 
g <- mkg(topology)

#chunk=mkSpurs(g,maxTotalPs=3,m1=90,pRows=T,doTights=T)$chunk
#chunk
#
#models=NULL
#mdlNames=rownames(chunk)
#lmdlNames=strsplit(mdlNames,split=".",fixed=TRUE)
#names(lmdlNames)<-mdlNames
#for (j in mdlNames) 
#  models[[j]]=mkModel(g,j,d,Kjparams=chunk[j,g$Z], pparams=chunk[j,c("p","m1"),drop=FALSE],indx=chunk[j,"indx"],
#      nParams=chunk[j,"nParams"])
#
#fmdl<-fitModel(models[["IIIJpm"]])
#

#ng=mkGrids(g,maxTotalPs=2,m1=90)
#ng$chunk
#
bigNms=c("localhost")
rackLen=c(12,12)
for (rack in 0:0) for (node in 0:rackLen[rack+1])  bigNms=c(bigNms,paste("compute",rack,node,sep="-"))
big=rep(4,length(bigNms))
names(big)<-bigNms
big
cpusPerHost=c("localhost" = 4,"compute-0-0"=4,"compute-0-1"=4,"compute-0-2"=4,
              "compute-0-3"=4,"compute-0-4"=4,"compute-0-5"=4,"compute-0-6"=4) # for tk2
if ((host=="tk1")|(host=="stdn")) cpusPerHost=cpusPerHost[1]
if ((host=="dck")|(host=="rnrClust")) cpusPerHost=cpusPerHost[1:2]
if (host=="ATP") cpusPerHost=c("localhost" = 1)
#if (host=="dck") cpusPerHost=c("localhost" = 3)
#if (host=="dck") cpusPerHost=big
tops=ems(d,g,pRows=T,doTights=T,m1=90,spurChunkSize=10,cpusPerHost=cpusPerHost,maxTotalPs=3,ptype="SOCK",topN=50,transform="none")
#tops=ems(d,g,pRows=T,doTights=T,m1=90,maxTotalPs=3,transform="none")
#tops=ems(d,g,m1=90,maxTotalPs=3,pRows=T,doTights=T,transform="none")
save(tops,file="case/results/Rt3tops")


