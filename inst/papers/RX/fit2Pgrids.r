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
print(host)
  topology <- list(
      heads=c("R1X0","R2X2","R4X4","R6X6"), # s-sites are already filled only in (j>1)-mer head nodes 
      sites=list(                    
          a=list(                                                              # a-site       thread #
              m=c("R1X1"),                                                 # monomer          1
              d=c("R2X3","R2X4"),                                          # dimer            2
              t=c("R4X5","R4X6","R4X7","R4X8"),                            # tetramer         3
              h=c("R6X7","R6X8","R6X9","R6X10", "R6X11", "R6X12")          # hexamer          4
          ),
          h=list( ## tails of a-site threads are heads of h-site threads       # h-site
              m=c("R1X2"),                                                 # monomer          5
              d=c("R2X5", "R2X6"),                                         # dimer            6
              t=c("R4X9", "R4X10","R4X11", "R4X12"),                       # tetramer         7
              h=c("R6X13", "R6X14", "R6X15","R6X16", "R6X17", "R6X18")     # hexamer          8
          )
      )
  )

g <- mkg(topology)
#dim(mkSpurs(g,maxTotalPs=4,minTotalPs=1,m1=90,p=1,forceM1=T,forceP=T,batchSize=6000)$chunk)
mkSpurs(g,maxTotalPs=3,minTotalPs=1,m1=90,batchSize=6000,p=1)$chunk[1:20,]

mkGrids(g,maxTotalPs=4,minTotalPs=1,m1=90,p=1,forceM1=T,forceP=T)$chunk


ng=mkGrids(g,maxTotalPs=4,minTotalPs=1,m1=90,p=1,forceM1=T,forceP=T)
chunk=ng$chunk
chunk
Kmapping=mkKd2Kj(g)
Keqs=ng$Keqs
keqs=ng$keqs
mdlNames=rownames(chunk)

lmdlNames=strsplit(mdlNames,split=".",fixed=TRUE)
names(lmdlNames)<-mdlNames

dd=subset(RNR,(year==2002)&(fg==1)&(X>0),select=c(R,X,m,year))
names(dd)[1:2]=paste(strsplit(g$id,split="")[[1]],"T",sep="") # e.g. c("RT","XT")
#dd=dd[-12,]  # outlier that raises things to cause the downturn
#dd=dd[-14,]  # causes downturn    
#dd=dd[-c(1,12),]  # kill speculative bias at 90 and kill high at 540
#dd=dd[-1,]  # remove bogus first data point since 90 is built into model anyway, or it should be estimated freely
dd

for (j in mdlNames) 
{
  mdl=mkModel(g,j,dd,Kdparams=chunk[j,2:(g$nZ+1),drop=FALSE], 
        Keq=Keqs[[j]],                Kd2KjLst=Kmapping,
        pparams=chunk[j,c("p","m1"),drop=FALSE],indx=chunk[j,"indx"],nParams=chunk[j,"nParams"],transform="none")
  print(mdl)
#  fmdl=fitModel(mdl)
}

