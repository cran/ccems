rm(list=ls(all=TRUE))  # clean up left overs from previous run, though shouldn't need to.
# This code is meant to be run on a ROCKS 5.1 cluster
if (.Platform$OS.type=="windows") home="/users/radivot" else home="/home/radivot"
source(paste(home,"/start.r",sep=""))  # define host="machineName" in this file
setwd(home)  # place where models and results will go
if (1) library(ccems) else { # if 0 source in the package 
  pkgNms=dir(paste(home,"/case/active/ccems/ccems/R",sep=""),full.names=TRUE)
  for (i in pkgNms) source(i)  # source in all of the R files of the package
  RNR=read.table(file=paste(home,"/case/active/rnr/datasets/RNR.txt",sep=""),header=TRUE)
}


topology <- list(  # sah run on tk2 24 cores
    heads=c("R1X0","R2X0","R4X0","R6X0"), # no requirement that s-sites are filled in oligomers
    sites=list(                    
        s=list(                                                              # s-site       thread #
            m=c("R1X1"),                                                 # monomer          1
            d=c("R2X1","R2X2"),                                          # dimer            2
            t=c("R4X1","R4X2","R4X3","R4X4"),                            # tetramer         3
            h=c("R6X1","R6X2","R6X3","R6X4", "R6X5", "R6X6")             # hexamer          4
        ),
        a=list(                                                              # a-site       thread #
            m=c("R1X2"),                                                 # monomer          5
            d=c("R2X3","R2X4"),                                          # dimer            6
            t=c("R4X5","R4X6","R4X7","R4X8"),                            # tetramer         7
            h=c("R6X7","R6X8","R6X9","R6X10", "R6X11", "R6X12")          # hexamer          8
        ),
        h=list( ## tails of a-site threads are heads of h-site threads       # h-site
            m=c("R1X3"),                                                 # monomer          9
            d=c("R2X5", "R2X6"),                                         # dimer            10
            t=c("R4X9", "R4X10","R4X11", "R4X12"),                       # tetramer         11
            h=c("R6X13", "R6X14", "R6X15","R6X16", "R6X17", "R6X18")     # hexamer          12
        )
    )
)

if(0) {  #sa option bas
    
topology <- list(
    heads=c("R1X0","R2X0","R4X0","R6X0"), # no requirement that s-sites are filled in oligomers
    sites=list(                    
        s=list(                                                              # s-site       thread #
            m=c("R1X1"),                                                 # monomer          1
            d=c("R2X1","R2X2"),                                          # dimer            2
            t=c("R4X1","R4X2","R4X3","R4X4"),                            # tetramer         3
            h=c("R6X1","R6X2","R6X3","R6X4", "R6X5", "R6X6")             # hexamer          4
        ),
        a=list(                                                              # a-site       thread #
            m=c("R1X2"),                                                 # monomer          5
            d=c("R2X3","R2X4"),                                          # dimer            6
            t=c("R4X5","R4X6","R4X7","R4X8"),                            # tetramer         7
            h=c("R6X7","R6X8","R6X9","R6X10", "R6X11", "R6X12")          # hexamer          8
        )
    )
)
}

if(1){  #ah option from ccems help
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
}

g <- mkg(topology)
dd=subset(RNR,(year==2002)&(fg==1)&(X>0),select=c(R,X,m,year))
names(dd)[1:2]=paste(strsplit(g$id,split="")[[1]],"T",sep="") # e.g. c("RT","XT")
cpusPerHost=c("localhost" = 4,"compute-0-0"=4,"compute-0-1"=4,"compute-0-2"=4,"compute-0-3"=4,"compute-0-4"=4) # for tk2
if ((host=="tk1")|(host=="stdn")) cpusPerHost=cpusPerHost[1]
if ((host=="dck")|(host=="rnrClust")) cpusPerHost=cpusPerHost[1:4]
tops=ems(dd,g,cpusPerHost=cpusPerHost,maxTotalPs=3,ptype="SOCK",topN=5,IC=100) 

