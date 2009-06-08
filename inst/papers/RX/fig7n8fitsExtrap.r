rm(list=ls(all=TRUE))  
if (.Platform$OS.type=="windows") home="/users/radivot" else home="/home/radivot"
setwd(home)  # place where models and results will go
if (0) library(ccems) else { # if 0 source in the package 
  pkgNms=dir(paste(home,"/case/active/ccems/ccems/R",sep=""),full.names=TRUE)
  for (i in pkgNms) source(i)  # source in all of the R files of the package
  RNR=read.table(file=paste(home,"/case/active/rnr/datasets/RNR.txt",sep=""),header=TRUE)
}



# next block makes sure that the right C code is in place from simulateData below
topology <- list(
    heads=c("R1X0","R2X2","R4X4","R6X6"), 
    sites=list(                # s-sites are already filled only in (j>1)-mers 
        a=list(  #a-site                                                    thread
            m=c("R1X1"),                                            # monomer   1
            d=c("R2X3","R2X4"),                                     # dimer     2
            t=c("R4X5","R4X6","R4X7","R4X8"),                       # tetramer  3
            h=c("R6X7","R6X8","R6X9","R6X10", "R6X11", "R6X12")     # hexamer   4
        ), # tails of a-site threads are heads of h-site threads
        h=list(   # h-site
            m=c("R1X2"),                                            # monomer   5
            d=c("R2X5", "R2X6"),                                    # dimer     6
            t=c("R4X9", "R4X10","R4X11", "R4X12"),                  # tetramer  7
            h=c("R6X13", "R6X14", "R6X15","R6X16", "R6X17", "R6X18")# hexamer   8
        )
    )
)
g=mkg(topology,TCC=TRUE) 
#g
load("case/results/RX1top100.RData") 
mkHTML(globalTopN)  # dump out to make sure it's what you think it is
#globalTopN[[1]]
d=globalTopN[[1]]$d # grab data from first model

# read in the top models
#load("case/results/RX2top100k112.RData") 
load("case/results/RX3top100nok.RData") 
#load("case/results/RX2top100.RData") 
#load("case/results/RX2top100k12.RData") 
mkHTML(globalTopN)  # dump out to make sure it's what you think it is
#globalTopN[[1]]
d=globalTopN[[1]]$d # grab data from first model
if (d$XT[1]<10) d=d[-1,]  # kill first point strictly to plot logs on x axis
#d 

lgx=log(d$XT)
upr=range(lgx)[2]
lwr=range(lgx)[1]
del=(upr-lwr)/50
fineX=exp(seq(lwr,upr,by=del))
newPnts <- data.frame(RT = rep(d$RT[1],length(fineX)), XT = fineX)
Top=3
imds=1:Top
MDS=NULL
for (i in 1:Top) 
  MDS[[i]] <- simulateData(globalTopN[[imds[i]]],predict=newPnts,typeYP="m")  
#MDS
if (!is.null(dev.list())) for (i in 2:max(dev.list())) dev.off(i);
windows(width = 8, height = 4,restoreConsole = TRUE) 
par(mfcol=c(1,2),mar=c(5,4,1,1)+.1)
plot(d$XT,d$m,type="p",pch=1, xlab="[ATP] (uM)", ylab="Mass (kDa)",log="x")
text(50,530,"A")
for (i in 1:Top)
  lines(MDS[[i]]$predict$XT,MDS[[i]]$predict$EY,type="l",lty=i)
getID <- function(x) { paste(row.names(x$params)[x$params$opt],collapse=".")}
mnms=sapply(globalTopN,getID)
mnms
names(globalTopN)<-mnms
#legend(1000,200,c("Best", "Second Best", "Third Best"),lty=1:3)
legend(350,270,mnms[imds],lty=1:5,cex=0.9,bty="n")
plot(globalTopN[[1]]$d$EY,globalTopN[[1]]$res,type="p",pch=1, xlab="Fitted Value", ylab="Residual")
text(100,50,"B")

#detach(globalTopN[[1]]$d)
#dev.copy2pMDS(file="results/fig5fit.pMDS")
lgx=log(d$XT)
upr=log(1e6)
lwr=range(lgx)[1]
del=(upr-lwr)/50
fineX=exp(seq(lwr,upr,by=del))
newPnts <- data.frame(RT = rep(d$RT[1],length(fineX)), XT = fineX)
#MDS[[2]]$d
imds=c(2,11,15); Top=length(imds)
MDS=NULL
for (i in 1:Top) 
  MDS[[i]] <- simulateData(globalTopN[[imds[i]]],predict=newPnts,typeYP="m")  

windows(width = 4, height = 4,restoreConsole = TRUE,ypos=75) 
par(mfcol=c(1,1),mar=c(5,4,1,1)+.1)
plot(MDS[[1]]$predict$XT,MDS[[1]]$predict$EY,type="l",lty=1, xlab="[ATP] (uM)", ylab="Mass (kDa)",log="x",ylim=c(0,540))
for (i in 2:Top)
  lines(MDS[[i]]$predict$XT,MDS[[i]]$predict$EY,type="l",lty=i)
legend(800,170,mnms[imds],lty=1:5,cex=0.9,bty="n")
MDS[[2]]

