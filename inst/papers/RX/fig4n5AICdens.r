if (.Platform$OS.type=="windows") home="/users/radivot" else home="/home/radivot"
source(paste(home,"/start.r",sep=""))  # define host="machineName" in this file
setwd(home)  # place where models and results will go
#  if (.Platform$OS.type=="windows") setwd("/users/radivot/case/active/rnr/RX") 
library(ccems)
#load("case/results/RXglobSOCKic100.RData")  
load("case/results/RX3nok.RData")  
#load("case/results/RXF3.RData")  
sAllModels[1]
length(sAllModels)


getnums<-function(model) {
  lets=c(LETTERS,letters)
  nums=c(paste(0:9),"_")
  allnms=row.names(model$sreport)
  rnms=strsplit(row.names(model$sreport)[1:(which(allnms=="p")-1)],split=NULL)
#  print(rnms)
  DF=data.frame(NULL)
  for (k in 1:length(rnms)) {
    react=rnms[[k]]
    tmpn=""
    vals=NULL
    nms=NULL
    tmp=""
    i=1
    while (i < length(react)) {
      if(react[i]%in%lets) nms=c(nms,react[i])
      i=i+1
      tmpn=""
      addOne=FALSE
      while (react[i]%in%nums) {
        if (react[i]=="_") {addOne=TRUE; i=i+2} else {
          tmpn=paste(tmpn,react[i],sep=""); 
          i=i+1}
      }
      vals=c(vals,ifelse(addOne,as.numeric(tmpn)+1,as.numeric(tmpn)))
    }
    DF=rbind(DF,vals)
    if(k==1) names(DF)<-nms
  }
#  print(DF)
  DF
}

tt=lapply(sAllModels,getnums)

#takeOut=which(sapply(lapply(tt,is.na),is.null)) # remove "nothing to fit" = IIIIIII...
#if (length(takeOut)>0) {
#  tt=tt[-takeOut]
#  sAllModels=sAllModels[-takeOut]
#}



for (i in 1:length(tt)) {
#  if (dim(tt[[i]])[1]==0) next
  sAllModels[[i]]$DF=tt[[i]]
  hypoth="noh"
  pows=NULL
  for (j in 1:dim(tt[[i]])[1]) {
#    if (tt[[i]][j,1]==1)  if (tt[[i]][j,2]==2) hypoth="h" 
    pows=c(pows,sum(tt[[i]][j,])-1)
    if (tt[[i]][j,1]==1)  if (tt[[i]][j,2]==3) hypoth="h"  #
    if (tt[[i]][j,1]>1)   if (tt[[i]][j,2]>2*tt[[i]][j,1])  hypoth="h" 
  }
  sAllModels[[i]]$maxPow=max(pows)
  sAllModels[[i]]$class=hypoth
}
#sAllModels

jmer=c(1,2,4,6)
for (i in 1:length(tt)) {
  maxRats=c(0,0,0,0) # mono, di, tet hex
  increasing=TRUE
  gt1p5="noh"  # was FALSE but to keep downstream code the same, map less than 4 filled to noh
  for (j in 1:dim(tt[[i]])[1]) {
    Ij=which(jmer==tt[[i]][j,1])
    maxRats[Ij]=max(maxRats[Ij],rt<-tt[[i]][j,2]/tt[[i]][j,1])
    if (rt>1.5) gt1p5="h"  # here having an h site is mapped to having more than 3 a-sites filled
  }
#  print(maxRats)
  maxRats=maxRats[maxRats>0]
#  print(maxRats)
  if (length(maxRats)>1)
    for (j in 1:(length(maxRats)-1)) if (maxRats[j+1]<maxRats[j]) increasing=FALSE
  sAllModels[[i]]$increasing=increasing
  sAllModels[[i]]$gt1p5=gt1p5
}

#sAllModels[1:10]

getIncreasing <- function(x) { x$increasing}
getAIC <- function(x) { x$sreport["AIC","final"]}
getSSE <- function(x) { x$sreport["SSE","final"]}
getClass <- function(x) {x$class}
getGt1p5 <- function(x) {x$gt1p5}
getMaxPow <- function(x) {x$maxPow}



getSuccess <- function(x) { x$sreport["cpu","confidenceInterval"]=="fit succeeded"}
getFailed <- function(x) { x$sreport["cpu","confidenceInterval"]=="fit failed"}
getSing <- function(x) { x$sreport["cpu","confidenceInterval"]=="hessian singular"}
getM1s <- function(x) { x$sreport["m1","confidenceInterval"]!="fixed"}
getPs <- function(x) { x$sreport["p","confidenceInterval"]!="fixed"}
getnOptParams <- function(x) {dim(x$ci)[1]}
infProps <- function(x) {sum(is.infinite(x$ci[,"upper"]))/dim(x$ci)[1]}

# this block gets rid of mass and p estimate models
Im1=sapply(sAllModels,getM1s)
Ip=sapply(sAllModels,getPs)
sum(Ip)
sum(Im1)
sum(Im1&Ip)
sum(!(Ip|Im1))
sAllModels=sAllModels[!(Ip|Im1)]


incOnly=TRUE
incOnly=FALSE

if (incOnly){
  Inc=sapply(sAllModels,getIncreasing)
  allMods=sAllModels[Inc]
  length(allMods)
} else  allMods=sAllModels  

length(allMods)
Ifail=sapply(allMods,getFailed)
#allMods[Ifail]
cat("Fit failed = ",sum(Ifail),"\n")
allMods=allMods[!Ifail]
length(allMods)
Isuc=sapply(allMods,getSuccess)
Ising=sapply(allMods,getSing)
cat("Fit succeeded = ", sum(Isuc),"   Hessian Singular = ", sum(Ising),
    "\nOut of a total of ",length(allMods),"\n")

aic=sapply(allMods,getAIC)
aicA=aic[Isuc]
fupMods=allMods[Isuc] # finite upper bound models
summary(as.factor(sapply(fupMods,infProps)))


lowMods=fupMods[aicA<180]
(s1<-summary(as.factor(sapply(lowMods,infProps))))
sum(s1[1])/sum(s1)
if (!incOnly){
  medMods=fupMods[(aicA>180)&(aicA<210)]
  (s2<-summary(as.factor(sapply(medMods,infProps))))
  sum(s2[1])/sum(s2)
}
hiMods=fupMods[aicA>210]
(s3<-summary(as.factor(sapply(hiMods,infProps))))
sum(s3[1])/sum(s3)

aicB=aic[Ising]
singMods=allMods[Ising] # finite upper bound models
summary(as.factor(sapply(singMods,infProps)))

#  *****************  plot clusters and their removal *******
if (!is.null(dev.list())) for (i in 2:max(dev.list())) dev.off(i);
windows(width = 8, height = 8,restoreConsole = TRUE) 
par(mfrow=c(2,2),mar=c(5,4,2,1)+.1)
plot(aicA,ylab="AIC",xlab="Rank",ylim=c(140,220),main="Non-Singular Hessians")
text(1,220,"A")
plot(aicB,ylab="AIC",xlab="Rank",ylim=c(140,220),main="Singular Hessians")
text(1,220,"B")
sse=sapply(allMods,getSSE)
minsse=min(sse)
minMods=allMods[sse<(2*minsse)]
#minMods[1:4]
length(minMods)
Isuc=sapply(minMods,getSuccess)
Ifail=sapply(minMods,getFailed)
Ising=sapply(minMods,getSing)
cat("Fit succeeded = ", sum(Isuc),"   Hessian Singular = ", sum(Ising),"   Fit failed = ",sum(Ifail),
    "\nOut of a total of ",length(Isuc),"\n")
aic=sapply(minMods,getAIC)
aicA=aic[Isuc]
aicB=aic[Ising]
plot(aicA,ylab="AIC",xlab="Rank",ylim=c(140,160))
text(1,160,"C")
plot(aicB,ylab="AIC",xlab="Rank",ylim=c(140,160))
text(1,160,"D")
par(mfrow=c(1,1))
#  *****************  DONE ploting AIC clusters and their removal *******



#h=sapply(minMods,getClass)
h=sapply(minMods,getGt1p5)
cat("No h site =", sum(h=="noh"),"With h site =",sum(h=="h"),"\n")
#aich=aic[h=="h"]
#aicnoh=aic[h=="noh"]
##aicnoh=aicnoh[-length(aicnoh)]
#ks.test(aich,aicnoh)
#
#hisaich=hist(aich,breaks=10,plot=FALSE)
#hisaicnoh=hist(aicnoh,breaks=10,plot=FALSE)
#
##pdf("case/results/fig4dens.pdf")
##png("case/results/fig4dens.png")
#windows(width = 4, height = 4,restoreConsole = TRUE) 
#par(mfrow=c(1,1),mar=c(5,4,2,1)+.1)
#plot(hisaicnoh$mids,hisaicnoh$density,xlab="AIC",ylab="density",type="l",lty=1,xlim=range(c(hisaicnoh$mids,hisaich$mids)))
#lines(hisaich$mids,hisaich$density,lty=2)
##legend(-33,.32,legend=c(paste("no h site (",sum(hisaicnoh$counts),")",sep=""),
##       paste("h site (",sum(hisaich$counts),")",sep="")),lty=1:2)
#legend(141,ifelse(incOnly,.38,.28),legend=c(paste("no h site (",sumnoh<-sum(hisaicnoh$counts),")",sep=""),
#        paste("h site (",sumh<-sum(hisaich$counts),")",sep="")),lty=1:2,bty="n",cex=0.9)
##dev.off()
Isuc=sapply(minMods,getSuccess)
Ising=sapply(minMods,getSing)
nOpts=sapply(minMods,getnOptParams)
windows(width = 8, height = 8,restoreConsole = TRUE) 
par(mfrow=c(2,2),mar=c(5,4,2,1)+.1)

aich=aic[h=="h"]
aicnoh=aic[h=="noh"]
ks.test(aich,aicnoh)
hisaich=hist(aich,breaks=10,plot=FALSE)
hisaicnoh=hist(aicnoh,breaks=10,plot=FALSE)
plot(hisaicnoh$mids,hisaicnoh$density,xlab="AIC",ylab="density",type="l",lty=1,xlim=range(c(hisaicnoh$mids,hisaich$mids)))
lines(hisaich$mids,hisaich$density,lty=2)
legend(148.5,ifelse(incOnly,.38,.28),legend=c(paste("no h site (",sumnoh<-sum(hisaicnoh$counts),")",sep=""),
        paste("h site (",sumh<-sum(hisaich$counts),")",sep="")),lty=1:2,bty="n",cex=0.9)
mtext("A",line=-1.1,side=3,adj=0.02)
(sumnoh30<-sum(aicnoh<=aic[30]))
sumh30<-sum(aich<=aic[30])
prop.test(c(sumnoh30,sumnoh),c(30,sumnoh+sumh))

aich=aic[(h=="h")&(nOpts<=2)]
aicnoh=aic[(h=="noh")&(nOpts<=2)]
ks.test(aich,aicnoh)
hisaich=hist(aich,breaks=10,plot=FALSE)
hisaicnoh=hist(aicnoh,breaks=10,plot=FALSE)
plot(hisaicnoh$mids,hisaicnoh$density,xlab="AIC",ylab="density",type="l",lty=1,xlim=range(c(hisaicnoh$mids,hisaich$mids)))
lines(hisaich$mids,hisaich$density,lty=2)
legend(146.5,ifelse(incOnly,.38,.34),legend=c(paste("no h site (",sumnoh<-sum(hisaicnoh$counts),")",sep=""),
        paste("h site (",sumh<-sum(hisaich$counts),")",sep="")),lty=1:2,bty="n",cex=0.9)
mtext("B",line=-1.1,side=3,adj=0.02)

aich=aic[(h=="h")&((nOpts==3)&Isuc)]
aicnoh=aic[(h=="noh")&((nOpts==3)&Isuc)]
ks.test(aich,aicnoh)
hisaich=hist(aich,breaks=10,plot=FALSE)
hisaicnoh=hist(aicnoh,breaks=10,plot=FALSE)
plot(hisaicnoh$mids,hisaicnoh$density,xlab="AIC",ylab="density",type="l",lty=1,xlim=range(c(hisaicnoh$mids,hisaich$mids)))
lines(hisaich$mids,hisaich$density,lty=2)
legend(149,ifelse(incOnly,.43,.3),legend=c(paste("no h site (",sumnoh<-sum(hisaicnoh$counts),")",sep=""),
        paste("h site (",sumh<-sum(hisaich$counts),")",sep="")),lty=1:2,bty="n",cex=0.9)
mtext("C",line=-1.1,side=3,adj=0.02)

aich=aic[(h=="h")&Ising]
aicnoh=aic[(h=="noh")&Ising]
ks.test(aich,aicnoh)
hisaich=hist(aich,breaks=10,plot=FALSE)
hisaicnoh=hist(aicnoh,breaks=10,plot=FALSE)
plot(hisaicnoh$mids,hisaicnoh$density,xlab="AIC",ylab="density",type="l",lty=1,xlim=range(c(hisaicnoh$mids,hisaich$mids)))
lines(hisaich$mids,hisaich$density,lty=2)
legend(150,ifelse(incOnly,.38,.5),legend=c(paste("no h site (",sumnoh<-sum(hisaicnoh$counts),")",sep=""),
        paste("h site (",sumh<-sum(hisaich$counts),")",sep="")),lty=1:2,bty="n",cex=0.9)
mtext("D",line=-1.1,side=3,adj=0.02)
