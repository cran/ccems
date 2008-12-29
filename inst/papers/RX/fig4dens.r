if (.Platform$OS.type=="windows") home="/users/radivot" else home="/home/radivot"
source(paste(home,"/start.r",sep=""))  # define host="machineName" in this file
setwd(home)  # place where models and results will go
#  if (.Platform$OS.type=="windows") setwd("/users/radivot/case/active/rnr/RX") 
library(ccems)
load("case/results/RXglobSOCKic100.RData")  

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

takeOut=which(sapply(lapply(tt,is.na),is.null)) # remove "nothing to fit" = IIIIIII...
if (length(takeOut)>0) {
  tt=tt[-takeOut]
  sAllModels=sAllModels[-takeOut]
}

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

getAIC <- function(x) { x$sreport["AIC","final"]}
getSSE <- function(x) { x$sreport["SSE","final"]}
getClass <- function(x) {x$class}
getMaxPow <- function(x) {x$maxPow}


aic=sapply(sAllModels,getAIC)
aicFt=aic[aic<0]
cat("Fitted =", length(aicFt),", out of a total of ",length(aic),"\n")
aic

#aic[1:10]
#sse=sapply(sAllModels,getSSE)
#sse[1:10]
#indx=sort.list(sse)
#sseAllModels=sAllModels[indx]
#sse=sapply(sseAllModels,getSSE)
#sse[1:10]

h=sapply(sAllModels,getClass)
cat("No h site =", sum(h=="noh"),"With h site =",sum(h=="h"),"\n")

aic=sapply(sAllModels,getAIC)
sse=sapply(sAllModels,getSSE)
minsse=min(sse)

aich=aic[(h=="h")&(sse<(2*minsse))]
aicnoh=aic[(h=="noh")&(sse<(2*minsse))]
ks.test(aich,aicnoh)

hisaich=hist(aich,breaks=10)
hisaicnoh=hist(aicnoh,breaks=10)

#pdf("case/results/fig4dens.pdf")
png("case/results/fig4dens.png")
plot(hisaicnoh$mids,hisaicnoh$density,xlab="AIC",ylab="density",type="l",lty=1)
lines(hisaich$mids,hisaich$density,lty=2)
legend(-33,.32,legend=c(paste("no h site (",sum(hisaicnoh$counts),")",sep=""),
       paste("h site (",sum(hisaich$counts),")",sep="")),lty=1:2)
dev.off()

