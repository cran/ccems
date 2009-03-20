rm(list=ls(all=TRUE))  # clean up left overs from previous run, though shouldn't need to.
if (.Platform$OS.type=="windows") home="/users/radivot" else home="/home/radivot"
setwd(home)  # place where models and results will go
if (1) library(ccems) else { # if 0 source in the package to save install time 
  pkgNms=dir(paste(home,"/case/active/ccems/ccems/R",sep=""),full.names=TRUE)
  for (i in pkgNms) source(i)  # source in all of the R files of the package
  TK1=read.table(file=paste(home,"/case/active/papers/TK1/data/TK1.txt",sep=""),header=TRUE)
}
# This block is straight from the paper
#library(ccems) # The TK1 data is included in ccems. 
topology <- list(  
    heads=c("E1S0"), # E1S0 = substrate free E
    sites=list(                    
        c=list(    # c for catalytic site  
            t=c("E1S1","E1S2","E1S3","E1S4")   
        ) # t for tetramer 
    ) 
) # TK1 is 25kDa = 25mg/umole, so 1 mg = .04 umoles 
g <-mkg(topology, activity=T,TCC=F)

redo1=TRUE
redo2=TRUE
doWLS=TRUE
doWLS=FALSE

dd=subset(TK1,(year==1993),select=c(E,S,v))

names(dd)[1:2]= c("ET","ST")#v was in µmoles/min/mg  
dd=transform(dd, ET=ET/4,v=ET*v/(.04*60))# now uM/sec

if (redo1) tops=ems(dd,g,maxTotalPs=8,kIC=10,topN=96) else { # ~9 min/1 cpu
  load("results/ES8K1k10Top96.RData")  # load it to save 9 minutes above
  tops=globalTopN
}
# this creates Table 1 in the application notes
# it also creates this TopN file in the results directory 


getKk <- function(x) {t(x$report[c(paste("E1S",0:3,"_S",sep=""),
              paste("kE1S",1:4,sep="")),"final",drop=FALSE])}
getAIC <- function(x) { x$report["AIC","final"]}
getSSE <- function(x) { x$report["SSE","final"]}
numClose<-function(tops) {
  aic=sapply(tops,getAIC)
  sse=sapply(tops,getSSE)
  minsse=min(sse)
  aich=aic[(sse<(1.2*minsse))]
  cat(length(aich),"models out of ",length(aic),"are within 20% of the SSE  min.\n")
}
numClose(tops)

if (doWLS) {
# the next chunk does weighted least squares
  sse=sapply(tops,getSSE)
  minsseName=names(sort(sse)[1])
  minsseName
  minsseModel=tops[[minsseName]]
  I=(minsseModel$d$EY<max(minsseModel$d$EY)/2)
  lv=var(minsseModel$res[I])
  hv=var(minsseModel$res[!I])
  lv
  hv
  wts=minsseModel$res  # just to initiate a vector of the right size 
  wts[I]=1/sqrt(lv)   # weights multiply residuals and thus get squared later
  wts[!I]=1/sqrt(hv)
  wts=wts/mean(wts) # normalize to keep SSE roughly where it was
  dd=data.frame(dd,weights=wts)
}



# now redo just the binaries
if (redo2) tops=ems(dd,g,maxTotalPs=8,doSpurs=FALSE,fullGrid=TRUE,kIC=10,topN=64) else {
  load(paste("results/ES8K1k10Top64",ifelse(is.null(dd$weights),"","W"),".RData",sep="")) 
  tops=globalTopN
}

numClose(tops)
Kk=lapply(tops,getKk)
nms=names(Kk)
rowList=data.frame(NULL)
for (j in nms) {
  rowList=rbind(rowList,Kk[[j]])
}
rownames(rowList)<-nms
aic=sapply(tops,getAIC)
sse=sapply(tops,getSSE)
eDelAIC=exp(-(aic-min(aic)))
wgts=eDelAIC/sum(eDelAIC)
print(sum(wgts))
df=data.frame(aic,sse,wgts,rowList)
M=as.matrix(rowList)
ma=exp(wgts%*%log(M)) # average in space of gibbs free energy changes
options(digits=2)
ma

plotMA<-function(ma,g,dd){
  nms=colnames(ma)
  vals=as.vector(ma)
  names(vals)<-nms
  Kmapping=mkKd2Kj(g)
  mdl=mkModel(g,"ma",Kdparams=vals[1:4], Kd2KjLst=Kmapping,kparams=vals[5:8])
  lgx=log(dd$ST)
  upr=range(lgx)[2]
  lwr=range(lgx)[1]
  del=(upr-lwr)/50
  fineX=exp(seq(lwr,upr,by=del))
  newPnts <- data.frame(ET = rep(dd$ET[1],length(fineX)), ST = fineX)
  df <- simulateData(mdl,predict=newPnts,typeYP="v")$predict
  par(mfrow=c(1,2))
  plot(dd$ST,dd$v,type="p", xlab="[dT] (uM)", ylab="v (uM/s)",main="Model Average")
  lines(df$ST,df$EY) 
  mdl <- simulateData(mdl)
  plot(mdl$d$EY,mdl$res,xlab="fitted value", ylab="residual")
  par(mfrow=c(1,1))
}

plotMA(ma,tops[[1]],dd)
#dev.copy2pdf(file="results/hillFit.pdf")
