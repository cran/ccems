# this file is scrap. Meaningful segments get shipped out with new file names.
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
#            t=c("E1S1","E1S2","E1S3","E1S4")   
      h=c("E1S1","E1S2","E1S3","E1S4","E1S5","E1S6")   
) # t for tetramer 
    )
)   # in transform below, TK1 is 25kDa => 25mg/umole
g <-mkg(topology, activity=TRUE,TCC=FALSE)
dim(mkGrids(g,fullGrid=T)$chunk)




## Warning: the next line clears all existing figures!!
if (!is.null(dev.list())) for (i in 2:max(dev.list())) dev.off(i);

if (.Platform$OS.type=="windows")  # now create a window for ccems fits 
  windows(width = 10, height = 2,restoreConsole = TRUE) else X11(width=10,height=2)
par(mfcol=c(1,5),mar=c(4,4,2,1)+.1)
load("case/results/outs")

options(stringsAsFactors = FALSE)
mymins=list(NULL)
doAICS=FALSE
modAve=TRUE
for (i in 1:5) {
  d=subset(TK1,index==i,select=c(E,S,k,frstAut,year)) # needed only to get titles in plot in next line
  with(outs[[i]]$df,plot(numP,sse,log="y",ylab="SSE",xlab="Number of Parameters",main=paste(d[1,"frstAut"],d[1,"year"])))
  outs[[i]]$df=data.frame(outs[[i]]$df,names=rownames(outs[[i]]$df))
  print(outs[[i]]$df[1:5,])
  if (i==5) outs[[i]]$df=subset(outs[[i]]$df,subset=numP<7)  # trim neg correction denoms for model averages below
  minmod<-outs[[i]]$df[which(outs[[i]]$df$sse==min(outs[[i]]$df$sse)),]
  if (doAICS) minmod<-outs[[i]]$df[which(outs[[i]]$df$aic==min(outs[[i]]$df$aic)),]
  mymins=rbind(mymins,minmod)
}

if (.Platform$OS.type=="windows")  # now create a window for ccems fits 
  windows(width = 12.5, height = 5,restoreConsole = TRUE,ypos=-25) else X11(width=12.5,height=5)
par(mfcol=c(2,5),mar=c(4,4,2,1)+.1)

KS=c("E1S0_S","E1S1_S","E1S2_S","E1S3_S")
kS=c("kE1S1","kE1S2","kE1S3","kE1S4")

for (i in 1:5) {
  d=subset(TK1,index==i,select=c(E,S,k,frstAut,year))
  plot(d$S,d$k,xlab="Total [dT]",log="xy", ylab="k (1/sec)", 
      main=paste(d[1,"frstAut"],d[1,"year"]))
  names(d)[1:2]= c("ET","ST")
  lgx=log(d$ST)
  upr=range(lgx)[2]
  lwr=range(lgx)[1]
  del=(upr-lwr)/50
  fineX=exp(seq(lwr,upr,by=del))
  predict <- data.frame(ET = rep(d$ET[1],length(fineX)), ST = fineX)
  Keq=NULL;keq=NULL
  nm=mymins[i,"names"]
  if (modAve) {  Keq=NULL; keq=NULL; nm="Model Average"}
  Kmapping=mkKd2Kj(g)
  Kic=mymins[i,KS,drop=FALSE]
  if (modAve) Kic=outs[[i]]$ma[1:4]
  names(Kic)<-KS
  kic=mymins[i,kS,drop=FALSE]
  if (modAve) kic=outs[[i]]$ma[,5:8]
  names(kic)<-kS
  mdl=mkModel(g,nm,d,Kdparams=Kic, Kd2KjLst=Kmapping,kparams=kic)
  df <- simulateData(mdl,predict=predict,typeYP="k")$predict  
  smdl <- simulateData(mdl,init=TRUE)  
  fmdl=smdl
  lines(df$ST,df$EY) 
  options(digits=2)
  Kstrn="         K = "
  for (j in KS) 
    Kstrn=paste(Kstrn,sprintf("%3.2f",Kic[j]),sep=ifelse(j!="E1S0_S",", ",""))
  Kstrn 
  kstrn="         k = "
  for (j in kS) 
    kstrn=paste(kstrn,sprintf("%3.2f",kic[j]),sep=ifelse(j!="kE1S1",", ",""))
  kstrn
  mtext(Kstrn,line=-2,side=1,font=1,cex=0.7)
  mtext(kstrn,line=-3,side=1,font=1,cex=0.7)
#  mtext(paste("N = ",length(d$k),", P = ",mymins[i,"numP"],sep=""),line=-4,side=1,font=1,cex=0.7)
    plot(fmdl$d$EY,fmdl$res,xlab="Fitted Value",ylab="Residual",main=fmdl$mid)
}



#  if (i==1) {nm="DDLM.DFFM";Keq=c(E1S1_S="E1S0_S"); keq=c(kE1S2="kE1S1")  }
#  if (i==2) {nm="DFLM.DDLM";keq=c(kE1S2="kE1S1")  }
#  if (i==3) {nm="DFLM.DDLM";keq=c(kE1S2="kE1S1")  }
#  if (i==4) {nm="DFLM.DDDM";keq=c(kE1S2="kE1S1",kE1S3="kE1S1")  }
#  if (i==5) {nm="DFFF.DFLM";Keq=c(E1S2_S="E1S1_S",E1S3_S="E1S1_S") }
#  if (doAICS) {  # more constraints in this case <=> fewer parameters
#    if (i==1) {nm="DDLL.DDDM";Keq=c(E1S3_S="E1S2_S",                E1S1_S="E1S0_S");keq=c(kE1S2="kE1S1",kE1S3="kE1S1")               }
#    if (i==2) {nm="DDDM.DDDD";Keq=c(                E1S2_S="E1S0_S",E1S1_S="E1S0_S");keq=c(kE1S2="kE1S1",kE1S3="kE1S1",kE1S4="kE1S1") }
#    if (i==3) {nm="DDDD.DFFF";Keq=c(E1S3_S="E1S0_S",E1S2_S="E1S0_S",E1S1_S="E1S0_S");keq=c(              kE1S3="kE1S2",kE1S4="kE1S2") }
#    if (i==4) {nm="DDDD.DDDD";Keq=c(E1S3_S="E1S0_S",E1S2_S="E1S0_S",E1S1_S="E1S0_S");keq=c(kE1S2="kE1S1",kE1S3="kE1S1",kE1S4="kE1S1") }
#    if (i==5) {nm="DDDD.DDDD";Keq=c(E1S3_S="E1S0_S",E1S2_S="E1S0_S",E1S1_S="E1S0_S");keq=c(kE1S2="kE1S1",kE1S3="kE1S1",kE1S4="kE1S1") }
#  }



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


