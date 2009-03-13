rm(list=ls(all=TRUE))  # clean up left overs from previous run, though shouldn't need to.
# This code reproduces the ccems Bioinformatics Applications Note Example
if (.Platform$OS.type=="windows") home="/users/radivot" else home="/home/radivot"
setwd(home)  # place where models and results will go
if (1) library(ccems) else { # if 0 source in the package 
  pkgNms=dir(paste(home,"/case/active/ccems/ccems/R",sep=""),full.names=TRUE)
  for (i in pkgNms) source(i)  # source in all of the R files of the package
  TK1=read.table(file=paste(home,"/case/active/papers/TK1/data/TK1.txt",sep=""),header=TRUE)
}
topology <- list(  
    heads=c("E1S0"), #one E is a tetramer
    sites=list(                    
        c=list(    # c-site = catyltic site  
            t=c("E1S1","E1S2","E1S3","E1S4")   
        )
    )
)
g <- mkg(topology,hubChar="E",activity=T,TCC=F)
dd=subset(TK1,(year==2000),select=c(E,S,v))
names(dd)[1:2]= c("ET","ST")
# the next line runs ~5 minutes on 1 cpu. It makes Table 1.
tops=ems(dd,g,maxTotalPs=5,kIC=30000) 

# this block shows that 25 out of 70 models are within 20% of the best SSE
load("results/ES5K1k30000.RData")  
getAIC <- function(x) { x$sreport["AIC","final"]}
getSSE <- function(x) { x$sreport["SSE","final"]}
aic=sapply(sAllModels,getAIC)
sse=sapply(sAllModels,getSSE)
minsse=min(sse)
aich=aic[(sse<(1.2*minsse))]
length(aich)  
length(aic)

# this time only fit the grid models, and fit all of them, including the full grid multiplied into k space  
tops=ems(dd,g,kIC=30000, maxTotalPs=8,showConstr=TRUE,doSpurs=FALSE,fullGrid=TRUE)
load("results/ES8K1k30000.RData")  
getAIC <- function(x) { x$sreport["AIC","final"]}
getSSE <- function(x) { x$sreport["SSE","final"]}
getKk <- function(x) { t(x$sreport[c(paste("E1S",0:3,"_S",sep=""),paste("kE1S",1:4,sep="")),"final",drop=FALSE])}
aic=sapply(sAllModels,getAIC)
sse=sapply(sAllModels,getSSE)
Kk=lapply(sAllModels,getKk)
Kk
nms=names(Kk)
rowList=data.frame(NULL)
  for (j in nms) {
    rowList=rbind(rowList,Kk[[j]])
  }
rownames(rowList)<-nms
rowList
eDelAIC=exp(-(aic-min(aic)))
wgts=eDelAIC/sum(eDelAIC)
df=data.frame(aic,sse,wgts,rowList)
df
N=4  # this gets 87% of the weight. N=5 gets 96%, but the 5th model has an outlier k1 value that skews the average, see df top.  
w=wgts[1:N]
sum(w)
w=w/sum(w)
M=as.matrix(rowList)[1:N,]
M
options(digits=2) # allow easier pasting into Word
w%*%M     # average in K and k space directly (not advisable, just to compare)
exp(w%*%log(M)) # average in space of gibbs free energy changes

# This is not in the Bioninformatics paper. It plots the top two fits/results.
load("results/EStop.RData") 
dd=subset(TK1,(year==2000),select=c(E,S,v))
names(dd)[1:2]= c("ET","ST")
attach(dd)
plot(ST,v,type="p",pch=1, xlab="[dT] (uM)", ylab="v")
lgx=log(ST)
upr=range(lgx)[2]
lwr=range(lgx)[1]
del=(upr-lwr)/100
fineX=exp(seq(lwr,upr,by=del))
newPnts <- data.frame(ET = rep(ET[1],length(fineX)), ST = fineX)
for (i in 1:3) {
  df <- simulateData(globalTopN[[i]],predict=newPnts,typeYP="v")$predict  
  lines(df$ST,df$EY,type="l",lty=i)
}
detach(dd)
dev.copy2pdf(file="results/ESfit.pdf")
#dev.off()
