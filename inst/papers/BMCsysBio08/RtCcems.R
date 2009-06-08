# This file is a manual version of autoRt.r, but not completely as only a short list is specified.
# This file is also a ccems version of the original Rt.r that went with BMC SB 2008, but not completely,
# as again, the list is shortened, but also because E-shaped graphs are not supported in ccems.
# The list of 10 models at the bottom is fitted using snow the old fashion way, i.e. without ems.
rm(list=ls(all=TRUE))  
if (.Platform$OS.type=="windows") home="/users/radivot" else home="/home/radivot"
setwd(home)  # place where models and results will go
if (1) library(ccems) else { # if 0 source in the package 
  pkgNms=dir(paste(home,"/case/active/ccems/ccems/R",sep=""),full.names=TRUE)
  for (i in pkgNms) source(i)  # source in all of the R files of the package
  RNR=read.table(file=paste(home,"/case/active/rnr/datasets/RNR.txt",sep=""),header=TRUE)
}
topology=list(  
    heads=c("R1t0","R2t0"),
    sites=list(
        s=list(
            m=c("R1t1"),
            d=c("R2t1","R2t2")
        )
    )
) 
g=mkg(topology) 
d1=subset(RNR,(year==2001)&(fg==1)&(G==0)&(t>0),select=c(R,t,m,year))
d2=subset(RNR,year==2006,select=c(R,t,m,year)) 
dd=rbind(d1,d2) 
names(d1)[1:2]=c("RT","tT")# d1 = just the Scott et al 2001 Rt data
names(dd)[1:2]=c("RT","tT")# dd also includes Rofougaran et al 2006
Kmappings=mkKd2Kj(g)
models=list(
    mkModel(g,"HDLL",d1,Kdparams=c(R2t0=150,R1t0_t=25,R2t0_t=.55, R2t1_t=.55),Keq=c(R2t1_t="R2t0_t"),Kd2Kj=Kmappings),
    mkModel(g,"IIIJp",d1,Kjparams=c(R2t0=Inf, R1t1=Inf, R2t1=Inf, R2t2=0),pparams=c(p=1)),
    mkModel(g,"IIIJ",d1,Kjparams=c(R2t0=Inf, R1t1=Inf, R2t1=Inf, R2t2=1))
)
fmodels=lapply(models,fitModel) 
pt=c(.1,1:20)
predict <- data.frame(RT = rep(7.6,length(pt)), tT = pt)
df1 <- simulateData(fmodels[[1]],predict=predict,typeYP="m")$predict  
df2 <- simulateData(fmodels[[2]],predict=predict,typeYP="m")$predict  
df3 <- simulateData(fmodels[[3]],predict=predict,typeYP="m")$predict  
plot(fmodels[[1]]$d$tT,fmodels[[1]]$d$m,cex.main=.9,xlab="Total  [dTTP]  (uM)",ylab="Average Mass  (kDa)",ylim=c(85,185))
lines(df1$tT,df1$EY,lty="solid") 
lines(df2$tT,df2$EY,lty="dotted")  
lines(df3$tT,df3$EY,lty="dashed") 

#redo with both datasets to show that 3M is now #1 (3Rp was for d1)
models=list(
    mkModel(g,"HDFF",dd,Kdparams=c(R2t0=150,R1t0_t=25,R2t0_t=.55, R2t1_t=.55),Keq=c(R2t1_t="R2t0_t"),Kd2Kj=Kmappings),
    mkModel(g,"IIIJp",dd,Kjparams=c(R2t0=Inf, R1t1=Inf, R2t1=Inf, R2t2=0),pparams=c(p=1)),
    mkModel(g,"IIIJ",dd,Kjparams=c(R2t0=Inf, R1t1=Inf, R2t1=Inf, R2t2=1))
)
fmodels=lapply(models,fitModel) # assign to avoid dumps

# now make the 10 models as a doubling of 5, i.e. as in Rt.r
models=list(
    mkModel(g,"HDDD",dd,Kdparams=c(R2t0=150,R1t0_t=25,R2t0_t=.55, R2t1_t=.55),Keq=c(R1t0_t="R2t0_t",R2t1_t="R2t0_t"),Kd2Kj=Kmappings),
    mkModel(g,"HDFF",dd,Kdparams=c(R2t0=150,R1t0_t=25,R2t0_t=.55, R2t1_t=.55),Keq=c(R2t1_t="R2t0_t"),Kd2Kj=Kmappings),
    mkModel(g,"IIIJ",dd,Kjparams=c(R2t0=Inf, R1t1=Inf, R2t1=Inf, R2t2=1)),
    mkModel(g,"IIJI",dd,Kjparams=c(R2t0=Inf, R1t1=Inf,   R2t1=1,   R2t2=Inf)),
    mkModel(g,"IIJJ",dd,Kjparams=c(R2t0=Inf, R1t1=Inf,   R2t1=1,   R2t2=1))
    )
pmodels=models
for (i in 1:length(pmodels)) {
  pmodels[[i]]$mid=paste(pmodels[[i]]$mid,"p",sep="")
  pmodels[[i]]$params["p","opt"]=TRUE
}
models=c(models,pmodels)

# Since the model list is already made, ccems cannot help us.
# Also, to keep in touch with how ems works, direct snow use is good exercise. 
library(snow)
# use both cpus on a laptop
strn=rep("localhost",times=2)
cl <- makeCluster(strn, type="SOCK", verbose=TRUE) 
print(clusterCall(cl, function() Sys.info()[c("nodename","machine")]))
clusterEvalQ(cl, require(odesolve))
clusterEvalQ(cl, require(ccems))  # functions already all there from global library(ccems) in parent script
fmodels=clusterApplyLB(cl,models,fitModel)
stopCluster(cl)
getShort<-function(model) list(mid=model$mid,sreport=model$report[model$report$confidenceInterval!="absent",])
lapply(fmodels,getShort)
