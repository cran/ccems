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
#d2=subset(RNR,year==2006,select=c(R,t,m,year)) 
#dd=rbind(d1,d2) 
names(d1)[1:2]=c("RT","tT")# d1 = just the Scott et al 2001 Rt data
#names(dd)[1:2]=c("RT","tT")# dd also includes Rofougaran et al 2006
dd=d1
getID <- function(x) { x$mid}


#Kmappings=mkKd2Kj(g)
#orig=mkModel(g,"orig",d1,Kdparams=c(R2t0=150,R1t0_t=25,R2t0_t=.55, R2t1_t=.55),Keq=c(R2t1_t="R2t0_t"),Kd2Kj=Kmappings)
load("case/results/Rt3tops")
mnms=sapply(tops,getID)
names(tops)<-mnms

if (!is.null(dev.list())) for (i in 2:max(dev.list())) dev.off(i);
windows(width = 6, height = 6,restoreConsole = TRUE) 
par(mfcol=c(1,1),mar=c(4,4,0,0)+.1)
pt=c(.1,1:7,7.6,8:20)
predict <- data.frame(RT = rep(7.6,length(pt)), tT = pt)
#df0 <- simulateData(orig,predict=predict,typeYP="m")$predict  
df1 <- simulateData(tops[[1]],predict=predict,typeYP="m")$predict  
df2 <- simulateData(tops[[2]],predict=predict,typeYP="m")$predict  
#df3 <- simulateData(tops[[3]],predict=predict,typeYP="m")$predict  
df4 <- simulateData(tops[["HDFF"]],predict=predict,typeYP="m")$predict  
plot(tops[[1]]$d$tT,tops[[1]]$d$m,cex.lab=1.2,cex.lab=1.2,xlab="Total  [dTTP]  (uM)",ylab="Average Mass  (kDa)",ylim=c(85,185))
lines(df1$tT,df1$EY,lty="solid",lwd=2) 
lines(df2$tT,df2$EY,lty="solid",lwd=1)  
#lines(df3$tT,df3$EY,lty="solid",lwd=0.5)  
lines(df4$tT,df4$EY,lty="dotted") 
legend(1,182,c(names(tops)[1:2],"HDFF"),lwd=c(2,1,1),lty=c(1,1,3))

