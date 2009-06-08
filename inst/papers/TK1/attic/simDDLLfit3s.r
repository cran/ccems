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
            t=c("E1S1","E1S2","E1S3","E1S4")   
        ) # t for tetramer 
    )
)   # in transform below, TK1 is 25kDa => 25mg/umole
g <-mkg(topology, activity=TRUE,TCC=FALSE)

KS=c("E1S0_S","E1S1_S","E1S2_S","E1S3_S")
kS=c("kE1S1","kE1S2","kE1S3","kE1S4")
if (!is.null(dev.list())) for (i in 2:max(dev.list())) dev.off(i);
#if (.Platform$OS.type=="windows")  # now create a window for ccems fits 
#  windows(width = 12, height = 2,restoreConsole = TRUE,ypos=100) else X11(width=10,height=6)
par(mfrow=c(1,1))
options(digits=3)
SEQ=c(0.01,0.02,0.05, .1, .2, .5, 1, 2, 5,10 )
fineX=rep(SEQ,each=3) 
predict <- data.frame(ET = rep(.0001,length(fineX)), ST = fineX)
Kmapping = mkKd2Kj(g)
#par(mfrow=c(1,5))#,mar=c(2,2,1,0)+.1,cex=1,oma=c(2,2,0,0))
modS="DDLL.DDDD"
Keq=c(E1S3_S="E1S2_S",E1S1_S="E1S0_S"); 
keq=c(kE1S2="kE1S1",kE1S3="kE1S1",kE1S4="kE1S1")
Kt = c(0.9, 0.9, 0.5, 0.5) 
kt = c(4,4,4,4)
names(Kt)<-KS
names(kt)<-kS
mdl = mkModel(g,modS,Kdparams=Kt, Kd2KjLst=Kmapping,Keq=Keq,kparams=kt,keq=keq)
df <- simulateData(mdl,predict=predict,typeYP="k")$predict  
names(df)[3]<-"k"
df=transform(df,k=rnorm(k,k,0.05*k))
#df=transform(df,weights=1000)
plot(df$ST,df$k,xlab="Total [dT]",log="xy", ylab="k (1/sec)",cex.main=.7, 
    pch=1,main=modS)
#lines(df$ST,df$k)
ng=mkGrids(g,minTotalPs=3,maxTotalPs=3,KIC=0.6,kIC=4,m1=-100)
chunk=ng$chunk
chunk

Keqs=ng$Keqs
keqs=ng$keqs
mdlNames=rownames(chunk)
lmdlNames=strsplit(mdlNames,split=".",fixed=TRUE)
names(lmdlNames)<-mdlNames
models=NULL
for (j in mdlNames) 
  models[[j]]=mkModel(g,j,df,Kdparams=chunk[j,2:(g$nZ+1),drop=FALSE], 
      Keq=Keqs[[lmdlNames[[j]][1]]],Kd2KjLst=Kmapping,
      pparams=chunk[j,c("p","m1"),drop=FALSE],indx=chunk[j,"indx"],nParams=chunk[j,"nParams"],
      kparams=chunk[j,paste("k",g$Z,sep="")],keq=keqs[[lmdlNames[[j]][2]]])
models
fmodels=lapply(models,fitModel)
print(modS)
for (j in mdlNames) 
  lines(fmodels[[j]]$d$ST,fmodels[[j]]$d$EY,lty=2,lwd=0.5)







