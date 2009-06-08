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
g <-mkg(topology, activity=TRUE,TCC=TRUE)
load("case/results/outs")
KS=c("E1S0_S","E1S1_S","E1S2_S","E1S3_S")
kS=c("kE1S1","kE1S2","kE1S3","kE1S4")
gk=rep(0,4) # global average
gK=rep(0,4)
ngk=rep(0,4) # global average
ngK=rep(0,4)
mngk=0
mngK=0
for (i in 1:5) {
  Kic=outs[[i]]$ma[1:4]
  kic=outs[[i]]$ma[,5:8]
  names(Kic)<-KS
  names(kic)<-kS
  gk=gk+kic
  gK=gK+Kic
  mngk=mngk+mean(kic)
  mngK=mngK+mean(Kic)
  ngk=ngk+kic/mean(kic)
  ngK=ngK+Kic/mean(Kic)
}
gk=gk/5
gK=gK/5
options(digits=2)
gk
gK
ngk=(mngk/5)*ngk/5  # average amplitudes and shapes separately
ngK=(mngK/5)*ngK/5
ngk
ngK

# simulate only in region where literature data exists
SEQS=c(0.05, .1, .2, .5, 1, 2, 5,10 )
SEQE=0.0001
fineX=rep(SEQS,each=1) 
predict <- data.frame(ET = rep(SEQE,length(fineX)), ST = fineX)
Kmapping=mkKd2Kj(g)
mdl=mkModel(g,"LA",Kdparams=ngK, Kd2KjLst=Kmapping,kparams=ngk) # literature average
df <- simulateData(mdl,predict=predict,typeYP="k")$predict  
names(df)[3]="k"
df   # this is now treated as the starting dataset

ng=mkGrids(g,minTotalPs=3,maxTotalPs=3,KIC=0.6,kIC=4)
chunk=ng$chunk[c("DDLL.DDDD","DDDD.DFFF"),]
Keqs=ng$Keqs
keqs=ng$keqs
mdlNames=rownames(chunk)
lmdlNames=strsplit(mdlNames,split=".",fixed=TRUE)
names(lmdlNames)<-mdlNames
  print(mdlNames)
  print(lmdlNames)
  models=NULL
  for (j in mdlNames) 
     models[[j]]=mkModel(g,j,df,Kdparams=chunk[j,2:(g$nZ+1),drop=FALSE], 
          Keq=Keqs[[lmdlNames[[j]][1]]],Kd2KjLst=Kmapping,
          pparams=chunk[j,"p",drop=FALSE],indx=chunk[j,"indx"],nParams=chunk[j,"nParams"],
          kparams=chunk[j,paste("k",g$Z,sep="")],keq=keqs[[lmdlNames[[j]][2]]])
models
fmodels=lapply(models,fitModel)
fmodels

SEQS=c(0.01,0.02,0.05, .1, .2, .5, 1, 2, 5,10 )
SEQE=c(0.001,0.01,0.1, 1)
fineX=rep(SEQS,each=1) 
predict <- data.frame(ET = rep(SEQE,each=length(fineX)), ST = fineX)
predict
df1 <- simulateData(fmodels[[1]],predict=predict,typeYP="k")$predict  
df2 <- simulateData(fmodels[[2]],predict=predict,typeYP="k")$predict  

if (!is.null(dev.list())) for (i in 2:max(dev.list())) dev.off(i);
windows(width = 4, height = 4,restoreConsole = TRUE)
library(rgl)

bdf=rbind(rbind(df1,df2))
n=dim(predict)[1]
bdf

plot3d(log10(bdf$ST),log10(bdf$ET),bdf$EY,col=c(rep("violet",n),rep("black ",n)),
type="s",radius=c(rep(.04,2*n)),ylab="[ET]",xlab="[ST]",zlab="k")
material3d(alpha=1)
surface3d(log10(SEQS),log10(SEQE),matrix(df1$EY,ncol=length(SEQE)),col="violet")
# the snipping tool in vista was used to capture the images in Figure 1


