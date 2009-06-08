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

kill3=TRUE
kill3=FALSE

kill1=FALSE
kill1=TRUE

#load("case/results/outs01")
#for (i in 1:5) print(outs[[i]]$ma)
#load("case/results/outs01G")
#for (i in 1:5) print(outs[[i]]$ma)


if ((!kill1)&(!kill3)) load("case/results/outs00")
if ((kill1)&(kill3))   load("case/results/outs31G")
if ((kill1)&(!kill3))  load("case/results/outs01G")



## Warning: the next line clears all existing figures!!
if (!is.null(dev.list())) for (i in 2:max(dev.list())) dev.off(i);


# this is the K trumps k section
if (.Platform$OS.type=="windows")  # now create a window for ccems fits 
  windows(width = 10, height = 2,restoreConsole = TRUE,ypos=250) else X11(width=10,height=2)
par(mfcol=c(1,5),mar=c(4,4,2,1)+.1)
x=1:4
cp=1
for (i in 1:5) {
  df=subset(outs[[i]]$df,subset=wgts>1e-6)
  n=dim(df)[1]
  d=subset(TK1,index==i,select=c(E,S,k,frstAut,year)) # needed only to get titles in plot in next line
  print(n)
  ak=NULL
  aK=NULL
  for (j in 1:n) {
    yK=as.numeric(df[j,5:8])
    yK=yK/mean(yK)
#    print(yK)
#    print(x)
    yk=as.numeric(df[j,9:12])
    yk=yk/mean(yk)
    lmK=lm(yK~x)
    lmk=lm(yk~x)
    aK=c(aK,coef(lmK)["x"])
    aK[abs(aK)>cp]=NA
    ak=c(ak,coef(lmk)["x"])
    ak[abs(ak)>cp]=NA
    #    print(lmK)
#    print(lmk)
  }
  df=data.frame(df,aK=aK,ak=ak)
  ur=dim(subset(df,subset=(ak>0)&(aK>0)))[1]
  lr=dim(subset(df,subset=(ak<0)&(aK>0)))[1]
  ll=dim(subset(df,subset=(ak<0)&(aK<0)))[1]
  ul=dim(subset(df,subset=(ak>0)&(aK<0)))[1]
  pur=sum(subset(df,subset=(ak>0)&(aK>0))$wgts)
  plr=sum(subset(df,subset=(ak<0)&(aK>0))$wgts)
  pll=sum(subset(df,subset=(ak<0)&(aK<0))$wgts)
  pul=sum(subset(df,subset=(ak>0)&(aK<0))$wgts)
  with(df,plot(aK,ak,ylim=c(-cp,cp),pch=19,cex=0.9,xlim=c(-cp,cp),
          ylab="k slope",xlab="K slope",main=paste(d[1,"frstAut"],d[1,"year"])))
  abline(h=0,v=0)
  text(.9,.9,ur)
  text(.9,-.9,lr)
  text(-.9,-.9,ll)
  text(-.9,.9,ul)
  print(c(pll,pul,pur,plr))
  #  with(outs[[i]]$df,plot(numP,aK,pch=1,ylab="K slopes",xlab="Number of Parameters",main=paste(d[1,"frstAut"],d[1,"year"])))
  print(df[1:4,])
}



if (.Platform$OS.type=="windows")  # now create a window for ccems fits 
  windows(width = 10, height = 2,restoreConsole = TRUE) else X11(width=10,height=2)
par(mfcol=c(1,5),mar=c(4,4,2,1)+.1)


doAICS=TRUE
modAve=TRUE
#TK1=subset(TK1,subset=(S>.09)&(S<3.1))
#TK1


options(stringsAsFactors = FALSE)
mymins=list(NULL)
for (i in 1:5) {
  d=subset(TK1,index==i,select=c(E,S,k,frstAut,year)) # needed only to get titles in plot in next line
  with(outs[[i]]$df,plot(numP,sse,log="y",ylab="SSE",xlab="Number of Parameters",main=paste(d[1,"frstAut"],d[1,"year"])))
  print(outs[[i]]$df[1:4,])
  outs[[i]]$df=data.frame(outs[[i]]$df,names=rownames(outs[[i]]$df))
  minmod<-outs[[i]]$df[which(outs[[i]]$df$sse==min(outs[[i]]$df$sse)),]
  if (doAICS) minmod<-outs[[i]]$df[which(outs[[i]]$df$aic==min(outs[[i]]$df$aic)),]
  mymins=rbind(mymins,minmod)
}
mymins


if (.Platform$OS.type=="windows")  # now create a window for ccems fits 
  windows(width = 12.5, height = 5,restoreConsole = TRUE,ypos=-25) else X11(width=12.5,height=5)
par(mfcol=c(2,5),mar=c(4,4,2,1)+.1)

for (i in 1:5) {
  d=subset(TK1,index==i,select=c(E,S,k,frstAut,year))
  if ((i==1)&(kill3)) d=d[-(14:16),]
  if ((i==3)&(kill1)) d=d[-1,]
  plot(d$S,d$k,xlab="Total [dT]",log="xy", ylab="k (1/sec)", 
      main=paste(d[1,"frstAut"],d[1,"year"]),cex.lab=1.1)
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
  options(digits=3)
  Kstrn="         K = "
  for (j in KS) 
    Kstrn=paste(Kstrn,sprintf("%3.2f",Kic[j]),sep=ifelse(j!="E1S0_S",", ",""))
  Kstrn 
  kstrn="         k = "
  for (j in kS) 
    kstrn=paste(kstrn,sprintf("%3.2f",kic[j]),sep=ifelse(j!="kE1S1",", ",""))
  kstrn
  mtext(Kstrn,line=-2,side=1,font=1,cex=0.7,adj=.4)
  mtext(kstrn,line=-3,side=1,font=1,cex=0.7,adj=.4)
#  mtext(paste("N = ",length(d$k),", P = ",mymins[i,"numP"],sep=""),line=-4,side=1,font=1,cex=0.7)
  plot(fmdl$d$EY,fmdl$tres,xlab="Fitted Value",ylab="Transformed Residual",main=fmdl$mid,cex.lab=1.1)
}





