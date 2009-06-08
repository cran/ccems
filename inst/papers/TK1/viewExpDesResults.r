mkMaTabs<-function(outs) {
  BG=NULL
  nms=names(outs)
  for (j in 1:length(nms) )
  {
    df=outs[[j]]$ma
    BG=rbind(BG,df)
  }
#  row.names(BG)<-NULL
  BG=data.frame(Dataset=1:length(nms),BG)
  names(BG)[2:9]<-c(paste("K",1:4,sep=""),paste("k",1:4,sep=""))
  BG
}

load("case/results/adaptOuts") 


if (!is.null(dev.list())) for (i in 2:max(dev.list())) dev.off(i);
windows(width = 10, height = 7.5,restoreConsole = TRUE,ypos=0) 
par(mfrow=c(3,4),mar=c(4,4,2,1)+.1)
par(oma=c(3,3,0,0))
#par(mar=c(0,0,0,0))

for (reps in 1:3) {
  load(paste("case/results/FFoutsR",reps,sep="")) 
  print(params<-outs[[length(outs)-1]])
  outs=outs[-((length(outs)-1):length(outs))]
  kt=params[5:8]
  Kt=params[1:4]
  names(kt)<-(kS<-paste("k",1:4,sep=""))
  names(Kt)<-(KS<-paste("K",1:4,sep=""))
  mdlNames=row.names(outs[[1]]$df)
  Wgts=rep(0,length(mdlNames))
  names(Wgts)<-mdlNames
  Firsts=Wgts
  for (i in 1:(length(outs))) 
  {
    Wgts=Wgts+outs[[i]]$df[mdlNames,"wgts"]
    name1=row.names(outs[[i]]$df)[1]
    Firsts[name1]=Firsts[name1]+1
  }
  Wgts=Wgts/length(outs)
  print(sort(Firsts, decreasing = TRUE))
  print(sort(Wgts, decreasing = TRUE))
  tabMAs=mkMaTabs(outs)
  print(summary(tabMAs))
  print(mns<-mean(tabMAs))
  plot((1:4)-.3,Kt,pch="-",cex=2,ylab="uM",xlab="Parameter",xaxt="n",ylim=c(0.1,1.2),xlim=c(0.5,4.5))
#  plot(Kt,pch="-",xlab="K",xaxt="n",type="n",yaxt="n",ylim=c(0.1,2))
  axis(1,at=1:4,label=KS) 
  mtext(paste(reps," Replicate",ifelse(reps>1,"s",""),sep=""),line=5.2,side=2,font=1,cex=0.8,adj=0.5)
  if (reps==1) mtext("Full Factorial",line=5.2,side=1,font=1,cex=0.8,adj=0.5)
  for (i in 1:4) points(jitter(rep(i,100),a=.2),tabMAs[,1+i],pch=19,cex=.6)
  points(1:4,mns[2:5],pch=19,cex=2)
  #  if (reps==1) axis(1) 
  plot.new()
  plot((1:4)-.3,kt,pch="-",cex=2,xlab="Parameter",ylab="Per Second",xaxt="n",ylim=c(1.2,5),xlim=c(0.5,4.5))
  for (i in 1:4) points(jitter(rep(i,100),a=.2),tabMAs[,5+i],pch=19,cex=.6)
  axis(1,at=1:4,label=kS) 
  points(1:4,mns[6:9],pch=19,cex=2)
  if (reps==1) mtext("Full Factorial",line=5.2,side=1,font=1,cex=0.8,adj=0.5)
  plot.new()
  #  plot(kt,pch="-",xlab="k",xaxt="n",type="n",yaxt="n",ylim=c(1,4.2))
#  axis(2) 
#  if (reps==1) axis(3) 
}


#threeOnly<-function(outs,gibbs=TRUE) {
#  KS=c("E1S0_S","E1S1_S","E1S2_S","E1S3_S")
#  kS=c("kE1S1","kE1S2","kE1S3","kE1S4")
#  for (i in 1:length(outs) ) {
#    outs[[i]]$df=outs[[i]]$df[outs[[i]]$df$numP==3,]
#    outs[[i]]$df$wgts=outs[[i]]$df$wgts/sum(outs[[i]]$df$wgts)
#    M=as.matrix(outs[[i]]$df[,c(KS,kS)])
#    if (gibbs)     ma=exp(outs[[i]]$df$wgts%*%log(M)) # average in space of gibbs free energy changes
#    if (!gibbs)     ma=outs[[i]]$df$wgts%*%M # average in straight params
#    outs[[i]]$ma=ma
#  } 
#  outs
#}
#outs=threeOnly(outs)   # comment to get lit ave of Fig. 4, leave in to get it for 3Ps only
#outs
#

