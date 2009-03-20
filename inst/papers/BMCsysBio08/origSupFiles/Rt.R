# it is assumed that odesolve, snow, and Rmpi or Rpvm have been installed 
rm(list=ls(all=TRUE))  # clean up left overs from previous run, though shouldn't need to.
setwd("/home/radivot/case/active/rnr/Rt/R")  # directory where the next two files live
load("RNR.RData")  # load RNR adata (all data) list object
source("fRt.r")    # R function definitions are in here 

# the next line generates and compiles C code in the models directory (which it creates if needed)
g=mkgObj("Rt",c("Rt","RR","RRt","RRtt")) # generic g=0 model. first arg is id, second is complexes
RtData2=adata[c("f1a01")]  # pull out just the Scott et al 2001 Rt data
RtData4=adata[c("f1a01","f1a06","f1b06")]  # pull out the Rt data of Rofougaran et al 2006 also
# the following two functions map Kd parameters into Kj parameters as shown in Table 1
Eshape<-function(x) c(x[1], 2*x[2],x[1]*x[3], 2*x[1]^2*x[4])   
nshape<-function(x) c(x[1], 2*x[2],x[2]*x[3], 2*x[2]*x[3]*x[4])

# Table 2
RtData=RtData2  # Scott et al only
models=list(
mkModelObj(RtData,g,"2E",Kdparams=c(R_t=25, R_R=75, RR_t=.55, RRt_t=.55),Keq=c(RRt_t="RR_t"),Kd2Kj=nshape),
mkModelObj(RtData,g,"3Rp",Kjparams=c(Rt=Inf,RR=Inf, RRt=Inf, RRtt=0),pparams=c(pRT=1)),
mkModelObj(RtData,g,"3M",Kjparams=c(Rt=Inf,RR=Inf, RRt=Inf, RRtt=1))
)
models=fitMS(models,"MS2")

# Figure 4 plots of Table 2 fits against data
options("warn"=-1)
if (!is.null(dev.list())) for (i in 2:max(dev.list())) dev.off(i); # clear off old plot windows
if (.Platform$OS.type=="windows") windows(width = 4, height = 4,restoreConsole = TRUE) else  
                 X11(width = 4, height = 4);   # open 4x4 plot window on either OS type 
par(mfrow=c(1,1),mar=c(4.1,4.1,0.5,0.5))
plot(models[[3]]$f1a01$d,cex.main=.9,xlab="Total  [dTTP]  (uM)",ylab="Average Mass  (kDa)",ylim=c(85,185))
lines(models[[3]]$f1a01$isd,lty="solid") 
lines(models[[3]]$f1a01$sd,lty="dotted")  
lines(models[[1]]$f1a01$sd,lty="dashed") 
lines(models[[2]]$f1a01$sd,lty="dotdash") 
options("warn"=1) 
dev.copy(pdf,file="results/fig4.pdf",width=4,height=4);dev.off()


# Table 4 
RtData=RtData4  # both datasets fitted jointly
models=list(
mkModelObj(RtData,g,"2E",Kdparams=c(R_t=25, R_R=75, RR_t=.55, RRt_t=.55),Keq=c(RRt_t="RR_t"),Kd2Kj=nshape),
mkModelObj(RtData,g,"3Rp",Kjparams=c(Rt=Inf,RR=Inf, RRt=Inf, RRtt=0),pparams=c(pRT=1)),
mkModelObj(RtData,g,"3M",Kjparams=c(Rt=Inf,RR=Inf, RRt=Inf, RRtt=1))
)
fitMS(models,"MS4")



#fits to whole model space of Table 1
RtData=RtData4 # using both datasets jointly
models=list(
mkModelObj(RtData,g,"2A",Kdparams=c(R_t=1,  R_R=1,  Rt_R=1, Rt_Rt=1),Keq=c(Rt_R="R_R",Rt_Rt="R_R"),Kd2Kj=Eshape),
mkModelObj(RtData,g,"2B",Kdparams=c(R_t=1,  R_R=1,  Rt_R=1, Rt_Rt=1),Keq=c(Rt_R="R_R"),Kd2Kj=Eshape),
mkModelObj(RtData,g,"2C",Kdparams=c(R_t=1,  R_R=1,  Rt_R=1, Rt_Rt=1),Keq=c(Rt_Rt="Rt_R"),Kd2Kj=Eshape),
mkModelObj(RtData,g,"2D",Kdparams=c(R_t=1,  R_R=1,  Rt_R=1, Rt_Rt=1),Keq=c(Rt_Rt="R_R"),Kd2Kj=Eshape),
mkModelObj(RtData,g,"2E",Kdparams=c(R_t=25, R_R=75, RR_t=.55, RRt_t=.55),Keq=c(RRt_t="RR_t"),Kd2Kj=nshape),
mkModelObj(RtData,g,"2F",Kdparams=c(R_t=1,  R_R=1,  Rt_R=1, Rt_Rt=1),Kd2Kj=Eshape),
mkModelObj(RtData,g,"2G",Kdparams=c(R_t=1,  R_R=Inf,Rt_R=1, Rt_Rt=1),Keq=c(Rt_Rt="Rt_R"),Kd2Kj=Eshape),
mkModelObj(RtData,g,"2I",Kdparams=c(R_t=1,  R_R=1,  Rt_R=Inf, Rt_Rt=1),Keq=c(Rt_Rt="R_R"),Kd2Kj=Eshape),
mkModelObj(RtData,g,"2K",Kdparams=c(R_t=1,  R_R=1,  Rt_R=1, Rt_Rt=Inf),Keq=c(Rt_R="R_R"),Kd2Kj=Eshape),
mkModelObj(RtData,g,"2M",Kdparams=c(R_t=Inf,R_R=75, RR_t=.55, RRt_t=.55),Keq=c(RRt_t="RR_t"),Kd2Kj=nshape),
mkModelObj(RtData,g,"2N",Kdparams=c(R_t=Inf,R_R=75, RR_t=.55, RRt_t=.55),Kd2Kj=nshape),
mkModelObj(RtData,g,"3A",Kjparams=c(Rt=1,  RR=1,   RRt=1,   RRtt=1)),
mkModelObj(RtData,g,"3B",Kjparams=c(Rt=1,  RR=Inf, RRt=1,   RRtt=1)),
mkModelObj(RtData,g,"3C",Kjparams=c(Rt=1,  RR=1,   RRt=Inf, RRtt=1)),
mkModelObj(RtData,g,"3D",Kjparams=c(Rt=1,  RR=1,   RRt=1,   RRtt=Inf)),
mkModelObj(RtData,g,"3E",Kjparams=c(Rt=1,  RR=Inf, RRt=1,   RRtt=1)),
mkModelObj(RtData,g,"3F",Kjparams=c(Rt=1,  RR=Inf, RRt=Inf, RRtt=1)),
mkModelObj(RtData,g,"3G",Kjparams=c(Rt=1,  RR=Inf, RRt=1,   RRtt=Inf)),
mkModelObj(RtData,g,"3H",Kjparams=c(Rt=1,  RR=1,   RRt=Inf, RRtt=Inf)),
mkModelObj(RtData,g,"3I",Kjparams=c(Rt=Inf,RR=Inf, RRt=1,   RRtt=1)),
mkModelObj(RtData,g,"3J",Kjparams=c(Rt=Inf,RR=1,   RRt=Inf, RRtt=1)),
mkModelObj(RtData,g,"3K",Kjparams=c(Rt=Inf,RR=1,   RRt=1,   RRtt=Inf)),
mkModelObj(RtData,g,"3L",Kjparams=c(Rt=1,  RR=Inf, RRt=Inf, RRtt=Inf)),
mkModelObj(RtData,g,"3M",Kjparams=c(Rt=Inf,RR=Inf, RRt=Inf, RRtt=1)),
mkModelObj(RtData,g,"3N",Kjparams=c(Rt=Inf,RR=Inf, RRt=1,   RRtt=Inf)),
mkModelObj(RtData,g,"3O",Kjparams=c(Rt=Inf,RR=1,   RRt=Inf, RRtt=Inf)),
mkModelObj(RtData,g,"3P",Kjparams=c(Rt=Inf,RR=Inf, RRt=Inf, RRtt=Inf)),
mkModelObj(RtData,g,"3Q",Kjparams=c(Rt=0,  RR=Inf, RRt=Inf, RRtt=Inf)),
mkModelObj(RtData,g,"3R",Kjparams=c(Rt=Inf,RR=Inf, RRt=Inf, RRtt=0)),
mkModelObj(RtData,g,"3S",Kjparams=c(Rt=Inf,RR=Inf, RRt=0,   RRtt=Inf)),
mkModelObj(RtData,g,"3T",Kjparams=c(Rt=Inf,RR=0,   RRt=Inf, RRtt=Inf))
)
pmodels=models
for (i in 1:length(pmodels)) {
  pmodels[[i]]$mid=paste(pmodels[[i]]$mid,"p",sep="")
  pmodels[[i]]$params["pRT","opt"]=TRUE
  }
models=c(models,pmodels)
smodels=fitMS(models,"M4s",cpus=4,ptype="PVM",smart=1) 
smodels=fitMS(models,"M4f",cpus=4,ptype="PVM",smart=0)

#Table 2 data only
RtData=RtData2 
models=list(
mkModelObj(RtData,g,"2A",Kdparams=c(R_t=1,  R_R=1,  Rt_R=1, Rt_Rt=1),Keq=c(Rt_R="R_R",Rt_Rt="R_R"),Kd2Kj=Eshape),
mkModelObj(RtData,g,"2B",Kdparams=c(R_t=1,  R_R=1,  Rt_R=1, Rt_Rt=1),Keq=c(Rt_R="R_R"),Kd2Kj=Eshape),
mkModelObj(RtData,g,"2C",Kdparams=c(R_t=1,  R_R=1,  Rt_R=1, Rt_Rt=1),Keq=c(Rt_Rt="Rt_R"),Kd2Kj=Eshape),
mkModelObj(RtData,g,"2D",Kdparams=c(R_t=1,  R_R=1,  Rt_R=1, Rt_Rt=1),Keq=c(Rt_Rt="R_R"),Kd2Kj=Eshape),
mkModelObj(RtData,g,"2E",Kdparams=c(R_t=25, R_R=75, RR_t=.55, RRt_t=.55),Keq=c(RRt_t="RR_t"),Kd2Kj=nshape),
mkModelObj(RtData,g,"2F",Kdparams=c(R_t=1,  R_R=1,  Rt_R=1, Rt_Rt=1),Kd2Kj=Eshape),
mkModelObj(RtData,g,"2G",Kdparams=c(R_t=1,  R_R=Inf,Rt_R=1, Rt_Rt=1),Keq=c(Rt_Rt="Rt_R"),Kd2Kj=Eshape),
mkModelObj(RtData,g,"2I",Kdparams=c(R_t=1,  R_R=1,  Rt_R=Inf, Rt_Rt=1),Keq=c(Rt_Rt="R_R"),Kd2Kj=Eshape),
mkModelObj(RtData,g,"2K",Kdparams=c(R_t=1,  R_R=1,  Rt_R=1, Rt_Rt=Inf),Keq=c(Rt_R="R_R"),Kd2Kj=Eshape),
mkModelObj(RtData,g,"2M",Kdparams=c(R_t=Inf,R_R=75, RR_t=.55, RRt_t=.55),Keq=c(RRt_t="RR_t"),Kd2Kj=nshape),
mkModelObj(RtData,g,"2N",Kdparams=c(R_t=Inf,R_R=75, RR_t=.55, RRt_t=.55),Kd2Kj=nshape),
mkModelObj(RtData,g,"3A",Kjparams=c(Rt=1,  RR=1,   RRt=1,   RRtt=1)),
mkModelObj(RtData,g,"3B",Kjparams=c(Rt=1,  RR=Inf, RRt=1,   RRtt=1)),
mkModelObj(RtData,g,"3C",Kjparams=c(Rt=1,  RR=1,   RRt=Inf, RRtt=1)),
mkModelObj(RtData,g,"3D",Kjparams=c(Rt=1,  RR=1,   RRt=1,   RRtt=Inf)),
mkModelObj(RtData,g,"3E",Kjparams=c(Rt=1,  RR=Inf, RRt=1,   RRtt=1)),
mkModelObj(RtData,g,"3F",Kjparams=c(Rt=1,  RR=Inf, RRt=Inf, RRtt=1)),
mkModelObj(RtData,g,"3G",Kjparams=c(Rt=1,  RR=Inf, RRt=1,   RRtt=Inf)),
mkModelObj(RtData,g,"3H",Kjparams=c(Rt=1,  RR=1,   RRt=Inf, RRtt=Inf)),
mkModelObj(RtData,g,"3I",Kjparams=c(Rt=Inf,RR=Inf, RRt=1,   RRtt=1)),
mkModelObj(RtData,g,"3J",Kjparams=c(Rt=Inf,RR=1,   RRt=Inf, RRtt=1)),
mkModelObj(RtData,g,"3K",Kjparams=c(Rt=Inf,RR=1,   RRt=1,   RRtt=Inf)),
mkModelObj(RtData,g,"3L",Kjparams=c(Rt=1,  RR=Inf, RRt=Inf, RRtt=Inf)),
mkModelObj(RtData,g,"3M",Kjparams=c(Rt=Inf,RR=Inf, RRt=Inf, RRtt=1)),
mkModelObj(RtData,g,"3N",Kjparams=c(Rt=Inf,RR=Inf, RRt=1,   RRtt=Inf)),
mkModelObj(RtData,g,"3O",Kjparams=c(Rt=Inf,RR=1,   RRt=Inf, RRtt=Inf)),
mkModelObj(RtData,g,"3P",Kjparams=c(Rt=Inf,RR=Inf, RRt=Inf, RRtt=Inf)),
mkModelObj(RtData,g,"3Q",Kjparams=c(Rt=0,  RR=Inf, RRt=Inf, RRtt=Inf)),
mkModelObj(RtData,g,"3R",Kjparams=c(Rt=Inf,RR=Inf, RRt=Inf, RRtt=0)),
mkModelObj(RtData,g,"3S",Kjparams=c(Rt=Inf,RR=Inf, RRt=0,   RRtt=Inf)),
mkModelObj(RtData,g,"3T",Kjparams=c(Rt=Inf,RR=0,   RRt=Inf, RRtt=Inf))
)
pmodels=models
for (i in 1:length(pmodels)) {
  pmodels[[i]]$mid=paste(pmodels[[i]]$mid,"p",sep="")
  pmodels[[i]]$params["pRT","opt"]=TRUE
  }
models=c(models,pmodels)
smodels=fitMS(models,"MS2s",cpus=4,ptype="PVM",smart=1) # works on centos
smodels=fitMS(models,"MS2f",cpus=4,ptype="PVM",smart=0) # works on centos

