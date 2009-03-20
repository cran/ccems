# fRt.r :  
mkgObj<-function(id,Z,monomerMass=90)  
{ 
 nZ=length(Z);
 reactantS=strsplit(Z,NULL);
 names(reactantS)<-Z; 
 atomS=unique(unlist(reactantS))
 print(atomS)
 print(id) 
 stateS=c(atomS,Z)
 nStateS=length(stateS)
 print(stateS)  
 WS=union("R",atomS)
 W=NULL
 for (nameWS in WS) { 
    w=sapply(gregexpr(nameWS,stateS),length); 
    w[regexpr(nameWS,stateS)<0]=0;
    W=cbind(W,w)
   }
 W=as.data.frame(W,row.names=stateS)
 names(W)<-WS
 nAtomS=length(atomS)
 initialStateTCC=rep(0,nAtomS) 
 names(initialStateTCC)<-atomS
# parmsTCC=c(rep(1,nAtomS),1,rep(1,nZ),rep(1,nZ))
# pnamesTCC=c(paste(atomS,"T",sep=""),"pRT",paste("Kj",Z,sep="_"),paste("Msk",Z,sep="_"))
 parmsTCC=c(rep(1,nAtomS),1,rep(1,nZ))
 pnamesTCC=c(paste(atomS,"T",sep=""),"pRT",paste("Kj",Z,sep="_"))
 names(parmsTCC)<-pnamesTCC
 gObj=list(id=id,Z=Z,nZ=nZ,
      atomS=atomS, nAtomS=nAtomS,
      stateS=stateS, nStateS=nStateS,
      reactantS=reactantS,
      W=W,monomerMass=monomerMass,
      wDir=sub("C:","",getwd()),
      sstime = 1e6,rtol=1e-5,atol=1e-7,
      parmsTCC=parmsTCC,
      initialStateTCC= initialStateTCC
 ); 
 gObj=mkZ(gObj)
library(odesolve); gObj=mkgC(gObj); testgC(gObj); 
 gObj
}

mkZ<-function(gObj) {
    if (sum(dir()=="models")==0) system("mkdir models")
    attach(gObj)
    setwd(paste(wDir,"models",sep="/"))
    ofile=file(fn<-paste(id,".r",sep=""),"wt")
    cat("fback<-function(y, p) {\n",file=ofile)
    for (i in 1:nAtomS)
      cat(sprintf("%10s = y[%d];\n",atomS[i],i),file=ofile)
    for (Zj in Z)
      cat(sprintf("%10s = %s/p[\"Kj_%s\"];\n",Zj,paste(reactantS[[Zj]],collapse="*"),Zj),file=ofile)
#      cat(sprintf("%10s = p[\"Msk_%s\"]*%s/p[\"Kj_%s\"];\n",Zj,Zj,paste(reactantS[[Zj]],collapse="*"),Zj),file=ofile)
    cat("c(",file=ofile)
    for (i in 1:nStateS)
       cat(sprintf("%s=%s%s",stateS[i],stateS[i],ifelse(i<length(stateS),",",")\n}\n\n")),file=ofile)
    close(ofile)
    source(paste(wDir,"/models/",fn,sep=""),local=TRUE)
    unlink(paste(wDir,"/models/",fn,sep=""))
    gObj$fback=fback
    setwd(wDir)
    detach(gObj)
    gObj
}


mkgC<-function(gObj) {
    attach(gObj)
    setwd(paste(wDir,"models",sep="/"))
    print("in mkgC")
    ofile=file(fn<-paste(id,".c",sep=""),"wt")
    cat(sprintf("/* compile with R CMD SHLIB %s.c or (on Windows) Rcmd SHLIB %s.c */\n",id,id),file=ofile)
    cat("#include <R.h> /* gives F77_CALL through R_ext/RS.h */\n\n",file=ofile)
    pnames=names(parmsTCC)
    cat(sprintf("static double parms[%d];\n",length(pnames)),file=ofile)
    for (i in 1:length(pnames))
      cat(sprintf("#define %25s  parms[%d]\n",pnames[i],i-1),file=ofile)
    cat(sprintf("void %s(void (* odeparms)(int *, double *))\n",id),file=ofile)
    cat(sprintf("{ int N=%d;\n",length(pnames)),file=ofile)   
    cat(" odeparms(&N, parms);",file=ofile); cat("}\n",file=ofile)
    cat("\nvoid myderivs(int *neq, double *curtime, double *statevec, double *ODERHSvec)\n",file=ofile)
    cat("{",file=ofile); print(stateS)
    cat(sprintf(" double %s;\n",paste(stateS,sep='',collapse=',')),file=ofile)
    print(atomS)
    for (i in 1:nAtomS) cat(sprintf("  %s = statevec[%d];",atomS[i],i-1),file=ofile)
    print(Z)
    cat("\n",file=ofile)
   for (Zj in Z) cat(sprintf("%10s = %s/Kj_%s;\n",Zj,paste(reactantS[[Zj]],collapse="*"),Zj),file=ofile)
   for (i in 1:nAtomS)
      {
      if (i==1) cat(sprintf("ODERHSvec[%d] = pRT*%sT ",i-1,atomS[i]),file=ofile) else
                cat(sprintf("ODERHSvec[%d] =     %sT ",i-1,atomS[i]),file=ofile)
      for (j in 1:nStateS) 
          if (j<=nAtomS) 
             {if (W[stateS[j],atomS[i]]==1) 
                cat(sprintf("-%s",stateS[j]),file=ofile)}
           else         
             if (W[stateS[j],atomS[i]]==0)
                        cat("              ",file=ofile)  else
             if (W[stateS[j],atomS[i]]==1)
                cat(sprintf(" - %s  ",stateS[j]),file=ofile) else
                cat(sprintf(" - %d*%s",W[stateS[j],atomS[i]],stateS[j]),file=ofile) 
#                cat(sprintf(" - Msk_%s*%s  ",stateS[j],stateS[j]),file=ofile) else
#                cat(sprintf(" - Msk_%s*%d*%s",stateS[j],W[stateS[j],atomS[i]],stateS[j]),file=ofile) 
      cat(";\n",file=ofile)
      }
   cat("}\n",file=ofile)
   close(ofile)
# note: cpu time savings of putting in blanks for multiplications by 0 and 1 was neglibible => something else is making this 2x slower than 
# hard coding without masks.  Need masks because unlike R, C code doesn't have Inf and doesn't handle x/Inf =0
    if (.Platform$OS.type=="windows") system(sprintf("Rcmd SHLIB %s.c\n",id)) else  system(sprintf("R CMD SHLIB %s.c\n",id))
#    system(sprintf("mv %s%s  bin\n",id,.Platform$dynlib.ext))
    unlink(sprintf("%s.o",id))
    unlink(sprintf("%s.d",id))
    gObj$code=readLines(fn)
setwd(wDir)
detach(gObj)
gObj
}

testgC<-function(gObj) {
    attach(gObj)
    times <- seq(0,10,1)
    print("************* start testgC **************************")
    dyn.load(strn<-paste(wDir,"/models/",id,.Platform$dynlib.ext,sep=""))
    out1 <- lsoda(initialStateTCC,times,"myderivs", parmsTCC, rtol=rtol,atol=atol, dllname=id)
    out=data.frame(out1)
    print(out)
    dyn.unload(strn)
    print("************* end testgC **************************")
    detach(gObj)
}



mkModelObj<-function(Data,g,mid,Kjparams=NULL,Kdparams=NULL,Keq=NULL,pparams=c(pRT=1i),Kd2Kj=NULL)  
{ 
# numbers on imaginary axis are fixed, i=inhibited from moving
# one of Kj or Kd should be NULL, the other then is being used
  if (is.null(Kjparams)&is.null(Kdparams)) {print("One of Kjparams or Kdparams must be specified."); return(0)}
  if (!is.null(Kjparams)&!is.null(Kdparams)) {print("Only one of Kjparams or Kdparams can be specified."); return(0)}
  if (!is.null(Kjparams)) Kj=TRUE else Kj=FALSE
  if(Kj) Kparams=Kjparams else Kparams=Kdparams
#  print(Kparams)
  vparams=c(Kparams,pparams) 
  nSysParams=length(Kparams);
  nParams=length(vparams)
  params = data.frame(initial=Mod(vparams),final=Mod(vparams),opt=(Re(vparams)>0),
                      constr=rep("none",nParams),stringsAsFactors=FALSE)
#  params = data.frame(initial=Mod(vparams),final=Mod(vparams),opt=(Re(vparams)>0),
#                      present=as.numeric(Re(vparams)!=Inf),constr=rep("none",nParams),stringsAsFactors=FALSE)
  if (!is.null(Keq)) {
    params[names(Keq),"constr"]=Keq
    params[names(Keq),"opt"]=FALSE 
    }
  params[Mod(vparams)==Inf,"opt"]=FALSE
  params[Mod(vparams)==0,"opt"]=FALSE
  params[Mod(vparams)==0,"initial"]=.0001
  params[Mod(vparams)==0,"final"]=.0001
  
  c(g,list(mid=mid,params=params,Kparams=Kparams,Kj=Kj,Kd2Kj=Kd2Kj,
#  AIC=list(initial=0,final=0),
#  SSE=list(initial=0,final=0),
  fitid="noID",fit="not fitted yet", 
  Data=Data
  )) 
}




lenData<-function(profiles)
{
nData=0; 
for (i in 1:length(profiles)) 
    nData=nData+dim(profiles[[i]]$d)[1]; 
nData
}


simulateData<-function(model,fine=FALSE,init=FALSE) {
# create sdata (smooth/fine for plotting)if fine=1 and edata for SSE evaluation otherwise
options("warn"=-1)
times <- c(0,4 * 10^c(-1,0,2,7))
SSE=0
attach(model)
# pull items from attachment to save typing, save them to real object model
nms=names(Kparams)
if (Kj) {model$parmsTCC[paste("Kj",Z,sep="_")]= params[nms,"final"]  
#         model$parmsTCC[paste("Msk",Z,sep="_")]=params[nms,"present"]
} else  {
  bound=row.names(params)[params[,"constr"]!="none"]
  free=params[bound,"constr"]
  model$params[bound,"final"]=params[free,"final"]
  model$parmsTCC[paste("Kj",Z,sep="_")]= Kd2Kj(model$params[nms,"final"])
#  model$parmsTCC[paste("Msk",Z,sep="_")]=params[nms,"present"]
  }
model$parmsTCC["pRT"]=pRT=params["pRT","final"]
#print(params)
#print(model$parmsTCC)
detach(model) # reattach to transmit the changes made to model 
attach(model)
#print(parmsTCC)
for (p in names(Data)) {
   d=Data[[p]]$d
   nonzeros=Data[[p]]$nonzeros
   dyn.load(strn<-paste(wDir,"/models/",id,.Platform$dynlib.ext,sep="")) 
   icNames=c(names(d)[1],names(nonzeros));
   RT=nonzeros["R"]
   xVarData=d[,1]
   if (fine) {
       upperx=1.1*max(xVarData);
       lowerx=0.5*min(xVarData);
       delx=upperx/40;
       xvals=seq(lowerx,upperx,delx);
       } else  {
       xvals=xVarData; 
      } # end if fine
     dfr=data.frame(xvals,rep(0,length(xvals)));                     names(dfr)=names(d) 
     ndf=length(xvals)
     SS=data.frame(matrix(rep(0,ndf*nStateS,),nrow=ndf));                   names(SS)<-stateS
     chk=data.frame(matrix(rep(0,ndf*2*length(icNames)),nrow=ndf));  names(chk)<-c(icNames,paste(icNames,"Q",sep=""))
     initialStateTCC[]=rep(0,nAtomS)
     parmsTCC[1:nAtomS]=0 # these are overwritten where they are nonzero in the current figure 
     for (ii in 1:ndf)
       { icValues=c(xvals[ii],nonzeros)
         model$parmsTCC[paste(icNames,"T",sep="")]=icValues # note that pRT is built into RT
         e=try(out1 <- lsoda(initialStateTCC,times,"myderivs", model$parmsTCC, rtol=rtol,atol=atol, dllname=id)) 
         out=data.frame(out1); names(out)<-c("time",atomS)
         nout=dim(out)[1]
         lastSS=abs(out1[nout,2:(nAtomS+1)])
         SS[ii,] = fback(out[nout,2:(nAtomS+1)],parmsTCC)         
         chk[ii,icNames]=icValues
       }  # end loop through rows of dataframe
     chk[,(length(icNames)+1):(2*length(icNames))]=as.matrix(SS)%*%as.matrix(W[,icNames])
     dfr[,2]=monomerMass*((as.matrix(SS)%*%as.matrix(W[,"R"]^2)) + (1-pRT)*RT)/RT # W here includes Free R
     if (fine) {
               model[[p]]$schk=chk # s for smooth
               model[[p]]$sSS=SS
               model[[p]]$sd=dfr
               if (init) model[[i]]$isd=dfr
     } else {  model[[p]]$echk=chk 
               model[[p]]$eSS=SS
               model[[p]]$d=d
               model[[p]]$ed=dfr
               if (init) model[[p]]$ied=dfr
               res=(d-dfr)  /mean(d[,2])
               SSE=SSE+t(res[,2])%*%(res[,2])
               model[[p]]$res=res
               model[[p]]$res[,1]=d[,1]
              } # end not fine
#     detach(Data[[p]])
 dyn.unload(strn)
   }  # end for loop on p
 N=model$nData=lenData(Data)
 P=sum(params[,"opt"]) 
 detach(model) 
 if (init) {
 model$SSE$initial=SSE; 
 model$AIC$initial=N*log(SSE/N)+2*P + 2*P*(P+1)/(N-P-1)
 }
 if (!fine) {
 model$SSE$final=SSE; 
 model$AIC$final=N*log(SSE/N)+2*P + 2*P*(P+1)/(N-P-1)
 }
# print(model$SSE)
# print(model$AIC)
options("warn"=1)
 return(model)
 }



fopt<-function(pars,model)
{
 pars=exp(pars)
 model$params[names(pars),"final"]=pars
 model=simulateData(model,fine=FALSE) 
 model$SSE$final   # for a given model, number of params is fixed, so goal is to minimize SSE
}


fitModel<-function(model) {
  model=simulateData(model,fine=FALSE,init=TRUE)  # computes the initial SSE 
  p0=model$params[model$params[,"opt"],"initial"]
  names(p0)<-row.names(model$params)[model$params[,"opt"]]
  if ((model$nOptParams<-length(p0))>0) {                    
    p0=lapply(p0,log)
#    print(p0)
    if (model$nOptParams>1) opt=optim(p0,fopt,hessian=T,model=model) else
        opt=optim(p0,fopt,method="BFGS",hessian=T,model=model);       
#  print("here")
    attach(model)
    nOptp=length(tmp<-intersect(c("pRT"),names(p0)) ) # zero or one
    model$fitid=paste(paste(model$id,nZ,sep=""),nOptParams-nOptp,nOptp,sep=".")
    sg<-sqrt(opt$value/(nData-nOptParams))
    if (det(opt$hessian)>0) 
        {sig=sg*sqrt(diag(solve(opt$hessian/2)));model$hess=TRUE} else
        {sig=Inf;model$hess=FALSE} 
    upper=signif(opt$par+1.96*sig,3)
    lower=signif(opt$par-1.96*sig,3)
    point=signif(opt$par,3)
    CI=cbind(lower,point,upper)
    CI=exp(CI); opar<-exp(opt$par)
    model$params[names(opar),"final"]=opar
    detach(model)
    model=simulateData(model,fine=FALSE)  # computes the final SSE 
    cat("\n")
    model$CI=CI
    model=simulateData(model,fine=TRUE) # create smooth curve for publication
  } else   {#print("length p0=0 => nothing to fit"); 
            model$fit="nothing to fit"; 
            model=simulateData(model,fine=FALSE,init=TRUE); }
 # if (is.null(model$AIC$final)) model$AIC$final=100
 # if (is.na(model$AIC$final)) model$AIC$final=100
  model
}

lfitModel<-function(model) {   # to allow use with lapply on a list of models
     strns=c("present","absent")[(model$params[!model$params$opt,"final"]==Inf)+1]
     strns[strns=="present"]="fixed"
     strns[strns!="fixed"]="absent"
     strns[model$params[!model$params$opt,"constr"]!="none"]="constrained"
     t0 = Sys.time();      
     e=try(model<-fitModel(model))
     model$cpu<-difftime(Sys.time(),t0,units="mins")[[1]]
     model$report=rbind(subset(model$params,select=c("initial","final"),subset=opt),
                           subset(model$params,select=c("initial","final"),subset=(!opt)),
                           SSE=model$SSE,
                           AIC=model$AIC,
                           cpu=c(0,model$cpu))
     if (class(e)=="try-error") 
      {try(detach(model)); 
       model$fit="fit failed"
       model$report=cbind(model$report,confidenceInterval=c(rep("fit failed",dim(model$report)[1]-length(strns)-3),strns,c(rep("",2),"fit failed")))
       model$report["cpu","final"]=model$cpu       } else 
      if (model$nOptParams>0)   {   
      if (model$hess)
     model$report=cbind(model$report,confidenceInterval=c(sprintf("(%4.3f,%4.3f)",model$CI[,"lower"],model$CI[,"upper"]),
     strns,rep("",2),"fit succeeded")) else
     model$report=cbind(model$report,confidenceInterval=c(sprintf("(%4.3f,%4.3f)",model$CI[,"lower"],model$CI[,"upper"]),
     strns,rep("",2),"hessian singular")) 
      }
     else {model$report=cbind(model$report,confidenceInterval=c(strns,rep("",2),"nothing to fit"));model$report["cpu","final"]=model$cpu       }
     print(model$mid)
     print(model$report)
cat('\n\n')
     model
}

fitMS<-function(models,msid="unMS",cpus=1,ptype="",smart=0) {
if (sum(dir()=="results")==0) system("mkdir results")
library(snow)
mysum<-function(x) sum(x$params$opt)
myXAIC<-function(x) {#print(x$AIC)
if(is.null(x$AIC)) {x=list(); x$AIC$final=100}
                     x$AIC$final[1]}
nump=sapply(models,mysum)
print("nump=")
print(nump)
if (smart==1) { best=1e10; fmodels=NULL} 
t0 = Sys.time();      
if (ptype!="") {
 cl <<- makeCluster(cpus,type = ptype)
 clusterEvalQ(cl, library(odesolve))
 #print(ls(.GlobalEnv))
 clusterExport(cl,ls(.GlobalEnv))
 if (smart==1) {
  for (i in min(nump):max(nump))  {
    fmods=clusterApplyLB(cl,models[nump==i],lfitModel) 
    aic=sapply(fmods,myXAIC)
  #  aic=as.numeric(aic)
  #  aic[is.null(aic)]=100
    print(aic)
    nbest=min(aic,na.rm=TRUE)   
    fmodels=c(fmodels,fmods)
    if (nbest>best) break else best=nbest
   }
  } else fmodels=clusterApplyLB(cl,models,lfitModel) # else brute force across whole model space 
  stopCluster(cl)
} else fmodels=lapply(models,lfitModel) # else ptype = "" => single processor 
print(totTime<-difftime(Sys.time(),t0,units="mins"))
saic=sort(aic<-sapply(fmodels,myXAIC))
iBest=NULL
for (i in 1:length(saic)) iBest=c(iBest,which(aic==saic[i]))
smodels=fmodels[unique(iBest)]
makeHTML(smodels,msid,totTime,smart,cpus)
makeLaTeX(smodels,msid,totTime)
tops=lapply(smodels,"[[","report")
save(tops,smodels,saic,totTime,file=paste("results/",msid,".RData",sep=""))
smodels
} # fitsmart


makeHTML<-function(fmodels,fname,totTime,smart,cpus) {
.HTML.file=file(paste("results/", fname,".htm",sep=""),"wt")
cat("<html><h1>Fitted Model Space of",fmodels[[1]]$id,"</h1>",file = .HTML.file)
cat("<html><h2>Total CPU time using ",cpus," cpu(s) = ",totTime," minutes.</h2>",file = .HTML.file,sep="")
cat("<html><h2>Model space search method:  ",c("brute force","semi-intelligent")[smart+1],".</h2>",file = .HTML.file,sep="")
repS=c("<TABLE BORDER=2 CELLPADDING=4>","<TR><TH>Model</TH><TH>Parameter</TH><TH>Initial Value</TH><TH>Optimal
Value</TH><TH>Confidence Interval</TH></TR>")
nMS=length(fmodels)
for (i in 1:nMS) {
 attach(fmodels[[i]])
 nRows=dim(report)[1]
 repS=c(repS,sprintf("<TR><TH>%s</TH><TD>%s</TD><TD>%4.3f</TD><TD>%4.3f</TD><TD>%s</TD></TR>",mid,row.names(report)[1],report[1,"initial"],
  report[1,"final"],report[1,"confidenceInterval"]))
for (j in 2:(nRows))
 repS=c(repS,sprintf("<TR><TD>&nbsp;</TD><TD>%s</TD><TD>%4.3f</TD><TD>%4.3f</TD><TD>%s</TD></TR>",row.names(report)[j],report[j,"initial"],report[j,"final"],report[j,"confidenceInterval"]))
 detach(fmodels[[i]])
 } # loop on i
repS=c(repS,"</TABLE>")
writeLines(repS,con=.HTML.file)
cat("\n</html>", append = TRUE,file = .HTML.file)
close(.HTML.file)
} # makeHTML

mysb<-function(x) gsub("_","\\\\_",x)

makeLaTeX<-function(fmodels,fname,totTime) {
f1=file(paste("results/", fname,".tex",sep=""),"wt")
repS= "\\documentclass[10pt]{article}"    
repS=c(repS,"\\begin{document}")
repS=c(repS,"\\begin{tabular}{|c|c|c|c|c|} \\hline")
repS=c(repS,"Model & Parameter & Initial Value  & Optimal Value & Confidence Interval \\\\ \\hline")
nMS=length(fmodels)
for (i in 1:nMS) {
 attach(fmodels[[i]])
 nRows=dim(report)[1]
 repS=c(repS,sprintf("%s & %s & %4.3f & %4.3f & %s\\\\ \\hline",mid,mysb(row.names(report)[1]),report[1,"initial"],
 report[1,"final"],report[1,"confidenceInterval"]))
for (j in 2:(nRows))
 {repS=c(repS,sprintf("& %s & %4.3f & %4.3f & %s\\\\ \\hline", mysb(row.names(report)[j]),report[j,"initial"],
  report[j,"final"], mysb(report[j,"confidenceInterval"])))
 }
 detach(fmodels[[i]])
 } # loop on i
repS=c(repS,"\\end{tabular}")
repS=c(repS,"\\end{document}")
writeLines(repS,con=f1)
close(f1)
} # makeLaTeX



