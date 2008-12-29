`mkg` <-
    function(strct,hubChar="R",monomerMass=90,TCC=FALSE)  
{ 
  lets=c(LETTERS,letters)
  nums=paste(0:9)
#    browser()
  mkgC <- function(g) {
#      setwd(paste(g$wDir,"models",sep="/"))
    if (sum(dir()=="models")==0) system("mkdir models")
    setwd("models")
#		print("in mkgC")
    ofile=file(fn<-paste(g$id,".c",sep=""),"wt")
    cat(sprintf("/* compile with R CMD SHLIB %s.c or (on Windows) Rcmd SHLIB %s.c */\n",g$id,g$id),file=ofile)
    cat("#include <R.h> /* gives F77_CALL through R_ext/RS.h */\n\n",file=ofile)
    pnames=names(g$parmsTCC)
    cat(sprintf("static double parms[%d];\n",length(pnames)),file=ofile)
    for (i in 1:length(pnames))
      cat(sprintf("#define %25s  parms[%d]\n",pnames[i],i-1),file=ofile)
    cat(sprintf("void %s(void (* odeparms)(int *, double *))\n",id),file=ofile)
    cat(sprintf("{ int N=%d;\n",length(pnames)),file=ofile)   
    cat(" odeparms(&N, parms);",file=ofile); cat("}\n",file=ofile)
    cat("\nvoid myderivs(int *neq, double *curtime, double *statevec, double *ODERHSvec)\n",file=ofile)
    cat("{",file=ofile); 
#		print(g$specieS)
    cat(sprintf(" double %s;\n",paste(g$specieS,sep='',collapse=',')),file=ofile)
#		print(g$atomS)
    for (i in 1:g$nAtomS) cat(sprintf("  %s = statevec[%d];",g$atomS[i],i-1),file=ofile)
    print(g$Z)
    cat("\n",file=ofile)
    for (Zj in g$Z) cat(sprintf("%10s = %s/Kj_%s;\n",Zj,paste(g$reactantS[[Zj]],collapse="*"),Zj),file=ofile)
    for (i in 1:g$nAtomS)
    {
      if (i==1)      cat(sprintf("ODERHSvec[%d] = p*%sT ",i-1,g$atomS[i]),file=ofile) else
        cat(sprintf("ODERHSvec[%d] =     %sT ",i-1,g$atomS[i]),file=ofile)
      for (j in 1:g$nSpecieS) 
        if (j<=g$nAtomS) 
        {if (g$W[g$specieS[j],g$atomS[i]]==1) 
            cat(sprintf("-%s",g$specieS[j]),file=ofile)}
        else         
        if (g$W[g$specieS[j],g$atomS[i]]==0)
          cat("              ",file=ofile)  else
        if (g$W[g$specieS[j],g$atomS[i]]==1)
          cat(sprintf(" - %s  ",g$specieS[j]),file=ofile) else
          cat(sprintf(" - %d*%s",g$W[g$specieS[j],g$atomS[i]],g$specieS[j]),file=ofile) 
      cat(";\n",file=ofile)
    }
    cat("}\n",file=ofile)
    close(ofile)
# note: cpu time savings of putting in blanks for multiplications by 0 and 1 was neglibible => something else is making this 2x slower than 
# hard coding without masks.  
    if (.Platform$OS.type=="windows") system(sprintf("Rcmd SHLIB %s.c\n",g$id)) else  system(sprintf("R CMD SHLIB %s.c\n",g$id))
#    system(sprintf("mv %s%s  bin\n",id,.Platform$dynlib.ext))
    unlink(sprintf("%s.o",g$id))
    unlink(sprintf("%s.d",g$id))
    g$code=readLines(fn)
#    setwd(g$wDir)
    setwd("..")
    g
  }
  
  mkZ <- function(g) {  # Adds the R function fback to the g object. 
#        if (sum(dir()=="models")==0) system("mkdir models")
#        setwd("models")
#        #    setwd(paste(g$wDir,"models",sep="/"))
#        ofile=file(fn<-paste(g$id,".r",sep=""),"wt")
#        cat("fback<-function(atoms, parmsTCC) {\n",file=ofile)
#        for (i in 1:g$nAtomS)
#            cat(sprintf("%10s = atoms[%d];\n",g$atomS[i],i),file=ofile)
#        for (Zj in g$Z)
#            cat(sprintf("%10s = %s/parmsTCC[\"Kj_%s\"];\n",Zj,paste(g$reactantS[[Zj]],collapse="*"),Zj),file=ofile)
#        cat("c(",file=ofile)
#        for (i in 1:g$nSpecieS)
#            cat(sprintf("%s=%s%s",g$specieS[i],g$specieS[i],ifelse(i<length(g$specieS),",",")\n}\n\n")),file=ofile)
#        close(ofile)
#        source(fn,local=TRUE)
#        unlink(fn)
#    source(paste(g$wDir,"/models/",fn,sep=""),local=TRUE)
#    unlink(paste(g$wDir,"/models/",fn,sep=""))
#        g$fback=fback
    strn="fback<-function(atoms, parmsTCC) {\n"
    for (i in 1:g$nAtomS)
      strn=paste(strn,sprintf("%10s = atoms[%d];\n",g$atomS[i],i),sep="")
    for (Zj in g$Z)
      strn=paste(strn,sprintf("%10s = %s/parmsTCC[\"Kj_%s\"];\n",Zj,paste(g$reactantS[[Zj]],collapse="*"),Zj),sep="")
    strn=paste(strn,"c(",sep="")
    for (i in 1:g$nSpecieS)
      strn=paste(strn,sprintf("%s=%s%s",g$specieS[i],g$specieS[i],ifelse(i<length(g$specieS),",",")\n}\n\n")),sep="")
    g$fback=eval(parse(text=strn))
    #		setwd(g$wDir)
#        setwd("..")
    g
  }
  
  mkRP <- function(g) {
#        if (sum(dir()=="models")==0) system("mkdir models")
#        setwd("models")
#        ofile=file(fn<-paste("rp",g$id,".r",sep=""),"wt")
#        cat("frp<-function(parmsTCC,kis) {\n",file=ofile)
#        cat("        E0 = parmsTCC[\"RT\"]\n         S = parmsTCC[\"ST\"]\n",file=ofile)
#        for (Zj in g$Z)
#            cat(sprintf("%10s = %s/parmsTCC[\"Kj_%s\"]\n",Zj,paste(strsplit(gsub("R","",Zj,fixed="T"),split="")[[1]],collapse="*"),Zj),file=ofile)
#        cat("#free R cancels from num and denom, so in here only to label complex assoc with term\n",file=ofile)
#        cat("denom=1",file=ofile)
#        for (i in 1:g$nZ)
#            cat(sprintf("+%s%s",g$Z[i],ifelse(i<g$nZ,"","\n")),file=ofile)
#        cat("num=E0*(",file=ofile)
#        for (i in 1:nZ) # look for a built in character within string counter to replace the 2 maker
#            cat(sprintf("%d*kis[\"k%s\"]*%s%s",length(which(strsplit(g$Z[i],split="")[[1]]=="S")),
#                            g$Z[i],g$Z[i],ifelse(i<g$nZ,"+",")\n")),file=ofile)
#        cat("EY=num/denom\nEY\n}\n\n",file=ofile)
#        close(ofile)
#        source(fn,local=TRUE)
#        unlink(fn)
#    source(paste(g$wDir,"/models/",fn,sep=""),local=TRUE)
#		unlink(paste(g$wDir,"/models/",fn,sep=""))
    strn="frp<-function(parmsTCC,kis) {\n"
    strn=paste(strn,"        E0 = parmsTCC[\"RT\"]\n         S = parmsTCC[\"ST\"]\n",paste="")
    for (Zj in g$Z)
      strn=paste(strn,sprintf("%10s = %s/parmsTCC[\"Kj_%s\"]\n",Zj,paste(strsplit(gsub("R","",Zj,fixed="T"),split="")[[1]],collapse="*"),Zj),sep="")
    strn=paste(strn,"denom=1",sep="")
    for (i in 1:g$nZ)
      strn=paste(strn,sprintf("+%s%s",g$Z[i],ifelse(i<g$nZ,"","\n")),sep="")
    strn=paste(strn,"num=E0*(",sep="")
    for (i in 1:nZ) # look for a built in character within string counter to replace the 2 maker
      strn=paste(strn,sprintf("%d*kis[\"k%s\"]*%s%s",length(which(strsplit(g$Z[i],split="")[[1]]=="S")),
              g$Z[i],g$Z[i],ifelse(i<g$nZ,"+",")\n")),sep="")
    strn=paste(strn,"EY=num/denom\nEY\n}\n\n",sep="")
    g$frp=eval(parse(text=strn))
#		setwd(g$wDir)
#        setwd("..")
    g
  }
  
  conv2long<-function(react) {
    tmp=""
    i=1
    while (i < length(react)) {
      if(react[i]%in%lets) curLet=react[i]
      i=i+1
      tmpn=""
      while (react[i]%in%nums) {tmpn=paste(tmpn,react[i],sep=""); i=i+1}
      tmp=paste(tmp,paste(rep(curLet,as.numeric(tmpn)),collapse=""),sep="")
    }
    tmp
  }
  
  getnums<-function(react) {
    vals=NULL
    nms=NULL
    tmp=""
    i=1
    while (i < length(react)) {
      if(react[i]%in%lets) nms=c(nms,react[i])
      i=i+1
      tmpn=""
      while (react[i]%in%nums) {tmpn=paste(tmpn,react[i],sep=""); i=i+1}
      vals=c(vals,as.numeric(tmpn))
    }
    names(vals)<-nms
    vals
  }
  
  oneLess<-function(react){
    i = length(react)
    tmpn=""
    while (react[i]%in%nums) {tmpn=paste(react[i],tmpn,sep=""); i=i-1}
    strn=paste(react,collapse="")
    basK=substr(strn, 1, i)
    basn=as.numeric(tmpn)-1
    endK=substr(strn, i,i)
    paste(basK,basn,"_",endK,sep="")
  }
  if (length(strct$heads)==length(strct$sites[[1]])) {
    # check to see if first head is root and if so, eliminate it
    head0=strsplit(strct$heads[1],NULL)[[1]];
#    print(head0)
    head0N=getnums(head0)
#    print(head0N)
    if ((head0N[1]==1)&(head0N[2]==0)) strct$heads=strct$heads[-1]
  }
  Z=unlist(strct)
  nZ=length(Z);
  reacts=strsplit(Z,NULL);
  names(reacts)<-Z
  atomS=union(hubChar,setdiff(unlist(reacts),paste(0:9,sep=""))) # union makes sure central/hub protein is first
  nAtomS=length(atomS)
  id=paste(atomS,collapse="")
#	print(id)
  reactantS=strsplit(sapply(reacts,conv2long),NULL)
  specieS=c(atomS,Z)
  nSpecieS=length(specieS)
  names(specieS)<-specieS
#	print(specieS)  
  reactaNts=sapply(reacts,getnums)
#	print(reactaNts)
  eye=diag(nAtomS)
  rownames(eye)=atomS
  colnames(eye)=atomS
  W=rbind(eye,t(reactaNts))
  W=as.data.frame(W)
#	print(W)
#	print(class(W))
  hdS=strct$heads
#	print(hdS)
  tmp=NULL
  for (i in 1:length(hdS)) tmp=c(tmp,which(Z==hdS[i])) # convert node names to indices
  hds=tmp
  print(hds)
  KdS=sapply(reacts,oneLess) # binary Kds 
  KdS[hds]=Z[hds] #of non-head nodes
  initialStateTCC=rep(0,nAtomS) 
  names(initialStateTCC)<-atomS
#  parmsTCC=c(rep(1,nAtomS),1,rep(100,nZ))  # this won't work because spurs matrix ICs overwrite this
  parmsTCC=c(rep(1,nAtomS),1,rep(1,nZ))
  pnamesTCC=c(paste(atomS,"T",sep=""),"p",paste("Kj",Z,sep="_"))
  names(parmsTCC)<-pnamesTCC
  gObj=list(id=id,
      Z=Z,nZ=nZ,
      atomS=atomS, nAtomS=nAtomS,
      specieS=specieS, nSpecieS=nSpecieS,
      hubChar=hubChar, 
      reactantS=reactantS,
      strct=strct,
      W=W,
      KdS=KdS,
      hdS=hdS,hds=hds,
      monomerMass=monomerMass,
#			wDir=sub("C:","",getwd()),
      TCC=TCC,
      sstime = 1e6,rtol=1e-5,atol=1e-7,
      parmsTCC=parmsTCC,
      initialStateTCC=initialStateTCC
  ); 
  gObj=mkZ(gObj)
  library(odesolve); 
  if (TCC) { gObj=mkgC(gObj); testgC(gObj)} else gObj=mkRP(gObj)
  gObj
}


