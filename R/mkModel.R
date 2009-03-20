`mkModel` <-
    function(g,mid,d=NULL,Kjparams=NULL,Kdparams=NULL,Keq=NULL,Kd2KjLst=NULL,pparams=c(p=1i),kparams=NULL,keq=NULL,
        tightLogic=TRUE,indx=NULL,nParams=NULL)  
{ 
# Numbers on the left half plane or imaginary axis are fixed to their (complex number) modulus. Thus, either purely imaginary numbers or 
# negative numbers can be used as initial parameter values as a conventient way to indicate that the parameter values are to be fixed. 
# 
#	One of Kj or Kd should be NULL. The other then should be used
#	print(Kdparams)
  Kjparams=unlist(Kjparams);pparams=unlist(pparams);kparams=unlist(kparams)
  if (is.null(Kjparams)&is.null(Kdparams)) {print("One of Kjparams or Kdparams must be specified."); return(0)}
  if (!is.null(Kjparams)&!is.null(Kdparams)) {print("Only one of Kjparams or Kdparams can be specified."); return(0)}
  if (!is.null(Kjparams)) Kj=TRUE else Kj=FALSE
  if(Kj) Kparams=Kjparams else Kparams=unlist(Kdparams)
#   print(Kparams)
#   print(pparams)
#   print(kparams)
  vparams=c(Kparams,pparams,kparams) 
#  print(vparams)
  nSysParams=length(Kparams);
  nParams=length(vparams)
#  print(vparams)
  
  oneMore<-function(binKdnm){
    KdChars=strsplit(binKdnm,split="")[[1]]
    i = length(KdChars)
    KdChars=KdChars[-c((i-1),i)]
    tmpn=""
    i=i-2
    nums=paste(0:9)
    while (KdChars[i]%in%nums) {tmpn=paste(KdChars[i],tmpn,sep=""); i=i-1}
    basK=substr(binKdnm, 1, i)
    basn=as.numeric(tmpn)+1
    paste(basK,basn,sep="")
  }
  
  params = data.frame(initial=Mod(vparams),final=Mod(vparams),opt=(Re(vparams)>0),  # note negative treated like imaginary here
      constr=rep("none",nParams),stringsAsFactors=FALSE)
#  print(params)
  if (!is.null(Keq)) {
    params[names(Keq),"constr"]=Keq
    params[names(Keq),"opt"]=FALSE 
  }
  if (!is.null(keq)) {
    params[names(keq),"constr"]=keq
    params[names(keq),"opt"]=FALSE 
  }
  params[Mod(vparams)==Inf,"opt"]=FALSE
  params[Mod(vparams)==0,"opt"]=FALSE
  
#  print(params)
  
# switch between these two for numeric versus logical
  if (tightLogic)
    params[(Mod(vparams)==0)&((1:nParams)<=nSysParams),c("initial","final")]=0  else 
  {
    eps=.0001
    params[(Mod(vparams)==0)&((1:nParams)<=nSysParams),c("initial","final")]=eps
    Kparams[(Mod(Kparams)==0)]=eps
  }
  
#	codeS=NULL  
  codeS="default"
  bin2dec <- function(x)  sum(x * 2^(rev(seq(along=x)) - 1))  # by Burt Gunter, uses c(1,0,1,1)=11 format
  if(!is.null(Kd2KjLst)) 
    if(length(Kd2KjLst)>1) { 
      midChars=strsplit(mid,split="")[[1]]  # model ID characters split into a character vector (first of list)
      KjIndices=which(midChars=="H")   # indices of all heads
      KjIndices=KjIndices[KjIndices>max(g$hds)]  # only want indices of inserted heads
      oldKdnms=row.names(params)[KjIndices]
      newKjnms=sapply(oldKdnms,oneMore)
      row.names(params)[KjIndices]<-newKjnms
      g$KdS[KjIndices]=newKjnms
      names(Kparams)[KjIndices]<-newKjnms
      nSites=length(g$threadsWithinSites) # rows in g$dfThreads
      nOligos=length(g$threadsWithinSites[[1]])  # cols in g$dfThreads
      pThreads=unlist(g$threadsWithinSites[-nSites]) # last filled site never inserted
      pntrs=NULL
#			print(pThreads)
      codes=rep(0,length(pThreads))
      for (i in 1:length(pThreads)){
#				print(g$threads[[i]])
#				print(g$threads[[i]]$nodes)
        pntrs=c(pntrs,max(g$threads[[i]]$nodes))
      }
      codes[midChars[pntrs]=="H"]=1
#      print(codes)
      currInfs=matrix(codes,byrow=TRUE,nrow=nSites-1)
#      print(currInfs)
      codes=apply(currInfs,2,bin2dec)
#      print(codes)
      codeS=paste(codes,collapse="")
      print(codeS)
#			ncode=paste(strsplit(codeS,split=NULL)[[1]][-1],collapse="") # this was to strip the monomer code
#			print(ncode)
#			codeS=ncode
    }
  
  typeYD=NULL
  posY=NULL
  posReactantsD=NULL
  if (!is.null(d)) {
    nms=names(d)
    if (length(out<-which(nms=="v"))>0) {typeYD="v";posY=out}
    if (length(out<-which(nms=="m"))>0) {typeYD="m";posY=out}
    posReactantsD=rep(0,g$nAtomS)
    for (i in 1:g$nAtomS) posReactantsD[i]=which(nms == paste(g$atomS[i],"T",sep=""))
  }
  c(g,list(mid=mid,params=params,Kparams=Kparams,Kj=Kj,codeS=codeS,Kd2Kj=Kd2KjLst[[codeS]],
          fitS="not fitted yet", typeYD=typeYD,nParams=nParams,
          posY=posY,posReactantsD=posReactantsD,d=d,indx=indx   # position of first reactant is position of central protein
      )) 
}

