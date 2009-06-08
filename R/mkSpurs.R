`mkSpurs` <-
    function(g,state=list(globMdlIndex=0,globCmbIndex=0,relCmbIndex=0,config=NULL),
        maxnKjPs=NULL,maxTotalPs=NULL,minTotalPs=NULL,batchSize=500,
        doTights=FALSE,atLeastOne=TRUE,atLeastOneOfEach=FALSE,KIC=1,kIC=1,m1=-90,p=-1,
        forceM1=FALSE,forceP=FALSE) {
  if (p>0) pRows=TRUE else pRows=FALSE
  #if (pRows&(batchSize%%100 !=0)) print("please make batchSize a multiple of 100")
  if (is.null(maxnKjPs)) maxnKjPs=g$nZ
  if (!g$activity&!is.null(maxTotalPs)) maxnKjPs=maxTotalPs
  # *************** cnfgs list made below *******************************	
  if (forceM1) maxnKjPs=maxnKjPs-1
  if (forceP) maxnKjPs=maxnKjPs-1
  ri=state$relCmbIndex
  nbP=length(state$config)
  if ((ri==0)&(nbP>0)) nbP=nbP+1  # i.e. last one ended on a dime at the end of a block 
  lastCompleted=1 # number of parameters
  cnfgs=list(NULL)
  i=0
  flag=0
  maxReached=FALSE
  if (nbP==0) {cnfgs[[1]]=NULL; i=1; nbP=1}
#	while(i<batchSize) {
  while((i<batchSize)&(nbP<=maxnKjPs)) {
    X=combn(g$nZ,nbP)
#		print(X)
    nbC<-dim(X)[2]
#		print(nbC)
#		print("in outer")
#		cat("bs1 = ",batchSize,"nbC = ",nbC,"ri = ",ri,"i = ",i,"\n")
    while ((ri<nbC)&(i<batchSize)) {
#			cat("inner: bs = ",batchSize,"nbC = ",nbC,"ri = ",ri,"i = ",i,"\n")
      ri=ri+1;i=i+1
      cnfgs[[i]]=X[,ri]; 
#			print(cnfgs)
    }
    if (ri==nbC){lastCompleted=nbP; nbP=nbP+1;ri=0} 
  }
#	print(i)
  if (nbP>maxnKjPs) maxReached=TRUE
  batchSize=i  # 
#	print(cnfgs)
  
  totN=2^g$nZ
  if(batchSize>totN) { print("Batch size cut to full model space");
    batchSize = totN}  
  if (state$globCmbIndex + batchSize > totN) batchSize = totN-state$globCmbIndex
  mbatchSize=batchSize # m for mat, see below; this is the default
  if (is.null(state$config)&(doTights)) 
  {
    if(batchSize<(g$nZ+1) ) 
    { 
      print("Batch size increased to all single edges");
      batchSize = g$nZ+1 
    }  # increase to first comb matrix
    mbatchSize = batchSize + g$nZ
  } # make space for tights
#  if (pRows&!forceP) mbatchSize=2*mbatchSize  # double all if p variable included
  if (pRows) mbatchSize=2*mbatchSize  # double all if p variable included
  cat("batchSize = ",batchSize,"\n")
  state$relCmbIndex=ri
  state$config=cnfgs[[length(cnfgs)]]
# *************** cnfgs list made above, use it below *******************************	
  mat=matrix(rep(Inf,mbatchSize*(g$nZ+3)),ncol=g$nZ+3)
#  mat[,g$nZ+3]=(state$globMdlIndex+1):(state$globMdlIndex+mbatchSize) 
  mat[,g$nZ+2]=-abs(p) # was default to -p = fixed p, then overwrite with +p if variable
#	print(mat)
  mi=1
  i=1 # reset to create the output version of this variable
  while (mi <= mbatchSize){
    nbP=length(cnfgs[[i]])
    mat[mi,1]=nbP
    if (nbP>0) {
      mat[mi,cnfgs[[i]]+1]=1;	
    }
    mi=mi+1
    if ((nbP==1)&(doTights)) {
      mat[mi,cnfgs[[i]]+1]=0;	mat[mi,1]=nbP-1; mi=mi+1
      if (pRows){mat[mi,cnfgs[[i]]+1]=0; mat[mi,g$nZ+2]=p; mat[mi,1]=nbP; mi=mi+1}
    }
    if (pRows){mat[mi,cnfgs[[i]]+1]=1; mat[mi,g$nZ+2]=p; mat[mi,1]=nbP+1; mi=mi+1}
    i=i+1
  }
  state$globCmbIndex=state$globCmbIndex+i-1  # i is the index through graph topologies and thus cnfgs of the combn matrix
  #	print(mat)
  KS=mat
  spurKs=mat
# nKs
  KS[KS==Inf]<-"I" ;KS[,2:(g$nZ+1)][mat[,2:(g$nZ+1)]==1]<-"J"
  KS[,g$nZ+2][mat[,g$nZ+2]>0]  = "p"
  KS[,g$nZ+2][mat[,g$nZ+2]<0] = ""
  KS[,g$nZ+3] = ""
  # KS=cbind(KS,rep("S",dim(KS)[1]))
#print(KS)
  rownames(spurKs)<-apply(KS[,-1,drop=FALSE],1,paste,collapse="")
  colnames(spurKs)<-c("nParams",g$Z,"p","indx")
  spurKs=data.frame(spurKs)
#  print(spurKs)
  if(atLeastOne) { #then model must have at least one term of maximum size oligo
    zor1=lapply(spurKs,"%in%",c(0,1)) 
    newW=g$W[-(1:(g$nAtomS-1)),g$hubChar]  # leave one atom row in as spacer for nParams column of chunk 
    colNums=which(newW==max(newW))
    Irows<-apply(as.matrix(as.data.frame(zor1[colNums])),1,sum)>0
    spurKs=spurKs[Irows,] # grab rows with at least one complex of maximum size 
    mbatchSize=mbatchSize-sum(!Irows)
    spurKs[,"indx"]=(state$globMdlIndex+1):(state$globMdlIndex+mbatchSize)  
    uniW=unique(newW)
    uniW=uniW[uniW<max(uniW)]
    uniW=uniW[uniW>1]
#    print(uniW)
  }
  cat("mbatchSize = ",mbatchSize,"\n")

  if(atLeastOneOfEach) { #then model must have at least one term of each oligo size
    for (i in uniW) {
      zor1=lapply(spurKs,"%in%",c(0,1)) 
#      newW=g$W[-(1:(g$nAtomS-1)),g$hubChar]  # leave one atom row in as spacer for nParams column of chunk 
      colNums=which(newW==i)          # i.e. leave in one atom row to make this index right
      Irows<-apply(as.matrix(as.data.frame(zor1[colNums])),1,sum)>0
      spurKs=spurKs[Irows,] # grab rows with at least one complex of maximum size 
      mbatchSize=mbatchSize-sum(!Irows)
      spurKs[,"indx"]=(state$globMdlIndex+1):(state$globMdlIndex+mbatchSize)  
    }
  }
  cat("mbatchSize = ",mbatchSize,"\n")
  
  
  ncols=dim(spurKs)[2]
  KCols=spurKs[,2:(g$nZ+1)]
  KCols[KCols==1]=KIC
  spurKs[,2:(g$nZ+1)]=KCols
#  mainCols=spurKs[,2:(ncols-2)]
#  mainCols[mainCols==1]=IC
#  spurKs[,2:(ncols-2)]=mainCols
  spurs=spurKs
  if (g$activity) {
    kCols=KCols
    kCols[KCols==KIC]=kIC
    kCols[KCols==Inf]=0
    names(kCols)<-paste("k",g$Z,sep="")
    spurs=cbind(spurKs[,1:(g$nZ+1)],kCols,spurKs[,(g$nZ+2):(g$nZ+3)])
    Mrnms=lapply(uni<-unique(spurs[,"nParams"]),mkSingleThread,g)
    names(Mrnms)<-paste("p",uni,sep="")
#    print(Mrnms)    
    condens<-function(x) apply(x,1,paste,collapse="")
    Crnms=lapply(Mrnms,condens)
    print(Crnms)
    lens=lapply(Crnms,length)
    print(lens)
    prepkPs<-function(x) apply(x,1,unique)
    countkPs<-function(x) sapply(x,length)
    unikPs=lapply(Mrnms,prepkPs)
    print(unikPs)
    nkPs=lapply(unikPs,countkPs)
    print(nkPs)
    rowList=NULL
    for (i in 1:dim(spurs)[1]){
      tmp=spurs[i,,drop=FALSE]
      pn=paste("p",tmp$nParams,sep="")
      knms=Crnms[[pn]];
      for (j in 1:length(knms)) {
        tmp0=tmp        
        Knm=rownames(tmp)
        knm=strsplit(Knm,split=NULL)[[1]]
        knmsj=strsplit(knms[j],split=NULL)[[1]]
        pntrs=which(knm=="J")
#        print("pntrs are")
#        print(pntrs)
#        print("knm is")
#        print(knm)
#        print("knmsj is")
#        print(knmsj)
        knm=rep("0",length(knm))
        if (length(pntrs)>0) for (kk in 1:length(pntrs)) knm[pntrs[kk]]=knmsj[kk]
        knm=paste(knm,collapse="")
        rownames(tmp0)<-paste(Knm,knm,sep=".")
        tmp0$nParams=tmp0$nParams+nkPs[[pn]][j]
        rowList=rbind(rowList,tmp0)
      }
    } 
    rowList=rowList[order(rowList$nParams),]
    rowList$indx=1:(dim(rowList)[1])
    print(rowList)
    spurs=rowList
  }
  if (!g$activity) keqs=NULL else keqs=mkEq(g,spurs,activity=TRUE)
  spurs=mkm(spurs,m1) # mkm (make mass column) is defined at the end of the mkGrids file
  # now trim it down to just those with up to maxTotalPs parameters
  if (!is.null(maxTotalPs)) spurs=spurs[spurs$nParams<=maxTotalPs,] 
  if (!is.null(minTotalPs)) spurs=spurs[spurs$nParams>=minTotalPs,] 
  if (forceM1) spurs=spurs[spurs$m1>0,] 
  if (forceP) spurs=spurs[spurs$p>0,] 
  state$globMdlIndex=state$globMdlIndex+dim(spurs)[1]
  chunk=spurs
  I=order(chunk$nParams)
  chunk=chunk[I,]
  chunk$indx=1:(dim(chunk)[1])
  list(chunk=chunk,state=state,maxReached=maxReached,lastCompleted=lastCompleted,keqs=keqs)
}

#incConfig<-function(state,N) {
#	nbP=length(state$config)
#	print(nbP)
#	if (nbP==0) state$config=1  
#	if (nbP==1) 
#		if(state$config<N) state$config=state$config+1 else state$config=c(1,2)  
#	if (nbP>1) 
#		if(state$config[nbP]<N) state$config[nbP]=state$config[nbP]+1 else 
#		if ((prod(diff(state$config))==1)& state$config[nbP]==N) state$config=1:(nbP+1) else
#		{ii=max(which(diff(state$config)>1));print(ii)
#			state$config[ii]=state$config[ii]+1
#		}
##	print(incPos(state$config,g$nZ,nbP))  	
#	state
#}
