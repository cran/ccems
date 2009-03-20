`mkGrids` <-
    function(g,maxTotalPs=NULL,pRows=FALSE,contig=TRUE,atLeastOne=TRUE,IC=1,kIC=1,fullGrid=FALSE) {
  # This function makes a chunk of grid graph models. It is less advanced than its counterpart mkSpurs because I still need a way
  # of going through this model space systematically with increasing numbers of parameters. Thus, it does 
  # not have state inputs and output. Note that it is implicit that R is always a thread head too
  
  
  mkCurtain<-function(mat0,g,iSite)
  # This function generates a curtain of j-mers for the current site.  
  # When there are more than two threads Infs come in from the left and right
  {
    threadNumbs=g$threadsWithinSites[[iSite]]
    mat=mat0
    jj=1
    for (patchLen in 2:length(threadNumbs)) {
      if (!contig) 
      {M=combn(threadNumbs,patchLen)
        for (k in 1:dim(M)[2]) {
          nodePositions=NULL
          currThreads=t(M[,k,drop=FALSE])
          for (kk in 1:length(currThreads)) nodePositions=c(nodePositions,g$threads[[currThreads[kk]]]$nodes)
          mat[jj,nodePositions]=g$threads[[min(currThreads)]]$let
          mat=rbind(mat,mat0)
          jj=jj+1
        } # k loop
      }  else # contiguous patch
        for (startThread in min(threadNumbs):(max(threadNumbs)-patchLen+1)) {
#         cat("\nstartThread=",startThread," patchLen=",patchLen,"\n")
          nodePositions=NULL
          for (kk in startThread:(startThread+patchLen-1)) nodePositions=c(nodePositions,g$threads[[kk]]$nodes)
          mat[jj,]=mat0 # reset in case not properly reset
          mat[jj,nodePositions]=g$threads[[startThread]]$let  # fill whole band with start thread label
          mat=rbind(mat,mat1<-mat[jj,])
          jj=jj+1
          #         print(mat)
          if (1) { #all combs of placing the I's
            outside=setdiff(threadNumbs,startThread:(startThread+patchLen-1))
            nOuts=length(outside)
#            print(outside)
#            cat("\nstartThread=",startThread," patchLen=",patchLen,"nOuts=",nOuts,"outside",outside,"\n")
            if (nOuts>1) 
              for (iOuts in 1:nOuts) { 
                X=combn(outside,iOuts)
#               print(X)
                for (kOuts in 1:dim(X)[2]) {
                  currOuts=t(X[,kOuts,drop=FALSE])
                  nodePositions=NULL
#                 maxes=NULL
                  for (kkOuts in 1:length(currOuts)) {
                    nodePositions=c(nodePositions,g$threads[[currOuts[kkOuts]]]$nodes)
#                   maxes=c(maxes,max(threads[[currOuts[kkOuts]]]$nodes))
                  }
                  mat[jj,nodePositions]="I"
                  if((iSite==1)&(max(currOuts)>1)) {
                    hdOuts=currOuts[currOuts!=1]-1
                    nHdOuts=length(hdOuts)
#                   cat("\nstartThread=",startThread," patchLen=",patchLen,"nHdOuts=",nHdOuts,"hdOuts",hdOuts,"\n")
                    matH=mat[jj,]
                    if (nHdOuts>1) 
                      for (iHdOuts in 1:nHdOuts) {
                        H=combn(hdOuts,iHdOuts)
#                       print(H)
                        for (kHdOuts in 1:dim(H)[2]) {
                          currHdOuts=t(H[,kHdOuts,drop=FALSE])
                          mat=rbind(mat,matH)
                          jj=jj+1 
                          mat[jj,currHdOuts]="I"
                        }
                      } else  {# only one hd taken out
                      mat=rbind(mat,matH)
                      jj=jj+1 
#                   mat[jj,maxes]="B"
                      mat[jj,hdOuts]="I"
                    }
                  }
                  if ((iOuts==nOuts)&(kOuts == dim(X)[2])) mat=rbind(mat,mat0) else mat=rbind(mat,mat1)
                  jj=jj+1
                } # kOuts loop
              } else if (nOuts==1) {
              mat[jj,g$threads[[outside]]$nodes]="I"
              if((iSite==1)&(outside>1)) {
#               mat[jj,max(threads[[outside]]$nodes)]="B"
                mat=rbind(mat,mat[jj,])
                jj=jj+1 
                mat[jj,outside-1]="I"
              }
              mat=rbind(mat,mat0)
              jj=jj+1
            }
          } else {  # do with I coming in from left and right
            leftGap=startThread-min(threadNumbs)
            rightGap=max(threadNumbs) - (startThread + patchLen -1)  
#           if (leftGap==0) mat=rbind(mat,mat0) else mat=rbind(mat,mat[jj,])
            if (leftGap>0) 
              for (iLeft in 1:leftGap) {
#               cat("\nstartThread=",startThread," patchLen=",patchLen,"leftGap=",leftGap,"iLeft=",iLeft,"\n")
                mat[jj,g$threads[[iLeft-1+min(threadNumbs)]]$nodes]="I"
#             matL=mat[jj,]
                mat=rbind(mat,matL<-mat[jj,])
                jj=jj+1
                if (rightGap>0) 
                  for (iRight in 1:rightGap) {
#                 cat("\nstartThread=",startThread," patchLen=",patchLen,"leftGap=",leftGap,"iLeft=",iLeft,"rightGap=",rightGap,"iRight=",iRight,"\n")
                    mat[jj,g$threads[[max(threadNumbs)-iRight+1]]$nodes]="I"
                    if (iRight==rightGap) mat=rbind(mat,matL) else mat=rbind(mat,mat[jj,])
                    jj=jj+1
                  }
              }
            if ((leftGap==0)&(rightGap>0)) {
              for (iRight in 1:rightGap) {
#             cat("\nstartThread=",startThread," patchLen=",patchLen,"leftGap is 0 rightGap=",rightGap,"iRight=",iRight,"\n")
                mat[jj,g$threads[[max(threadNumbs)-iRight+1]]$nodes]="I"
                if (iRight==rightGap) mat=rbind(mat,mat0) else mat=rbind(mat,mat[jj,])
                jj=jj+1
              }
            }
          }
        } # loop on startThread
    } # loop on patchLen 
    mat[jj,]=mat0  # last row is always mat0, even if all block sizes are 1 => redundant with spur graph
    list(mat=mat,jj=jj)
  }
#  END mkCurtain function definition
  
  
  spurKOThreads<-function(mat,mat0,g,iSite,jj) {
    # This function takes the full block model and generates spur like knock downs of threads (like nodes)
  threadNumbs=g$threadsWithinSites[[iSite]]
    for (iBspurs in 1:length(threadNumbs)) { 
      S=combn(threadNumbs,iBspurs)
#     print(S)
      for (kSp in 1:dim(S)[2]) {
        mat=rbind(mat,mat0) # for each new model/row mat0 is like the full model (i.e. we start with it and delete).
        jj=jj+1
        currSpurs=t(S[,kSp,drop=FALSE])
        nodePositions=NULL
        for (kkSp in 1:length(currSpurs)) nodePositions=c(nodePositions,g$threads[[currSpurs[kkSp]]]$nodes)
        mat[jj,nodePositions]="I"  # knock down threads like nodes
        if((iSite==1)&(max(currSpurs)>1)) {
          hdOuts=currSpurs[currSpurs!=1]-1
          nHdOuts=length(hdOuts)
#          cat("\nIn Spurs: startThread=",startThread," patchLen=",patchLen,"nHdOuts=",nHdOuts,"hdOuts",hdOuts,"\n")
#          cat("\nIn Spurs: nHdOuts=",nHdOuts,"hdOuts",hdOuts,"\n")
          matHS=mat[jj,]
          if (nHdOuts>1) 
            for (iHdOuts in 1:nHdOuts) {
              H=combn(hdOuts,iHdOuts)
              print(H)
              for (kHdOuts in 1:dim(H)[2]) {
                currHdOuts=t(H[,kHdOuts,drop=FALSE])
                mat=rbind(mat,matHS)
                jj=jj+1 
                mat[jj,currHdOuts]="I"
              }
            } else  {# only one hd taken out
            mat=rbind(mat,matHS)
            jj=jj+1 
#                   mat[jj,maxes]="B"
            mat[jj,hdOuts]="I"
          }
        }
      } # kSp loop
    }
#    cat("jj is",jj,"\n")
    mat
  }  # END spurKOThreads
  
  
  mkBigMats<-function(mats,g)
  { # This function forms the product model space of the site spaces (matrices) of the list mats
    mylen<-function(x) dim(x)[1]
#    nSites=length(g$strct$sites)
    lens=sapply(mats,mylen)
# cat("\n Number of rows in mats matrices are:",lens,"\n")
    Bmats=matrix(nrow=prod(lens),ncol=g$nZ)  # this is bigMats
    mkRows<-function(lens,d) {
#   for (i in 1:lens[d]){
      if (d>0) B=Recall(lens,d-1) else {
#     print("in the bottom call")
        return(matrix(nrow=prod(lens),ncol=length(lens)))
      }
      B[,d]=rep(1:lens[d],times=ifelse(d==1,1,prod(lens[1:(d-1)])), each=ifelse(d==length(lens),1,prod(lens[(d+1):length(lens)])))
#     print(B)
      B        # forgot what this was, but may have more to do with Bridges than being Big
    }
#   lens=c(2,3,4)
    rowsI=mkRows(lens,length(lens))
#    print(rowsI)  
# print(Bmats)
    for (i in 1:prod(lens))
      for (iSite in 1:g$nSites){
        Bmats[i,g$nodesWithinSites[[iSite]]]=mats[[iSite]][rowsI[i,iSite],g$nodesWithinSites[[iSite]]]
        Bmats[i,g$hds]=mats[[1]][rowsI[i,1],g$hds]
      }
    n=dim(Bmats)[1]
    m=dim(Bmats)[2]
#    bytes=object.size(Bmats)
# cat("\nSize of ",n,"x",m," big matrix object in bytes is:",bytes, "or ",bytes/(n*m) ,"bytes per element\n")
# print(g$dfThreads)  
#    print(Bmats)
    for (iRow in 1:n)
      for (iOligo in 1:dim(g$dfThreads)[2]){
        if (g$nSites>=2)
          for (iSite in g$nSites:2){
#         print(Bmats[i,threads[[dfThreads[iSite-1,iOligo]]]$nodes[1]])
            if ((Bmats[iRow,g$threads[[g$dfThreads[iSite,iOligo]]]$nodes[1]]!="I")&
                (Bmats[iRow,g$threads[[g$dfThreads[iSite-1,iOligo]]]$nodes[1]]=="I")){
              Bmats[iRow,max(g$threads[[g$dfThreads[iSite-1,iOligo]]]$nodes)]="H"
#           print("in loop")
            } 
          }
      }
#    print(Bmats)
    colnames(Bmats)<-g$KdS
    rownames(Bmats)<-apply(Bmats,1,paste,collapse="")
    sumb<-function(x) sum(x=="H")
    xx=sapply(lapply(apply(Bmats,1,unique),setdiff,c("H","I")),length)
    bb=apply(Bmats,1,sumb)
    xx=xx+bb
    nBmats=Bmats
    nBmats[]=1
    nBmats[Bmats=="I"]=Inf
    Bmats=nBmats
    Bmats=cbind(nParams=xx,as.data.frame(Bmats),p = -1)
    if (pRows) {pBmats=Bmats; 
      rownames(pBmats)=paste(rownames(Bmats),"p",sep=""); 
      pBmats[,"p"]=1; 
      pBmats[,"nParams"]=pBmats[,"nParams"]+1; 
      Bmats=rbind(Bmats,pBmats)
    }
    Bmats=Bmats[order(Bmats$nParams),]
    Bmats=cbind(Bmats,indx=1:dim(Bmats)[1])
    Bmats[2:(g$nZ+1)]=lapply(Bmats[2:(g$nZ+1)],as.numeric)
    Bmats
  }   # END mkBigMats function definition
  
  mkk<-function(Kchunk,kchunk,g,kIC){
    ncols=length(Kchunk)
    meat=kchunk[,2:(g$nZ+1)]
    meat[]=kIC
    names(meat)<-paste("k",g$Z,sep="")
    rnms1=rownames(Kchunk)
    rnms2=rownames(kchunk)
    #    newBmats=cbind(2*Bmats[1,1,drop=F],Bmats[1,2:(g$nZ+1),drop=F],meat[1,,drop=F],Bmats[1,(g$nZ+2):ncols,drop=F])
    newBmats=NULL
#    print(newBmats)
#    break
    brnms=NULL
    for (i in 1:length(rnms1))
      for (j in 1:length(rnms2)) {
#        for (j in ifelse(i==1,2,1):length(rnms)) {
        newBmats=rbind(newBmats,cbind(Kchunk[i,1,drop=F]+kchunk[j,1,drop=F],
                Kchunk[i,2:(g$nZ+1),drop=F],meat[j,,drop=F],Kchunk[i,(g$nZ+2):ncols,drop=F]))
        brnms=c(brnms,paste(rnms1[i],rnms2[j],sep=".") )
      }
#    print(brnms)
    row.names(newBmats)<-brnms
    newBmats=newBmats[order(newBmats$nParams),]
    newBmats
  }
  
  # ***************** FUNCTION DEFINITIONS ABOVE **************************
  
  options(stringsAsFactors = FALSE)
#  
#  print(mkSingleThread(3,g))
#  print(mkSingleThread(2,g))
#  
  mat=matrix(ncol=g$nZ)
  mat[1,1:length(g$hds)]="H"
  for (j in 1:length(g$threads)) mat[1,g$threads[[j]]$nodes]=g$threads[[j]]$let
#  print(mat)
  mat0=mat
  mats=list(NULL)
  for (iSite in 1:g$nSites) { # go through sites (curtain rods are within sites)
    if(g$singleThread) {
      mat=mkSingleThread(length(g$nodesWithinSites[[1]]),g)
    } else { # multiple threads
      mj=mkCurtain(mat0,g,iSite)
      mat=mj$mat
      jj =mj$jj
      mat=spurKOThreads(mat,mat0,g,iSite,jj)
    }
#    print(mat)
    mats[[iSite]]=mat
  } # loop on iSite
#  print("mats is")  
#  print(mats)
  Bmats=mkBigMats(mats,g)
#  break
#  if (!is.null(maxnPs)) Bmats=Bmats[Bmats$nParams<=maxnPs,]
  chunk=Bmats
  Keqs=mkEq(g,chunk)
  keqs=mkEq(g,chunk,activity=TRUE)
  if (!fullGrid){ # if full grid not wanted, remove it
    pNulls=sapply(Keqs,is.null)
# the next two lines remove grids that are really spurs
    chunk=chunk[-which(pNulls),] 
    Keqs=Keqs[-which(pNulls)]
    keqs=keqs[-which(pNulls)]
  }
  if(atLeastOne) { #then also remove models without at least one maximum size oligo
    zor1=lapply(chunk,"%in%",c(0,1)) 
    newW=g$W[-(1:(g$nAtomS-1)),g$hubChar]  # leave one atom row in as spacer for nParams column of chunk 
    colNums=which(newW==max(newW))
    Irows<-apply(as.matrix(as.data.frame(zor1[colNums])),1,sum)>0
    chunk=chunk[Irows,] # grab rows with at least one complex of maximum size 
    Keqs=Keqs[Irows]  # shortening in these next two lines may not be critical
    keqs=keqs[Irows]
  }
#  print(chunk)
#  print(Bmats)
  if (g$activity) chunk=mkk(chunk,Bmats,g,kIC) 
#  print(chunk)
  if (!is.null(maxTotalPs)) chunk=chunk[chunk$nParams<=maxTotalPs,]
  chunk$indx=1:(dim(chunk)[1])
# new block to handle user defined ICs on K params (for readability only bases are entered in chunk)
  KCols=chunk[,2:(g$nZ+1)]
  KCols[KCols==1]=IC
  chunk[,2:(g$nZ+1)]=KCols
  if (!g$activity) keqs=NULL 
  list(chunk=chunk,Keqs=Keqs,keqs=keqs)
}

