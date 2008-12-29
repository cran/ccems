`mkGrids` <-
		function(g,maxnPs=NULL,pRows=FALSE,contig=TRUE,atLeastOne=TRUE,IC=1) {
	# This function makes a chunk of grid graph models. It is less advanced than its counterpart mkSpurs because I still need a way
	# of going through this model space systematically with increasing numbers of parameters. Thus, it does not have state inputs and outputs. 
	# Note that it is implicit that R is always a thread head too
	
	lstEq <- function(g,grids) {
# This function converts a dataframe of K equality models into a list structure that maps the trackers onto the independent parameter estimate. }
		nms=colnames(grids)
		nms=nms[c(-1,-(length(nms)-1),-length(nms))] # eliminate num params and p param, leave backbones in as place holders
#	print(nms)
		rnms<-rownames(grids)
		KeqsG=rep(list(NULL),length(rnms))
		for (iRow in 1:length(rnms)) {
			rChars=strsplit(rnms[iRow],split="")[[1]]
			constraintsVec=NULL
			for (curChar in g$usedLets) {
				pntrs=which(rChars==curChar)
				cnts=length(pntrs)
				if (cnts >1) {
					tmp=rep(nms[pntrs[1]],cnts-1)
					names(tmp)<-nms[pntrs[2:cnts]]
					constraintsVec=c(constraintsVec,tmp)
				}
			}
			if(!is.null(constraintsVec)) KeqsG[[iRow]]=constraintsVec
		}
		names(KeqsG)<-rnms
		KeqsG
	}
	
	options(stringsAsFactors = FALSE)
	hds=g$hds
	nHds=length(hds)
	jj=nHds # jj points to nodes
	iThread=1 # count up all the blocks and index them in order
	threads=list(NULL)
	iSite=1
	mylets=c("D","E","F","L","M","N","Q","V","Y",paste(2:9),"e","f","m","n","q","v","w","y") # leave I and J for inf and free as before
#	print(mylets)
	threadsWithinSites=list(NULL)
	nodesWithinSites=list(NULL)
	usedLets=character(0)
	nSites=length(g$strct$sites)
	dfThreads=data.frame(NULL)
	for (iSite in 1:nSites) {
		currSite=g$strct$sites[[iSite]]
		threadsWithinSites[[iSite]]=as.numeric(NULL)
		nodesWithinSites[[iSite]]=as.numeric(NULL)
		for (iOligo in 1:length(currSite)) {
			currOligo=currSite[[iOligo]]
			threadSize=length(currOligo)
			nodes=(jj+1):(jj+threadSize)
			names(nodes)<-currOligo
			threads[[iThread]]=list(site=iSite,let=mylets[iThread],nodes=nodes)
			usedLets=c(usedLets,mylets[iThread])
			threadsWithinSites[[iSite]]=c(threadsWithinSites[[iSite]],iThread)
			nodesWithinSites[[iSite]]=c(nodesWithinSites[[iSite]],nodes)
			iThread=iThread+1
			jj=jj+threadSize
		} # iOligo loop through oligos, within sites, i.e. things that can be equal
		dfThreads=rbind(dfThreads,threadsWithinSites[[iSite]])
		iSite=iSite+1
	} # currSite loop over things that cannot be equal
#	print("nodesWithinSites") 
#	print(nodesWithinSites)
#	print("threadsWithinSites") 
#	print(threadsWithinSites)
#	print("dfThreads dataframe is") 
	names(dfThreads)<-names(g$strct$sites[[1]])
	rownames(dfThreads)<-names(g$strct$sites)
#	print(dfThreads)
#	break
	nsites=iSite-1
	names(threads)<-paste("t",1:(nThreads<-iThread-1),sep="")
#	print("threads are") 
#		print(threads)
	mat=matrix(ncol=g$nZ)
#	mat[1,]="J"  # holder: these should all be replaced
	mat[1,1:nHds]="H"
	for (j in 1:length(threads)) mat[1,threads[[j]]$nodes]=threads[[j]]$let
#	print(mat)
	mat0=mat
	mats=list(NULL)
	for (iSite in 1:nSites) { # go through sites
		threadNumbs=threadsWithinSites[[iSite]]
#		flag=1
		mat=mat0
		jj=1
#		while ((length(threadNumbs)>2)|((length(threadNumbs)==2)&(flag==1))) {
#			flag=0
		for (patchLen in 2:length(threadNumbs)) {
			if (!contig) 
			{M=combn(threadNumbs,patchLen)
				for (k in 1:dim(M)[2]) {
					nodePositions=NULL
					currThreads=t(M[,k,drop=FALSE])
					for (kk in 1:length(currThreads)) nodePositions=c(nodePositions,threads[[currThreads[kk]]]$nodes)
					mat[jj,nodePositions]=threads[[min(currThreads)]]$let
					mat=rbind(mat,mat0)
					jj=jj+1
				} # k loop
			}  else # contiguous patch
				for (startThread in min(threadNumbs):(max(threadNumbs)-patchLen+1)) {
#					cat("\nstartThread=",startThread," patchLen=",patchLen,"\n")
					nodePositions=NULL
					for (kk in startThread:(startThread+patchLen-1)) nodePositions=c(nodePositions,threads[[kk]]$nodes)
					mat[jj,]=mat0 # reset in case not properly reset
					mat[jj,nodePositions]=threads[[startThread]]$let
					mat=rbind(mat,mat1<-mat[jj,])
					jj=jj+1
					#					print(mat)
					if (1) { #all combs of placing the I's
						outside=setdiff(threadNumbs,startThread:(startThread+patchLen-1))
						nOuts=length(outside)
#						print(outside)
#						cat("\nstartThread=",startThread," patchLen=",patchLen,"nOuts=",nOuts,"outside",outside,"\n")
						if (nOuts>1) 
							for (iOuts in 1:nOuts) { 
								X=combn(outside,iOuts)
#								print(X)
								for (kOuts in 1:dim(X)[2]) {
									currOuts=t(X[,kOuts,drop=FALSE])
									nodePositions=NULL
#									maxes=NULL
									for (kkOuts in 1:length(currOuts)) {
										nodePositions=c(nodePositions,threads[[currOuts[kkOuts]]]$nodes)
#										maxes=c(maxes,max(threads[[currOuts[kkOuts]]]$nodes))
									}
									mat[jj,nodePositions]="I"
									if((iSite==1)&(max(currOuts)>1)) {
										hdOuts=currOuts[currOuts!=1]-1
										nHdOuts=length(hdOuts)
#										cat("\nstartThread=",startThread," patchLen=",patchLen,"nHdOuts=",nHdOuts,"hdOuts",hdOuts,"\n")
										matH=mat[jj,]
										if (nHdOuts>1) 
											for (iHdOuts in 1:nHdOuts) {
												H=combn(hdOuts,iHdOuts)
#												print(H)
												for (kHdOuts in 1:dim(H)[2]) {
													currHdOuts=t(H[,kHdOuts,drop=FALSE])
													mat=rbind(mat,matH)
													jj=jj+1 
													mat[jj,currHdOuts]="I"
												}
											} else  {# only one hd taken out
											mat=rbind(mat,matH)
											jj=jj+1 
#										mat[jj,maxes]="B"
											mat[jj,hdOuts]="I"
										}
									}
									if ((iOuts==nOuts)&(kOuts == dim(X)[2])) mat=rbind(mat,mat0) else mat=rbind(mat,mat1)
									jj=jj+1
								} # kOuts loop
							}	else if (nOuts==1) {
							mat[jj,threads[[outside]]$nodes]="I"
							if((iSite==1)&(outside>1)) {
#								mat[jj,max(threads[[outside]]$nodes)]="B"
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
#						if (leftGap==0) mat=rbind(mat,mat0) else mat=rbind(mat,mat[jj,])
						if (leftGap>0) 
							for (iLeft in 1:leftGap) {
#								cat("\nstartThread=",startThread," patchLen=",patchLen,"leftGap=",leftGap,"iLeft=",iLeft,"\n")
								mat[jj,threads[[iLeft-1+min(threadNumbs)]]$nodes]="I"
#							matL=mat[jj,]
								mat=rbind(mat,matL<-mat[jj,])
								jj=jj+1
								if (rightGap>0) 
									for (iRight in 1:rightGap) {
#									cat("\nstartThread=",startThread," patchLen=",patchLen,"leftGap=",leftGap,"iLeft=",iLeft,"rightGap=",rightGap,"iRight=",iRight,"\n")
										mat[jj,threads[[max(threadNumbs)-iRight+1]]$nodes]="I"
										if (iRight==rightGap) mat=rbind(mat,matL) else mat=rbind(mat,mat[jj,])
										jj=jj+1
									}
							}
						if ((leftGap==0)&(rightGap>0)) {
							for (iRight in 1:rightGap) {
#							cat("\nstartThread=",startThread," patchLen=",patchLen,"leftGap is 0 rightGap=",rightGap,"iRight=",iRight,"\n")
								mat[jj,threads[[max(threadNumbs)-iRight+1]]$nodes]="I"
								if (iRight==rightGap) mat=rbind(mat,mat0) else mat=rbind(mat,mat[jj,])
								jj=jj+1
							}
						}
						
					}
				} # loop on startThread
		} # loop on patchLen
		mat[jj,]=mat0  # last row is always mat0, even if all block sizes are 1 => redundant with spur graph
		# make block spurs next
		for (iBspurs in 1:length(threadNumbs)) { 
			S=combn(threadNumbs,iBspurs)
#			print(S)
			for (kSp in 1:dim(S)[2]) {
				mat=rbind(mat,mat0)
				jj=jj+1
				currSpurs=t(S[,kSp,drop=FALSE])
				nodePositions=NULL
				for (kkSp in 1:length(currSpurs)) nodePositions=c(nodePositions,threads[[currSpurs[kkSp]]]$nodes)
				mat[jj,nodePositions]="I"
				if((iSite==1)&(max(currSpurs)>1)) {
					hdOuts=currSpurs[currSpurs!=1]-1
					nHdOuts=length(hdOuts)
#					cat("\nIn Spurs: startThread=",startThread," patchLen=",patchLen,"nHdOuts=",nHdOuts,"hdOuts",hdOuts,"\n")
					matHS=mat[jj,]
					if (nHdOuts>1) 
						for (iHdOuts in 1:nHdOuts) {
							H=combn(hdOuts,iHdOuts)
#							print(H)
							for (kHdOuts in 1:dim(H)[2]) {
								currHdOuts=t(H[,kHdOuts,drop=FALSE])
								mat=rbind(mat,matHS)
								jj=jj+1 
								mat[jj,currHdOuts]="I"
							}
						} else  {# only one hd taken out
						mat=rbind(mat,matHS)
						jj=jj+1 
#										mat[jj,maxes]="B"
						mat[jj,hdOuts]="I"
					}
				}
			} # kSp loop
		}
		mats[[iSite]]=mat
	} # loop on i/sites
#				if(iLeft==leftGap) mat=rbind(mat,mat0) else mat=rbind(mat,mat[jj,])
	# below is a temporary fix: the J are either free or inf, thus only a doubling of models
#	print(jj)
#	print(mats)
#	print(dim(mats[[1]]))
#	print(dim(mats[[2]]))
#	break
	
	mylen<-function(x) dim(x)[1]
	lens=sapply(mats,mylen)
#	cat("\n Number of rows in mats matrices are:",lens,"\n")
	Bmats=matrix(nrow=prod(lens),ncol=g$nZ)
	
	
	mkRows<-function(lens,d) {
#		for (i in 1:lens[d]){
		if (d>0) B=Recall(lens,d-1) else {
#			print("in the bottom call")
			return(matrix(nrow=prod(lens),ncol=length(lens)))
		}
		B[,d]=rep(1:lens[d],times=ifelse(d==1,1,prod(lens[1:(d-1)])), each=ifelse(d==length(lens),1,prod(lens[(d+1):length(lens)])))
#			print(B)
		B
	}
#		lens=c(2,3,4)
	rowsI=mkRows(lens,length(lens))
#	print(rowsI)	
#	print(Bmats)
	
	for (i in 1:prod(lens))
		for (iSite in 1:nSites){
			Bmats[i,nodesWithinSites[[iSite]]]=mats[[iSite]][rowsI[i,iSite],nodesWithinSites[[iSite]]]
			Bmats[i,hds]=mats[[1]][rowsI[i,1],hds]
		}
#	print(Bmats)
	n=dim(Bmats)[1]
	m=dim(Bmats)[2]
	bytes=object.size(Bmats)
#	cat("\nSize of ",n,"x",m," big matrix object in bytes is:",bytes, "or ",bytes/(n*m) ,"bytes per element\n")
#	print(dfThreads)	
	for (iRow in 1:n)
		for (iOligo in 1:dim(dfThreads)[2]){
			if (nSites>=2)
				for (iSite in nSites:2){
#					print(Bmats[i,threads[[dfThreads[iSite-1,iOligo]]]$nodes[1]])
					if ((Bmats[iRow,threads[[dfThreads[iSite,iOligo]]]$nodes[1]]!="I")&(Bmats[iRow,threads[[dfThreads[iSite-1,iOligo]]]$nodes[1]]=="I")){
						Bmats[iRow,max(threads[[dfThreads[iSite-1,iOligo]]]$nodes)]="H"
#						print("in loop")
					} 
				}
		}
	
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
	if (!is.null(maxnPs)) Bmats=Bmats[Bmats$nParams<=maxnPs,]
	g$dfThreads=dfThreads
	g$threads=threads
	g$threadsWithinSites=threadsWithinSites
	g$nodesWithinSites=nodesWithinSites
	g$usedLets=usedLets
	
	chunk=Bmats
	Keqs=lstEq(g,chunk)
	pNulls=sapply(Keqs,is.null)
# the next two lines remove grids that are really spurs
	chunk=chunk[-which(pNulls),] 
	Keqs=Keqs[-which(pNulls)]
  if(atLeastOne) { #then also remove models without at least one maximum size oligo
    zor1=lapply(chunk,"%in%",c(0,1)) 
    newW=g$W[-(1:(g$nAtomS-1)),"R"]  # leave one atom row in as spacer for nParams column of chunk 
    colNums=which(newW==max(newW))
    Irows<-apply(as.matrix(as.data.frame(zor1[colNums])),1,sum)>0
    chunk=chunk[Irows,] # grab rows with at least one complex of maximum size 
    Keqs=Keqs[Irows]
  }
	chunk$indx=1:(dim(chunk)[1])
  # new block to handle user defined ICs on K params (for readability only bases are entered in chunk)
  ncols=dim(chunk)[2]
  mainCols=chunk[,2:(ncols-2)]
  mainCols[mainCols==1]=IC
  chunk[,2:(ncols-2)]=mainCols
  
	list(g=g,chunk=chunk,Keqs=Keqs)
}

