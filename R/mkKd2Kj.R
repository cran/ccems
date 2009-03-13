`mkKd2Kj` <-function(g) {
  # the number of these mappings is 2^(number of threads in the site)
  bin2dec <- function(x)  sum(x * 2^(rev(seq(along=x)) - 1))  # by burt gunter, better since uses e.g. c(1,0,1,1)=11
#    bin2dec1 <- function(x) {  by marc schwartz, for strings or one number, e.g. 110010
#       x <- as.character(as.numeric(x))
#       b <- as.numeric(unlist(strsplit(x, "")))
#       pow <- 2 ^ ((length(b) - 1):0)
#       sum(pow[b == 1])
#     }
# at most 4 sites are allowed and thus codes are form 0 to 7 for the Inf status of the first three sites filled
  hds=g$hds
  nHds=length(hds)
  nodeStrn=NULL
  for (i in 1:nHds) nodeStrn[i]=paste("x[",i,"]", sep="")
  jj=nHds+1
  for (iSite in 1:dim(g$dfThreads)[1]) {
    for (iOligo in 1:dim(g$dfThreads)[2]) {
      ithread=g$dfThreads[iSite,iOligo]
#			print(ithread)
      currNodes=g$threads[[ithread]]$nodes
      hdPowers=g$W[names(currNodes)[1],]-c(0,1) # find seed head powers for this strip
      hd=paste(g$atomS,hdPowers,sep="",collapse="") # make things like "R4X4"
#			print(hd)
      # next is true if  head is root R
      if ((hdPowers[1]==1)&(hdPowers[2]==0)) thread="1" else  thread=paste("x[",which(g$Z==hd),"]", sep="")
#			cat("\nInitial thread is ",thread,"\n")
      threadSize=length(currNodes)
      ifactors=threadSize:1
      ifactors=ifactors/rev(ifactors); 
#			print(ifactors)
      for (kk in 1:threadSize) {
        thread=paste(thread,"*x[",jj,"]*",ifactors[kk],sep="")
        nodeStrn[jj]=thread
        jj=jj+1
      } # kk loop within blocks
    } # iOligo loop through oligos, within sites, i.e. things that can be equal
  } # iSite loop 
  strn0=nodeStrn
  nshape<-eval(parse(text=paste("function(x) c(\n",paste(nodeStrn,collapse=",\n"),"\n)",sep="")))
  shpLst=list(nshape)
#	print(shpLst)  # this is the complete function, now set parts to Inf and insert bridges
  nSites=length(g$threadsWithinSites) # rows in g$dfThreads
  nOligos=length(g$threadsWithinSites[[1]])  # cols in g$dfThreads
  pThreads=unlist(g$threadsWithinSites[-nSites]) # last filled site never inserted
  #  print(pThreads)  # e.g. 1:8 for (s, a, h)
  if (length(pThreads)>0) {
    codes0=rep(0,length(pThreads))
    codes=codes0
    codeS=paste(rep(0,nOligos),collapse="")  # all zeros already in as seed
    jj=1
    for (i in 1:length(pThreads)){
      M=combn(pThreads,i)
      for (k in 1:dim(M)[2]) {
        currThreads=t(M[,k,drop=FALSE])
        codes=codes0
        codes[currThreads]=1
#        cat("\n\njj =",jj,"\n"); jj=jj+1
#                print(codes)
        currInfs=matrix(codes,byrow=TRUE,nrow=nSites-1)
#                print(currInfs)
        codes=apply(currInfs,2,bin2dec)
#        print(codes)
        codeS=c(codeS,paste(codes,collapse=""))
#                print(codeS)
        strn=strn0 # start with template base of function
        for (jcol in 1:dim(currInfs)[2]) {
          nodePositions=NULL
          for (irow in 1:dim(currInfs)[1]) 
            if(currInfs[irow,jcol]==1) 
              nodePositions = c(nodePositions, g$threads[[g$dfThreads[irow,jcol]]]$nodes)
          strn[nodePositions] = Inf
          if (!is.null(nodePositions)) 
            strn[max(nodePositions)] = paste("x[",max(nodePositions),"]", sep="")
          if (bin2dec(currInfs[,jcol]) == 5)  { # special case of two I's both needing insertions
            nodePositions=g$threads[[g$dfThreads[1,jcol]]]$nodes
            cat("\nIn Kd2Kj with a 101 pattern with first site nodes:",nodePositions,"\n")
            strn[max(nodePositions)] = paste("x[",max(nodePositions),"]", sep="")
          }
        }
        nshape<-eval(parse(text=paste("function(x) c(\n",paste(strn,collapse=",\n"),"\n)",sep="")))
        shpLst=c(shpLst,nshape)
      }
    }
    names(shpLst)<-codeS
    trimMonomer=FALSE
    if(trimMonomer) {
      shpLst=shpLst[sapply(strsplit(codeS,split=NULL),"[",1)=="1"]
      (nms<-names(shpLst))
      (nmsL<-strsplit(nms,split=NULL))
      (nmsL2<-lapply(nmsL,"[",-1))
      (nms2<-sapply(nmsL2,paste,collapse=""))
      names(shpLst)<-nms2
    }
  } else
    names(shpLst)<-"default"
  shpLst
}
