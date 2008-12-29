`ems` <-
    function(d,g,cpusPerHost=c("localhost" = 1),ptype="",
        chunkParams=list(size=500,n=1,maxnPs=2,extend2maxP=TRUE),smart=FALSE,
        pRows=FALSE,doTights=FALSE,doGrids=TRUE,doSpurs=TRUE,topN=5,showConstr=FALSE,atLeastOne=TRUE,IC=1) {
  
  if(chunkParams$maxnPs<2) doGrids=FALSE
  
  mysum<-function(x) sum(x$params$opt)  # sum the numer of parameters fitted/optimized
  
  # get short objects so all fits can be stored and AIC histograms plotted
  getShort<-function(model) list(mid=model$mid,indx=model$indx,ci=model$CI,nParams=model$nParams,Kj=model$Kj,
        sreport=model$report[model$report$confidenceInterval!="absent",])
  # for use with lapply on model lists
  
  myXAIC <- function(x) {#cat("x$AIC are ",x$AIC[[1]],x$AIC[[2]],"\n")
    # if(is.null(x$AIC)|is.nan(x$AIC)) {x=list(); x$AIC$final=100}
    if(!is.numeric(x$AIC$final)) {x=list();
      x$AIC$final=1e6}
    x$AIC$final[1]}
  
  grabTopN <- function(fmodels,topN=10) {# grab the best topN models from the model list
# print(fmodels)
#   aic<-as.numeric(sapply(fmodels,myXAIC))
#   fmodels=fmodels[-which(aic==1e6)]  # do not allow failures into top 10
    aic<-as.numeric(sapply(fmodels,myXAIC))
# print("ingrab10")
#   print(aic) 
    saic=sort(aic)
    iBest=NULL
    for (i in 1:length(saic)) iBest=c(iBest,which(aic==saic[i]))
    smodels=fmodels[unique(iBest)]
    nMS=length(smodels)
    smodels[[1]]$nMS=nMS
    if (nMS>topN) smodels=smodels[1:topN]  # only look at and save the top 10 models
    smodels
  }
  
  
  printAICs <- function(models) { # for monitoring progress 
    for (i in 1:length(models))
#cat("The final AIC of the",i, "th model with index ",models[[i]]$indx," is ",models[[i]]$AIC[[2]],"\n")
      cat(sprintf("%3d Model %3d; nbp=%2d; id=%20s; AIC=%8.4f\n",
              i,models[[i]]$indx,models[[i]]$nOptParams,models[[i]]$mid, models[[i]]$AIC[[2]]))
  }
  
  
  cpus=sum(cpusPerHost)
  hosts=names(cpusPerHost)
  slaves=hosts[-1]
  if (ptype=="") {cpus=1}

#  # Blocked commented out for windows checks
#  if (ptype=="RMPI") {  
#    if (!is.loaded("mpi_initialize")) { 
#      require("Rmpi")  # see rpvm comments below, same applies here
#      mpi.spawn.Rslaves(nslaves=(cpus-1)) 
#    }
#    cat("RMPI univ size is ",mpi.universe.size(),"\n\n")
#    print(mpi.remote.exec(paste("I am",mpi.comm.rank(),"of",mpi.comm.size()))) 
#    mpi.bcast.cmd(library(odesolve)) 
#    mpi.bcast.cmd(library(ccems)) 
#  }
  
  
  if ((ptype=="PVM")|(ptype=="NWS")|(ptype=="SOCK")|(ptype=="MPI")) 
  { require(snow)
    
#  ############## next Block commented to allow check in windows     
#    setDefaultClusterOptions(outfile="/home/radivot/debug.t")
#     if (ptype=="PVM") {  
#            tryCatch(require("rpvm",quietly=TRUE,character.only = TRUE), # from Chambers 2008 page 76
#                    warning = function(w) FALSE,  # tried to get through vista package check without error or warning 
#                    error = function(w) FALSE)  #    but this didn't help. Needs pvm in suggests, not enhances.       
#      require(rpvm)
#      .PVM.start.pvmd()
#      print(slaves)
#      .PVM.addhosts(slaves)
#    }
#    if ((ptype=="PVM")|(ptype=="MPI")) {  
#      cl <- makeCluster(ifelse(ptype=="PVM",cpus,cpus-1), type=ptype, verbose=TRUE)   # different env on the compute nodes
#    }

    
    if ((ptype=="SOCK")|(ptype=="NWS")) {  
      strn=rep(hosts,times=cpusPerHost)
 			cl <- makeCluster(strn, type=ptype, verbose=TRUE) # this one can work either way
    }
    print(clusterCall(cl, function() Sys.info()[c("nodename","machine")]))
    clusterEvalQ(cl, require(odesolve))
    clusterEvalQ(cl, require(ccems))  # functions already all there from global library(ccems) in parent script
#    clusterExport(cl,ls(.GlobalEnv))
  }
  t0 = Sys.time();      
  if (sum(dir()=="results")==0) system("mkdir results")
  globalTopN=NULL
  pCnts=rep(0,6)
  gBestAics=rep(1e6,topN)
  cat("\n nChunks*chunkSize is",chunkParams$n,"*",chunkParams$size," = ",chunkParams$n*chunkParams$size,"\n")
  cat("maxnPs is ",chunkParams$maxnPs,"\n\n")
  i=ifelse(doGrids,0,1) 
  done=FALSE
  allModels=NULL  # shortened objects but all of those fitted
  maxReached=FALSE
  lastCompleted=0
  print(i)
  while (!done) 
  { # loop on chunks
    if (i==1) state=list(globMdlIndex=ifelse(doGrids,dim(chunk)[1],0),globCmbIndex=0,relCmbIndex=0,config=NULL)
    if (i==0) {
      ng=mkGrids(g,chunkParams$maxnPs,pRows=pRows,IC=IC) # ng for n-shaped grid
      g=ng$g
      chunk=ng$chunk
      Keqs=ng$Keqs
      Kmapping=mkKd2Kj(g)
    } else  {
      sp=mkSpurs(g,state,chunkParams$maxnPs,batchSize=chunkParams$size,pRows=pRows,doTights=doTights,IC=IC)
      chunk=sp$chunk
      state=sp$state
#			print(sp)
      maxReached=sp$maxReached
      lastCompleted=sp$lastCompleted
    }
    pCnts=pCnts+hist(chunk$nParams,breaks=seq(-.5,40.5,1))$counts[1:6]
    models=NULL
    mdlNames=rownames(chunk)
#		print(Kmapping)
    for (j in mdlNames){ 
#			print(j);print(i) 
      if (i>0)  models[[j]]=mkModel(g,j,d,Kjparams=chunk[j,g$Z],pparams=chunk[j,"p",drop=FALSE],
                                         indx=chunk[j,"indx"],nParams=chunk[j,"nParams"])  else 
        models[[j]]=mkModel(g,j,d,Kdparams=chunk[j,2:(g$nZ+1),drop=FALSE],
            Keq=Keqs[[j]],Kd2KjLst=Kmapping,pparams=chunk[j,"p",drop=FALSE],indx=chunk[j,"indx"],nParams=chunk[j,"nParams"])
    }
    msid=paste(g$id,chunkParams$n,"p",i,ptype,sep="")
#		print(msid)
#		print(models)
    nMS=length(models)
    nump=sapply(models,mysum);
    curBestAics=rep(1e6,topN)
    # uncomment next line to bring back MPI without snow
##    if (ptype=="RMPI")  fmodels=mpi.applyLB(models,fitModel) 
    if ((ptype=="PVM")|(ptype=="NWS")|(ptype=="SOCK")|(ptype=="MPI")) fmodels=clusterApplyLB(cl,models,fitModel)
    if (ptype=="")     fmodels=lapply(models,fitModel) # single processor
    allModels=c(allModels,lapply(fmodels,getShort))
    aic=sapply(fmodels,myXAIC)   # fmodels = fitted models
    tmp=tapply(aic,nump,min)
    curBestAics[sort(unique(nump))+1]=tmp      
    print(totTime<-difftime(Sys.time(),t0,units="mins")) # print("about to call grab10 on fmodels")
    smodels=grabTopN(fmodels,topN=topN) # cat("\n length of smodels after grabTop10 is ",length(smodels),"\n")
    smodels[[1]]$msid=msid                # smodels = sorted fitted models
    smodels[[1]]$totTime=as.numeric(totTime)
    smodels[[1]]$cpus=cpus
    smodels[[1]]$smart=smart
    smodels[[1]]$hists=pCnts
    smodels[[1]]$bestAics=curBestAics
    smodels[[1]]$ptype=ptype
    printAICs(smodels)
    print(msid)
    mkHTML(smodels,showConstr=showConstr)
    if (.Platform$OS.type!="windows")	system("r2w") # this puts the file out where it can be grabbed on the internet
    gBestAics= apply(cbind(gBestAics,curBestAics),1,min)
    bestAicsCompleted = gBestAics[1:(lastCompleted+1)]      
#		print(gBestAics)
    print(bestAicsCompleted)
    globalTopN=c(globalTopN,smodels)
    globalTopN=grabTopN(globalTopN,topN=topN)
    i=i+1
#		cat("i=",i,"\n")
#		print(sp$maxReached)
    if ((i > chunkParams$n)&(!chunkParams$extend2maxP)) done=TRUE
    if ((!doSpurs)&(i==1)) done=TRUE
    if (maxReached) done=TRUE
    if (lastCompleted>0) if ( smart & (bestAicsCompleted[lastCompleted+1]>bestAicsCompleted[lastCompleted]) )  done = TRUE
  }  # loop through chunks indexed by i
  globalTopN[[1]]$msid=paste(g$id,"glob",ptype,sep="")
  globalTopN[[1]]$nMS=ifelse(!doSpurs,dim(chunk)[1],state$globMdlIndex)  # total over all chunks
  globalTopN[[1]]$hists=pCnts
  globalTopN[[1]]$bestAics=gBestAics
  globalTopN[[1]]$ptype=ptype
  print(totTime<-difftime(Sys.time(),t0,units="mins")) 
  globalTopN[[1]]$totTime=as.numeric(totTime)
  save(globalTopN,file=paste("results/",globalTopN[[1]]$msid,"topN.RData",sep=""))
  getAIC <- function(x) { x$sreport["AIC","final"]}
  aic=sapply(allModels,getAIC)
  indx=sort.list(aic)
  sAllModels=allModels[indx]
  save(sAllModels,file=paste("results/",globalTopN[[1]]$msid,"ic",IC,".RData",sep=""))
  aicFt=aic[aic!=1e6]  # next two lines are FYI
  cat("Fitted =", length(aicFt),", out of a total of ",length(aic),"\n")
  mkHTML(globalTopN,showConstr=showConstr)
# uncomment next two to bring these approaches back in to play
##  if (ptype=="RMPI") mpi.close.Rslaves()
##  if (ptype=="PVM") {.PVM.delhosts(slaves); .PVM.halt()}
  if ((ptype=="PVM")|(ptype=="NWS")|(ptype=="SOCK")|(ptype=="MPI")) stopCluster(cl)
  printAICs(globalTopN)
  globalTopN
} 

