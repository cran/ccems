`simulateData` <-
    function(g,init=FALSE,predict=NULL,typeYP=NULL) {
# create sdata (smooth/fine for plotting)if fine=1 and edata for SSE evaluation otherwise
# print(g$parmsTCC)
  doPredict=FALSE
  if (!is.null(predict)) doPredict=TRUE
  if (!doPredict) { 
    posReactants=g$posReactantsD
    dfr=g$d 
    typeY=g$typeYD
  } else {
    nms=names(predict)
    posReactants=rep(0,g$nAtomS)
#    print(posReactants)
#    print(g$free)
#    print(ifelse(g$free,"F","T"))
#    print( paste(g$atomS[2],ifelse(g$free,"F","T"),sep=""))
#   print(nms)
    posReactants[1]=which(nms == paste(g$hubChar,"T",sep=""))
    for (i in 2:g$nAtomS) posReactants[i]=which(nms == paste(g$atomS[i],ifelse(g$free,"F","T"),sep=""))
#    for (i in 1:g$nAtomS) posReactants[i]=which(nms == paste(g$atomS[i],"T",sep=""))
    dfr=predict
    typeY=typeYP
  }  
  options("warn"= -1)
  times <- c(0,4 * 10^c(-1,0,2,7))
#  SSE=0
  
  bound=row.names(g$params)[g$params[,"constr"]!="none"]  # bound = follower, free = leader
  free=g$params[bound,"constr"] # these are the free params to which the bound ones are bound
  g$params[bound,"final"]=g$params[free,"final"]
  nms=names(g$Kparams)
  if (g$Kj) g$parmsTCC[paste("Kj",g$Z,sep="_")]= g$params[nms,"final"]  else 
    g$parmsTCC[paste("Kj",g$Z,sep="_")]= g$Kd2Kj(g$params[nms,"final"])
  g$parmsTCC["p"]=p=g$params["p","final"]
  m1=g$params["m1","final"]
  if (typeY=="k") {# we have additional output linkage g parameters k
    g$kis=g$params[g$kS,"final"]
#    g$kis=g$params[(length(g$Kparams)+3):dim(g$params)[1],"final"]
#    names(g$kis)<-paste("k",g$Z,sep="")  # from RP version
    names(g$kis)<-g$kS
  }  
  
  boxcox<-function(yin,lam=0.5) {
    gm=prod(yin)^(1/length(yin))
    if (lam==0) yout=log(yin)/gm  else
      yout=(yin^lam-1)/(lam*gm^(lam-1))
    yout
  }
  
  logit<-function(x) log(x/(1-x))  # p= e^(AX)/(1+e^(AX)) => p/(1-p) = e^AX
  
  
#    attach(g)
# pull items from attachment to save typing, save them to real object g
#    bound=row.names(params)[params[,"constr"]!="none"]  # bound = follower, free = leader
#    free=params[bound,"constr"] # these are the free params to which the bound ones are bound
#    g$params[bound,"final"]=params[free,"final"]
#    nms=names(Kparams)
#    if (Kj) g$parmsTCC[paste("Kj",Z,sep="_")]= params[nms,"final"]  else 
#        g$parmsTCC[paste("Kj",Z,sep="_")]= Kd2Kj(g$params[nms,"final"])
#    g$parmsTCC["p"]=p=params["p","final"]
#    if (typeY=="v") {# we have additional output linkage g parameters k
#        g$kis=params[(length(Kparams)+2):dim(params)[1],"final"]
#        names(g$kis)<-paste("k",g$Z,sep="")  # from RP version
#    }  
  # print(kis)
#print(g$parmsTCC)
  icNames=names(dfr)[posReactants]; # g is hard coded in reactant/atom order given in atomS
  ndf=dim(dfr)[1] 
  SS=data.frame(matrix(rep(0,ndf*g$nSpecieS,),nrow=ndf));                   
  names(SS)<-g$specieS
  chk=data.frame(matrix(rep(0,ndf*2*length(g$atomS)),nrow=ndf));  
  names(chk)<-c(g$atomS,paste(g$atomS,"Q",sep=""))
  hubQ=paste(g$hubChar,"Q",sep="")
  if (typeY=="P") {
    uniW=unique(g$W[,g$hubChar])
    uniW=uniW[uniW>1]  # e.g. 2,4,6
    EoutNames=paste("EP",uniW,sep="")
    outNames=paste("P",uniW,sep="")
  }
#  p=g$parmsTCC["p"]
# print(SS)
# print(Kparams)
  tights=FALSE
# if (any(g$parmsTCC[(nAtomS+2):(2*nAtomS+1)]==0)) tights=TRUE
  if (any(g$Kparams==0)) tights=TRUE
# cat("\ntights is ",tights,"\n\n")
  if ((tights)&(g$TCC)) {   ###  this is written strictly for the Rt system since I'm not sure that I want to pursue this concept in general
    posZero=which(g$Kparams==0)
    sinkZ=g$Z[posZero]
#  cat("\n posZero is ",posZero," Z is ",sinkZ,"\n\n")
    tightAtomS=unique(g$reactantS[[sinkZ]])
    tightAtomST=paste(unique(g$reactantS[[sinkZ]]),"T",sep="")
    for (ii in 1:ndf)
    { 
      tots=as.numeric(dfr[ii,tightAtomST])
      tots[1]=p*tots[1]
      if (g$W[sinkZ,tightAtomS[1]]*tots[1]>=g$W[sinkZ,tightAtomS[2]]*tots[2]) {
        R=tots[1]-(g$W[sinkZ,tightAtomS[1]]/g$W[sinkZ,tightAtomS[2]])*tots[2]
        t=0} else   {
        t=tots[2]-(g$W[sinkZ,tightAtomS[2]]/g$W[sinkZ,tightAtomS[1]])*tots[1]
        R=0} 
      cmplx=min(tots/g$W[sinkZ,])
      SS[ii,tightAtomS] = c(R,t)          
      SS[ii,sinkZ] = cmplx          
      if (typeY=="m") dfr[ii,"EY"]=m1*((as.matrix(SS[ii,])%*%as.matrix(g$W[,1]^2)) + 
              (1-p)*dfr[ii,posReactants[1]])/dfr[ii,posReactants[1]] 
#         print(as.matrix(SS[ii,]))
#         print(t(as.matrix(kis))) 
#      if (typeY=="v") dfr[ii,"EY"]=as.matrix(SS[ii,(g$nAtomS+1):g$nSpecieS])%*%as.matrix(g$kis*g$W[g$Z,"S"]) 
      #  if ((typeY=="v")&(!TCC)) # ligands >> hub  => must have step function response, e.g. as M-M g slope is k/K and thus infinity  
      #   cat("\ntots are ",tots[1]," and ",tots[2],"ii is ",ii,"\n\n")
    }
# *****************************************
  } else {  # use ODEs instead of logic
# *****************************************
    g$parmsTCC[1:g$nAtomS]=0 # these are overwritten where they are nonzero in the current figure ####### WHERE is this used? 
#        dyn.load(strn<-paste(wDir,"/gs/",id,.Platform$dynlib.ext,sep="")) 
    if (g$TCC&!g$free) dyn.load(strn<-paste("models/",g$id,.Platform$dynlib.ext,sep="")) 
    g$initialStateTCC[]=rep(0,ifelse(g$free,1,g$nAtomS))
    for (ii in 1:ndf) {
      g$parmsTCC[paste(g$hubChar,"T",sep="")]=as.numeric(dfr[ii,posReactants[1]])
      g$parmsTCC[paste(g$atomS[-1],ifelse(g$free,"F","T"),sep="")]=as.numeric(dfr[ii,posReactants[-1]])
      if (g$TCC) {
#            cat("IC is",g$initialStateTCC,"  and parmsTCC is ",g$parmsTCC,"\n")
#        if (!g$free) ic=g$initialStateTCC else ic=g$initialStateTCC[1]
#        e=try(out1 <- lsoda(ic,times,"myderivs", g$parmsTCC, rtol=g$rtol,atol=g$atol, dllname=g$id)) 
      if (!g$free){ # move this down three lines to test ODE single polynom
        e=try(out1 <- lsoda(g$initialStateTCC,times,"myderivs", g$parmsTCC, rtol=g$rtol,atol=g$atol, dllname=g$id)) 
        out=data.frame(out1); if (g$free) names(out)<-c("time",g$atomS[1]) else names(out)<-c("time",g$atomS)
        nout=dim(out)[1]
          SS[ii,] = g$fback(out[nout,2:(g$nAtomS+1)],g$parmsTCC)       
      }    else {
#        freeRa=c(out[nout,2],g$parmsTCC[2:g$nAtomS]) # and uncomment this to try ODE single polynom
        freeRa=c(g$fpoly(g$parmsTCC),g$parmsTCC[2:g$nAtomS])   # and comment this instead
#        print(freeRa)
        SS[ii,] = g$fback(freeRa,g$parmsTCC)     
      }
        chk[ii,g$atomS]=dfr[ii,posReactants]
        chk[ii,(length(g$atomS)+1):(2*length(g$atomS))]=as.matrix(SS[ii,])%*%as.matrix(g$W[,g$atomS])
        chk[ii,hubQ]=chk[ii,hubQ]/p  # correct to expectation
        if (typeY=="m") dfr[ii,"EY"]=m1*((as.matrix(SS[ii,])%*%as.matrix(g$W[,g$hubChar]^2)) + 
                (1-p)*dfr[ii,posReactants[1]])/dfr[ii,posReactants[1]] 
        if (typeY=="k") dfr[ii,"EY"]=(1/(dfr[ii,posReactants[1]]*max(g$W[,"S"])))*
              as.matrix(SS[ii,(g$nAtomS+1):g$nSpecieS])%*%as.matrix(g$kis*g$W[g$Z,"S",drop=FALSE])
        if (typeY=="P") {
          dfr[ii,"EP1"]=(sum(SS[ii,g$W[,g$hubChar]==1])+(1-p)*dfr[ii,posReactants[1]])/dfr[ii,posReactants[1]] 
#          print(EoutNames)
          for (kk in 1:length(EoutNames)) 
            dfr[ii,EoutNames[kk]]=uniW[kk]*sum(SS[ii,g$W[,g$hubChar]==uniW[kk]])/dfr[ii,posReactants[1]] 
#          print(dfr)
        }  
      }
      if ((typeY=="k")&(!g$TCC)) {
#              print(g$kis)
#              print(g$parmsTCC)
        dfr[ii,"EY"]=g$frp(g$parmsTCC,g$kis)
#              print(dfr)
      }
    }  # end loop through rows of dataframe
    if (g$TCC&!g$free) {dyn.unload(strn)}
    if (g$TCC) {g$echk=chk}
    } # end block of ODE solutions (i.e. infinitely tight handled by K=.0001)
  if (g$TCC) g$eSS=SS
  if (!doPredict) {
    g$d=dfr
#  if (typeY=="k") res=(dfr[,g$posY]-dfr[,"EY"])/sqrt(dfr[,"EY"])
   if (typeY!="P") {
     res=dfr[,g$posY]-dfr[,"EY"]
     g$res=res  # these are always straight
     if( g$transform=="none") tres=res
     if( g$transform=="sqrt")  tres=sqrt(dfr[,g$posY])-sqrt(dfr[,"EY"])
     if( g$transform=="log")  tres=log(dfr[,g$posY])-log(dfr[,"EY"])
     if( g$transform=="boxCox")  tres=boxcox(dfr[,g$posY],lam=g$lam)-boxcox(dfr[,"EY"],lam=g$lam)
     if( g$transform=="relResid") tres=(dfr[,g$posY]-dfr[,"EY"])/dfr[,g$posY]
     if (!is.null(g$d$weights)) tres=res*g$d$weights # leave this option in in case means given with different sample numbers
     g$tres=tres  
     SSE=t(tres)%*%tres  # sum of squared transformed residuals
   } else
   { SSE=0
     for (kk in 1:length(outNames)) {
#       tres=logit(dfr[,outNames[kk]])-logit(dfr[,EoutNames[kk]])
       tres=dfr[,outNames[kk]]-dfr[,EoutNames[kk]]
       SSE=SSE+t(tres)%*%tres
     }
   }
#    SSE=SSE+t(res)%*%res
    N=g$nData=ndf
    if (typeY=="P") N=N*length(uniW)
    P=sum(g$params[,"opt"]) + 1 # include the variance
#    phi=SSE/(N-P-1)  # excluse the variance in this calculation
#    mTwoLLM= N*log(2*pi) +             N                 + N*log(SSE/N)
#    mTwoLLk= N*log(2*pi) + sum(log(phi*dfr[,"EY"]))    + SSE/phi
#    cat("phi=",phi,"; -2LLM=",mTwoLLM,"; -2LLk=",mTwoLLk,"\n")
    #    P=sum(g$params[,"opt"]) + ifelse(is.null(g$d$weights),1,length(unique(g$d$weights))) # one param for each weight
    aic=N*log(SSE/N)+2*P + 2*P*(P+1)/(N-P-1) + N*log(2*pi) + N
    if (P>=N-1) aic=Inf
    if (init) {
      g$SSE$initial=SSE; 
      g$AIC$initial = aic
    } # last two terms just to keep it all consistent
    g$SSE$final=SSE; 
    g$AIC$final =  aic
  }    else  {
    if (g$TCC&!g$free){
      keepers=(abs(chk[,g$hub]-chk[,paste(g$hub,"Q",sep="")])< 0.01*chk[,g$hub])&!is.nan(chk[,paste(g$hub,"Q",sep="")]) # it should be rare to not be a keeper
      g$predict=dfr[keepers,]
    } else g$predict=dfr
  } 
  options("warn"=1)
#    print(g)
  return(g)
}

