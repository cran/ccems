`fitModel` <-
    function(model) {   # to allow use with lapply on a list of models
  
  fitModelSub <- function(model) { # subroutine of fitModel. This function existed well before the parent fitmodel function and
    # there may be a way to integrate them together to have smoother code, i.e. their relationship was 
    # not designed but rather driven by convenience given that fitModelSub already existed. 
#    print(getwd())
    fopt <- function(pars,model) {
      pars=exp(pars)
      spt=which(names(pars)=="p")
      if (length(spt)==1) pars[spt]=pars[spt]/(1+pars[spt])
#      if ("p")
      model$params[names(pars),"final"]=pars  
#			if (model$TCC) model=simulateData(model) else model=simulateDataRP(model) 
      model=simulateData(model) 
#      print(model$SSE$final)
      model$SSE$final   # for a given model, number of params is fixed, so goal is to minimize SSE
    }
    
    model=simulateData(model,init=TRUE)  # computes the initial SSE 
    p0=model$params[model$params[,"opt"],"initial"]
    names(p0)<-row.names(model$params)[model$params[,"opt"]]
#    print(p0)
    p0nms=names(p0)
    if ((model$nOptParams<-length(p0))>0) {                    

      for ( j in 1:model$nOptParams)         
        if ((p0nms[j]!="p") & (p0nms[j]!="m1") & (length(grep("_",p0nms[j]))==0) & (length(grep("k",p0nms[j]))==0)) {
          biRcts=sum(model$W[p0nms[j],])-1 # number of binary reactions 
          p0[j]=p0[j]^biRcts 
        }
      
      np0=sapply(p0,log)
#      print(np0)
      spt=which(names(np0)=="p")
#      print(spt)
      if (length(spt)==1) np0[spt] = log(.9/(1-.9)) # if optimized assume it starts at 0.9 instead of 1 to avoid Inf
      p0=np0
#      print("here0")
#      print(p0)                                                        # maxit default=500 for Nelder-Meade (default)
      if (model$nOptParams>1) opt<-optim(p0,fopt,hessian=TRUE,model=model,control=list(maxit=5000)) else
        opt<-optim(p0,fopt,method="BFGS",hessian=TRUE,model=model);       
#  print("here")
#    attach(model)
#			nOptp=length(tmp<-intersect(c("p"),names(p0)) ) # zero or one
#			model$fitid=paste(paste(model$id,model$nZ,sep=""),model$nOptParams-nOptp,nOptp,sep=".") # try life without this
      sg<-sqrt(opt$value/(model$nData-model$nOptParams))
#      print("sg="); print(sg)
      if (det(opt$hessian)>0) 
      {sig=sg*sqrt(diag(solve(opt$hessian/2)));model$hess=TRUE} else
      {sig=Inf;model$hess=FALSE} 
      upper=signif(opt$par+1.96*sig,3)
      lower=signif(opt$par-1.96*sig,3)
      point=signif(opt$par,3)
      CI=cbind(lower,point,upper)
      CI=exp(CI); opar<-exp(opt$par)
      if (length(spt)==1) {
        CI[spt,]=CI[spt,]/(1+CI[spt,])
        opar[spt]=opar[spt]/(1+opar[spt])
      }
        model$params[names(opar),"final"]=opar
#    detach(model)
      model=simulateData(model)  # computes the final SSE 
      cat("\n")
#      print(CI)
#      print(opar)
      model$CI=CI
#     model=simulateData(model,fine=TRUE) # create smooth curve for publication
    } else   {#print("length p0=0 => nothing to fit"); 
      model$fitS="nothing to fit"; 
    } 
    # model=simulateData(model,fine=FALSE,init=TRUE); }
    # if (is.null(model$AIC$final)) model$AIC$final=100
    # if (is.na(model$AIC$final)) model$AIC$final=100
    model
  }   # end fitModelSub
  # now start actual function fitModels
  # This stuff was initially driven by a need to use lapply on lists of models.  It grew from there
  # and now is primarily concerned with making the report ready for the html output, e.g. time taken to fit. 
  
  strns=c("present","absent")[(model$params[!model$params$opt,"final"]==Inf)+1]
  strns[strns=="present"]="fixed"
  strns[strns!="fixed"]="absent"
  strns[model$params[!model$params$opt,"constr"]!="none"]="constrained"
  t0 = Sys.time();      
  e=try(model<-fitModelSub(model))
  model$cpu<-difftime(Sys.time(),t0,units="mins")[[1]]
  if (is.null(model$SSE)) model$SSE=list(initial=1e6,final=1e6)
  if (is.null(model$AIC)) {model$AIC=list(initial=1e6,final=1e6); model$nOptParams=0}
#      if (is.null(model$nOptParams)) model$nOptParams=0
  model$report=rbind(model$params[model$params$opt,c("initial","final")],
      model$params[!model$params$opt,c("initial","final")],
#	model$report=rbind(subset(model$params,select=c("initial","final"),subset=model$params$opt),
#			subset(model$params,select=c("initial","final"),subset=(!model$params$opt)),
      SSE=model$SSE,
      AIC=model$AIC,
      cpu=c(0,model$cpu))
  if (class(e)=="try-error") 
  {try(detach(model)); 
    model$fitS="fit failed"
    model$report=cbind(model$report,confidenceInterval=c(rep("fit failed",dim(model$report)[1]-length(strns)-3),strns,c(rep("",2),"fit failed")))
    model$report["cpu","final"]=model$cpu       } else 
  if (model$nOptParams>0)    {
    cistrn=character(0) # CI string
    pestrn=character(0) # point estimate string
    for ( j in 1:model$nOptParams) {        
      biRcts=sum(model$W[row.names(model$report)[j],])-1 # number of binary reactions 
#						 biRcts=nchar(row.names(model$report)[j])-1
      if ((row.names(model$report)[j]=="p") | (row.names(model$report)[j]=="m1") |(length(grep("_",row.names(model$report)[j]))>0) 
            | (length(grep("k",row.names(model$report)[j]))>0)) {
        cistrn=c(cistrn,sprintf("(%4.3f, %4.3f)",model$CI[j,"lower"],model$CI[j,"upper"])) 
        pestrn=c(pestrn,sprintf("%4.3f",model$CI[j,"point"])) 
      }
      else {
        cistrn=c(cistrn,sprintf("(%4.3f^%d, %4.3f^%d)",model$CI[j,"lower"]^(1/biRcts),biRcts,model$CI[j,"upper"]^(1/biRcts),biRcts))
        pestrn=c(pestrn,sprintf("%4.3f^%d",model$CI[j,"point"]^(1/biRcts),biRcts))
      }
    } # for loop on j to fill cistrn
    model$report=cbind(model$report, 
        confidenceInterval=c(cistrn,strns,rep("",2),ifelse(model$hess,"fit succeeded","hessian singular")), 
        pointEstimate=c(pestrn,rep("",length(strns)+3)) ) 
#		model$report=cbind(model$report,confidenceInterval=c(cistrn,strns,rep("",2),ifelse(model$hess,"fit succeeded","hessian singular"))) 
  } else 
  {model$report=cbind(model$report,confidenceInterval=c(strns,rep("",2),"nothing to fit"));model$report["cpu","final"]=model$cpu       }
  print(model$mid)
  print(model$report)
  cat('\n\n')
  model
}

