mkEq <- function(g,chunk,activity=FALSE) {
# This function converts a dataframe of K equality models into a list structure that
#  maps the trackers onto the independent parameter estimate. }
  nms=colnames(chunk)
  nms=nms[c(-1,-(length(nms)-1),-length(nms))] # eliminate num params and p param, leave backbones in as place holders
#  nms=nms[c(-1,-(length(nms)-2),-(length(nms)-1),-length(nms))] # eliminate num params and p param, leave backbones in as place holders
# print(nms)
  rnms<-rownames(chunk)
  KeqsG=rep(list(NULL),length(rnms))
  for (iRow in 1:length(rnms)) {
    rChars=strsplit(rnms[iRow],split="")[[1]]
    pntrs=which(rChars==".")
#    print(pntrs)
    if (length(pntrs)>0) rChars=rChars[-(1:pntrs)]
#    print(rChars)
    constraintsVec=NULL
    for (curChar in g$usedLets) {
      pntrs=which(rChars==curChar)
      cnts=length(pntrs)
      if (cnts >1) {
        if(!activity) {
          tmp=rep(nms[pntrs[1]],cnts-1)
          names(tmp)<-nms[pntrs[2:cnts]]
        } else {
          tmp=rep(paste("k",g$Z[pntrs[1]],sep=""),cnts-1)
          names(tmp)<-paste("k",g$Z[pntrs[2:cnts]],sep="")
        }
        constraintsVec=c(constraintsVec,tmp)
      }
    }
    if(!is.null(constraintsVec)) KeqsG[[iRow]]=constraintsVec
  }
  names(KeqsG)<-rnms
  KeqsG
}

