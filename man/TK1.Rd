\name{TK1}
\alias{TK1}
\docType{data}
\title{ Thymidine Kinase 1 Data}
\description{
  Human thymidine kinase 1 (i.e. cytosolic) data.   
}
\usage{data(TK1)}
\format{
  A data frame with the following columns.
  \describe{
    \item{\code{E}}{The total concentration of TK1 enzyme. }
    \item{\code{S}}{The total concentration of dT (the nucleoside substrate).}
    \item{\code{X}}{The total concentration of ATP.}
    \item{\code{v}}{The velocity of the kinase reaction.}
    \item{\code{fg}}{The figure number in the original paper.}
    \item{\code{year}}{The year of the original paper.}
    \item{\code{jrnl}}{The journal of the data source: 
    \code{PEP} is Protein Expression and Purification, 
    \code{EJB} is  European Journal of Biochemistry, 
    \code{BBRC} is Biochem Biophys Res Commun, 
    \code{JBC} is Journal of Biological Chemistry.
    }
    \item{\code{vol}}{The volume number of the article.}
    \item{\code{page}}{The page number of the article.}
    \item{\code{frstAut}}{The first author of the article. This is a factor.  }
    \item{\code{index}}{The articles indexed as 1 through 5.  }
  }
}

\details{All concentrations are in micromolar.   }
\source{ The figure number, year, first author, journal, volume and page number are all included in the dataframe. }
\examples{

## Note that two windows devices will end up exactly on top of each other. 
## Please move device 3 below device 2 to compare their residual plots. 
library(ccems)
## Warning: the next line clears all existing figures!!
if (!is.null(dev.list())) for (i in 2:max(dev.list())) dev.off(i);
if (.Platform$OS.type=="windows") 
 windows(width = 9, height = 4,restoreConsole = TRUE) else X11(width=9,height=4)
par(mfcol=c(2,5))
for (i in 1:5) {
d=subset(TK1,index==i,select=c(E,S,v,frstAut,year))
## for Birringer, 30 uL vessel and 306 pmoles/L= .000306 pmoles/uL 
##   and thus .000306*30= .00918 pmoles/30 uL, i.e. 
if (i==1) d=transform(d,v=v/(60*.00918))# v now in pmoles/sec/pmole = 1/sec
if (i>1) d=transform(d,v=v/(.04*60))# v now in umole/s/umole = 1/sec
plot(d$S,d$v,xlab="Total [dT]",log="xy", ylab="reaction velocity (1/sec)", 
         main=paste(d[1,"frstAut"],d[1,"year"]))
## Note that the specific activity is ~16 fold higher in non-Birringer data
## In the Hill fits V is the rate at which enzyme bound substrate becomes product  
hillda<-nls(v~V*alpha*S^h/(1+alpha*S^h),d,start=list(V=5,alpha=1,h=1))
print(hillda)
## print(summary(hillda))
S50=1/(hillda$m$getPars()["alpha"]^(1/hillda$m$getPars()["h"]))
print(S50)
lgx=log(d$S)
upr=range(lgx)[2]
lwr=range(lgx)[1]
del=(upr-lwr)/50
fineX=exp(seq(lwr,upr,by=del))
lines(fineX,predict(hillda,list(S=fineX)),col="blue",lwd=1)
plot(hillda$m$fitted(),hillda$m$resid(),xlab="fitted value",
                                        ylab="residual",main="Hill")
## Note that variance increases with the mean in non-Birringer data
## and that the 2000 and 1993 Hill fits are poor at low v (and [S])
}

if (.Platform$OS.type=="windows")  # now create a window for ccems fits 
 windows(width = 9, height = 4,restoreConsole = TRUE) else X11(width=9,height=4)
library(ccems)
topology <- list(  
    heads=c("E1S0"), # E1S0 = substrate free E
    sites=list(                    
        c=list(    # c for catalytic site  
            t=c("E1S1","E1S2","E1S3","E1S4")   
        ) # t for tetramer 
    )
)   # in transform below, TK1 is 25kDa => 25mg/umole
g <-mkg(topology, activity=TRUE,TCC=FALSE)

getKk <- function(x) {t(x$report[c(paste("E1S",0:3,"_S",sep=""),
                                 paste("kE1S",1:4,sep="")),"final",drop=FALSE])}
getAIC <- function(x) { x$report["AIC","final"]}
getSSE <- function(x) { x$report["SSE","final"]}

outs=list(NULL)
par(mfcol=c(2,5))
for (i in 1:5) {
d=subset(TK1,index==i,select=c(E,S,v,frstAut,year))
if (i==1) d=transform(d,E=E/4,v=E*v/(60*.00918))# v now in uM/sec
if (i>1) d=transform(d,E=E/4,v=E*v/(.04*60))# v now in uM/sec
plot(d$S,d$v,xlab="Total [dT]",log="xy", ylab="reaction velocity (1/sec)", 
         main=paste(d[1,"frstAut"],d[1,"year"]))
names(d)[1:2]= c("ET","ST")
tops=ems(d,g,maxTotalPs=3,doSpurs=FALSE)# takes ~15 sec for each dataset
lgx=log(d$ST)
upr=range(lgx)[2]
lwr=range(lgx)[1]
del=(upr-lwr)/50
fineX=exp(seq(lwr,upr,by=del))
predict <- data.frame(ET = rep(d$ET[1],length(fineX)), ST = fineX)
df <- simulateData(tops[[1]],predict=predict,typeYP="v")$predict  
lines(df$ST,df$EY) 
Kk=lapply(tops,getKk)
nms=names(Kk)
rowList=data.frame(NULL)
for (j in nms) {
    rowList=rbind(rowList,Kk[[j]])
  }
rownames(rowList)<-nms
aic=sapply(tops,getAIC)
sse=sapply(tops,getSSE)
eDelAIC=exp(-(aic-min(aic)))
wgts=eDelAIC/sum(eDelAIC)
print(sum(wgts))
df=data.frame(aic,sse,wgts,rowList)
M=as.matrix(rowList)
ma=exp(wgts\%*\%log(M)) # average in space of gibbs free energy changes
dataID=paste(d[1,"frstAut"],d[1,"year"],sep="")
outs[[dataID]]$df=df
outs[[dataID]]$ma=ma
plot(tops[[1]]$d$EY,tops[[1]]$res,xlab="fitted value",
                           ylab="residual",main=tops[[1]]$mid)
## Note that the 2000 and 1993 fits are now improved 
}
outs=outs[-1] # remove leading NULL
print(outs)   # compare model averages across datasets
par(mfrow=c(1,1))

}
\keyword{datasets}
