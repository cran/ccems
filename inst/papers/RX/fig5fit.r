if (.Platform$OS.type=="windows") home="/users/radivot" else home="/home/radivot"
source(paste(home,"/start.r",sep=""))  # define host="machineName" in this file
setwd(home)  # place where models and results will go
library(ccems)
load("case/results/RXglobSOCKtopN.RData") 
mkHTML(globalTopN)

globalTopN[[1]]
n=length(globalTopN)
attach(globalTopN[[1]]$d)
pdf("case/results/fig5fit.pdf")
plot(XT,m,type="p",pch=1, xlab="[ATP] (uM)", ylab="Mass (kDa)")
lgx=log(XT)
upr=range(lgx)[2]
lwr=range(lgx)[1]
del=(upr-lwr)/200
fineX=exp(seq(lwr,upr,by=del))
newPnts <- data.frame(RT = rep(RT[1],length(fineX)), XT = fineX)
newPnts
for (i in 1:n) {
  df <- simulateData(globalTopN[[1]],predict=newPnts,typeYP="m")$predict  
  lines(df$XT,df$EY,type="l",lty=i)
}
#  legend(2000,200,c("Best", "Second Best", "Third Best"),lty=1:n, col = 1:n)
detach(globalTopN[[1]]$d)
dev.off()
