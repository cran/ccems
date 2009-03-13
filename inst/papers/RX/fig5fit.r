if (.Platform$OS.type=="windows") home="/users/radivot" else home="/home/radivot"
source(paste(home,"/start.r",sep=""))  # define host="machineName" in this file
setwd(home)  # place where models and results will go
library(ccems)

# next block makes sure that the right C code is in place from simulateData below
topology <- list(
    heads=c("R1X0","R2X2","R4X4","R6X6"), 
    sites=list(                # s-sites are already filled only in (j>1)-mers 
        a=list(  #a-site                                                    thread
            m=c("R1X1"),                                            # monomer   1
            d=c("R2X3","R2X4"),                                     # dimer     2
            t=c("R4X5","R4X6","R4X7","R4X8"),                       # tetramer  3
            h=c("R6X7","R6X8","R6X9","R6X10", "R6X11", "R6X12")     # hexamer   4
        ), # tails of a-site threads are heads of h-site threads
        h=list(   # h-site
            m=c("R1X2"),                                            # monomer   5
            d=c("R2X5", "R2X6"),                                    # dimer     6
            t=c("R4X9", "R4X10","R4X11", "R4X12"),                  # tetramer  7
            h=c("R6X13", "R6X14", "R6X15","R6X16", "R6X17", "R6X18")# hexamer   8
        )
    )
)
g=mkg(topology,TCC=TRUE) 

# read in the top 5 models
load("case/results/RXglobSOCKtopN.RData") 
mkHTML(globalTopN)  # dump out to make sure it's what you think it is
#globalTopN[[1]]
attach(globalTopN[[1]]$d) # grab data from first model
plot(XT,m,type="p",pch=1, xlab="[ATP] (uM)", ylab="Mass (kDa)")
lgx=log(XT)
upr=range(lgx)[2]
lwr=range(lgx)[1]
del=(upr-lwr)/50
fineX=exp(seq(lwr,upr,by=del))
newPnts <- data.frame(RT = rep(RT[1],length(fineX)), XT = fineX)
newPnts
# plot top 3
for (i in 1:3) {
  df <- simulateData(globalTopN[[i]],predict=newPnts,typeYP="m")$predict  
  lines(df$XT,df$EY,type="l",lty=i)
}
legend(2000,200,c("Best", "Second Best", "Third Best"),lty=1:3)
detach(globalTopN[[1]]$d)
dev.copy2pdf(file="results/fig5fit.pdf")
