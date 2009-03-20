# this short script generates and fits the Rt model space automatically 
setwd("/users/radivot")
library(ccems)
topology=list(  
    heads=c("R1t0","R2t0"),
    sites=list(
        s=list(
            m=c("R1t1"),
            d=c("R2t1","R2t2")
        )
    )
) 
g=mkg(topology) 
d1=subset(RNR,(year==2001)&(fg==1)&(G==0)&(t>0),select=c(R,t,m,year))
d2=subset(RNR,year==2006,select=c(R,t,m,year)) 
dd=rbind(d1,d2) 
names(dd)[1:2]=c("RT","tT")
# the next line takes <13 minutes 
tops=ems(dd,g,pRows=T,doTights=T)
