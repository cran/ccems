rm(list=ls(all=TRUE))  # clean up left overs from previous run, though shouldn't need to.
if (.Platform$OS.type=="windows") home="/users/radivot" else home="/home/radivot"
setwd(home)  # place where models and results will go
source(paste(home,"/start.r",sep=""))  # define host="machineName" in this file
if (1) library(ccems) else { # if 0 source in the package to save install time 
  pkgNms=dir(paste(home,"/case/active/ccems/ccems/R",sep=""),full.names=TRUE)
  for (i in pkgNms) source(i)  # source in all of the R files of the package
  load(paste(home,"/case/active/papers/TK1/data/TK1.rda",sep=""))
}

topology <- list(  
    heads=c("E1S0"), # E1S0 = substrate free E
    sites=list(                    
        c=list(    # c for catalytic site  
            t=c("E1S1","E1S2","E1S3","E1S4")   
        ) # t for tetramer 
    )
)   # in transform below, TK1 is 25kDa => 25mg/umole
g <-mkg(topology, activity=TRUE,TCC=FALSE)

KS=c("E1S0_S","E1S1_S","E1S2_S","E1S3_S")
kS=c("kE1S1","kE1S2","kE1S3","kE1S4")

load("case/results/twouts")
nu2=twouts
for (i in 1:length(nu2))
  for (j in 1:length(nu2[[i]])) {
    nu2[[i]][[j]]$df=nu2[[i]][[j]]$df[1:5,]
  }
nu2

# the line up here is perfect, also at the level of SSEs and AICs
# below, though SSEs and AICs part, the model probs track very well
for (j in 1:length(nu2[["none"]])) {
  print(nu2[["none"]][[j]]$df)
  print(nu2[["none"]][[j]]$ma)
  print(nu2[["lam1"]][[j]]$df)
  print(nu2[["lam1"]][[j]]$ma)
}

# in this comparison the 1993 data yields only slight differences
for (j in names(nu2[["none"]])) {
  print(j)
  print(nu2[["sqrt"]][[j]]$df)
  print(nu2[["sqrt"]][[j]]$ma)
  print(nu2[["lam0p5"]][[j]]$df)
  print(nu2[["lam0p5"]][[j]]$ma)
}

# in this comparison not the substantial difference in Li et al. 
for (j in names(nu2[["none"]])) {
  print(j)
  print(nu2[["log"]][[j]]$df)
  print(nu2[["log"]][[j]]$ma)
  print(nu2[["lam0"]][[j]]$df)
  print(nu2[["lam0"]][[j]]$ma)
}

# this shows that most model averages are qualitatively consistent across lambda values.
# the only exception is in 1993 with a log transform differing from the other two lams
for (j in names(nu2[["none"]])) {
  print(j)
  print(nu2[["lam0"]][[j]]$df)
  print(nu2[["lam0"]][[j]]$ma)
  print(nu2[["lam0p5"]][[j]]$df)
  print(nu2[["lam0p5"]][[j]]$ma)
  print(nu2[["lam1"]][[j]]$df)
  print(nu2[["lam1"]][[j]]$ma)
}

#  note here that relative residuals are a bit like the log transform in that they perturb the 1993
# dataset much more than the others. 
for (j in names(nu2[["lam0p5"]])) {
  print(j)
  print(nu2[["relResid"]][[j]]$df)
  print(nu2[["relResid"]][[j]]$ma)
  print(nu2[["lam0p5"]][[j]]$df)
  print(nu2[["lam0p5"]][[j]]$ma)
}
