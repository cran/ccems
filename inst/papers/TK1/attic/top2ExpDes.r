rm(list=ls(all=TRUE))  # clean up left overs from previous run, though shouldn't need to.
if (.Platform$OS.type=="windows") home="/users/radivot" else home="/home/radivot"
setwd(home)  # place where models and results will go
source(paste(home,"/start.r",sep=""))  # define host="machineName" in this file
if (0) library(ccems) else { # if 0 source in the package to save install time 
  pkgNms=dir(paste(home,"/case/active/ccems/ccems/R",sep=""),full.names=TRUE)
  for (i in pkgNms) source(i)  # source in all of the R files of the package
  load(paste(home,"/case/active/papers/TK1/data/TK1.rda",sep=""))
}

# This was chased for a while as an experimental design strategy in an earlier version of the MS. 
# The idea was to assume one of the top 2 models and fit the other, then assume the other and fit the one.
# This didn't make it into the final paper. The lack of fits seemed too small. 

#for (i in c(.02,.03,.04))
  for (i in c(.07,.06,.04))
    print(power.t.test(n = 100, delta = i, sd = NULL, sig.level = 0.05,
    power = .80,
    type = "one.sample",
    alternative = "two.sided",
    strict = FALSE))


topology <- list(  
    heads=c("E1S0"), # E1S0 = substrate free E
    sites=list(                    
        c=list(    # c for catalytic site  
            t=c("E1S1","E1S2","E1S3","E1S4")   
        ) # t for tetramer 
    )
)   # in transform below, TK1 is 25kDa => 25mg/umole
g <-mkg(topology, activity=TRUE,TCC=FALSE,doSqrt=TRUE)


KS=c("E1S0_S","E1S1_S","E1S2_S","E1S3_S")
kS=c("kE1S1","kE1S2","kE1S3","kE1S4")

options(digits=2)

fineX=seq(.1,2,by=.1)
predict <- data.frame(ET = rep(.0001,length(fineX)), ST = fineX)
mdls=list(NULL)

{nm="DDDM.DDDD";Keq=c(                E1S2_S="E1S0_S",E1S1_S="E1S0_S");keq=c(kE1S2="kE1S1",kE1S3="kE1S1",kE1S4="kE1S1")
  Kic=c( 0.8,0.8,0.8,0.2);                               kic=c(4,4,4,4)}
Kmapping = mkKd2Kj(g)
names(Kic)<-KS
names(kic)<-kS
mdls[[1]] = mkModel(g,nm,Kdparams=Kic, Kd2KjLst=Kmapping,Keq=Keq,kparams=kic,keq=keq)
df <- simulateData(mdls[[1]],predict=predict,typeYP="k")$predict  
names(df)[3]<-"k"
mdls[[1]] = mkModel(g,nm,df,Kdparams=Kic, Kd2KjLst=Kmapping,Keq=Keq,kparams=kic,keq=keq)

{nm="DFFF.DFFF";Keq=c(E1S3_S="E1S1_S",E1S2_S="E1S1_S"                );keq=c(              kE1S3="kE1S2",kE1S4="kE1S2")
#  Kic= c(.63, .8,.8,.8);                               kic=c(2.7,4.8,4.8,4.8)}
  Kic= c(.9, .55,.55,.55);                               kic=c(2,4,4,4)}
Kmapping = mkKd2Kj(g)
names(Kic)<-KS
names(kic)<-kS
mdls[[2]] = mkModel(g,nm,df,Kdparams=Kic, Kd2KjLst=Kmapping,Keq=Keq,kparams=kic,keq=keq)


{nm="DFFF.DFFF";Keq=c(E1S3_S="E1S1_S",E1S2_S="E1S1_S"                );keq=c(              kE1S3="kE1S2",kE1S4="kE1S2")
  Kic= c(.9, .55,.55,.55);                               kic=c(2,4,4,4)}
Kmapping = mkKd2Kj(g)
names(Kic)<-KS
names(kic)<-kS
mdls[[3]] = mkModel(g,nm,Kdparams=Kic, Kd2KjLst=Kmapping,Keq=Keq,kparams=kic,keq=keq)
df <- simulateData(mdls[[2]],predict=predict,typeYP="k")$predict  
names(df)[3]<-"k"
mdls[[3]] = mkModel(g,nm,df,Kdparams=Kic, Kd2KjLst=Kmapping,Keq=Keq,kparams=kic,keq=keq)

{nm="DDDM.DDDD";Keq=c(                E1S2_S="E1S0_S",E1S1_S="E1S0_S");keq=c(kE1S2="kE1S1",kE1S3="kE1S1",kE1S4="kE1S1")
                Kic=c( 0.8,0.8,0.8,0.2);                               kic=c(4,4,4,4)}
Kmapping = mkKd2Kj(g)
names(Kic)<-KS
names(kic)<-kS
mdls[[4]] = mkModel(g,nm,df,Kdparams=Kic, Kd2KjLst=Kmapping,Keq=Keq,kparams=kic,keq=keq)

fmdls=lapply(mdls,fitModel)

bdf=data.frame(fmdls[[2]]$d,fmdls[[4]]$d)
bdf
write.table(bdf,file="TK1table1.txt",sep="\t")



