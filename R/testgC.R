`testgC` <-
function(g) {
    times <- seq(0,10,1)
    print("************* start testgC **************************")
#    dyn.load(strn<-paste(g$wDir,"/models/",g$id,.Platform$dynlib.ext,sep=""))
    dyn.load(strn<-paste("models/",g$id,.Platform$dynlib.ext,sep=""))
    out1 <- lsoda(g$initialStateTCC,times,"myderivs", g$parmsTCC, rtol=g$rtol,atol=g$atol, dllname=g$id)
    out=data.frame(out1)
    print(out)
    dyn.unload(strn)
    print("************* end testgC **************************")
}

