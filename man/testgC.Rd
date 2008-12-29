\name{testgC}
\alias{testgC}
\title{ Test g Object TCC C Compilation }
\description{ This function runs \code{lsoda} using the compiled C code of the total concentration constraint (TCC) 
              system equations that should have been created by \code{mkg}. }
\usage{ testgC(g) }
\arguments{
  \item{g}{ The output of \code{mkg}. Since \code{mkg} already calls this function, it has no independent use. }
}
\details{
  Components of \code{g} used by this function are: \code{g$wDir} which is the location of the \code{models} directory that holds the compiled code; 
  \code{g$parmsTCC} which includes the total concentrations, the fraction of active protein p, and the spur graph complete dissociation
  constant parameter values Kj; and the \code{lsoda} ODE integration
  tolerance parameters \code{g$rtol} and \code{g$atol}.
}
\value{ No value returned. Ten time points are printed to the screen. The free concentrations should show signs of convergence toward a steady state. }
\author{ Tom Radivoyevitch }
\seealso{ \code{\link[odesolve]{lsoda}} }
\examples{
library(ccems)
topology <- list(  
        heads=c("R1t0","R2t0"),  
        sites=list(       
                s=list(                     # s-site    thread #
                        m=c("R1t1"),        # monomer      1
                        d=c("R2t1","R2t2")  # dimer        2
                )
        )
) 
g=mkg(topology,TCC=TRUE) 
testgC(g) # redundantly, this checks it again  ... there should never be a reason to call this function.
}
\keyword{ internal }% __ONLY ONE__ keyword per line
