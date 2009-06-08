\name{simulateData}
\alias{simulateData}
\title{ Simulate Data }
\description{
  This function generates expected values of responses at the total concentrations 
  of the dataframe \code{g$d}  
  or at points specified in \code{predict} if \code{predict} is not \code{NULL}.
  }
\usage{ simulateData(g, init = FALSE, predict = NULL, typeYP = NULL) }
\arguments{
  \item{g}{ A specific model/hypothesis created by \code{mkModel}. }
  \item{init}{ This is \code{TRUE} only in first calls to this function by the parameter estimate optimization algorithm. 
               When \code{TRUE} the initial AIC value is set. }
  \item{predict}{ A dataframe of total concentrations of the reactants 
  at which response predictions are desired.  }
  \item{typeYP}{ The type of output data desired for the predictions. 
  Options are "m" and "v" for average mass and reaction velocity, respectively.
  The option "P" (proportions) for mass distribution data is still in development.   }
}
\details{
  This function is the workhorse core of the nonlinear least squares 
  algorithm, so its speed is 
  critical which is why it uses compiled C code when \code{g$TCC = TRUE}. 
  In addition to model fitting, this function can also be used to predict 
  system response surfaces over grids of
  physiologically relevant total concentrations of the reactants. 
  It is assumed that the model used to formulate predictions is based on data, and that
  the output type of the predictions is the same as the output type used to build the model. 
}
\value{ 
The input model object augmented to include the following fields if \code{predict} is \code{NULL}.
  \item{echk}{A matrix that checks the TCC solver and \code{g$fback}. Matrix column names 
  that end in Q should match their sans-Q counterparts. }
  \item{eSS}{The expected steady state concentrations of species (complexes and free reactants). 
  For each row of the data dataframe there is a row in this
             matrix. Its contents are the TCC solver solution (free reactant 
             expected concentrations) and the result of applying 
             \code{g$fback} to them to create expected complex concentrations.  }
  \item{res}{The residuals of the fit.}
  \item{nData}{The number of data points (i.e. rows) in the data dataframe \code{g$d}.}
  \item{SSE}{The initial and final sum of squared errors (i.e. residual sum of squares).}
  \item{AIC}{The initial and final Akaike Information Criterion values, corrected for small samples. S
     ince nonlinear least squares is used  
        \code{AIC = N*log(SSE/N)+2*P + 2*P*(P+1)/(N-P-1) + N*log(2*pi) + N } where 
        \code{N = nData} and \code{P} is the 
        number of estimated parameters (including the variance).}
If \code{predict} is not \code{NULL} the value is the input \code{g} augmented only by the field 
   \item{predict}{The input argument \code{predict} with an additional 
         expected system response column named "EY".}
}
\note{The function \code{fitModel} augments the input model object by the same six fields above because 
      it calls this function iteratively. 

  Measurements are often made at total concentrations that are substantially 
  higher than physiological values due to 
  signal-to-noise limitations. Thus, predictions in physiologically relevant 
  (and thus important) regions tend to be weak. 
       
       This work was supported by the National Cancer Institute (K25CA104791). 
       }
\author{ Tom Radivoyevitch (txr24@case.edu) }
\references{ Radivoyevitch, T. (2008) Equilibrium model selection: dTTP induced R1 dimerization. 
             \emph{BMC Systems Biology} \bold{2}, 15. }
\seealso{The experimental design example \code{expDesign} in the \code{docs} directory. }
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
g <- mkg(topology,TCC=TRUE) 
d=subset(RNR,(year==2001)&(fg==1)&(t>0)&(G==0),select=c(R,t,m,year))
names(d)[1:2] <-c("RT","tT") 
mdl=mkModel(g,"IIIJ",d,Kjparams=c(R2t0=Inf, R1t1=Inf,R2t1=Inf, R2t2=1),
                     pparams=c(p=1))
fmdl <- fitModel(mdl)
pt=c(.1,1:20)
predict <- data.frame(RT = rep(7.6,length(pt)), tT = pt)
df <- simulateData(fmdl,predict=predict,typeYP="m")$predict  
plot(d$tT,d$m,type="p",  xlab="[dTTP] (uM)", ylab="Weight averaged R1 mass", 
     main="Scott et al. Biochemistry, 2001, Fig. 1 (DLS data)")
lines(df$tT,df$EY) 

}
\keyword{math}
