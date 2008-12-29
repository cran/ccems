\name{fitModel}
\alias{fitModel}
\title{Fit Model}
\description{
  This function fits a model/hypothesis created by \code{mkModel}. 
  It is typically passed to \code{lapply}, \code{clusterApplyLB} or \code{mpi.applyLB} to 
  fit a list of such model objects (a model space) and this is typically done within \code{ems}. 
}
\usage{ fitModel(model) }
\arguments{
  \item{model}{  The output list of \code{mkModel}. }
}
\details{ The main output of this function is the \code{report} component of its value (see below) which is also echoed to the screen during 
           computations.}  
\value{
The input argument \code{model} extended to include the following fields:  
  \item{echk}{A matrix that checks the TCC solver and \code{model$fback}. Matrix column names that end in Q should match their sans-Q counterparts. }
  \item{eSS}{The expected steady state concentrations of species (complexes and free reactants). For each row of the data dataframe there is a row in this
             matrix. Its contents are the TCC solver solution (free reactant expected concentrations) and the result of applying 
             \code{model$fback} to them to create expected complex concentrations.  }
  \item{res}{The residuals of the fit.}
  \item{nData}{The number of data points (i.e. rows) in the data dataframe \code{model$d}.}
  \item{SSE}{The initial and final sum of squared errors (i.e. residual sum of squares).}
  \item{AIC}{The initial and final Akaike Information Criterion values, corrected for small samples. Since nonlinear least squares is used  
        \code{AIC = N*log(SSE/N)+2*P + 2*P*(P+1)/(N-P-1) + N*log(2*pi) + N } where \code{N = nData} and \code{P} is the 
        number of estimated parameters (including the variance).}
  \item{nOptParams}{The number of optimized parameters, i.e. the length of the parameter vector sent to \code{optim}. }
  \item{hess}{This is \code{TRUE} if the determinant of the Hessian of the log-likelihood evaluated at the optimum is greater than zero, 
              i.e. if the hessian can be inverted to create a parameter estimate covariance matrix.  }
  \item{CI}{Confidence intervals. Unlike those in \code{model$report} these are numeric rather than strings and these are 
            not expressed as concentrations raised to integer powers (in cases of complete dissociation constants).}
  \item{cpu}{The amount of computing time (in minutes) taken to fit the model. }
  \item{report}{An extension of \code{model$params} to include parameter point estimates and confidence intervals (see \code{CI} above). The \code{final}
         column holds numerics and the \code{pointEstimate} column holds strings of the same numbers expressed as powers in cases 
         of complete dissociation constants.}
}
\note{This work was supported by the National Cancer Institute (K25CA104791). }

\author{ Tom Radivoyevitch (txr24@case.edu) }
\references{ Radivoyevitch, T. (2008) Equilibrium model selection: dTTP induced R1 dimerization. \emph{BMC Systems Biology} \bold{2}, 15. }
\seealso{\code{\link{mkModel}},\code{\link{ems}},\code{\link{ccems}}}
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
data(RNR)
d1 <- subset(RNR,(year==2001)&(fg==1)&(G==0)&(t>0),select=c(R,t,m,year))
d2 <- subset(RNR,year==2006,select=c(R,t,m,year)) 
dRt <- rbind(d1,d2)
names(dRt)[1:2] <- paste(strsplit(g$id,split="")[[1]],"T",sep="") # e.g. to form "RT"
rownames(dRt) <- 1:dim(dRt)[1] # lose big number row names of parent dataframe

\dontrun{
models <- list(
       mkModel(g,"IIJJ",dRt,Kjparams=c(R2t0=Inf, R1t1=Inf,R2t1=1,   R2t2=1)), 
       mkModel(g,"IIIJ",dRt,Kjparams=c(R2t0=Inf, R1t1=Inf,R2t1=Inf, R2t2=1))
       )
fmodels <- lapply(models,fitModel) # a list of two models are fitted in series here on a single processor 
}
# Note that fitModel always delivers a summary of the fit to the screen as a byproduct. If no assignment is made, the returned large
# fitted list of models makes this summary time consuming to find (i.e. too much scrolling up). 
}
\keyword{ models }% __ONLY ONE__ keyword per line
