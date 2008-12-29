\name{mkHTML}
\alias{mkHTML}
\title{ Make HTML File }
\description{ Writes the best models to an html summary file in the \code{results} directory. }
\usage{ mkHTML(smodels,showConstr=FALSE) }
\arguments{
  \item{smodels}{ A list of fitted model objects, preferably already sorted by best AIC.  }
  \item{showConstr}{ This is \code{TRUE} if constrained (fixed or tracking) parameters are to be included 
                      in the html output file. This is used for error checking.   }
}
\details{ This function is already called internally by \code{ems} after each chunk of models is fitted. Thus, the user should not need to call it. }
\value{No value. The output is an html file in the \code{results} directory. }

\author{ Tom Radivoyevitch (txr24@case.edu) }
\references{ Radivoyevitch, T. (2008) Equilibrium model selection: dTTP induced R1 dimerization. \emph{BMC Systems Biology} \bold{2}, 15. }
\seealso{\code{\link{mkModel}}, \code{\link{ems}} }
\keyword{ internal }
