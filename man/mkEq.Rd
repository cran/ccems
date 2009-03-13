\name{mkEq}
\alias{mkEq}
\title{ Make Equality Constraints Lists }
\description{ This function converts a dataframe/chunk of models with either k or K equalities
 into a list structure that maps the trackers onto the independent parameter estimate. }
\usage{ mkEq(g,chunk,activity=FALSE) }
\arguments{
  \item{g}{ The output of \code{mkg}. }
  \item{chunk}{ The \code{chunk} component of the output of \code{mkGrids} and \code{mkSpurs}.  }
  \item{activity}{ This is \code{TRUE} if the desired output is \code{keq} and \code{FALSE} if it is \code{Keq}.   }
}
\details{ This function is already called internally by \code{mkGrids} and \code{mkSpurs}. 
               Thus, the user should not need to call it. Note that the input \code{activity} can be false
               even if the dataset is activity data (i.e. \code{activity} in \code{g} is \code{TRUE}) since 
               even then this function may be called to create K equalities rather than k equalities. }
\value{Lists of either k or K equality constraints  }

\author{ Tom Radivoyevitch (txr24@case.edu) }
\references{ Radivoyevitch, T. (2008) Equilibrium model 
            selection: dTTP induced R1 dimerization. \emph{BMC Systems Biology} \bold{2}, 15. }
\seealso{\code{\link{mkGrids}}, \code{\link{mkSpurs}} }
\keyword{ internal }
