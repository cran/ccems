\name{mkSingleThread}
\alias{mkSingleThread}
\title{ Make Equality Constraints Lists }
\description{ This function creates k or K equality matrices for single threads.} 
\usage{ mkSingleThread(nNodes,g) }
\arguments{
  \item{nNodes}{ The number of nodes in the thread.  }
  \item{g}{ The output of \code{mkg}. This is used only to grab the \code{mylets} field. }
}
\details{ This function is called internally by \code{mkGrids} and \code{mkSpurs}. 
          Thus, the user should not need to call it. }
\value{Matrix of characters of contiguous blocks of K or k equality constraints within a thread.  }

\author{ Tom Radivoyevitch (txr24@case.edu) }
\seealso{\code{\link{mkGrids}}, \code{\link{mkSpurs}} }
\keyword{ internal }
