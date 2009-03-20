\name{mkGrids}
\alias{mkGrids}
\title{ Make Grid Model Space }
\description{ This function takes \code{g} from \code{mkg} and maps it into a 
              dataframe of grid model definitions. It also augments and returns the input list \code{g}. }
\usage{ mkGrids (g,maxTotalPs=NULL,pRows=FALSE,contig=TRUE, atLeastOne=TRUE, 
                                 IC=1,kIC=1,fullGrid=FALSE) }
\arguments{
  \item{g}{ The output of \code{mkg}. }
  \item{maxTotalPs}{ The maximum number of parameters of models in the model space. 
        If \code{NULL} all models are found (but see \code{fulGrid} below). }
  \item{pRows}{ Set to \code{TRUE} if additional rows are to be added to the output dataframe to include models 
               with \code{p} freely estimated. }
  \item{contig}{ Set to \code{TRUE} to allow Kd equalities of threads only in contiguous runs. Tentatively, 
  this should always be \code{TRUE}.}
  \item{atLeastOne}{ Leave \code{TRUE} if only models with at least one complex of maximal size are 
            to be considered. Set \code{FALSE} if there is no
                     prior knowledge supportive of the assertion that the largest oligomer must be in the model.}
  \item{IC}{ The initial condition of all K parameters optimized. The default is \code{IC=1} (in uM). }
  \item{kIC}{ The initial condition of all k parameters optimized. The default is \code{kIC=1} (in 1/seconds per occupied active site). }
   \item{fullGrid}{ Set \code{TRUE} if full binary K models are wanted, e.g. to obtain grid 
   reparameterizations of full
   spur models to form model averages across common binary K parameters. }
}
\details{In a run of equal threads, the first head node of the run is the leader (optimized parameter) and the remaining nodes in the run are followers 
(i.e. constrained to track the leader in parameter estimate optimizations); 
a default in \code{ems} is that only leader estimates are reported in html
    in the \code{results} subdirectory.
	In contrast to its counterpart \code{\link{mkSpurs}}, a means of traversing
	the grid model space incrementally with increasing numbers of parameters remains to be 
	found and implemented, i.e. 
	\code{mkGrids} does not have \code{state} inputs and outputs and the whole space is found in one batch.  
	When this drawback 
	limits research due to too much memory usage, 
	attempts will be made to identify an approach that, similar to what has been implemented for 
	\code{mkSpurs}, specifies
	chunks of grid models, fits them, and then specifies the next chunk 
	based on knowledge of where the previous chunk stopped. 
	 Note that setting 
	\code{maxTotalPs} to smaller values will not help this foreseen memory problem 
	as the entire grid model space dataframe is first generated 
	and then, only later, truncated to \code{maxTotalPs}.
}

\value{A list with components 
  \item{chunk}{The entire K equality model space requested by the arguments. 
  In this dataframe each row specifies a model. If the \code{activity} field of \code{g} is \code{TRUE} 
  this dataframe includes k columns. The row names encode the equality constraints. 
  Therein . separates K models from k models, I stands for infinity, 
  J stands for freely estimated (in spur components), and other letters are the same when 
  parameters that correspond to their positions equal each other.
  Unmatched other letters are freely estimated and thus just like J's but in
  k and binary K components of the model names. 
  }
  \item{Keqs}{A list of K equality constraints indexed by 
  model names where each element is a vector of character strings whose names are followers and
             values are leaders. }
 \item{keqs}{Similar to \code{Keqs} but for activities parameter constraints. This is \code{NULL} if 
    \code{activity} in \code{g} is \code{TRUE}. }
}

\note{This work was supported by the National Cancer Institute (K25CA104791).}
\author{ Tom Radivoyevitch (txr24@case.edu) }
\references{ Radivoyevitch, T. (2008) Equilibrium model selection: dTTP induced R1 dimerization. 
            \emph{BMC Systems Biology} \bold{2}, 15. }
\seealso{\code{\link{mkSpurs}},  \code{\link{ccems}} }
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
gridL <- mkGrids(g) 
print(gridL$chunk)
topology <- list(
    heads=c("R1X0","R2X2","R4X4","R6X6"), 
    sites=list(                # s-sites are already filled only in (j>1)-mers 
        a=list(  #a-site                                                    thread
            m=c("R1X1"),                                            # monomer   1
            d=c("R2X3","R2X4"),                                     # dimer     2
            t=c("R4X5","R4X6","R4X7","R4X8"),                       # tetramer  3
            h=c("R6X7","R6X8","R6X9","R6X10", "R6X11", "R6X12")     # hexamer   4
        ), # tails of a-site threads are heads of h-site threads
        h=list(   # h-site
            m=c("R1X2"),                                            # monomer   5
            d=c("R2X5", "R2X6"),                                    # dimer     6
            t=c("R4X9", "R4X10","R4X11", "R4X12"),                  # tetramer  7
            h=c("R6X13", "R6X14", "R6X15","R6X16", "R6X17", "R6X18")# hexamer   8
        )
    )
)


g <- mkg(topology,TCC=TRUE) 
gridL <- mkGrids(g,maxTotalPs=2) 
print(gridL$chunk)

gridL <- mkGrids(g,maxTotalPs=4) 
# the next line should be run separately since its output is large.
# print(gridL$chunk)
}
\keyword{ models }
