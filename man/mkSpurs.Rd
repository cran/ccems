\name{mkSpurs}
\alias{mkSpurs}
\title{ Make Spur Model Space }
\description{
  This function takes \code{g} from \code{mkg} and maps it into a dataframe of spur graph model definitions.
}
\usage{ mkSpurs(g,state=list(globMdlIndex=0,globCmbIndex=0,relCmbIndex=0,config=NULL), maxnPs=NULL, 
          batchSize=500, pRows=FALSE, doTights=FALSE,atLeastOne=TRUE,IC=1) }
\arguments{
  \item{g}{ The generic model output list of \code{mkg}. }
  \item{state}{ The current state of model space generation. Its components are: 
            \code{globMdlIndex} which is the global model space index 
            (when \code{mkSpurs} is called within \code{ems} spurs begin where grids end); 
            \code{globCmbIndex} which is the global spur model index 
            (this index is used to trim the last \code{batchSize}, if needed, to avoid attempts
              to extend the spur space size beyond its upper limit of \code{2^g$nZ}); and 
            \code{relCmbIndex} which is the most critical component (this is the column number of
             the current matrix output of \code{combn}, i.e. it is a relative index.  
             The spur space \code{chunk} returned by this function begins just after this column.); 
            \code{config} is a vector of the integer positions in \code{g$Z} of the 
            last model's finite K (i.e. the last model of the previous
            chunk). \code{config}  is 
            the \code{relCmbIndex} column of the current \code{combn} matrix. 
            Its length is the current number of K parameters in the model unless \code{relCmbIndex=0}, 
            in which case the number of 
            parameters is one more than this 
            (in this case the end of the last batch coincides with the end of a \code{combn} matrix). }
  
  \item{maxnPs}{ The maximum number of parameters of models in the model space. Full chunks 
  are created and then trimmed, so this cannot solve 
                  "out of memory" problems. }
  \item{batchSize}{ This is the number of K infinity models fitted per batch. Chunk sizes 
                    are bigger than this if \code{doTights} is \code{TRUE} and/or  
                    if \code{pRows} is \code{TRUE}. Values less than ~1000 are recommended. }
  \item{pRows}{ Set to \code{TRUE} if additional rows are to be added to the output dataframe 
                to include models with \code{p} freely estimated. }
  \item{doTights}{ This should be \code{TRUE} if infinitely tight binding models 
                   of single edge spur graphs are to be created.}
   \item{atLeastOne}{ Leave \code{TRUE} if only models with at least one complex of maximal size are to be considered. Set \code{FALSE} if there is no
                     prior knowledge supportive of the assertion that the largest oligomer must be in the model.}
   \item{IC}{ The initial condition of all K parameters optimized. The default is \code{IC=1}. }
                           
}
\details{ This function is complicated by the fact that one readily runs out of memory 
with 29 complexes and thus roughly 500,000,000 spur models (in this case 8 GB RAM allows at most 16 bytes per model!).  
Thus, chunks of the spur model space must be created, fitted and summarized in sizes small enough to fit into memory. 
As both an input and an output, \code{state} links successive calls to this function. 
It keeps track of where we are in the spur graph model space and it allows searches through 
the low parameter number models without first defining all of the higher
parameter number models (and thus consuming all of the RAM in the process). }

\value{A list with components 
  \item{chunk}{A dataframe where each row is a spur model. }
  \item{state}{The state, defined in the same way as the input argument \code{state} since the output \code{state} of one call is
               the input \code{state} of the next call.}
  \item{maxReached}{This is \code{TRUE} if the maximum number of parameter has been reached.}
  \item{lastCompleted}{When all of the models with  \eqn{j} K parameters have been specified, \code{lastCompleted} equals \eqn{j}.  
                       This is useful in \code{ems} when \code{smart = TRUE} as it defines how far into the model space AICs should be tested 
                       for increases due to overparameterization.    }
}

\note{This work was supported by the National Cancer Institute (K25CA104791).}
\author{ Tom Radivoyevitch (txr24@case.edu) }
\references{ Radivoyevitch, T. (2008) Equilibrium model selection: dTTP induced R1 dimerization. \emph{BMC Systems Biology} \bold{2}, 15. }
\seealso{\code{\link{mkGrids}}, \code{\link{ccems}}, \code{\link{combn}} }
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
mkSpurs(g,pRows=TRUE,doTights=TRUE)
}
\keyword{ models}
