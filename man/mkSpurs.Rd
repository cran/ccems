\name{mkSpurs}
\alias{mkSpurs}
\title{ Make Spur Model Space }
\description{
  This function takes \code{g} from \code{mkg} and maps it into a dataframe of spur graph model definitions.
}
\usage{mkSpurs(g,state=list(globMdlIndex=0,globCmbIndex=0,relCmbIndex=0,
         config=NULL), maxnKjPs=NULL, maxTotalPs=NULL,minTotalPs=NULL, 
         batchSize=500,doTights=FALSE,
         atLeastOne=TRUE,atLeastOneOfEach=FALSE,
         KIC=1,kIC=1,m1=-90,p=-1,forceM1=FALSE,forceP=FALSE) }

\arguments{
  \item{g}{ The generic model output list of \code{mkg}. }
  \item{state}{ The current state of model space generation. This is a list with several components. 
            \code{globMdlIndex} is the global model space index. 
            When \code{mkSpurs} is called within \code{ems} the \code{globMdlIndex} of 
            the spurs begins where it ends for grids. 
            \code{globCmbIndex} is the global spur model index. 
            This index is used to trim the last \code{batchSize}, if needed, to avoid attempts
              to extend the spur space size beyond its upper limit of \code{2^g$nZ}. 
            \code{relCmbIndex} is the most critical component of the \code{state}. 
            This is the column number of
             the current matrix output of \code{combn}, i.e. it is a relative index.  
            The spur space \code{chunk} returned by \code{mkSpurs} begins just after this column. 
            \code{config} is a vector of the integer positions in \code{g$Z} of the 
            last model's finite K (i.e. the last model of the previous
            chunk). \code{config}  is 
            the \code{relCmbIndex}th column of the current \code{combn} matrix. 
            Its length is the current number of K parameters in the model unless \code{relCmbIndex=0}, 
            in which case the number of 
            parameters is one more than this 
            (in this case the end of the last batch coincides with the end of a \code{combn} matrix). }
  
  \item{maxnKjPs}{ The maximum number of Kj parameters of models in the model space. Full chunks 
  are created and then trimmed, so descreases in the value of this option will not solve 
                  "out of memory" problems. }
  \item{maxTotalPs}{ This is the maximum number of freely estimated k or K parameters.  }
  \item{minTotalPs}{ The minimum number of parameters of models in the model space.
        If \code{NULL} no minimum is imposed. }
  \item{batchSize}{ This is the number of K infinity models fitted per batch. Chunk sizes 
                    are bigger than this if \code{doTights} is \code{TRUE} and/or  
                    if \code{pRows} is \code{TRUE}, and/or if activity parameter constraints
                    split models further. Values less than ~1000 are recommended for 
                    quad core 8 GB motherboards. }
  \item{doTights}{ This should be \code{TRUE} if infinitely tight binding models 
                   of single edge spur graphs are to be created.}
   \item{atLeastOne}{ Leave \code{TRUE} if only models with at least one complex 
                    of maximal size are to be considered. Set \code{FALSE} if there is no
                     prior knowledge supportive of the assertion that the largest oligomer must be in the model.}
   \item{atLeastOneOfEach}{ Set \code{TRUE} if only models with at least one complex 
                    of each oligomer size are to be considered. This is useful when the data are multivariate
                    proportions (i.e. mass distribution data) and each j-mer is clearly present. }
  \item{KIC}{ The initial condition of all K parameters optimized. The default is \code{IC=1} (in uM). }
  \item{kIC}{ The initial condition of all k parameters optimized. The default is \code{kIC=1} (in 1/seconds per occupied active site). }
  \item{m1}{The hub protein's monomer mass in kDa. 
               The default is 90 for the big (R1) subunit of ribonucleotide reductase (RNR).
               This only matters if the data is mass data. Negative numbers imply fixed values and
                positive numbers imply starting values to be fitted to the data.}
  \item{p}{ Probability that hub can oligomerize, i.e. is not damaged. 
      Set to a positive value if additional rows are to be added to the output dataframe to include models 
               with \code{p} freely estimated. Set negative to hold fixed. Value is the initial or fixed value.}
  \item{forceM1}{ Set \code{TRUE} to force all models to estimate M1, i.e. to not generate models with M1 fixed. } 
  \item{forceP}{ Set \code{TRUE} to force all models to estimate p, i.e. to not generate models with fixed p. } 
                
}
\details{ This function is complicated by the fact that one readily runs out of memory 
with 29 complexes and thus roughly 500,000,000 spur models (in this case 8 GB RAM 
allows at most 16 bytes per model!).  
Thus, chunks of the spur model space must be created, fitted and summarized in sizes 
small enough to fit into memory. 
As both an input and an output, \code{state} links successive calls to this function. 
It keeps track of where we are in the spur graph model space and it allows searches through 
the low parameter number models without first defining all of the higher
parameter number models (and thus consuming all of the RAM in the process). }

\value{A list with components 
  \item{chunk}{A dataframe where each row is a spur model. }
  \item{state}{The state, defined in the same way as the input argument \code{state} 
                since the output \code{state} of one call is
               the input \code{state} of the next call.}
  \item{maxReached}{This is \code{TRUE} if the maximum number of parameter has been reached.}
  \item{lastCompleted}{When all of the models with  \eqn{j} K parameters have 
  been specified, \code{lastCompleted} equals \eqn{j}.  
 This is useful in \code{ems} when \code{smart = TRUE} as it defines how far into the 
 model space the fitting process has gone so far in terms of numbers of model parameters. }
}

\note{This work was supported by the National Cancer Institute (K25CA104791).}
\author{ Tom Radivoyevitch (txr24@case.edu) }
\references{ Radivoyevitch, T. (2008) Equilibrium model selection: dTTP 
   induced R1 dimerization. \emph{BMC Systems Biology} \bold{2}, 15. }
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
g <- mkg(topology) 
mkSpurs(g,p=0.95,doTights=TRUE)
}
\keyword{ models}
