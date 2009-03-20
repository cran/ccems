\name{ems}
\alias{ems}
\title{Equilibrium Model Selection}
\description{ This is the main automation function of this package. It generates a space of  
                 combinatorially complex equilibrium models and fits them to data.}
\usage{ems(d, g, cpusPerHost=c("localhost" = 1), ptype="",
          spurChunkSize=1000, nSpurChunks=1,maxTotalPs=5,extend2maxP=TRUE, 
          smart=FALSE,pRows=FALSE,doTights=FALSE,doGrids=TRUE,
          doSpurs=TRUE,topN=5,showConstr=FALSE,atLeastOne=TRUE,
          IC=1,kIC=1,fullGrid=FALSE)}
\arguments{
  \item{d}{ The data as a dataframe.}
  \item{g}{ The list output of \code{mkg}. }
  \item{cpusPerHost}{ This is an integer vector where names are host names and values are their cpu numbers. }
  \item{ptype}{ Parallelization type: \code{""} for single cpus; \code{"SOCK"} and  \code{"NWS"} (networkspaces) 
  for \code{snow} options.   
%  \item{ptype}{ Parallelization type: \code{""} for single cpus; \code{"PVM", "SOCK", "MPI"} and  \code{"NWS"} (networkspaces) 
%  for \code{snow} options; and 
%\code{"RMPI"} for pure \code{Rmpi} using \code{mpi.applyLB} rather than \code{clusterApplyLB}. 
%Batch mode should be used for \code{"RMPI"}/\code{"MPI"} using e.g. 
%\code{ mympi ( ) { orterun -n 1 -bynode -machinefile ~/machines R CMD BATCH $1 rout & } }
%where machines is a file that contains the cluster's machine names (\code{cpusPreHost} is not used in this case). 
%For \code{"PVM"} batch mode should be used via
%\code{ myr ( ) { R CMD BATCH --no-save --no-restore $1 rout & } } since
%R exits (and PVM halts) at the end of \code{ems} anyway. \code{"SOCK"} and  \code{"NWS"} are recommended since they can be used
%interactively on clusters, and since they are easy to set up on both Windows and Linux. \bold{Warning:} as of ROCKS 5.1 and R 2.8.0,  
%\code{"RMPI"}/\code{"MPI"} causes the program to hang.   


}

  \item{spurChunkSize}{ The \code{batchSize} of spur model chunks, see \code{\link{mkSpurs}}}
  \item{nSpurChunks}{The number of spur model chunks requested 
               (this may increase internally if \code{extend2maxP} = \code{TRUE} or \code{smart=TRUE}).}
  \item{maxTotalPs}{ The maximum number of parameters of models that will 
               be fitted (internally, larger models may be generated but not fitted).}
  \item{extend2maxP}{ This logical is \code{TRUE} if \code{nSpurChunks} should be extended 
               (if needed) to reach \code{maxTotalPs}. }
  \item{smart}{ Set to \code{TRUE} to stop when models with \code{lastCompleted} parameters (see \code{\link{mkSpurs}}) 
                have an AIC that is bigger than that of the \code{lastCompleted-1} parameter models. }
  \item{pRows}{ Set to \code{TRUE} if models with estimated inactive protein fractions \code{p} are wanted in the model space, 
                  else \code{p=1} will be fixed for all models generated. }
  \item{doTights}{ Set  to \code{TRUE} if spur models with infinitely tight binding single edges (with K=0) are wanted in the model space.}
  \item{doGrids}{ Leave \code{TRUE} (the default) if grid models are wanted, set to \code{FALSE} if not (e.g. if only spur models are wanted). }
  \item{doSpurs}{ Leave \code{TRUE} if the spur model space is wanted, set to \code{FALSE} if not (e.g. if only grid models are wanted). }
  \item{topN}{ The number of best models of the current batch of models that will be carried over to compete with the next batch; such carryovers 
               are needed to allow fits of model spaces that are too large to reside in memory at one time. This number  is also the number of best models
               summarized in html in the \code{results} folder after fitting each batch.}
  \item{showConstr}{ Set to \code{TRUE} if constrained (fixed and tracking) parameters are to be included in the html report in \code{results}.}
  \item{atLeastOne}{ Leave \code{TRUE} if only models with at least one complex of maximal size are to be considered. Set \code{FALSE} if there is no
                     prior knowledge supportive of the assertion that the largest oligomer must be in the model.}
  \item{IC}{ The initial condition of all K parameters optimized. The default is \code{IC=1} (in uM). }
  \item{kIC}{ The initial condition of all k parameters optimized. The default is \code{kIC=1} (in 1/seconds per occupied active site). }
  \item{fullGrid}{ Set \code{TRUE} if a full binary K model is wanted, else grids that are 
                   equivalent to spurs are eliminated from the model space. }
 }
\details{ This is the highest level function in \code{ccems}. The other functions serve this function, though they may also be used to fit individual 
          models manually.   }
\value{ A list of the \code{topN} best (lowest AIC) models. This should be assigned to a variable 
to avoid large screen dumps. 
 An html report, the topN fitted models, and a brief summary of all fitted models, are saved to 
       \code{results} and are the main outputs and use of this function. }
\author{ Tom Radivoyevitch (txr24@case.edu) }
\references{ Radivoyevitch, T. (2009) Automated model generation and analysis methods 
              for combinatorially complex biochemical equilibriums. (In preparation) }
\note{ Spur and grid graph models have network topologies that either radiate 
       from the hub or can be overlaid on a city block lay out, respectively. 
       Though head node spur graph edges can be superimposed in curtain rods (see \code{\link{ccems}}) 
       to give these graphs a grid appearance, curtain rods are really sets 
       of nested arches. Thus curtains could be called spur-grid hybrid K equality graphs or simply hybrids
       (i.e. a term that is more tolerant than grid). Another option is to tolerate spur 
       edges to head nodes in a 
       broadened definition of the term grid. Advantages include an emphasis on parallel edges and thus 
       equality aspects of the graph (compared to the term hybrid), more compactness 
       (compared to the term K equality) and usage inertia. 
       Readers are thus asked to accept this broadened definition 
       of the term grid, i.e. to allow head node spur edges in grid graphs. 
       
       This work was supported by the National Cancer Institute (K25CA104791).
       }
\seealso{\code{\link{ccems}}, \code{\link{mkg}} }
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
dd <- rbind(d1,d2)
names(dd)[1:2] <- c("RT","tT")
rownames(dd) <- 1:dim(dd)[1] # lose big number row names of parent dataframe
# the call above ends sooner if maxTotalPs is reached
\dontrun{ 
top <- ems(dd,g,maxTotalPs=1)  # this takes roughly one minute }
}
\keyword{ models }  % note: internal makes the html file not show up in 00index.html
