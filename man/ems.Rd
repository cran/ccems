\name{ems}
\alias{ems}
\title{Equilibrium Model Selection}
\description{ This is the main automation function of this package. It generates a space of  
                 combinatorially complex equilibrium models and fits them to data.}
\usage{ems(d, g, cpusPerHost=c("localhost" = 1), ptype="",
          spurChunkSize=1000, nSpurChunks=1,
          maxTotalPs=5,minTotalPs=NULL,extend2maxP=TRUE, 
          smart=FALSE,doTights=FALSE,doGrids=TRUE,
          doSpurs=TRUE,topN=10,showConstr=FALSE,
          atLeastOne=TRUE,atLeastOneOfEach=FALSE,
          KIC=1,kIC=1,fullGrid=FALSE,
        transform=c("boxCox","relResid","none","sqrt","log"),lam=0.5,
        m1=-90,p=-1,forceM1=FALSE,forceP=FALSE)}
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
   \item{minTotalPs}{ The minimum number of parameters of models in the model space.
        If \code{NULL} no minimum is imposed. }
  \item{extend2maxP}{ This logical is \code{TRUE} if \code{nSpurChunks} should be extended 
               (if needed) to reach \code{maxTotalPs}. }
  \item{smart}{ Set to \code{TRUE} to stop when models with \code{lastCompleted} parameters (see \code{\link{mkSpurs}}) 
                have an AIC that is bigger than that of the \code{lastCompleted-1} parameter models. }
  \item{doTights}{ Set  to \code{TRUE} if spur models with infinitely tight binding single edges (with K=0) are wanted in the model space.}
  \item{doGrids}{ Leave \code{TRUE} (the default) if grid models are wanted, set to \code{FALSE} if not (e.g. if only spur models are wanted). }
  \item{doSpurs}{ Leave \code{TRUE} if the spur model space is wanted, set to \code{FALSE} if not (e.g. if only grid models are wanted). }
  \item{topN}{ The number of best models of the current batch of models that will be carried 
               over to compete with the next batch; such carryovers 
               are needed to allow fits of model spaces that are too large to reside in memory at one time. 
               This number  is also the number of best models
               summarized in html in the \code{results} folder after fitting each batch.}
  \item{showConstr}{ Set to \code{TRUE} if constrained (fixed and tracking) parameters are to be included in the html report in \code{results}.}
  \item{atLeastOne}{ Leave \code{TRUE} if only models with at least one complex of maximal size are to be considered. Set \code{FALSE} if there is no
                     prior knowledge supportive of the assertion that the largest oligomer must be in the model.}
  \item{atLeastOneOfEach}{ Set \code{TRUE} if only models with at least one complex 
                    of each oligomer size are to be considered. This is useful when the data are multivariate
                    proportions (i.e. mass distribution data) and each j-mer is clearly present. }
  \item{KIC}{ The initial condition of all K parameters optimized. The default is \code{KIC=1} (in uM). }
  \item{kIC}{ The initial condition of all k parameters optimized. The default is \code{kIC=1} (in 1/seconds per occupied active site). }
  \item{fullGrid}{ Set \code{TRUE} if a full binary K model is wanted, else grids that are 
                   equivalent to spurs are eliminated from the model space. }
  \item{transform}{If not \code{"none"} data and model are transformed before forming residuals. 
                   This is used to stabilize enzyme activity variances. Other options are \code{"boxCox"} for Box-Cox
                   transformations, in which case \code{lam} below is used as lambda, 
                   \code{"relResid"} to divide the residuals by the data, 
                   and square root and natural log transformations using \code{"sqrt"} and \code{"log"}, respectively.}
 \item{lam}{ The lambda parameter of the Box-Cox transformation, if used. }
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
g <- mkg(topology) 
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
