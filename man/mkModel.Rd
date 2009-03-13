\name{mkModel}
\alias{mkModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Make Specific Model }
\description{ This function takes a generic model input list \code{g} and augments it to include data and a specific hypothesis. The hypotheses come in the
              form of claims that certain Kj are so large that the data cannot discriminate them from being infinity, that certain Kd are so close in value
              that the data cannot distinguish them from being equal, and that the protein proportion that is active is so 
              close to 1 that the data cannot discriminate it from being 1. 
              Hypotheses regarding enzyme activity parameters k being essentially 
              equal to others on an average per occupied catalytic site basis can also be made.  }

\usage{ mkModel(g,mid,d=NULL,Kjparams=NULL,Kdparams=NULL,Keq=NULL,
        Kd2KjLst=NULL, pparams=c(p=1i), kparams=NULL,
        keq=NULL, tightLogic=TRUE, indx=NULL, nParams=NULL) }
\arguments{
  \item{g}{ The output of \code{mkg}.  }
  \item{mid}{ The name of the specific hypothesis/model. The convention is that I stands for infinity, 
    J for a freely estimated spur graph edge, H for a freely estimated grid head node, 
    and D, E, F, etc. (i.e. other characters) for grid Kd that are equal
    to each other. In the latter case the same letter in two different positions 
    indicates equality between them; positions in model name strings are the binary reaction product
    node positions in \code{g$Z}, see \code{\link{mkg}}). When k constraints exist they preceded by the 
    K name (a period separates them) and they follow  
    single thread curtain conventions (k constraints 
    for larger curtains remain to be worked out).
    The model name, which becomes \code{g$mid},
    should not be confused with the name of the biochemical system \code{g$id}.    }
  \item{d}{ The data as a dataframe.  }
  \item{Kjparams}{ If the hypothesis is a spur model, this is a numeric 
  vector of its initial complete dissociation constant parameter values.     }
  \item{Kdparams}{ If the hypothesis is a generalized grid model these 
           are the initial dissociation constant parameter values: 
                    head node spur edges are distinguished from thread 
                    edges by not having \code{"_"} in their names.   }
  \item{Keq}{ This character vector specifies which Kd parameters are equal to each other. Names are followers 
               and values leaders in the sense of parameters constrained to track each other.  }
  \item{Kd2KjLst}{ This is a list of functions (see \code{\link{mkKd2Kj}}) that maps generalized 
                     grid Kd parameters into full spur model Kj parameters. 
                   An appropriate component of this list is assigned to \code{g$Kd2Kj}.
                    Such functions are needed so that one generic full spur graph
                   model, typically compiled in C, can be used by all of the specific hypotheses/models of the model space. }
  \item{pparams}{ This is the fraction of hub protein that is active. The pure imaginary default of \code{1i} keeps it fixed at 1 (minus 1 also achieves this). 
                 An initial value of +1 indicates that it is to be optimized.     }
  \item{kparams}{ These are the enzyme activity parameter initial values, 
  if the data is reaction rate data, else it should be \code{NULL}.   }
  \item{keq}{ These are the equality constraints (if any) that are being placed on the enzyme activity parameters.  }
  \item{tightLogic}{ If \code{tightLogic} is true, instead of taking Kj to .001 as an approximation of infinitely tight binding, 
                              logic is used to model Kj = 0 exactly. }
  \item{indx}{ This is an integer index of the model.  The current
               hypothesis is the \code{indx}th model of the model space.  }
  \item{nParams}{ The number of model parameters, i.e. the first column of a \code{chunk} data frame. }
}
\details{
Infinite initial Kj parameters remain fixed at \code{Inf} and are passed from R to C properly to eliminate corresponding polynomial terms in the 
total concentration constraints.
}
\value{
The input object \code{g} augmented to include the arguments \code{d}, \code{mid} and \code{indx} and the following:  
  \item{params}{A dataframe specification of the parameter's initial values and whether they are optimized or fixed or 
  constrained to track others. Final value placeholders are initialized to initial values. }
  \item{Kparams}{This is either \code{Kjparams} or \code{Kdparams}. 
  It is the one of the two which is not \code{NULL}. In cases of hybrids \code{Kdparams} is used.}
  \item{codeS}{This is the name of the component of \code{Kd2KjLst} that is relevant to the current hypothesis. It 
     is a string of digits whose binary representation indicates which threads are infinite.  }
  \item{Kd2Kj}{This is the component of \code{Kd2KjLst} that is 
               relevant to the current hypothesis, i.e. \code{Kd2Kj = Kd2KjLst[[codeS]]}. }
%  \item{fitid}{Initialized to "noID", this will later become a collapsed string of the biochemical system name, the
%               number of full model complexes, and, separated by periods, the number of non-p and p optimized 
%               parameters.   }
  \item{fitS}{A string indicator of the status of the model fitting. It is initialized here to \code{"not fitted yet"}. }
  \item{typeY}{The type of output. This is \code{m} for average mass and \code{v} for reaction velocity. }
  \item{posY}{This is the column number of the output measurement in the data dataframe \code{d}. }
  \item{posReactantsD}{These are the column numbers of the total reactant 
  concentrations (system inputs) in the data dataframe \code{d}.}
The value returned by this function is a model object that is ready to be fitted by \code{fitModel}.  
}
\note{This work was supported by the National Cancer Institute (K25CA104791).}
\author{Tom Radivoyevitch (txr24@case.edu) }
\references{ Radivoyevitch, T. (2008) Equilibrium model selection: dTTP induced R1 dimerization. \emph{BMC Systems Biology} \bold{2}, 15. }
\seealso{\code{\link{mkg}}, \code{\link{mkKd2Kj}} and \code{\link{fitModel}}. }
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
dRt <- subset(RNR,(year==2001)&(fg==1)&(G==0)&(t>0),select=c(R,t,m,year))
names(dRt)[1:2] <- c("RT","tT") 

## first a simple spur model
mkModel(g,"IIIJ",dRt,Kjparams=c(R2t0=Inf, R1t1=Inf, R2t1=Inf, R2t2=1))
Kmapping=mkKd2Kj(g)
mdl=mkModel(g,"HDFF",dRt,Kdparams=c(R2t0=1, R0t0_t=1, R2t0_t=1, R2t1_t=1), 
                  Keq=c(R2t1_t="R2t0_t"),Kd2KjLst=Kmapping) 
fmdl <- fitModel(mdl)

## or mkGrids can be used to achieve the same thing as follows
gridL <- mkGrids(g,maxTotalPs=3)
chunk <- gridL$chunk
Keqs <- gridL$Keqs
mdl <- mkModel(g,"HDFF",dRt,Kdparams=chunk["HDFF",2:(g$nZ+1)], Keq=Keqs[["HDFF"]], Kd2KjLst=Kmapping,
       pparams=chunk["HDFF","p",drop=FALSE],indx=chunk["HDFF","indx"])
print(mdl)
print(chunk)
print(Keqs)

}
\keyword{ models }% __ONLY ONE__ keyword per line
