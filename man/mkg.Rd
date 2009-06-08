\name{mkg}
\alias{mkg}
\title{Make Generic Model}
\description{
  This function maps network topology information into a generic full spur graph model.   
  If \code{TCC} is \code{TRUE} (default) it also automatically generates 
  and compiles total concentration constraint C code and stores it in the \code{models} subdirectory.
}
\usage{ mkg(strct,TCC=TRUE,activity=FALSE,free=FALSE) }
\arguments{
  \item{strct}{ The thread topology of the equilibrium network 
  (see \code{\link{ccems}}). The slots of this list structure are:
      \code{heads} which is a character vector of the head nodes of threads 
      of the first site in \code{sites}; and
      \code{sites} which is a list of sites where each site is a list of 
      character vectors of non-head thread nodes. 
                           Naming conventions are given below under \code{Z}. 
                           }
  \item{TCC}{This is \code{TRUE} if total concentration constraints (TCCs) are to be used. If so, model C code is automatically generated
  and compiled. If \code{TCC} is \code{FALSE} a rational polynomial model is automatically generated as an R function component of the output.  
  Use of rational polynomial models assumes that the approximations [Xfree] = [XT] are reasonable for all non-hub species X (i.e. ligands).}
  \item{activity}{If \code{TRUE} data fitted is enzyme activity data. In this case the substrate variable should be named S. }
  \item{free}{Set \code{TRUE} if non-hub reactants are to be treated as free concentrations.  }
}


\details{This function is typically the first \code{ccems} function called in a script. 
  It creates a list of objects generic to the entire model space, i.e. a model space kernel.
  }
\value{
 A list comprised of the inputs and the following: 
  \item{id }{The biochemical equilibrium system ID. This is the set of single character reactant representations collapsed into one string. }
  \item{hubChar }{The first character of the first name in \code{strct}. As of ver 1.02 this is no longer a passed option. }
  \item{Z }{The names of the hub protein complexes expressed as single character reactant symbols followed by the number of copies of it in the complex.
  For example, R2t1 is a dimer of R with one bound ligand t, R2t2 is saturated dimer, and R2t0 is ligand free dimer. 
  The choice of the symbol Z for hub complexes derives from its common use as a complex number. }
  \item{nZ}{The length of \code{Z}. An \code{n} in front of a name often implies a length. }
  \item{atomS }{A character vector of atom/reactant names where the term atom refers to reactant molecules being indivisible 
              in the system of interest. A capital \code{S} at the end of a name implies a string.}
  \item{nAtomS }{The length of \code{atomS}.}
  \item{specieS }{A character vector equal to \code{c(atomS,Z)}, i.e. the names of all of the chemical species. }
  \item{nSpecieS }{The length of \code{specieS}.}
  \item{reactantS }{A list with component names equal to complex names and values equal 
  to vectors of the \code{atomS} that comprise them.}
  \item{W }{A dataframe of copy numbers/weights of atoms (column names) in each species (row names). }
  \item{KdS }{Subscript names of the dissociation constants, one for each element of \code{Z}. 
             Binary K sit in positions of their products
            and use \code{"_"} to separate their reactants; spur edges (e.g. head nodes)
            take their names from \code{Z}. 
 }
  \item{hdS }{The names of the head nodes.}
  \item{hds }{The positions of the head nodes within the Z vector of complex nodes.}
%  \item{wDir }{The working directory from which \code{models} and \code{results} are subdirectories.}
  \item{sstime }{The amount of ODE integration time used to reach steady state when solving TCCs (preset to 1e6). Since
                     integration uses variable step sizes which rapidly become large near steady state, overkill is OK.}
  \item{rtol }{An \code{lsoda} integration relative error tolerance parameter: \code{rtol}=1e-5.  }
  \item{atol }{An \code{lsoda} integration absolute error tolerance parameter: \code{atol}=1e-7.  }
  \item{parmsTCC}{The parameters of the total concentration constraints passed to the .dll (or .so) 
       when using \code{lsoda}. These include the total concentrations (i.e. system inputs). }
  \item{initialStateTCC}{The initial conditions (default zero) used to solve the TCC ODEs.} 
  \item{dfThread}{A dataframe of the thread within site and oligo structure.}
  \item{threads}{A nested list of the thread structure with threads at the top of the list and their contents and memberships below.} 
  \item{threadsWithinSites}{A list  of threads within sites, i.e. the list indices are sites.} 
  \item{nodesWithinSites}{A list of non-head nodes within sites.} 
  \item{usedLets}{A vector of the single characters used to label Kd equivalent threads. Currently, if there is 
  more than one thread, entire threads are either completely equal or all free to be different. If the curtain
  contains a single thread, then that thread can have patches of contiguous blocks that are equal. 
  Thus, a curtain with n threads has n used letters, but a curtain with one thread has as many used 
  letters as there are nodes in that one thread. }   
  \item{fback}{A function that maps a vector of free concentrations into 
  a vector of complex concentrations given the current set of complete 
      dissociation constant estimates of the optimization algorithm.}
  \item{code }{If \code{TCC = TRUE} this is the C code of the TCC ODE right hand side.  
                   This code is automatically compiled for use by \code{lsoda}.}
  The \code{TCC} flag, which passes unchanged from input to output, is used by \code{simulateData} to determine how the expected
  response surface is to be generated.  
}
\references{ Radivoyevitch, T. (2008) Equilibrium model selection: dTTP induced R1 dimerization. 
                      \emph{BMC Systems Biology} \bold{2}, 15. }
\author{Tom Radivoyevitch}
\note{
             Kd and Kj are generic grid and spur edge names, respectively.  The name \code{KdS} above is thus a grid-like name. That it
             is used to describe hybrid models that include head node spur edges is consistent with their 
             allowance in a generalized definition of grid graphs.
             
             In the description of \code{KdS} above "reactants" is used
             generically to mean a complex or a ligand, but by far, it refers most often to 
             purified substances of known controlled amounts (i.e. the experimentally 
             manipulated variables or the reactants that are initially present). This common usage is 
             synonymous with the use of "atoms" above, which emphasizes their indivisible building block nature
             in the non-covalent binding equilibriums of interest. 
                     
             This work was supported by the National Cancer Institute (K25CA104791). 
 }
\seealso{\code{\link{ccems}},\code{\link{ems}} }
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
mkg(topology) 
}
\keyword{ models }
