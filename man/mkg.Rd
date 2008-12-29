\name{mkg}
\alias{mkg}
\title{Make Generic Model}
\description{
  This function mapsss network topology information into a generic full spur graph model.   
  If \code{TCC} is \code{TRUE} it also automatically generates and compiles total concentration constraint C code and stores it in the \code{models} subdirectory.
}
\usage{ mkg(strct,hubChar="R",monomerMass=90,TCC=FALSE) }
\arguments{
  \item{strct}{ The thread topology of the equilibrium network 
  (see \code{\link{ccems}}). The slots of this list structure are:
      \code{heads} which is a character vector of the head nodes of threads 
      of the first site in \code{sites}; and
      \code{sites} which is a list of sites where each site is a list of 
      character vectors of non-head thread nodes. 
                           Naming conventions are given below under \code{Z}. 
                           }
  \item{hubChar}{ The hub protein's single character representation. The default is R for the big (R1) subunit of ribonucleotide reductase (RNR).  }
  \item{monomerMass}{The hub protein's monomer mass in kDa. The default is 90 for the R1 subunit of RNR. }
  \item{TCC}{This is \code{TRUE} if total concentration constraints (TCCs) are to be used. If so, model C code is automatically generated
  and compiled. If \code{TCC} is \code{FALSE} a rational polynomial model is automatically generated as an R function component of the output.  
  Use of rational polynomial models assumes that the approximations [Xfree] = [XT] are reasonable for all non-hub species X (i.e. ligands).}
}
\details{This function is typically the first \code{ccems} function called in a script. 
  It creates a list of objects generic to the entire model space, i.e. a model space kernel.
}
\value{
 A list comprised of the inputs and the following: 
  \item{id }{The biochemical equilibrium system ID. This is the set of single character reactant representations collapsed into one string. }
  \item{Z }{The names of the hub protein complexes expressed as single character reactant symbols followed by the number of copies of it in the complex.
  For example, R2t1 is a dimer of R with one bound ligand t, R2t2 is saturated dimer, and R2t0 is ligand free dimer. 
  The choice of the symbol Z for hub complexes derives from its common use as a complex number. }
  \item{nZ}{The length of \code{Z}. An \code{n} in front of a name often implies a length. }
  \item{atomS }{A character vector of atom/reactant names where the term atom refers to reactant molecules being indivisible 
              in the system of interest. A capital \code{S} at the end of a name implies a string.}
  \item{nAtomS }{The length of \code{atomS}.}
  \item{specieS }{A character vector equal to \code{c(atomS,Z)}, i.e. the names of all of the chemical species. }
  \item{nSpecieS }{The length of specieS.}
  \item{reactantS }{A list with component names equal to complex names and values equal to vectors of the atomS that comprise them.}
  \item{W }{A dataframe of copy numbers/weights of atoms (column names) in each species (row names). }
  \item{KdS }{Subscript names of the dissociation constants, one for each element of \code{Z}. Binary K sit in positions of their products
            and use \code{"_"} to separate their reactants; head spur edges take their names from \code{Z}. 
 }
  \item{hdS }{The names of the head nodes.}
  \item{hds }{The positions of the head nodes within the Z vector of complexes nodes.}
  \item{wDir }{The working directory from which \code{models} and \code{results} are subdirectories.}
  \item{sstime }{The amount of ODE integration time used to reach steady state when solving TCCs (preset to 1e6). Since
                     integration uses variable step sizes which rapidly become large near steady state, overkill is OK.}
  \item{rtol }{An \code{lsoda} integration relative error tolerance parameter: \code{rtol}=1e-5.  }
  \item{atol }{An \code{lsoda} integration absolute error tolerance parameter: \code{atol}=1e-7.  }
  \item{parmsTCC}{The parameters of the total concentration constraints passed to the .dll (or .so) file when using \code{lsoda}. 
                     These include the total concentrations, i.e. the system inputs. }
  \item{initialStateTCC}{The initial conditions (i.e. zero free concentrations) used to solve the TCC ODEs.} 
  \item{fback }{A function that maps a vector of free concentrations into 
  a vector of complex concentrations given the current dissociation constant estimates of the optimization algorithm.}
  \item{code }{C code of the TCC ODE right hand side.  This code is automatically compiled for use by \code{lsoda}.}
  The \code{TCC} flag, which passes unchanged from input to output, is used by \code{simulateData} to determine how the expected
  response surface is to be generated.  
}
\references{ Radivoyevitch, T. (2008) Equilibrium model selection: dTTP induced R1 dimerization. \emph{BMC Systems Biology} \bold{2}, 15. }
\author{Tom Radivoyevitch}
\note{
             Kd and Kj are generic grid and spur edge names, respectively.  The name \code{KdS} above is thus a grid-like name. That it
             is used to describe hybrid models that include head node spur edges is consistent with their 
             allowance in a generalized definition of grid graphs.
             
             In the description of \code{KdS} above "reactants" is used
             generically to mean a complex and a ligand. In most other places it refers to 
             purified substances of known controlled amounts (i.e. the experimentally 
             manipulated variables or the reactants that are initially present). This common usage is 
             synonymous with the use of "atoms" above, which emphasizes their indivisible building block nature
             in the non-covalent binding equilibrium systems of interest. 
                     
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
mkg(topology,TCC=TRUE) 
}
\keyword{ models }
