\name{mkKd2Kj}
\alias{mkKd2Kj}
\title{ Make Kd2Kj Mappings }
\description{
  This function takes the \code{g}  output  of \code{mkg} and converts it into a list of functions 
  that maps Kd (grid edge) values into Kj (spur edge) values.  
  Such functions are needed so that one generic full spur graph
  model can be used by all of the specific hypotheses of the model space.
}
\usage{ mkKd2Kj(g) }
\arguments{
  \item{g}{ The output \code{g} of \code{mkg}. }
}
\details{Suppose a hub protein ligand has two binding sites (call them a- and h-sites)
		that are filled in (a,h) order. 
		The names of the Kd2Kj functions produced by calls to this function are then strings of 0's and 1's where
       a '1' implies that the a-site thread is infinite while its corresponding (within the same j-mer)
        h-site thread is finite  
       and a '0' implies that the a-site thread is finite.
The relevance of this is that a '1' implies that a  
       finite tail Kj parameter must be inserted into an infinite a-site 
       thread to serve as the head of a downstream finite h-site thread; a '0'
       implies that no action need be taken. Thus, for monomers, dimers, tetramers and hexamers,    
'0000' is used when no bridge spur edge insertions are needed and '1111' is used 
       when all h-site threads are finite and all
       a-site threads are infinite. 
       
       More sites can also be handled. Suppose we have three ordered (s, a, h) across 
       j-mers (m, d, t, h).  The digits then follow binary codings of the s- and a-site threads being infinity: 
       '0' is niether (and thus no insertions), 
       '1' is an finite a-site thread, '2' is an infinite s-site thread, and '3' (binary 11) 
       is that both are infinite (this case is similar to  
       '1' in that only one a-site insertion is needed). If there were 4 binding 
       sites the digits would range from 
       '0' to '7' (binary 111).    
          }
          
\value{A list of functions that map vectors of Kd's into vectors of Kj's where, within 
stoichiometric integer factors 
in the numerators and denominators, an output Kj value is the product of 
the input Kd's that sit between the output Kj node and the root node 
(free hub). Here Kd's whose products are not heads
are binary and those which are heads are already spur edges, 
so they pass through an identity mapping.  Thus, the Kd inputs are of
hybrid/generalized grid graphs (i.e. hanging thread network 
topologies) and the Kj outputs are of 
pure spur full model graphs (i.e. the generic model output 
of \code{mkg} that is used by the entire model space).    
        }

\note{ Two functions with two different names can be identical. 
This happens when the monomer a-site thread is the same whether it is finite as a 
finite thread or finite because it was asked to
               support the subsequent monomer (e.g. h-site) thread;
               in longer threads these two scenarios yield different models: one where
               the thread nodes are all finite and one where only the tail node is finite while the rest are infinite.  

This work was supported by the National Cancer Institute (K25CA104791).
}

\author{Tom Radivoyevitch (txr24@case.edu) }
\references{ 
Radivoyevitch, T. (2008) Equilibrium model selection: dTTP induced R1 dimerization. \emph{BMC Systems Biology} \bold{2}, 15. 

}
\seealso{\code{\link{mkGrids}}, \code{\link{ccems}} }
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
g <- mkg(topology,TCC=TRUE) # generic model
mkKd2Kj(g) # when there is only one function its name is "default"

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

g <- mkg(topology) 
mkKd2Kj(g)  

\dontrun{
topology <- list(
    heads=c("R1X0","R2X0","R4X0","R6X0"), 
    # no requirement that s-sites are filled in oligomers
    sites=list(                    
        s=list(           # s-site                                           thread
            m=c("R1X1"),                                           # monomer    1
            d=c("R2X1","R2X2"),                                    # dimer      2
            t=c("R4X1","R4X2","R4X3","R4X4"),                      # tetramer   3
            h=c("R6X1","R6X2","R6X3","R6X4", "R6X5", "R6X6")       # hexamer    4
        ),
        a=list(          # a-site       
            m=c("R1X2"),                                           # monomer    5
            d=c("R2X3","R2X4"),                                    # dimer      6
            t=c("R4X5","R4X6","R4X7","R4X8"),                      # tetramer   7
            h=c("R6X7","R6X8","R6X9","R6X10", "R6X11", "R6X12")    # hexamer    8
        ),  # tails of a-site threads are heads of h-site threads
        h=list( #       # h-site
            m=c("R1X3"),                                            # monomer   9
            d=c("R2X5", "R2X6"),                                    # dimer     10
            t=c("R4X9", "R4X10","R4X11", "R4X12"),                  # tetramer  11
            h=c("R6X13", "R6X14", "R6X15","R6X16", "R6X17", "R6X18")# hexamer   12
        )
    )
)
g <- mkg(topology,TCC=TRUE) 
#gridL <- mkGrids(g,maxTotalPs=3)   
mkKd2Kj(g)     # this does not take much time to run but it will fill 
#          the buffer with many function definitions.
}


}
\keyword{ models}
