\name{TK1}
\alias{TK1}
\docType{data}
\title{ Thymidine Kinase 1 Data}
\description{
  Human thymidine kinase 1 (i.e. cytosolic) data.   
}
\usage{data(TK1)}
\format{
  A data frame with the following columns.
  \describe{
    \item{\code{E}}{The total concentration of TK1 enzyme. }
    \item{\code{S}}{The total concentration of dT (the nucleoside substrate).}
    \item{\code{X}}{The total concentration of ATP.}
    \item{\code{v}}{The velocity of the kinase reaction in 1/seconds.}
    \item{\code{fg}}{The figure number in the original paper.}
    \item{\code{year}}{The year of the original paper.}
    \item{\code{jrnl}}{The journal of the data source: 
    \code{PEP} is Protein Expression and Purification, 
    \code{EJB} is  European Journal of Biochemistry, 
    \code{BBRC} is Biochem Biophys Res Commun, 
    \code{JBC} is Journal of Biological Chemistry.
    }
    \item{\code{vol}}{The volume number of the article.}
    \item{\code{page}}{The page number of the article.}
    \item{\code{frstAut}}{The first author of the article. This is a factor.  }
    \item{\code{index}}{The articles indexed as 1 through 5.  }
  }
}

\details{All concentrations are in micromolar. 
  Here R stands for root (i.e. hub protein). 
}
\source{
  The figure number, year, first author, journal, volume and page number are all included in the dataframe.   
}
\examples{
library(ccems)
TK1[1:20,]
d=subset(TK1,year==2006,select=c(E,S,v))
plot(d$S,d$v,xlab="Total [dT]", ylab="reaction velocity", 
         main="Birringer et al. PEP 2006 Fig. 4a")
}
\keyword{datasets}
