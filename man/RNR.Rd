\name{RNR}
\alias{RNR}
\docType{data}
\title{ Ribonucleotide Reductase Data}
\description{
  This is a compilation of ribonucleotide reductase (RNR) literature data into one data frame. 
}
\usage{data(RNR)}
\format{
  A data frame with 546 observations on the following 23 variables.
  \describe{
    \item{\code{R}}{The total concentration of R1 (big R for the big subunit of RNR).}
    \item{\code{r}}{The total concentration of R2 dimer (small r for the small subunit of RNR).}
    \item{\code{U}}{The total concentration of UDP.}
    \item{\code{C}}{The total concentration of CDP.}
    \item{\code{G}}{The total concentration of GDP.}
    \item{\code{A}}{The total concentration of ADP.}
    \item{\code{t}}{The total concentration of dTTP.}
    \item{\code{tf}}{The concentration of free dTTP.}
    \item{\code{g}}{The total concentration of dGTP.}
    \item{\code{gf}}{The concentration of free dGTP.}
    \item{\code{a}}{The total concentration of dATP.}
    \item{\code{af}}{The concentration of free dATP.}
    \item{\code{X}}{The total concentration of ATP.}
    \item{\code{Mg}}{The total concentration of magnesium.}
    \item{\code{v}}{The velocity of the reduction reaction in 1/seconds; v in Scott et al should be multiplied by .00151 to create these units, 
                    save figure 7 (GDP reduction) which is in percentages of the maximum velocity.}
    \item{\code{m}}{The weight averaged average mass.}
    \item{\code{n}}{The number of ligands bound. Exception is Ingemarson et al. 1996 where it is a possibly negative BIAcore experiment parameter.}
    \item{\code{fg}}{The figure number in the original paper.}
    \item{\code{year}}{The year of the original paper.}
    \item{\code{jrnl}}{The journal of the data source: \code{Bioc} is Biochemistry, \code{JBC} is Journal of Biological Chemistry.}
    \item{\code{vol}}{The volume number of the article.}
    \item{\code{page}}{The page number of the article.}
    \item{\code{frstAut}}{The first author of the article. This is a 
    factor with levels \code{IngemarsonR}, \code{KashlanOB}, \code{Rofougaran}, \code{sanath}, and \code{ScottCP}.}
    \item{\code{organism}}{An element of \code{c("human","mouse","yeast","ecoli")}. 
                 This dataset is predominantly mouse RNR data. }
  }
}

\details{All concentrations are in micromolar. 
  To some extent, capital letters are used for ribonucleotides and small letters for deoxyribonucleotides. Free concentrations, where given,
  were model dependent and were thus converted back to total concentrations, i.e. to make them model independent. 
  Values of zero that lead a sequence of measurements of increasing concentrations were set to 1e-2 (essentially zero) to allow [x]>0 subsetting.
}
\source{
  The figure number, year, first author, journal, volume and page number are all included in the dataframe.  The main reference is
  Kashlan OB, Scott CP, Lear JD, Cooperman BS: A comprehensive model for the allosteric regulation of mammalian ribonucleotide 
  reductase. Functional consequences of ATP- and dATP-induced oligomerization of the large subunit. Biochemistry 2002, 41(2):462-474. 
  The first seven rows, by first author Sanath (Sanath Wijerathna of the Dealwis Lab at CWRU), are unpublished. 
}
\references{ Radivoyevitch, T. (2008) Equilibrium model selection: dTTP induced R1 dimerization. \emph{BMC Systems Biology} \bold{2}, 15. }
\examples{
library(ccems)
data(RNR)
RNR[1:40,]
d=subset(RNR,(year==2002)&(fg==1)&(X>0),select=c(R,X,m,year))
plot(d$X,d$m,xlab="Total [ATP]", ylab="Weight averaged R1 mass", main="Kashlan et al. Biochemistry, 2002, Fig. 1 (DLS data)")
}
\keyword{datasets}
