\name{aaStruct}
\alias{aaStruct}
\title{
Positional Amino Acid Composition Calculations
}
\description{
Calculates the probability of each amino acid residue at each position within
a peptide. A sub-function of \code{vMotif} and \code{vComp}.
}
\usage{
aaStruct(x, y, sigWeight = TRUE)
}
\arguments{
  \item{x}{
A \code{data.frame}, containg the peptides to be calculated
}
  \item{y}{
Object containing the signal set of interest for the defined peptides
in arguament \code{x}
}
  \item{sigWeight}{
\code{Logical} which determines if signal is incorporated into weight calculations
}
}
\details{
A sub - function of \code{vMotif} and \code{vComp}
}
\author{
Cody Moore}
\seealso{
\code{vComp, vMotif}
}
\keyword{ ~composition }
\keyword{ ~subfunction }
