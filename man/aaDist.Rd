\name{aaDist}
\alias{aaDist}
\title{
Position Independent Amino Acid Distributions
}
\description{
Generates Position Independent Amino Acid Ditributions within
VDAP data sets
}
\usage{
aaDist(x, plotName = NULL, linker = TRUE)
}
\arguments{
  \item{x}{
An R object, usually a \code{data.frame} generally created by the function \code{FLoad()}
}
  \item{plotName}{
A plot title may be entered here surrounded by "quotations" or a \code{class(character)} object
}
\item{linker}{
\code{Logical} determining if a 3 residue linker "GSG" is present or not. If \code{linker = TRUE},
the "GSG" linker portion of each peptide will be excluded from distribution calculations. Default is
FALSE.
}
}
\details{
Uses both \pkg{stringr} and \pkg{ggplot2} for peptide calculations and plotting
}
\value{
\code{aaDist} will return a \code{data.frame} that contains a table with the amino acid
distribution over the entire array object. A \code{ggplot2} object will also be displayed
with the same information as the histogram.
}
\author{
Cody Moore
}
\examples{
protEx <- data.frame(Peptides = c("PWRGPWARVGSG","GYNRVGQGSG","PWRGPWARVGSG","GYNRVGQGSG","GSG"))

## Plot example with GSG linker ##

aaDistEx <- aaDist(protEx,"aaDistEx Plot",linker = TRUE)
}
\keyword{ ~distributions }
\keyword{ ~composition }
