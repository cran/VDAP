\name{vSep}
\alias{vSep}
\title{
Select Peptides of a Particular Length/Charge Combination
}
\description{
Select Peptides that have a specified length/charge combination, a subfunction
for \code{lcScan}, and all methods of \code{LCMotif} and \code{LcComp}}
\usage{
vSep(File, Length = NULL, Charge = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{File}{
An R object, usually a \code{data.frame} generally created by the function \code{FLoad()}
}
  \item{Length}{
An \code{integer} value, specifies the desired length to subset.
}
  \item{Charge}{
An \code{integer} value, specified the desired charge to subset.
}
}
\value{
Returns a \code{data.frame} with peptides of the selected Length/Charge
combination.
}
\author{
Cody Moore
}
\examples{
protExChargeSep <- data.frame(Peptides = c("PWRGPWARVGSG","GYNRVGQGSG","PWRGPWARVGSG"),
Length = c(12,10,12), Charge = c(2,1,2))

## Length/Charge Combination ##

hitSelEx <- vSep(protExChargeSep,10,1)

## Charge only ##

hitSelEx <- vSep(protExChargeSep,Charge = 1)

## Length Only ##

hitSelEx <- vSep(protExChargeSep,Length = 12)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~Separation }
\keyword{ ~filtering }% __ONLY ONE__ keyword per line
