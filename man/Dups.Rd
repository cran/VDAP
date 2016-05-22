\name{Dups}
\alias{Dups}
\title{
Average duplicated peptides from a VDAP dataset
}
\description{
Looks for duplicate peptides in the first column of the dataset,
averages the signal of duplicates and replaces them with a single row.
A subfunction of \code{vFormat}
}
\usage{
Dups(x)
}
\arguments{
  \item{x}{
An R object, generally a \code{data.frame} with peptides in column 1,
followed by signal values at various concentrations.
}
}
\value{
Returns a \code{data.frame} without duplicated peptides. Duplicate
entries display the mean of the signal at each concentration
}
\author{
Cody Moore
}
\note{
Duplicated peptide entries will generally be at the top of the dataset
}
\examples{

protExDups <- data.frame(Peptides = c("PWRGPWARVGSG","GYNRVGQGSG","PWRGPWARVGSG"),
C_6uM = c(65011.48,47462.24,24778), C_3uM = c(62637.81,31899.85,21313.67),
C_1.5uM = c(57893.22,25911.35,10397.99))

exDups <- Dups(protExDups)
}

\keyword{ ~duplicate }
\keyword{ ~formatting }
