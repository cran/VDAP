\name{KdA}
\alias{KdA}
\title{
Peptide Dissociation Rate Constant (Kd) Calculations
}
\description{
Calculates the Kd of each peptide using a non-linear
single site specific binding model. A sub-function of \code{vFormat}
}
\usage{
KdA(x, y, z)
}
\arguments{
  \item{x}{
An R object, generally a \code{data.frame}, containing peptides in the
first column
}
  \item{y}{
The concentrations of each column used for Kd calculations, separated by commas.
The order must match the relative position of the columns.
}
  \item{z}{
The columns used for Kd calculations, expressed as a sequence.
Ex: Columns 2 through 4 = 2:4
}
}
\author{
Cody Moore
}
\note{
Uses the R package:
\href{https://cran.r-project.org/web/packages/drc/drc.pdf}{drc}
created by Christian Ritz and Jens C. Strebig
}
\examples{
protEx <- data.frame(Peptides = c("PWRGPWARVGSG","GYNRVGQGSG","PWRGPWARVGSG"),
C_6uM = c(65011.48,47462.24,24778), C_3uM = c(62637.81,31899.85,21313.67),
C_1.5uM = c(57893.22,25911.35,10397.99))

exKdA <- KdA(protEx,c(6,3,1.5),2:4)
}

\keyword{ ~kinetics }
\keyword{ ~formatting }
