\name{Attrib}
\alias{Attrib}
\title{
Calculate Peptide Length and Charge Attributes
}
\description{
Calculates the length and charge of peptides in the first column of a
given dataset. A sub-function of \code{vFormat}
}
\usage{
Attrib(x)
}
\arguments{
  \item{x}{
An R object, generally a \code{data.frame}, containing peptides in the
first column
}
}
\value{
Returns a \code{data.frame} of 3 columns, starting with Peptide, the peptide's length,
followed by charge.
}
\author{Cody Moore
}
\note{
Uses the R Package:
\href{https://cran.r-project.org/web/packages/stringr/stringr.pdf}{stringr}
created by Hadley Wickham
}

\examples{

protEx <- data.frame(Peptides = c("PWRGPWARVGSG","GYNRVGQGSG","PNGYRSGVKGSG"),
C_6uM = c(65011.48,47462.24,24778), C_3uM = c(62637.81,31899.85,21313.67),
C_1.5uM = c(57893.22,25911.35,10397.99))

attribEx <- Attrib(protEx)

}
\keyword{ ~attributes }
\keyword{ ~formatting }
