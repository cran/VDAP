\name{vFormat}
\alias{vFormat}
\title{
Length/Charge/Kd Peptide Calculations and File Assembly
}
\description{
Calculates the length, charge, and dissociation rate constant (Kd)
for each peptide and assembles the file into a universal format for
subsequent VDAP Functions.
}
\usage{
vFormat(x,Kd = FALSE,Concs,Cols)
}
\arguments{
  \item{x}{
An R object, usually a \code{data.frame} generally created by the function \code{FLoad()}
}
  \item{Kd}{Toggle to specify if dissociation rate constants (Kd) values should be calculated. If
\code{Kd = FALSE}, the nonlinear regression package \code{drc} will not be used.
}
  \item{Concs}{
The concentrations of each column used for Kd calculations, separated by commas.
The order must match the relative position of the columns.
}
  \item{Cols}{
The columns used for Kd calculations, expressed as a sequence.
Ex: Columns 2 through 4 = 2:4
}
}
\details{
The order of concentrations should not matter, as long as they are identical between
the \code{Concs} and \code{Cols} arguments. However, the columns must all be adjacent.
}
\value{
A \code{data.frame} will be returned with the Length, charge, and Kd if \code{Kd = TRUE} characteristics
placed in columns 2 - 4, followed by the signal at each concentration from the x argument.
This is followed by quality values such as std.error, p-value, and t-value from the Kd
of each peptide. Peptides will remain in column 1.
}
\author{
Cody Moore
}
\note{
Uses the R Package:
\href{https://cran.r-project.org/web/packages/stringr/stringr.pdf}{stringr}
created by Hadley Wickham
and
\href{https://cran.r-project.org/web/packages/drc/drc.pdf}{drc}
created by Christian Ritz and Jens C. Strebig
}
\seealso{
\code{\link{Dups}}, \code{\link{Attrib}}, \code{\link{KdA}}.
}
\examples{

## vFormat on example data set ##

protEx <- data.frame(Peptides = c("PWRGPWARVGSG","GYNRVGQGSG","PNGYRSGVKGSG"),
C_6uM = c(65011.48,47462.24,24778), C_3uM = c(62637.81,31899.85,21313.67),
C_1.5uM = c(57893.22,25911.35,10397.99))

## Preformatted protEx ##

      #Peptides    C_6uM    C_3uM  C_1.5uM
#1 PWRGPWARVGSG 65011.48 62637.81 57893.22
#2   GYNRVGQGSG 47462.24 31899.85 25911.35
#3 PNGYRSGVKGSG 24778.00 21313.67 10397.99


formatEx <- vFormat(protEx,Kd = TRUE, c(6,3,1.5), 2:4)

## Formatted output ##

       #Peptide Length Charge        Kd    C_6uM    C_3uM  C_1.5uM    Std..Dev   t.value    p.value
#1 PWRGPWARVGSG     12      2 0.2572361 65011.48 62637.81 57893.22 0.008441968 30.471112 0.02088507
#2   GYNRVGQGSG     10      1 2.8239730 47462.24 31899.85 25911.35 1.619385359  1.743855 0.33146423
#3 PNGYRSGVKGSG     12      2 3.3911868 24778.00 21313.67 10397.99 2.522251940  1.344508 0.40711826

}

\keyword{ ~VDAP }
\keyword{ ~Formatting }
