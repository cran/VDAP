\name{hitSel}
\alias{hitSel}
\title{
Signal Based Hits Selection for VDAP}
\description{
Filters the dataset based upon signal from the specified columns.
Can be normalized to the average signal of any given peptide at the
given concentration. Works for multiple RFU signal inputs or a single Kd
input.
}
\usage{
hitSel(File, AvgSet, CutOff, Kd = FALSE)
}
\arguments{
  \item{File}{
An R object, usually a \code{data.frame} generally created by the function \code{FLoad()}
}
  \item{AvgSet}{
An integer sequence, defines the columns that contain the concentration data to be used
for hits selection. A given peptide will have to qualify as a hit at all given
concentration columns to be considered a true peptide hit. Ex: Hits based
upon 3 concentrations in columns 5 through 8 = 5:8. If \code{Kd = TRUE}, then a single
column with the calculated Kd values (generally column 4 created by vFormat) should be entered.
}
  \item{CutOff}{
A character string that defines the peptide to to normalize to. Hits must be 5 times
higher in signal than the given peptide to be returned as hits. Normally "GSG".If
\code{Kd = TRUE}, hits will be defined as peptides that have a calculated
Kd less than one half of the \code{Cutoff} peptide
}
\item{Kd}{ Toggle that determines if hits will be selected by RFU signal or Kd values.
If \code{Kd = TRUE}, hits will be defined as peptides that have a calculated
Kd less than one half of the \code{Cutoff} peptide
}
}
\value{
A \code{data.frame} will be returned only with the peptides that are hits in the given
context. (Hits must have Avg signal 5 times greater than the average signal of the peptide
specified in the argument \code{Cutoff}. Or one fifth (0.2) the \code{Cutoff} Kd value if
\code{Kd = TRUE})
}
\author{
Cody Moore
}
\examples{

protEx.hitSel <- data.frame(Peptides = c("PWRGPWARVGSG","GYNRVGQGSG","PNGYRSGVKGSG","GSG"),
Kd = c(0.2572361,2.8239730,3.3911868,281.3058),C_6uM = c(65011.48,47462.24,24778,2613.03),
C_3uM = c(62637.81,31899.85,21313.67,1161.216),C_1.5uM = c(57893.22,25911.35,10397.99,630.4025))

## Hits selection by RFU signal ##

hitSelRFU <- hitSel(protEx.hitSel,3:5,"GSG",Kd = FALSE)

## Hits selection by calculated Kd ##

hitSelKd <- hitSel(protEx.hitSel,2,"GSG",Kd = TRUE)
}

\keyword{ ~hits }
\keyword{ ~Kd }
