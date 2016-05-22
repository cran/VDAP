\name{vMotif}
\alias{vMotif.lc}
\alias{vMotif.l}
\alias{vMotif}
\title{
Generate Signal Weighted Amino Acid Heat Maps by Position
}
\description{
Generate signal weighted amino acid composition maps by postion at
specific length or length/charge combinations. Weights are compared to the global
distribution of peptides at the particular length or length/charge.

Weights are centered at 1, meaning that there is no change in probability
or signal from the global set. Weights above 1 indicate higher probability
at the given position and/or signal while weights below 1 indicate lower
probability at the given position and/or signal.

When \code{Kd = TRUE}, weighting by Kd instead of signal is performed. Weights are
generated using (1/Kd) since lower Kd values generally indicate higher affinity interactions,
and would correlate with higher signal.
}
\usage{
vMotif.lc(Prot, ProtG, Length, Charge, SigCol, Kd = FALSE)

vMotif.l(Prot, ProtG, Length, SigCol, Kd = FALSE)
}
\arguments{
  \item{Prot}{
An \code{R} object, generally a \code{data.frame}. Contains peptides that are considered "hits"
or selected peptides with their length,charge, and signal/Kd attributes.
}
  \item{ProtG}{
An \code{R} object, generally a \code{data.frame}. Contains the set of
peptides from which the argument \code{Prot} were selected with their
corresponding length, charge, and signal information.
}
  \item{Charge}{
An \code{integer} value, indicating the desired charge to analyze
}
  \item{Length}{
An \code{integer} value, indicating the desired peptide length to analyze
}
  \item{SigCol}{
An \code{Integer} value, indicating the column that contains the desired signal data
at a given concentration
}
\item{Kd}{
An \code{logical} value, indicating if weights should be generated using signal or Kd data. Effects
signal weighting behavior. If \code{Kd = TRUE}, weights are generated using 1/\code{SigCol}.
}
}
\value{
Returns a \code{data.frame} that shows weights for each amino acid at each position
within the peptide of the selected length. Also output a positional heatmap using the
package \pkg{ggplot2}
}
\author{
Cody Moore
}
\seealso{
\code{vComp}
}
\examples{
protEx.Motif <- data.frame(Peptides = c("PWRGPWARVGSG","GYNRVGQGSG","PNGYRSGVKGSG","GSG"),
Length = c(12,10,12,3),Charge = c(2,1,2,0),Kd = c(0.2572361,2.8239730,3.3911868,281.3058),
C_6uM = c(65011.48,47462.24,24778,2613.03),C_6uM2 = c(62637.81,20723.85,21313.67,2300.216))

## vMotif Length/Charge and Length Signal Examples ##

vMotif.lcEx <- vMotif.lc(protEx.Motif,protEx.Motif, 12,2,5,Kd = FALSE)

vMotif.lEx <- vMotif.l(protEx.Motif,protEx.Motif, Length = 12,SigCol = 5,Kd = FALSE)

## vMotif Length/Charge Kd Example ##

vMotif.lcEx <- vMotif.lc(protEx.Motif,protEx.Motif, Length = 12,Charge = 2, SigCol = 5,Kd = TRUE)

}
\keyword{~Motif}

