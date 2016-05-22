\name{vComp}
\alias{vComp.lc}
\alias{vComp.l}
\alias{vComp}
\title{
Amino Acid Disbutions by Position at Various Length/Charge}
\description{
Generates the probability of each amino acid to appear in each position within a peptide
of a specific length or length/charge combination. Can either be the raw probability or the
ratio between the probabilities of 2 peptide sets.

Weights are centered at 1, meaning that there is no change in probability or signal from the
global set. Weights above 1 indicate higher probability at the given position while
weights below 1 indicate lower probability at the given position.
}
\usage{
vComp.lc(Prot, ProtG, Length, Charge)

vComp.l(Prot, ProtG, Length)

}
\arguments{
  \item{Prot}{
An \code{R} object, generally a \code{data.frame}. Contains peptides that are considered "hits"
or selected peptides with their length,charge, and signal information.
}
  \item{ProtG}{
An \code{R} object, generally a \code{data.frame}. Contains the set of
peptides from which the argument \code{Prot} were selected with their
corresponding length, charge, and signal information.
}
  \item{Length}{
An \code{integer} value, indicating the desired peptide length to analyze
}
  \item{Charge}{
An \code{integer} value, indicating the desired charge to analyze
}
}
\details{
If raw probabilities are desired, the same object can be loaded into both the
\code{Prot} and \code{ProtG} arguments.
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
\code{vMotif}, \code{genPep}
}
\examples{
protEx.Motif <- data.frame(Peptides = c("PWRGPWARVGSG","GYNRVGQGSG","PNGYRSGVKGSG","GSG"),
Length = c(12,10,12,3),Charge = c(2,1,2,0),Kd = c(0.2572361,2.8239730,3.3911868,281.3058),
C_6uM = c(65011.48,47462.24,24778,2613.03),C_6uM2 = c(62637.81,20723.85,21313.67,2300.216))

## Length/Charge Example ##

vComp.lcEx <- vComp.lc(protEx.Motif,protEx.Motif, 12,2)

## Length Example ##

vComp.lEx <- vComp.l(protEx.Motif,protEx.Motif, 12)
}

\keyword{ ~Composition }
\keyword{ ~Comp }
