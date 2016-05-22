\name{resSep}
\alias{resSep}
\title{
Select Peptides with the Specified Amino Acid Residue(s) at an Indicated Position
}
\description{
Allows the experimenter to subset peptide data based on a selected amino acid residue
or sequence a specified position(s). Requires the experimenter to select the residue(s)
and position(s) of interest at a given length or length/charge combination.
}
\usage{
resSep(File,Length,Charge = NULL,Pos,Res)
}
\arguments{
  \item{File}{
An object, generally a \code{data.frame}, the \code{vFormat} object with peptide and signal data.
}
\item{Length}{
An \code{integer}, the desired length of the peptides to separate.
}
\item{Charge}{
An \code{integer}, the desired charge of the peptides to separate.
Defaults to \code{Charge = NULL}, which carries out length separation only. 
}
\item{Pos}{
An \code{integer} or \code{sequence}, the position(s) to check for the 
residue(s) of interest.
}
\item{Res}{
A \code{character} input. The residue(s) to check for at the given position(s). The lengths of
the arguments \code{Pos} and \code{Res} must match. Multiple residues are entered as a single
character string. Ex: \code{Res = "RA"}.
}
}
\details{
The lengths of the arguments \code{Pos} and \code{Res} must match.

Sequence Positions are read from right to left. 

Ex: The residue "R" in 5-mer sequence "RSGSG" is at position 5. 

When typing in a sequence of interest, it will be in reverse with regard to the displayed sequence.

Ex: Sequence "SR" at positions 4:5 in the 5-mer"RSGSG"
}
\value{
A \code{data.frame} of the same format as the argument \code{File} containing only peptides
that contain the specified residue(s) at the indicated position(s).
}
\author{
Cody Moore
}
\seealso{
\code{vSep}
}
\examples{
## Example data.frame ##

protEx.resSep <- data.frame(Peptides = c("PWRGPWARVGSG","GYNRVGQGSG","PNGYRSGVKGSG","GSG"),
Length = c(12,10,12,3),Charge = c(2,1,2,0),Kd = c(0.2572361,2.8239730,3.3911868,281.3058),
C_6uM = c(65011.48,47462.24,24778,2613.03),C_6uM2 = c(62637.81,20723.85,21313.67,2300.216))

## Single Residue Separation ##

resSepEx1 <- resSep(protEx.resSep,12,2,5,"R")

## Positional Sequence Separation ##

resSepEx2 <- resSep(protEx.resSep,12,2,5:6,c("RA"))
}
\keyword{ ~Separation }
\keyword{ ~Residue }
