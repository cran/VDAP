\name{QCKd}
\alias{QCKd}
\title{
Quality Control of Peptides Based on Reproducibility and Kd
}
\description{
Filter out peptides based on reproducibility between replicate
concentrations and relative dissociation constants (Kd). Peptides
must have a signal ratio between 0.5 and 2.0. A second reference file
may be loaded with the same peptides referenced against another sample.
Peptides are then compared based upon relative Kd value which must be at
least one log10 apart.
}
\usage{
QCKd(File1, File2 = NULL, Kd = FALSE, QC = TRUE, ColSet1 = NULL,
ColSet2 = NULL, ColSet3 = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{File1}{
An R object, usually a \code{data.frame} generally created by the function \code{FLoad()}
}
  \item{File2}{
An R object, usually a \code{data.frame} generally created by the function \code{FLoad()}
}
  \item{Kd}{
A \code{logical} value, if \code{Kd = TRUE} then peptides will be filtered by Kd against
the argument \code{File2}
}
  \item{QC}{
A \code{logical} value, if \code{QC = TRUE} then peptides will be filtered by ratios of
signal between replicates. Ratios must be between 0.5 to 2.0 to remain in the dataset.
}
\item{ColSet1}{
A sequence value, represents the two columns that are replicates at a single concentration.
Peptides must fit QC criteria in all given \code{ColSets} to remain in the dataset.
ColSets may be omitted if less than three concentrations are to be compared. Ex: 2:3
}
  \item{ColSet2}{
A sequence value, represents the two columns that are replicates at a single concentration.
Peptides must fit QC criteria in all given \code{ColSets} to remain in the dataset.
ColSets may be omitted if less than three concentrations are to be compared. Ex: 6:7}
  \item{ColSet3}{
A sequence value, represents the two columns that are replicates at a single concentration.
Peptides must fit QC criteria in all given \code{ColSets} to remain in the dataset.
ColSets may be omitted if less than three concentrations are to be compared. Ex: 4:5}
}
\details{
Either the \code{QC} or \code{Kd} filter may be applied by itself of both simultaneously.
}
\value{
A \code{data.frame} will be returned with peptides filtered out that do not meet the
given criteria for either the \code{QC} or \code{Kd} filters.
}
\author{
Cody Moore}
\examples{

protEx.QCKd <- data.frame(Peptides = c("PWRGPWARVGSG","GYNRVGQGSG","PNGYRSGVKGSG","GSG"),
Length = c(12,10,12,3),Charge = c(2,1,2,0),Kd = c(0.2572361,2.8239730,3.3911868,281.3058),
C_6uM = c(65011.48,47462.24,24778,2613.03),C_6uM2 = c(62637.81,20723.85,21313.67,2300.216))

## All peptides filtered out due to same Kd value between files ##

QCKdEx <- QCKd(protEx.QCKd, protEx.QCKd,Kd = TRUE, QC = TRUE, ColSet1 = 5:6)

## QC control only ##

QCKdEx <- QCKd(protEx.QCKd, QC = TRUE, ColSet1 = 5:6)
}
\keyword{ ~Kd }
\keyword{ ~Replicates }
