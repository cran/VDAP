\name{lcScan}
\alias{lcScan}
\title{
Signal or Kd Distributions separated by Length/Charge attributes
}
\description{
Calculates the mean with standard error, and population peptides at each length/charge combination
within a VDAP dataset. If the argument \code{Glob = TRUE}, average signals will be compared against
a global set of peptides and p - values will be calculated for hypoethesis testing. \code{lcScan}
will also return a plot for visualization of signal, population, and hypothesis testing.
}
\usage{
lcScan(File,Glob = NULL, Conc = 5, Kd = FALSE)
}
\arguments{
  \item{File}{
An R object, usually a \code{data.frame} generally created by the function \code{FLoad()}
}
  \item{Glob}{
A second \code{data.frame} with the global set of peptides. If the original \code{File} argument contains
peptides hits, \code{Glob} should contain the dataset before hits were filtered out.
}
  \item{Conc}{
The column contianing the concentration or Kd data to be analyzed, an integer. Default is column 5
which is generally the highest concentration average according to the default formatting function
\code{vFormat}

Ex: Column 1 = 1
}
  \item{Kd}{
Toggle to calculate by a defined signal column or by calculated Kd values, effects final plot behavior
and labels. If \code{Kd = TRUE}, then the arguement \code{Conc} should be set to 4 if the file was
formatted by the default formatting function \code{vFormat}.
}
}
\value{
A \code{data.frame} will  be returned with columns for the mean, standard error,
and population of peptides at each length/charge combination that can be exported
for further analysis. Also uitilizes \pkg{ggplot2} and \pkg{reshape2} to create
a heat map plot that shows the signal distribution with corresponding populations that
can be exported.
}
\references{
Plot generation utilizes
\href{https://cran.r-project.org/package=ggplot2}{ggplot2}
created by Hadley Wickham [aut, cre] and Winston Chang [aut] and
\href{https://cran.r-project.org/package=reshape2}{reshape2}
created by Hadley Wickham
}
\author{
Cody Moore
}
\examples{

protEx.lcScan <- data.frame(Peptides = c("PWRGPWARVGSG","GYNRVGQGSG","PNGYRSGVKGSG","GSG"),
Length = c(12,10,12,3),Charge = c(2,1,2,0),Kd = c(0.2572361,2.8239730,3.3911868,281.3058),
C_6uM = c(65011.48,47462.24,24778,2613.03),C_6uM2 = c(62637.81,20723.85,21313.67,2300.216))

## Signal length/charge Analysis ##

lcScanEx <- lcScan(protEx.lcScan)

## Kd length/charge Analysis ##

lcScanEx <- lcScan(protEx.lcScan, Conc = 4, Kd = TRUE)
}
\keyword{ ~scan }
\keyword{ ~plot }
