\name{QCon}
\alias{QCon}
\title{
Subsetting for VDAP function QCKd
}
\description{
A sub - function of QCKd, subsets data for replicate control
}
\usage{
QCon(File1,ColSet)
}
\arguments{
  \item{File1}{
Input File.
}
\item{ColSet}{
ColSet (Same as QCKd)
}
}
\author{
Cody Moore
}
\seealso{
\code{QCKd}
}
\examples{
## The function is currently defined as
function(File1,ColSet){

    Sig <- File1[,min(ColSet)]                  ##Column Calls
    Sig2 <- File1[,max(ColSet)]
    FVari1 <- File1[Sig/Sig2 > 0.5 & Sig/Sig2 < 2.0,]
    FVari1 <- na.omit(FVari1)
    return(FVari1)
    }
}
\keyword{~subset}
\keyword{~sub-function}
