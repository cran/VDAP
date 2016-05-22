  ###############################Library Dependencies###########################################
library(stringr)
library(drc)
library(ggplot2)

################################Averaging Duplicates (Dups)##################################

Dups <- function(x) # x = File
{
  DuP <- x[duplicated(x[,1]),]

  if (nrow(DuP) == 0){
     Pavg <- x  # Added handling for no replicates##0.0.2##
     return(Pavg)
  }
  else{

  DuPep <- DuP[!duplicated(DuP[,1]),]  # Trimming of Duplicated List, leaves 1 iteration of each

  PepRep <- data.frame(DuPep[,1])  # Duplicated Peptides
  PepGlob <- data.frame(x[,1])  # Global Peps
  PepGSig <- x[,2:ncol(x)]  # Global Signal

  ScanAvg <- function (x)  # Duplicated set against Global set, averages the signal between replicates
  {
    a <- PepGlob == x
    b <- PepGSig[a,]
    c <- apply(b,2,mean)
    return(c)
  }

  Final <- apply(PepRep,1,ScanAvg)  # Application of function
  Finall <- t(Final)  # transpose from rows to columns

  Yupp <- data.frame(DuPep[,1],Finall)  # Attachment of peptide and length/Charge Columns
  names(Yupp) <- names(x)  # Matching of col names (DuPep$Peptide->Peptide)

  FProt <- rbind(Yupp,x)  # Binding to global set

  Pavg <- FProt[!duplicated(FProt[,1]),]  # Final List with averaged sig of duplicated peptides

  return(Pavg)
  }
}

##########################Attributes Columns (Attrib)######################################

Attrib <- function(x) # x = File
{
  Peps <- data.frame(as.character(x[,1]))
  colnames(Peps) <- "Peptide"

  Length <- data.frame(apply(Peps,1,nchar)) # Length Column
  colnames(Length) <- "Length"

  Peps <- cbind(Peps,Length)  #Binding

  Peps$Charge <- (str_count(Peps[,1],"R")+str_count(Peps[,1],"K")-str_count(Peps[,1],"D")-str_count(Peps[,1],"E"))
  return(Peps)
}

############################################Kd Analysis (KdA)##################################

KdA <- function(x,y,z) # x = File  y = c([]'s), z = Corresponding Columns
{
  Rest <- function()
  {
    Test1 <- data.frame(matrix(c(0,0,0,0),nrow=1))  # Function for tryCatch error
    Test[i,] <- Test1  # row of NA's and binding
  }

  Test <- data.frame(matrix(0,ncol=4,nrow=length(x[,1])))
  Concs <- data.frame(y)      ##x-values [concentrations]

  for(i in 1:length(x[,1]))
  {
    Sigs <- data.frame(Sigs=t(x[i,z]),row.names=NULL)  # y-values RFU signals by row
    Struct <- cbind(Sigs,Concs)

    tryCatch({
      model1 <- drm(Struct,fct=MM.2(fixed = c(NA, NA), names = c("max", "Kd"))) #Works #drc package MM regression
      Test1 <- data.frame(matrix(coef(summary(model1))["Kd:(Intercept)",],nrow=1)) #Extracting Kd and other stats from summary statistics
      Test[i,] <- Test1},error = function(f){Rest()})
  }
  names(Test) <- c("Kd","Std. Dev","t-value","p-value")
  return(Test)
}

#########################Assembled Function (Draformat)#########################

vFormat<-function(x,Kd = FALSE,Concs,Cols) # Length/Charge/Kd and File format
{

  a <- Dups(x)
  b <- Attrib(a)

  if (Kd == TRUE){

    c <- KdA(a,Concs,Cols)
    DaBusiness <- data.frame(b,c[,1],a[,2:ncol(a)],c[,2:4],row.names=NULL)
    colnames(DaBusiness)[4] <- "Kd"
    return(DaBusiness)
  }

  DaBusiness <- data.frame(b,a[,2:ncol(a)],row.names = NULL)
  return(DaBusiness)
}

################################################################################



