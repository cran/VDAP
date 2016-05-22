#########################################QCKd (Quality control and Kd)#######################################################

##################################################HitSel################################################

hitSel <- function(File,AvgSet,CutOff,Kd = FALSE)
{
  if(Kd == TRUE){
    Norm <- File[File[,1] == CutOff,]
    Final <- File[File[,AvgSet] < (0.2*Norm[,AvgSet]) & File[,AvgSet] != 0,]
    return(Final)
  }

  Norm <- File[File[,1] == CutOff,]
  Final <- File[File[,min(AvgSet)] > (5*mean(Norm[,min(AvgSet)])),]

  for (i in min(AvgSet):max(AvgSet))
  {
    Selection <- File[File[,i] > (5*mean(Norm[,i])),]
    Combined <- rbind(Selection,Final)
    Final <- Combined[duplicated(Combined),]
  }
  return(Final)
}


#####################################################QC###############################################

QCon <- function(File1,ColSet)
  {
    Sig <- File1[,min(ColSet)]                  ##Column Calls
    Sig2 <- File1[,max(ColSet)]
    FVari1 <- File1[Sig/Sig2 > 0.5 & Sig/Sig2 < 2.0,]
    FVari1 <- na.omit(FVari1)
    return(FVari1)
  }

##########################################################QCKd########################################

QCKd<-function(File1,File2 = NULL,Kd = FALSE,QC = TRUE,ColSet1 = NULL,ColSet2 = NULL,ColSet3 = NULL)
{
  if (Kd == TRUE)
  {
    Kds <- File1[abs(log10(File1[,4])-log10(File2[,4])) > 1,]
    File1 <- Kds
  }

  if (QC == TRUE){
    QCon <- function(File1,ColSet){

      Sig <- File1[,min(ColSet)]                  ##Column Calls
      Sig2 <- File1[,max(ColSet)]
      FVari1 <- File1[Sig/Sig2 > 0.5 & Sig/Sig2 < 2.0,]
      FVari1 <- na.omit(FVari1)
      return(FVari1)
    }

    if(is.null(ColSet2) == TRUE & is.null(ColSet3) == TRUE){
      CSet1 <- QCon(File1,ColSet1)
      return(CSet1)

    } else if(is.null(ColSet3) == TRUE){

        CSet1 <- QCon(File1,ColSet1)
        CSet2 <- QCon(File1,ColSet2)
        CSet <- rbind(CSet1,CSet2)
        CSet <- CSet[duplicated(CSet),]
        return(CSet)
      } else{

          CSet1 <- QCon(File1,ColSet1)
          CSet2 <- QCon(File1,ColSet2)
          CSet3 <- QCon(File1,ColSet2)
          CSet <- rbind(CSet1,CSet2)
          CSet <- CSet[duplicated(CSet),]
          CSet <- rbind(CSet,CSet3)
          CSet <- CSet[duplicated(CSet),]
          return(CSet)
        }
  } else{

      return(File1)
    }
}
