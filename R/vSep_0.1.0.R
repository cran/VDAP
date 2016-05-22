##Master's Separator 0.0.2

vSep <- function (File,Length = NULL,Charge = NULL){  # L/C Separation (File, Length, Charge)

  if(is.null(Charge) == TRUE & is.null(Length) == FALSE){
    L <- File[File[,2] == Length,]
    print("Length Separation \n")
    return(L)
  }

  if (is.null(Length) == TRUE & is.null(Charge) == FALSE){
    C <- File[File[,3] == Charge,]
    print("Charge Separation \n")
    return(C)
  }

  LC <- File[File[,2] == Length & File[,3] == Charge,]
  print("Length/Charge Separation \n")
  return(LC)
}
