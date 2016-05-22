resSep <- function(File,Length,Charge = NULL,Pos,Res){
  
  Peps <- vSep(File,Length,Charge)  # vSep function 
    
  if (length(Pos) == 1){
    
    AAs <- as.character(Peps[,1])  # ensuring character type data
    AAs <- strsplit(AAs,"")  # Splitting of Sequences into Matrix
    AAs <- do.call(rbind,AAs)  # Converting from list
  
    Result <- Peps[AAs[,ncol(AAs) - (Pos - 1)] == Res,]  # Conditional Selection
    return(Result)
  
  } else if (length(Pos > 1)){
    
    AAs <- as.character(Peps[,1])  # ensuring character type data
    AAs <- strsplit(AAs,"")  # Splitting of Sequences into Matrix
    AAs <- do.call(rbind,AAs)  # Converting from list
    AAs <- AAs[,(ncol(AAs) - (min(Pos) - 1)):(ncol(AAs) - (max(Pos) - 1))]  # Isolate desired positions
    
    paste_res <- function(x){  # Paste implementation with defined parameters
      
      pasted <- paste(x,sep = "",collapse = "")
      return(pasted)
    }
    
    AAs <- apply(AAs,1,paste_res)  # apply paste implementation
    
    Result <- Peps[AAs == Res,]  # Conditional Selection
    return(Result)
  }
}
