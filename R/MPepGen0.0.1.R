#######################genPep0.0.1##################################
genPep <- function(Struct,draw){
  
  pepStruct <- data.frame(matrix(0,ncol = ncol(Struct),nrow = draw))
  
  genRes <- function(x){
    coeffs <- x
    coeffs[is.nan(coeffs)] <- 0
    peps <- c("A","R","N","D","Q","E","G","H","L","K",
              "F","P","S","W","Y","V","C","I","M","T")
    sel <- sample(peps,draw,replace = TRUE, prob = coeffs^2)
    return(sel)
  }
  
  for (i in 1:ncol(Struct)){
    pepStruct[,i] <- genRes(Struct[,i])
  }
  return(pepStruct)
  }
