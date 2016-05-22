aaDist <- function(x,plotName = NULL, linker = TRUE)
{
  Peps <- x[,1]
  Long <- paste(Peps,collapse = "")
  Sep <- strsplit(Long,"")
  Tab <- data.frame(table(Sep))

  if (linker == TRUE){

    Tab[Tab == "G",2] <- Tab[Tab == "G",2] - (2*length(x[,1]))
    Tab[Tab == "S",2] <- Tab[Tab == "S",2] - (length(x[,1]))
  }

  Tab[,2] <- Tab[,2]/length(Sep[[1]])
  localenv <- environment()

  SigPlot <- ggplot(data = Tab,aes(x = Tab[,1], y = Tab[,2]),environment = localenv)+
    geom_bar(aes(fill = Tab[,1]),stat = "identity")+
    labs(fill = "Amino Acid", x = "Amino Acid",y = paste("Frequency, n = ",length(x[,1])),
         title = paste(plotName,"AA Distribution"))

  print(SigPlot)

  return(Tab)
}
