##########################################################################
##ggplot p-value reporting######Master's Alengths File####################
##2.0 hypothesis testing validated########Component Functions#############
###########################Function Start#################################

lcScan <- function(File, Glob = NULL, Conc = 5, Kd = FALSE){  #Start Here

  SeM <- function(x){  # Standard error of the mean Calc

    (sd(x)/sqrt(length(x)))
  }

  Lens <- min(File[,2]):max(File[,2])

  Cgen <- function(x){  # L/C Comb matrix gen

    a <- matrix(min(File[,3]):max(File[,3]),ncol = 1)
    b <- cbind(x,a)
    return(b)
  }

  LCmat <- lapply(Lens,Cgen)
  LCmat <- do.call(rbind,LCmat)  # Len/Char Matrix Assembly
  colnames(LCmat) <- c("Length","Charge")

  AlChargeSep <- function (x,y,z){  # L/C Separation Fx Specific

    LC<-x[x[,2] == y & x[,3] == z,Conc]
    return(LC)
  }

  if (is.null(Glob) == FALSE){

    Calcs <- function(x){  # Attribute Calculations

      Peps <- AlChargeSep(File,x[1],x[2])
      nDist <- AlChargeSep(Glob,x[1],x[2])

      if (length(Peps) > 1 & length(Glob) > 1){
        Difs <- t.test(Peps,nDist,paired = FALSE, equal.var = FALSE)
        Dif <- Difs$p.value  # must be numeric

      }else{
        Dif<-NA
      }
        if(Kd == TRUE){  # new control for Kd addition
        Mets <- matrix(c(mean(Peps),SeM(Peps),length(Peps), Dif),nrow = 1)
        return(Mets)

        }else if (Kd == FALSE){
        Mets <- matrix(c(mean(Peps),SeM(Peps),length(Peps), Dif),nrow = 1)
        return(Mets)
      }
    }

    LCalcs <- apply(LCmat,1,Calcs)
    LCalcs <- t(LCalcs)  # Attributes Maxtrix Assembly
    colnames(LCalcs) <- c("Mean","SeM","n","p-value")

    Final <- data.frame(cbind(LCmat,LCalcs)) # Final Matrix Assembly

    sig <- Final[Final[,6] < 0.05 & is.na(Final[,6]) == FALSE,]
    nsig <- Final[Final[,6] >= 0.05 | is.na(Final[,6]) == TRUE,]

    print(ggplot() + geom_raster(aes(x = Final[,1], y = Final[,2], fill = Final[,3]))+  # Create filled Rects
      geom_point(colour = "white", aes(x = nsig[,1], y = nsig[,2], size = nsig[,5]))+  #  non-significant pep sets
      geom_point(colour = "black", aes(x = sig[,1], y = sig[,2], size = sig[,5]))+  # significant pep sets
      geom_rect(size = 0.5,fill = NA,colour = "black",
                aes(xmin = Final[,1]-0.5,xmax = Final[,1]+0.5,ymin = Final[,2]-0.5,ymax = Final[,2]+0.5))+  # grid lines

      scale_x_continuous(breaks = min(Final[,1]):max(Final[,1]))+  # Scale and fill def
      scale_y_continuous(breaks = min(Final[,2]):max(Final[,2]))+
      scale_fill_distiller(ifelse(Kd == TRUE,"Kd","Signal"), palette = "YlOrRd", direction = ifelse(Kd ==TRUE,-1,1))+

      ggtitle("lcScan Comparative Analysis")+  # text options and theme
      labs(x = "Length", y = "Charge", size = "n")+
      theme_dark())

    return(Final)

  } else{

    Calcs <- function(x){  # Attribute Calculations

      Peps <- AlChargeSep(File,x[1],x[2])

      if (Kd == TRUE){  # new control for Kd addition
        Mets <- matrix(c(mean(Peps),SeM(Peps),length(Peps)),nrow = 1)
        return(Mets)

        }else if (Kd == FALSE){
        Mets <- matrix(c(mean(Peps),SeM(Peps),length(Peps)),nrow = 1)
        return(Mets)
      }
    }

    LCalcs <- apply(LCmat,1,Calcs)
    LCalcs <- t(LCalcs)  # Attributes Maxtrix Assembly
    colnames(LCalcs) <- c("Mean","SeM","n")

    Final <- data.frame(cbind(LCmat,LCalcs)) # Final Matrix Assembly

    print(ggplot() + geom_raster(aes(x = Final[,1], y = Final[,2], fill = Final[,3]))+  # Colored squares for heatmap
            geom_point(colour = "black", aes(x = Final[,1], y = Final[,2], size = Final[,5]))+  # Population rep
            geom_rect(size = 0.5,fill = NA,colour = "black",
                      aes(xmin = Final[,1]-0.5,xmax = Final[,1]+0.5,ymin = Final[,2]-0.5,ymax = Final[,2]+0.5))+  # Gridlines

            scale_x_continuous(breaks = min(Final[,1]):max(Final[,1]))+  # Scale defs and Fill
            scale_y_continuous(breaks = min(Final[,2]):max(Final[,2]))+
            scale_fill_distiller(ifelse(Kd == TRUE, "Kd","Signal"), palette = "YlOrRd", direction = ifelse(Kd == TRUE,-1,1))+

            ggtitle("lcScan Raw Analysis")+  # plot text and theme
            labs(x = "Length", y = "Charge", size = "n")+
            theme_dark())

    return(Final)
  }
}

############################Function End#########################




