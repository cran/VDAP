###0.0.2 Row.names for residues#####Motif and plot handling##################################################
###############################################AA Probability Matrix Start#####################################

## Unite AAstruct and AAstructnS functions using conditional switch == less code ##
aaStruct <- function (x,y, sigWeight = TRUE){ # x = Charge Sep object y = charge sep object$signal set

  Peps <- as.character(x[,1])
  Peps <- strsplit(Peps,"")  # Splitting of Seq's into Matrix
  Peps <- do.call(rbind,Peps)

  A <- Peps == "A"
  R <- Peps == "R"
  N <- Peps == "N"
  D <- Peps == "D"
  Q <- Peps == "Q"
  E <- Peps == "E"
  G <- Peps == "G"
  H <- Peps == "H"
  L <- Peps == "L"
  K <- Peps == "K"
  Ef <- Peps == "F"      #Amino Acid Logicals
  P <- Peps == "P"
  S <- Peps == "S"
  W <- Peps == "W"
  Y <- Peps == "Y"
  V <- Peps == "V"
  C <- Peps == "C"
  I <- Peps == "I"
  M <- Peps == "M"
  Tee <- Peps == "T"

  ProbTot <- length(x[,1])

  if (sigWeight == FALSE){

    cat("sigWeight = FALSE \n")

    A2 <- (colSums(A)/ProbTot)
    R2 <- (colSums(R)/ProbTot)  # R = Boolean Matrix, RawFile = Signal Intensity
    N2 <- (colSums(N)/ProbTot)  # (colSums(x)/ProbTot)-> Prob ratio of occurence
    D2 <- (colSums(D)/ProbTot)  # colSums(y)/colSums(x)-> Average Signal intensity of AA @ pos
    Q2 <- (colSums(Q)/ProbTot)
    E2 <- (colSums(E)/ProbTot)
    G2 <- (colSums(G)/ProbTot)
    H2 <- (colSums(H)/ProbTot)  # Average Signal each AA
    L2 <- (colSums(L)/ProbTot)
    K2 <- (colSums(K)/ProbTot)
    Ef2 <- (colSums(Ef)/ProbTot)
    P2 <- (colSums(P)/ProbTot)
    S2 <- (colSums(S)/ProbTot)
    W2 <- (colSums(W)/ProbTot)
    Y2 <- (colSums(Y)/ProbTot)
    V2 <- (colSums(V)/ProbTot)
    C2 <- (colSums(C)/ProbTot)
    I2 <- (colSums(I)/ProbTot)
    M2 <- (colSums(M)/ProbTot)
    Tee2 <- (colSums(Tee)/ProbTot)

    CPeps <- rbind(A2,R2,N2,D2,Q2,E2,G2,H2,L2,K2,Ef2,P2,S2,W2,Y2,V2,C2,I2,M2,Tee2)  # bringing sums together
    return(CPeps)

  } else if (sigWeight == TRUE){

  cat("sigWeight default = TRUE \n")

  RawFile <- matrix(0,nrow(Peps),ncol(Peps))  # empty vector for loop, sized by pep frag obj

  Sigs <- y

  for (i in 1:ncol(RawFile)){  # Creation of signal by row matrix (Signal Intensities)
    RawFile[,i] = y
  }

  A1 <- A*RawFile
  R1 <- R*RawFile
  N1 <- N*RawFile
  D1 <- D*RawFile
  Q1 <- Q*RawFile
  E1 <- E*RawFile
  G1 <- G*RawFile
  H1 <- H*RawFile  # Signal Assignment to AA position
  L1 <- L*RawFile
  K1 <- K*RawFile
  Ef1 <- Ef*RawFile
  P1 <- P*RawFile
  S1 <- S*RawFile
  W1 <- W*RawFile
  Y1 <- Y*RawFile
  V1 <- V*RawFile
  C1 <- C*RawFile
  I1 <- I*RawFile
  M1 <- M*RawFile
  Tee1 <- Tee*RawFile

  A2 <- (colSums(A)/ProbTot)*colSums(A1)/colSums(A)
  R2 <- (colSums(R)/ProbTot)*colSums(R1)/colSums(R)  # R = Boolean Matrix, RawFile = Signal Intensity
  N2 <- (colSums(N)/ProbTot)*colSums(N1)/colSums(N)  # (colSums(x)/ProbTot)-> Prob ratio of occurence
  D2 <- (colSums(D)/ProbTot)*colSums(D1)/colSums(D)  # colSums(y)/colSums(x)-> Average Signal intensity of AA @ pos
  Q2 <- (colSums(Q)/ProbTot)*colSums(Q1)/colSums(Q)
  E2 <- (colSums(E)/ProbTot)*colSums(E1)/colSums(E)
  G2 <- (colSums(G)/ProbTot)*colSums(G1)/colSums(G)
  H2 <- (colSums(H)/ProbTot)*colSums(H1)/colSums(H)  # Average Signal each AA
  L2 <- (colSums(L)/ProbTot)*colSums(L1)/colSums(L)
  K2 <- (colSums(K)/ProbTot)*colSums(K1)/colSums(K)
  Ef2 <- (colSums(Ef)/ProbTot)*colSums(Ef1)/colSums(Ef)
  P2 <- (colSums(P)/ProbTot)*colSums(P1)/colSums(P)
  S2 <- (colSums(S)/ProbTot)*colSums(S1)/colSums(S)
  W2 <- (colSums(W)/ProbTot)*colSums(W1)/colSums(W)
  Y2 <- (colSums(Y)/ProbTot)*colSums(Y1)/colSums(Y)
  V2 <- (colSums(V)/ProbTot)*colSums(V1)/colSums(V)
  C2 <- (colSums(C)/ProbTot)*colSums(C1)/colSums(C)
  I2 <- (colSums(I)/ProbTot)*colSums(I1)/colSums(I)
  M2 <- (colSums(M)/ProbTot)*colSums(M1)/colSums(M)
  Tee2 <- (colSums(Tee)/ProbTot)*colSums(Tee1)/colSums(Tee)

  CPeps <- rbind(A2,R2,N2,D2,Q2,E2,G2,H2,L2,K2,Ef2,P2,S2,W2,Y2,V2,C2,I2,M2,Tee2)  # bringing sums together

  CPepsNorm <- CPeps/(mean(y)) # Normalizing to average signal of Length/Charge Set
  return(CPepsNorm)

  } else {

    cat("Error in sigWeight, select either TRUE or FALSE")
  }
}

#############################################Length Charge and Prob Matrix Combined#######################################################

vMotif.lc <- function(Prot,ProtG,Length,Charge,SigCol,Kd = FALSE)  # AA Prob Matrix Input = Charge State, Length, Pepstruct$Signal Column

{
  Peps <- vSep(Prot,Length,Charge) # Creation of Hits Matrix

  if (Kd == TRUE){  # control for added Kd calcs

    a <- aaStruct(Peps,1/Peps[,SigCol])

  }else if (Kd == FALSE){

    a <- aaStruct(Peps,Peps[,SigCol])
  }

  Peps <- vSep(ProtG,Length,Charge)  # Creation of Global Matrix, normalize to all Peps of target length

  if (Kd == TRUE){

    b <- aaStruct(Peps,1/Peps[,SigCol])

  }else if (Kd == FALSE){

    b <- aaStruct(Peps,Peps[,SigCol])
  }

  Finalstruct <- (a/b)  # Hits/Global

  Final <- data.frame(Finalstruct, row.names = c("A","R","N","D","Q","E","G","H","L","K",
                                               "F","P","S","W","Y","V","C","I","M","T"))
  colnames(Final) <- Length:1

  aaMap <- melt(as.matrix(Final))
  colnames(aaMap) <- c("AA","Position","Weight")

  aaPlot <- ggplot(data = aaMap) + geom_raster(aes(x = aaMap[,2], y = aaMap[,1], fill = aaMap[,3],
                                               alpha = ifelse(is.nan(aaMap[,3]) == TRUE, 0, 1)))+
          scale_fill_distiller(palette = "Spectral")+
          geom_tile(size = 0.5, fill = NA, colour = "black", aes(x = aaMap[,2], y = aaMap[,1]))+
          scale_x_continuous(breaks = min(aaMap[,2]):max(aaMap[,2]))+

          theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
                panel.background = element_rect(fill = "white"), axis.text = element_text(size=12),
                axis.title = element_text(size = 14), plot.title = element_text(size = 17, face = "bold"))+

          labs(title = paste("vMotif.lc Analysis Length",Length, "Charge",Charge), x = "Position", y = "AA", fill = "Weight")+
          guides(alpha = FALSE)

  print(aaPlot)

  return(Final)
}


vMotif.l <- function(Prot,ProtG,Length,SigCol,Kd = FALSE)  # AA Prob Matrix Input = Charge State, Length, Pepstruct$Signal Column

{
  Peps <- vSep(Prot,Length)  # All Charge states are together for hits

  if (Kd == TRUE){  # Control for added Kd calcs

    a <- aaStruct(Peps,1/Peps[,SigCol])

  }else if (Kd == FALSE){

    a <- aaStruct(Peps,Peps[,SigCol])
  }

  Peps <- vSep(ProtG,Length)

  if (Kd == TRUE){

    b <- aaStruct(Peps,1/Peps[,SigCol])

  }else if (Kd == FALSE){

    b <- aaStruct(Peps,Peps[,SigCol])
  }

  Finalstruct <- (a/b)

  Final <- data.frame(Finalstruct, row.names = c("A","R","N","D","Q","E","G","H","L","K",
                                                      "F","P","S","W","Y","V","C","I","M","T"))
  colnames(Final) <- Length:1

  aaMap <- melt(as.matrix(Final))
  colnames(aaMap) <- c("AA","Position","Weight")

  aaPlot <- ggplot(data = aaMap) + geom_raster(aes(x = aaMap[,2], y = aaMap[,1], fill = aaMap[,3],
                                               alpha = ifelse(is.nan(aaMap[,3]) == TRUE, 0, 1)))+
          scale_fill_distiller(palette = "Spectral")+
          geom_tile(size = 0.5, fill = NA, colour = "black", aes(x = aaMap[,2], y = aaMap[,1]))+
          scale_x_continuous(breaks = min(aaMap[,2]):max(aaMap[,2]))+

          theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
                panel.background = element_rect(fill = "white"), axis.text = element_text(size = 12),
                axis.title = element_text(size = 14), plot.title = element_text(size = 17, face = "bold"))+

          labs(title = paste("vMotif.l Analysis Length", Length), x = "Position", y = "AA", fill = "Weight")+

          guides(alpha = FALSE)

  print(aaPlot)

  return(Final)
}

############################################Full Function####################################################

vComp.l <- function(Prot,ProtG,Length)   # No Signal Term in Prob Matrix and Length Sep Only

{
  Peps <- vSep(Prot,Length)
  a <- aaStruct(Peps,sigWeight = FALSE)

  Peps <- vSep(ProtG,Length)
  b <- aaStruct(Peps,sigWeight = FALSE)

  Finalstruct <- (a/b)

  Final <- data.frame(Finalstruct, row.names = c("A","R","N","D","Q","E","G","H","L","K",
                                               "F","P","S","W","Y","V","C","I","M","T"))
  colnames(Final) <- Length:1

  aaMap <- melt(as.matrix(Final))
  colnames(aaMap) <- c("AA","Position","Weight")

  aaPlot <- ggplot(data = aaMap) + geom_raster(aes(x = aaMap[,2], y = aaMap[,1], fill = aaMap[,3],
                                     alpha = ifelse(is.nan(aaMap[,3]) == TRUE, 0, 1)))+
    scale_fill_distiller(palette = "Spectral")+
    geom_tile(size = 0.5, fill = NA, colour = "black", aes(x = aaMap[,2], y = aaMap[,1]))+
    scale_x_continuous(breaks = min(aaMap[,2]):max(aaMap[,2]))+

      theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
          panel.background = element_rect(fill = "white"), axis.text = element_text(size=12),
          axis.title = element_text(size=14), plot.title = element_text(size = 17, face = "bold"))+

      labs(title = paste("vComp.l Analysis Length", Length), x = "Position", y = "AA", fill = "Weight")+
    guides(alpha = FALSE)

  print(aaPlot)

  return(Final)
}

vComp.lc <- function(Prot,ProtG,Length,Charge)  # No Signal Term but Length+Charge Separation Input = Charge, Length

{
  Peps <- vSep(Prot,Length,Charge)
  a <- aaStruct(Peps,sigWeight = FALSE)

  Peps <- vSep(ProtG,Length,Charge)
  b <- aaStruct(Peps,sigWeight = FALSE)

  Finalstruct <- (a/b)

  Final <- data.frame(Finalstruct, row.names = c("A","R","N","D","Q","E","G","H","L","K",
                                               "F","P","S","W","Y","V","C","I","M","T"))
  colnames(Final) <- Length:1

  aaMap <- melt(as.matrix(Final))
  colnames(aaMap) <- c("AA","Position","Weight")

  aaPlot <- ggplot(data = aaMap) + geom_raster(aes(x = aaMap[,2], y = aaMap[,1], fill = aaMap[,3],
                                               alpha = ifelse(is.nan(aaMap[,3]) == TRUE, 0, 1)))+
          scale_fill_distiller(palette = "Spectral")+
          geom_tile(size = 0.5, fill = NA, colour = "black", aes(x = aaMap[,2], y = aaMap[,1]))+
          scale_x_continuous(breaks = min(aaMap[,2]):max(aaMap[,2]))+

          theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
                panel.background = element_rect(fill = "white"), axis.text = element_text(size=12),
                axis.title = element_text(size=14), plot.title = element_text(size = 17, face = "bold"))+

          labs(title =paste("vComp.lc Analysis Length", Length,"Charge", Charge), x = "Position", y = "AA", fill = "Weight")+
          guides(alpha = FALSE)

  print(aaPlot)

  return(Final)
}

#########################################################################################################


