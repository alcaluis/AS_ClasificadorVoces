# PGM 

PGM <- function(signal) {
  
  melfcc <- melfcc(signal,
                  minfreq = 50,
                  maxfreq = 4000,
                  sr = signal@samp.rate, 
                  hoptime = 0.01,
                  wintime = 0.025, usecmp = TRUE) 
  
  coefs <- gsignal::idct(melfcc)
  
  medias_filas <- apply(coefs, 1, mean)
  covarianzas <- cov(t(coefs))
  diagonal <- diag(covarianzas)
  
  c1 <- (medias_filas - min(medias_filas)) / (max(medias_filas)- min(medias_filas))
  c2 <- (diagonal - min(diagonal)) / (max(diagonal)- min(diagonal))

  # Crear un data frame con medias y diagonas de varianzas
  meanvar_carac <- data.frame(
    medias = c1,
    varianzas = c2)
  
  return(meanvar_carac)
  
  
}



# PS

PS <- function(signal, sampling_rate, span=NULL) {

  # Calcular el espectro de potencia usando la función spectrum
  if (is.null(span)) {
    spec <- spectrum(signal, method = "pgram", log = "no", plot = FALSE)
  } else {
    spec <- spectrum(signal, method = "pgram", log = "no", plot = FALSE, spans = span)
  }  
  
  # Crear un data frame con las frecuencias y la potencia
  espectro_potencia <- data.frame(
    Frecuencia = spec$freq * sampling_rate,
    Potencia = spec$spec
  )
  
  return(espectro_potencia)
}




