# STE
STE <- function(wave, wlen, ovlp = 0) {
  tm <- 1 / wave@samp.rate
  N <- length(wave@left)
  t <- seq(0, (N-1)*tm, length.out=N)
  n <- length(wave)
  step10 <- seq(1, n-wlen, wlen - (ovlp * wlen/100))
  m10 <- length(step10)
  
  STE <- numeric(m10)
  time.frame <- numeric(m10)
  for (i in 1:m10) {
    frame.wave <- wave[step10[i]:(step10[i]+wlen)]
    STE[i] <- mean(frame.wave@left * frame.wave@left)
    time.frame[i] <- mean(t[step10[i]: (step10[i] + wlen)])
  }
  # NOTA: Devuelve normalizado
  return(cbind(time.frame, STE/max(STE)))
}



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







# POWER SPECTRUM ??

PS_2 <- function(signal, sampling_rate) {

    n <- length(signal)
    fs <- sampling_rate
    
    # Calcular la FFT (Transformada Rápida de Fourier)
    fft_result <- fft(signal)
    
    # Calcular el espectro de potencia (módulo al cuadrado de la FFT)
  
    
    xk <- (1/(n^2)) * abs(fft_result)^2
    
    # Calcular el espectro de potencia de una cara
    xk <- xk[1:(n/2+1)]
    
    # Duplicar los valores (excepto DC y Nyquist) para el espectro de una cara
    xk[2:(length(xk)-1)] <- 2 * xk[2:(length(xk)-1)]
    
    
    # Calcular las frecuencias
    frequencies <- (0:(n-1)) * (fs / n)
    
  
  # Crear un data frame con las frecuencias y la potencia
  espectro_potencia <- data.frame(
    Frecuencia = frequencies,
    Potencia = xk
  )
  
  return(espectro_potencia)
}