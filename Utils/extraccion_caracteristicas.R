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


# PS

PS <- function(signal, sampling_rate, span=NULL) {
  # Duración de la señal en tiempo
  longitud_senal <- length(signal) / sampling_rate
  
  # Calcular el espectro de potencia usando la función spectrum
if (is.null(span)) {
  spec <- spectrum(signal, method = "pgram", log = "dB", plot = FALSE)
} else {
  spec <- spectrum(signal, method = "pgram", log = "dB", plot = FALSE, spans = span)
}  

# Crear un data frame con las frecuencias y la potencia
espectro_potencia <- data.frame(
  Frecuencia = spec$freq * sampling_rate,
  Potencia = spec$spec
)
  
  return(espectro_potencia)
}


