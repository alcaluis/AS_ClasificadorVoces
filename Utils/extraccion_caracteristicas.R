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

# HNR
HNR <- function(wave) {
  
  library(wavelets)
  
  audio_signal <- db@left
  
  # Realizamos una descomposición wavelet de la señal
  wt <- dwt(as.numeric(audio_signal), filter = 'la8')
  
  # Interpretamos la Aprox. como parte harmónica y el Detalle como ruido
  harmonic <- wt@V[[1]]
  noise <- wt@W[[1]]
  
  # Calculamos la energía de cada componente
  harmonic_energy <- sum(abs(harmonic))
  noise_energy <- sum(abs(noise))
  
  # Habitualmente se representa la proporción de ambas usando un log base 10
  hnr_db <- 10 * log10(harmonic_energy / noise_energy)
  return(hnr_db)
}

# CE

CE <- function(wave) {
  
  sampling_rate <- wave@samp.rate
  audio_signal <- wave@left
  
  # Evitamos que la función meanspec devuelva gráficas o listas de números
  suppressMessages({
    # Obtenemos la amplitud media de la distribución freuencial
    # (notamos que también se podría usar spec)
    spectral_wave <- meanspec(audio_signal, f = sampling_rate, wl = 1024, wn = "hanning", plot=FALSE)
  })
  
  # Calculamos propiedades estadísticas del espectro de frecuencias
  spectral_prop <- specprop(spectral_wave, f=sampling_rate)
  
  # Devolvemos el centroide espectral presente en el specprop
  return(spectral_prop$cent)
}