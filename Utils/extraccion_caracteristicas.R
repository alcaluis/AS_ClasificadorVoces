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

# F0
FO <- function(signal, sampling_rate, window_size = 2) {
  # Duración de la señal en tiempo
  longitud_senal <- length(signal) / sampling_rate
  # Número de ventanas
  num_ventanas <- round(longitud_senal / window_size)
  
  # Vector para almacenar frecuencias fundamentales
  frecuencia_fundamental <- numeric(num_ventanas)
  
  for (i in 1:num_ventanas) {
    # Índices de inicio y fin para cada ventana
    inicio <- ((i - 1) * window_size * sampling_rate) + 1
    final <- min(i * window_size * sampling_rate, length(signal))
    
    # Extraer el segmento de la señal
    segment <- signal[inicio:final]
    n <- length(segment)  # Tamaño del segmento
    
    # Transformada de Fourier del segmento
    segment_fft <- fft(segment)
    mod_fft <- Mod(segment_fft)  # Módulo de la FFT
    
    # Crear vector de frecuencias
    freqs <- (0:(n - 1)) * (sampling_rate / n)
    
    # Tomar solo la mitad positiva del espectro
    half_n <- floor(n / 2)
    mod_fft <- mod_fft[1:half_n]
    freqs <- freqs[1:half_n]
    
    # Índice de la frecuencia máxima
    idx_max <- which.max(mod_fft)
    frecuencia_fundamental[i] <- freqs[idx_max]
  }
  # Crear pesos para la media ponderada
  pesos <- rep(1, num_ventanas)
  # Eliminamos la aportación de los 4 últimos segundos, pues en ellos las personas ya suelen haber terminado de hablar
  pesos[(length(pesos) - 1):length(pesos)] <- 0
  
  # Calcular la media ponderada
  media_ponderada <- sum(frecuencia_fundamental * pesos, na.rm = TRUE) / sum(pesos, na.rm = TRUE)
  
  return(media_ponderada)
}

# ZCR
ZCR_FUN<-function(signal,sampling_rate,window_size=2){
  # Duración de la señal en tiempo
  longitud_senal <- length(signal) / sampling_rate
  # Número de ventanas
  num_ventanas <- round(longitud_senal / window_size)
  vector_zcr<-zcr(signal,sampling_rate,wl=length(signal)/num_ventanas,ovlp = 50,plot = FALSE)
  # Crear pesos para la media ponderada
  pesos <- rep(1, num_ventanas)
  return(mean(vector_zcr))
}
