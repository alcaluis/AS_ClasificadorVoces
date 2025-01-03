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
ZCR_FUN <- function(signal, sampling_rate, window_size=2) {
  # Duración de la señal en tiempo
  longitud_senal <- length(signal) / sampling_rate
  # Número de ventanas
  num_ventanas <- round(longitud_senal / window_size)
  vector_zcr <- zcr(signal, sampling_rate, 
                    wl = length(signal)/num_ventanas,
                    ovlp = 50, plot = FALSE)
  # Crear pesos para la media ponderada
  pesos <- rep(1, num_ventanas)
  return(mean(vector_zcr))
}

# Función para calcular LPCC
calc_lpcc <- function(signal, frame_size, overlap, order) {
  
  # Preprocesamiento: segmentar la señal en frames
  frames <- list()
  step_size <- frame_size - overlap
  num_frames <- floor((length(signal) - overlap) / step_size)
  
  for (i in 1:num_frames) {
    start <- (i - 1) * step_size + 1
    end <- start + frame_size - 1
    frames[[i]] <- signal[start:end]
  }
  
  # Resultados de LPCC para cada frame
  lpcc_result <- list()
  
  for (frame in frames) {
    # Calcular los coeficientes LPC
    lpc_coeffs <- lpc(frame, order = order)
    
    # Inicializar el vector LPCC
    lpcc <- numeric(order)
    lpcc[1] <- lpc_coeffs[1]  # El primer LPCC es igual al primer LPC
    
    # Calcular los LPCC utilizando la transformación cepstral
    for (n in 2:order) {
      lpcc[n] <- lpc_coeffs[n]
      for (k in 1:(n - 1)) {
        lpcc[n] <- lpcc[n] - (k / n) * lpcc[k] * lpc_coeffs[n - k]
      }
    }
    
    # Guardar los coeficientes LPCC para el frame
    lpcc_result[[length(lpcc_result) + 1]] <- lpcc
  }
  
  return(lpcc_result)
}

# PGM 
PGM <- function(signal, coef=12) {
  melfcc <- melfcc(signal,
                  minfreq = 50,
                  maxfreq = 4000,
                  sr = signal@samp.rate, 
                  hoptime = 0.01,
                  numcep = coef,
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

# HNR
HNR <- function(wave) {
  
  audio_signal <- wave@left
  
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

