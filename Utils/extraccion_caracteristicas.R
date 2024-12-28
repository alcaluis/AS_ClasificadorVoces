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




# LPCC
LPCC <- function(wave, wlen, order, ovlp = 0) {
  
  N <- length(wave@left)          
  step <- wlen - (ovlp * wlen / 100) 
  indices <- seq(1, N - wlen, by = step) 
  m <- length(indices)           
  
  LPCC_features <- list()        
  
  for (i in 1:m) {
    # Seleccionar ventana de audio
    start_idx <- indices[i]
    end_idx <- min(start_idx + wlen - 1, N)
    frame <- wave@left[start_idx:end_idx]
    
    # Calcular los LPC
    lpc_coeffs <- lpc(frame, order = order)  # Coeficientes LPC
    
    # Transformar LPC a LPCC
    lpcc <- numeric(order)
    for (n in 1:order) {
      if (n <= length(lpc_coeffs)) {
        lpcc[n] <- lpc_coeffs[n]
      }
      if (n > length(lpc_coeffs)) {
        lpcc[n] <- 0
        for (k in 1:(n - 1)) {
          lpcc[n] <- lpcc[n] + (k / n) * lpcc[k] * lpc_coeffs[n - k]
        }
      }
    }
    
    LPCC_features[[i]] <- lpcc  # Guardar LPCC de esta ventana
  }
  
  return(LPCC_features)
}


LPCC2 <- function(wave, wlen, order, ovlp = 0) {
  
  N <- length(wave@left)          
  step <- wlen - (ovlp * wlen / 100) 
  indices <- seq(1, N - wlen, by = step) 
  m <- length(indices)           
  
  LPCC_features <- list()        
  
  # Inicializar un vector para acumular los LPCCs promedio
  lpcc_promedio <- numeric(order)  # 12 LPCCs por ventana
  
  for (i in 1:m) {
    # Seleccionar ventana de audio
    start_idx <- indices[i]
    end_idx <- min(start_idx + wlen - 1, N)
    frame <- wave@left[start_idx:end_idx]
    
    # Calcular los LPC
    lpc_coeffs <- lpc(frame, order = order)  # Coeficientes LPC
    
    # Transformar LPC a LPCC
    lpcc <- numeric(order)
    for (n in 1:order) {
      if (n <= length(lpc_coeffs)) {
        lpcc[n] <- lpc_coeffs[n]
      }
      if (n > length(lpc_coeffs)) {
        lpcc[n] <- 0
        for (k in 1:(n - 1)) {
          lpcc[n] <- lpcc[n] + (k / n) * lpcc[k] * lpc_coeffs[n - k]
        }
      }
    }
    
    # Acumular los LPCCs de esta ventana para calcular el promedio
    LPCC_features[[i]] <- lpcc  # Guardar LPCC de esta ventana
    lpcc_promedio <- lpcc_promedio + lpcc  # Sumar al promedio
  }
  
  # Calcular el promedio de LPCCs de todas las ventanas
  lpcc_promedio <- lpcc_promedio / m
  
  # Devolver tanto los LPCCs de cada ventana como el promedio
  return(lpcc_promedio)
}




# library(pracma)

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