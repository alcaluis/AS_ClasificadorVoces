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
    lpc_coeffs <- lpc(frame, order = order)$a  # Coeficientes LPC
    
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

