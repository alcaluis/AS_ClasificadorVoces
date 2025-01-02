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

# SH
SH <- function(signal){
  
  shimmer <- 20*log10(sum(abs(diff(signal))/abs(signal[-length(signal)]))/(length(signal)-1))
  
  return(shimmer)
}

PLP <- function(wave, wlen = 512, ovlp = 50, side = "l") {
  library(seewave)
  library(tuneR)
  library(phonTools)
  library(pracma)
  library(signal)
  fs <- wave@samp.rate
  # PASO 1: STFT, window name, see ftwindow (by default "hanning") <- DOCUMENTACION spectro
  if(side == 'l'){
    stft_result <- spectro(wave = wave@left, f = fs, wl = wlen, ovlp = ovlp, plot = F)
  }else{
    stft_result <- spectro(wave = wave@right, f = fs, wl = wlen, ovlp = ovlp, plot = F)
    
  }
  
  # PASO 2: PSD (Power Spectral Density)
  psd <- abs(stft_result$amp)^2
  
  # PASO 3: BARK FILTER-BANK
  freqs <- stft_result$freq
  bark_freqs <- hz2bark(freqs)
  
  num_bark_bands <- 24
  bark_band_edges <- seq(0, max(bark_freqs), length.out = num_bark_bands + 1)
  filterbank_weights <- matrix(0, nrow = num_bark_bands, ncol = length(freqs))
  
  # FILTRO TRIANGULAR
  for (j in 1:num_bark_bands) {
    lower_edge <- bark_band_edges[j]
    center <- bark_band_edges[j + 1]
    upper_edge <- ifelse(j + 2 <= length(bark_band_edges), bark_band_edges[j + 2], center + (center - lower_edge))  # Handle last band
    
    filterbank_weights[j, ] <- pmax(
      0,
      pmin(
        (bark_freqs - lower_edge) / (center - lower_edge),
        (upper_edge - bark_freqs) / (upper_edge - center)
      )
    )
  }
  
  bark_outputs <- matrix(0, nrow = num_bark_bands, ncol = ncol(psd))
  for (j in 1:num_bark_bands) {
    bark_outputs[j, ] <- colSums(filterbank_weights[j, ] * psd)
  }
  
  # PASO 4: EQUAL LOUDNESS PRE-EMPHASIS
  equal_loudness <- function(f) {
    numerator <- (f^2 + 1.44e6) * f^4
    denominator <- (f^2 + 1.6e5)^2 * (f^2 + 9.61e6)
    return(numerator / denominator)
  }
  equal_loudness_weights <- equal_loudness(freqs)
  
  bark_band_centers <- (bark_band_edges[-1] + bark_band_edges[-length(bark_band_edges)]) / 2
  bark_band_centers_hz <- bark2hz(bark_band_centers)
  bark_equal_loudness_weights <- interp1(freqs, equal_loudness_weights, bark_band_centers_hz, method = "linear")
  
  bark_outputs_pre_emphasis <- sweep(bark_outputs, 1, bark_equal_loudness_weights, `*`)
  
  # PASO 5: ^0.33 SIMULAR EL OÍDO HUMANO
  power_law_exponent <- 0.33
  
  perceptual_coefficients <- bark_outputs_pre_emphasis^power_law_exponent
  
  # PASO 6: Calcula los LPC (coeficientes lineales predictivos) para cada columna de la matriz de coeficientes perceptuales.
  lpc_order <- 12 
  num_columns <- ncol(perceptual_coefficients)
  lpc_matrix <- matrix(0, nrow = lpc_order + 2, ncol = num_columns) 
  
  for (i in 1:num_columns) {
    lpc_result <- lpc(perceptual_coefficients[, i], p = lpc_order)
    #print(lpc_result)
    lpc_matrix[, i] <- lpc_result  # Extract LPC coefficients
  }
  #print(lpc_matrix[1:10,1:10])
  #PASO 7: Cálculo de los coeficientes cepstrales
  cepstral_matrix <- lpc2cep(lpc_matrix, nrow(lpc_matrix))
  
  return(cepstral_matrix)
}
