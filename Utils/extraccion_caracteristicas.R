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
