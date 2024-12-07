carga_audios <- function(path) {
  library(tuneR)
  p_audios <- list.files(path = path, recursive=FALSE,
                         pattern = '\\.wav', full.names = TRUE)

  audios <- c()
  genero <- c()
  pista <- c()
  locutor <- c()
  for (p_audio in p_audios) {
    cat("Audio:", p_audio)
    # Metadata
    f_audio <- strsplit(p_audio, "/")[[1]][4]
    f_audio <- gsub(".wav", "", f_audio)
    f_audio <- strsplit(f_audio, "_")[[1]]
    genero <- append(genero, f_audio[1])
    pista <- append(pista, f_audio[2])
    locutor <- append(locutor, f_audio[3])
    
    # Leer audio
    audios <- append(audios, readWave(p_audio))
  }
  metadata <- data.frame(genero, pista, locutor)
  return(c(audios, metadata))
}

normalizacion <- function(audios) {
  audios_normalizados <- c()
  # WIP
  # Duracion
  # Amplitud
  return(audios_normalizados)
}

limpieza_ruido <- function(audios) {
  audios_s_ruido <- c()
  # WIP
  return(audios_s_ruido)
}

