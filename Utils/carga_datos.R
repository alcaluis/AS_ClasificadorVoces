carga_audios <- function(path) {
  library(tuneR)
  p_audios <- list.files(path = path, recursive=FALSE,
                         pattern = '\\.wav', full.names = TRUE)

  audios <- c()
  genero <- c()
  pista <- c()
  locutor <- c()
  id_audio <- c()
  id = 1
  for (p_audio in p_audios) {
    # Metadata
    f_audio <- strsplit(p_audio, "/")[[1]][4]
    f_audio <- gsub(".wav", "", f_audio)
    f_audio <- strsplit(f_audio, "_")[[1]]
    genero <- append(genero, f_audio[1])
    pista <- append(pista, f_audio[2])
    locutor <- append(locutor, f_audio[3])
    
    # Leer audio
    audios <- append(audios, readWave(p_audio))
    id_audio <- append(id_audio, id)
    id = id + 1
  }
  metadata <- data.frame(genero, pista, locutor, p_audios, id_audio)
  return(list(audios, metadata))
}

normalizacion <- function(audios, duracion=10, norm_amplitud=TRUE) {
  # Partiendo de que todas las frecuencias de muestreo son
  # iguales.
  for (id_audio in 1:length(audios)) {
    # Normalizar duracion
    nueva_len <- audios[[id_audio]]@samp.rate * duracion
    audios[[id_audio]]@left <- audios[[id_audio]]@left[1:nueva_len]
    audios[[id_audio]]@right <- audios[[id_audio]]@right[1:nueva_len]
    
    # Normalizar amplitud
    if (!norm_amplitud) {
      next
    }
    
    # Se asume canal MONO
    min_amp <- min(audios[[id_audio]]@left)
    max_amp <- max(audios[[id_audio]]@left)
    audios[[id_audio]]@left <- 2 * (audios[[id_audio]]@left - min_amp) /
                                   (max_amp - min_amp) - 1
    audios[[id_audio]]@right <- 2 * (audios[[id_audio]]@right - min_amp) /
                                    (max_amp - min_amp) - 1
    
  }
  
  return(audios)
}

limpieza_ruido <- function(audios, niveles = 3) {
  # WIP
  library(wavelets)
  
  for (id_audio in 1:length(audios)) {
    dwt_sn <- dwt(as.numeric(audios[[id_audio]]@left),
                  n.levels = niveles)
    for (id in 1:niveles) {
      dwt_sn@W[[id]] <- matrix(0,
                               nrow = nrow(dwt_sn@W[[id]]),
                               ncol = ncol(dwt_sn@W[[id]]))
    }

    idwt_sn <- idwt(dwt_sn)
    audios[[id_audio]]@left <- as.integer(idwt_sn)
    audios[[id_audio]]@right <- as.integer(idwt_sn)
  }
  
  return(audios)
}

# Función para preprocesar señales de audio
limpieza_señales <- function(audios) {

  # Cargar la biblioteca necesaria
  library(signal)
  
  # Inicializar lista de resultados
  audios_filtrados <- list()
  
  for (id_audio in 1:length(audios)) {
    # Validar que el elemento es un objeto Wave
    audio <- audios[[id_audio]]

    
    # Sacamos la frecuencia de muestreo y las señales
    fs <- audio@samp.rate
    signal_left <- audio@left
    signal_right <- audio@right
    
    # Corregimos la línea base
    signal_left_DC <- signal_left - mean(signal_left)
    signal_right_DC <- signal_right - mean(signal_right)
    
    # Diseño del filtro
    nyquist <- fs / 2
    filtro <- butter(4, c(60, 400) / nyquist, type = "pass")  # Orden 4, pasa banda
    
    # Aplicar el filtro a las señales
    signal_left_filtered <- filtfilt(filtro, signal_left_DC)
    signal_right_filtered <- filtfilt(filtro, signal_right_DC)
    
    # Guardar el audio procesado
    audio@left <- signal_left_filtered
    audio@right <- signal_right_filtered
    audios_filtrados[[id_audio]] <- audio
  }
  
  return(audios_filtrados)
}
