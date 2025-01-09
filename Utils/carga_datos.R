carga_audios <- function(path) {
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

  metadata <- data.frame(
    genero, pista, locutor,
    p_audios, id_audio)
  
  metadata <- metadata %>%
    mutate(locutor = factor(locutor, levels=c("01", "02", "03", 
                                              "05", "06", "04",
                                              "07", "08", "09",
                                              "10"))) %>%
    mutate(pista = factor(pista))
             
  return(list(audios, metadata))
}

normalizacion <- function(audios, duracion=10, norm_amplitud=TRUE) {
  # Partiendo de que todas las frecuencias de muestreo son
  # iguales.
  
  maximum <- max(audios[[1]]@left)
  minimum <- min(audios[[1]]@left)
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
    if (maximum < max(audios[[id_audio]]@left)) {
      maximum <- max(audios[[id_audio]]@left)
    }
    if (minimum < min(audios[[id_audio]]@left)) {
      minimum <- min(audios[[id_audio]]@left)
    }
  }
  
  if (norm_amplitud) {
    for (id_audio in 1:length(audios)) {
      audios[[id_audio]]@left <- 2 * (audios[[id_audio]]@left - minimum) /
        (maximum - minimum) - 1
      audios[[id_audio]]@right <- 2 * (audios[[id_audio]]@right - minimum) /
        (maximum - minimum) - 1
    }
  }
  
  return(audios)
}

limpieza_ruido <- function(audios, niveles = 3) {
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
limpieza_senales <- function(audios) {
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

carga_formants <- function(path, audios) {
  df_formants <- data.frame(f1_avg = numeric(70), f1_sd = numeric(70),
                            f2_avg = numeric(70), f2_sd = numeric(70),
                            f3_avg = numeric(70), f3_sd = numeric(70))
  
  fichero_praat <- readLines(path)
  
  # Eliminamos comillas
  fichero_praat <- gsub("\"", "", fichero_praat)
  
  # Saltamos primeras lineas
  f_praat <- fichero_praat[-c(1:4)]
  
  # Por cada pista... seguiremos el índice
  ids_pista <- grep("Table", f_praat)
  
  id_audio <- 1
  for (pista in ids_pista) {
    # Saltar cabecera datos
    lin <- pista + 12
    
    # 10 segundos de datos, hay por lo general más
    num_v <- numeric()
    df_pista <- data.frame(f1=rep(c(num_v, NA), 2000),
                           f2=rep(c(num_v, NA), 2000),
                           f3=rep(c(num_v, NA), 2000))
    
    ventanas <- as.numeric(f_praat[(lin-1)]) - 1
    for (i in 0:ventanas) {
      # Instante del cálculo
      t <- as.numeric(f_praat[lin + i * 9 + 1])
      
      # ¿Se habla en este instante?
      hablado = TRUE
      
      # Formantes
      # Número de formantes encontrados para este instante
      n_for <- f_praat[lin + i * 9 + 2]
      for (f in 1:min(c(n_for, 3))) {
        df_pista[[f]][i] <- as.numeric(f_praat[lin + i * 9 + 2 + f])
      }
    }
    
    # Guardar valores de la pista
    df_formants[id_audio, ] <- c(mean(df_pista[[1]], na.rm=TRUE), sd(df_pista[[1]], na.rm=TRUE),
                                 mean(df_pista[[2]], na.rm=TRUE), sd(df_pista[[2]], na.rm=TRUE),
                                 mean(df_pista[[3]], na.rm=TRUE), sd(df_pista[[3]], na.rm=TRUE))
    id_audio <- id_audio + 1
  }
  
  return(cbind("id_audio"=rep(1:length(audios)), df_formants))
}

