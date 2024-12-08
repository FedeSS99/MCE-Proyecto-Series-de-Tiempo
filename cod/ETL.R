##################
# Carga de datos #
##################

# Datos proporcionados para el proyecto
datosGA1 <- read.csv("./data/GA1.csv")
datosIVF <- read.csv("./data/IVF_SIAP.csv")

# Datos auxiliares 
dat_idf <- read.csv("./data/IDEFF_oct24.csv")
dat_idf <- dat_idf %>% group_by(YEAR, LEY, CONCEPTO, TIPO) %>% summarise(across(c(colnames(dat_idf)[7:18]), sum))

# Series de tiempo relativas a la delincuencia
aux <- dat_idf %>% filter(TIPO == 'CONTRA EL AMBIENTE Y LA GESTION AMBIENTAL')
for(r in 2015:2023){
  aux_1 <- aux %>% filter(YEAR == r)
  aux_1 <- t(aux_1[5:16])
  
  if (r == 2015){
    dat_idf1 <- aux_1
  }else{
    dat_idf1 <- rbind(dat_idf1, aux_1)
  }
}

aux <- dat_idf %>% filter(TIPO == 'VIAS DE COMUNICACION Y CORRESPONDENCIA')
for(r in 2015:2023){
  aux_2 <- aux %>% filter(YEAR == r)
  aux_2 <- t(aux_2[5:16])
  
  if (r == 2015){
    dat_idf2 <- aux_2
  }else{
    dat_idf2 <- rbind(dat_idf2, aux_2)
  }
}

# Series de tiempo relativas al clima
dat_pre <- data.frame()
dat_tem <- data.frame()

for (r in 2015:2023){
  
  if (r<2020){
    aux_1 <- data.frame(t(read_excel(paste("D:/Postgrado/Tercer_Semestre/Series_Tiempo/Proyecto/data/Precipitacion/", r, "Precip.xls", sep = ""), skip = 1, sheet = 1)[33,]))[2:13,1]  
  }else{
    aux_1 <- data.frame(t(read_excel(paste("D:/Postgrado/Tercer_Semestre/Series_Tiempo/Proyecto/data/Precipitacion/", r, "Precip.xlsx", sep = ""), skip = 1, sheet = 1)[33,]))[2:13,1]
  }
  aux_2 <- data.frame(t(read_excel(paste("D:/Postgrado/Tercer_Semestre/Series_Tiempo/Proyecto/data/Temperatura_Promedio_Excel/", r, "Tmed.xlsx", sep = ""), skip = 1, sheet = 1)[33,]))[2:13,1]
  
  aux_1 <- data.frame(aux_1)
  aux_2 <- data.frame(aux_2)
  
  if (r == 2015){
    dat_pre <- aux_1
    dat_tem <- aux_2
  }else{
    dat_pre <- rbind(dat_pre, aux_1)
    dat_tem <- rbind(dat_tem, aux_2)
  }
  
  
}

# Otras series de tiempo auxiliares
dat_die <- t(data.frame(read_excel("./data/Precios_promedio_diarios_y_mensuales_en_estaciones_de_servicio.xlsx", sheet = 5, skip = 4))[1, 14:95])

# intervalo de tiempo que comparten ene2015-dic2023
GA1 <- ts(datosGA1[265:372, 2], start = 2015, frequency = 12)
IVF <- ts(datosIVF[1:108, 2], start = 2015, frequency = 12)
IDF1 <- ts(dat_idf1, start = c(2015, 1), frequency=12)
IDF2 <- ts(dat_idf2, start = c(2015, 1), frequency=12)
PREC <- ts(dat_pre, start = c(2015, 1), frequency=12)
TEMP <- ts(dat_tem, start = c(2015, 1), frequency=12)
DIES <- ts(dat_die, start = c(2017, 1, 1), frequency = 12)
