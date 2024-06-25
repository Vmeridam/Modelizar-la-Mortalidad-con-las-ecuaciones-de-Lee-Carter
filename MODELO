
###### Importar librerias

library(fpp2)
library(tseries)
library(quantmod)
library(stats)
library("goftest")

###### FUNCIONES PROPIAS ######

diagnosys_phase = function(modelo,time_serie_data) {
  
  problems = 0 #Numero de errores en el diagnostico
  
  # COEFICIENTES NO SIGNIFICATIVOS
  
  coef_no_significativos = 0
  
  tryCatch( {
    
    for (i in 1:length(modelo$coef)) {
      
      if (qnorm(c(0.025,0.975),0,(modelo$var.coef^0.5)[i,i])[1] < modelo$coef[i] &  qnorm(c(0.025,0.975),0,(modelo$var.coef^0.5)[i,i])[2] > modelo$coef[i] ) {
        
        problems = problems + 1 
        
        coef_no_significativos = coef_no_significativos + 1
        
      }
      
    }
    
  }, error = function(e) {
    
    problems = 1000
    return(problems)
    
  }
  
  )
  
  # RESIDUOS INCORRELACIONADOS
  
  residuos_incorrelacionados = 0
  
  res <- residuals(modelo) #Residuals
  ggacf_tsd = ggAcf(time_serie_data)
  
  incorrelacion_pvalue = Box.test(res, lag = max(ggacf_tsd$data$lag), type = c("Ljung-Box"))
  
  residuos_incorrelacionados = incorrelacion_pvalue$p.value
  
  # if (incorrelacion$p.value < 0.05) {
  #   
  #   problems = problems + 1
  #   
  #   residuos_incorrelacionados = residuos_incorrelacionados + 1
  #   
  # }
  
  ### NORMALIDAD EN LOS RESIDUOS
  
  normalidad_residuos = 0
  
  CvM_normal = cvm.test(res,"pnorm",0,sd(res),estimated = TRUE) # H0: Normalidad
  ad_normal = ad.test(res,"pnorm",0,sd(res),estimated = TRUE)
  jb = jarque.bera.test(res)
  
  normalidad_residuos = CvM_normal$p.value
  
  if (CvM_normal$p.value < 0.05) {
    
    problems = problems + 1
    #normalidad_residuos = normalidad_residuos + 1
    
  }
  
  # Media marginal igual a cero - Esperanza nula
  
  media_marginal_cero = "SI"
  
  if (qnorm(c(0.025,0.975),0,sd(res)/sqrt(length(res)))[1] > mean(res) | qnorm(c(0.025,0.975),0,sd(res)/sqrt(length(res)))[2] < mean(res)) {
    
    problems = problems + 1
    media_marginal_cero = " NO "
    
  }
  
  
  possible_solutions_table <- matrix(
    c(problems, coef_no_significativos, residuos_incorrelacionados,
      normalidad_residuos, media_marginal_cero),
    ncol = 5
  )
  
  colnames(possible_solutions_table) <- c(
    "Diagnosys_Problems", "Coef_no_sign.",
    "Incorrelacion p.value", "Normalidad p.value", "Esperanza nula"
  )
  
  return(possible_solutions_table)
  
} 

Time_serie_VaR_TVaR_next_point = function(modelo_Arima, VaR_percentile) {
  
  fcast <- forecast(modelo_Arima, h = 1)
  
  VaR_99 = qnorm(VaR_percentile,fcast$mean[1],sqrt(modelo_Arima$sigma2))
  
  muestra_MC = rnorm(1000000,fcast$mean[1], sqrt(modelo_Arima$sigma2))
  
  TVaR_99 = mean(muestra_MC[muestra_MC>VaR_99])
  
  matrix_data = matrix(c(VaR_99, TVaR_99), ncol = 2)
  
  colnames(matrix_data) = c("VaR", "TVaR")
  
  return(matrix_data)
  
  
  
  
}

ARIMA_Summary <- function(time_serie_data) {
  
  count_problems_vector <- c()
  ar_vector <- c()
  difference_vector <- c()
  ma_vector <- c()
  BIC_vector <- c()
  AIC_vector <- c()
  coef_no_significativos <- c()
  residuos_incorrelacionados <- c()
  normalidad_residuos <- c()
  media_marginal_cero <- c()
  var_epsilon = c()
  
  contador <- 1
  
  for (i in 0:5) { # Ar
    for (j in 0:5) { # Ma
      for (k in 0:4) { # difference
        
        if (j == 0 & i == 0) {
          next
        } else {
          tryCatch({
            fit_ <- Arima(time_serie_data, c(i, k, j), include.constant = FALSE)
            count_problems_vector[contador] <- diagnosys_phase(fit_, time_serie_data)[1]
            ar_vector[contador] <- i
            ma_vector[contador] <- j
            difference_vector[contador] <- k
            BIC_vector[contador] <- fit_$bic
            AIC_vector[contador] <- fit_$aic
            coef_no_significativos[contador] <- diagnosys_phase(fit_, time_serie_data)[2]
            residuos_incorrelacionados[contador] <- diagnosys_phase(fit_, time_serie_data)[3]
            normalidad_residuos[contador] <- diagnosys_phase(fit_, time_serie_data)[4]
            media_marginal_cero[contador] <- diagnosys_phase(fit_, time_serie_data)[5]
            var_epsilon[contador] = fit_$sigma2 #error
            
            contador <- contador + 1
          }, error = function(e) {
            # Manejar el error aquí
            #next
            #cat("Error ajustando modelo ARIMA con parámetros:", i, k, j, "\n")
          })
        }
      }
    }
  }
  
  possible_solutions_table <- matrix(
    c(ar_vector, difference_vector, ma_vector, AIC_vector, BIC_vector,
      count_problems_vector, coef_no_significativos, residuos_incorrelacionados,
      normalidad_residuos, media_marginal_cero,var_epsilon),
    ncol = 11
  )
  colnames(possible_solutions_table) <- c(
    "p", "d", "q", "AIC", "BIC", "Diagnosys_Problems", "Coef_no_signi",
    "Incorrelacion", "Normalidad", "Media_marginal_cero", "Var_Epsilon - Sigma^2"
  )
  
  return(possible_solutions_table)
}

###### Datos

Mortality_South_Korea = read.csv2("https://raw.githubusercontent.com/UC3M-student/Proyecto_Series_Temporales/main/Mortality_South_Korea.csv", skip = 1)

personas_en_cartera = matrix(c(902, 1241, 1471,978,1621,566, 882, 1035, 234), ncol = 9)
colnames(personas_en_cartera) = c("45","46","47","48","49","50","51","52", "53" )

suma_asegurada_individual = 500000

###### Time Series - Serie temporal de cada edad

y_45<-ts(as.numeric(Mortality_South_Korea$X45))
y_46<-ts(as.numeric(Mortality_South_Korea$X46))
y_47<-ts(as.numeric(Mortality_South_Korea$X47))
y_48<-ts(as.numeric(Mortality_South_Korea$X48))
y_49<-ts(as.numeric(Mortality_South_Korea$X49))
y_50<-ts(as.numeric(Mortality_South_Korea$X50))
y_51<-ts(as.numeric(Mortality_South_Korea$X51))
y_52<-ts(as.numeric(Mortality_South_Korea$X52))
y_53<-ts(as.numeric(Mortality_South_Korea$X53))


######################################## ESTUDIO EDAD 45 ############################################

#### Idea preeliminar sobre la posible ARIMA objetivo ####
auto.arima(y_45)

#### Comprobacion estacionariedad ####

adf.test(y_45) #NO ESTACIONARIO
kpss.test(y_45) #NO ESTACIONARIO
pp.test(y_45) #NO ESTACIONARIO

# Diferimos un periodo y volvemos a comprobar la estacionariedad en media

adf.test(diff(y_45, lag = 1)) #ESTACIONARIO
kpss.test(diff(y_45, lag = 1)) # NO ESTACIONARIO
pp.test(diff(y_45, lag = 1)) #ESTACIONARIO

# La continua NO estacionariedad del KPSS puede indicar heterocedasticidad. 
# No lo resolvemos porque por un lado esta fuera del alcance de la asignatura y, por otro,
# No queremos diferenciar por segunda vez porque doblaríamos los errores y empeora la estimacion
# CONCLUSION: Hacemos el estudio con una unica diferencia

y_45_diff = diff(y_45)

#### Estudio preeliminar Box-Jenkins sobre autocorrelacion parcial y simple ####
ggAcf(y_45_diff) 
ggPacf(y_45_diff)  

# Se infieren las siguientes posibilidades: ARIMA(1,1,0) - ARIMA(2,1,0) - ARIMA(0,1,1) -ARIMA(3,1,0)

#### Estudiamos los modelo mas verosimiles segun el investigador ####

# Si estacionario 
y45_modelo1 <- Arima(y_45, c(1,1,0))
y45_modelo2 <- Arima(y_45, c(2,1,0))
y45_modelo3 <- Arima(y_45, c(0,1,1))
y45_modelo4 <- Arima(y_45, c(0,1,2))

## Siguiente FASE: mejor AIC/BIC y RMSE
y45_modelo1 <- Arima(y_45, c(1,1,0))

y45_AIC_BIC_models = matrix(c(y45_modelo1$aic, y45_modelo2$aic,y45_modelo3$aic,y45_modelo4$aic,
                     y45_modelo1$bic,y45_modelo2$bic,y45_modelo3$bic,y45_modelo4$bic), ncol = 2, nrow = 4)

rownames(y45_AIC_BIC_models) = c("y45_modelo1","y45_modelo2","y45_modelo3","y45_modelo4")
colnames(y45_AIC_BIC_models) = c("AIC", "BIC")

y45_AIC_BIC_models = y45_AIC_BIC_models[order(y45_AIC_BIC_models[,1]),]
y45_AIC_BIC_models

## Los modelos ARIMA(1,1,0) y ARIMA(2,1,0) pasan a la siguiente fase por mejores propiedades del AIC/BIC

#### FASE SIGNIFICATIVIDAD Y DIAGNOSIS (residuos son ruido blanco) de los candidatos PREDILECTOS 

y45_modelo1_diagnosis = diagnosys_phase(y45_modelo1, y_45) 

#Problemas en la diagnosis: 0

y45_modelo2_diagnosis = diagnosys_phase(y45_modelo2, y_45) 

# Problemas en la diagnosis: 1 ->  Coeficiente NO significativo por el segundo coef. del AR(2)


## Por ultimo, nos queda contrastar la estabilidad en varianza - Homocedasticidad de los residuos

y45_diff = diff(y_45)

# ARIMA(1,1,0)
y45_res_modelo1 <- residuals(y45_modelo1)
y45_df <-data.frame(cbind(y45_res_modelo1^2,lag(y45_diff,-1)))
colnames(y45_df) <- c('residuos','yt_1')
fitaux<-lm(y45_df$residuos~y45_df$yt_1)
df<-2 
bp.statistic<-(length(y45_res_modelo1)-df)*summary(fitaux)$r.squared
bp.pvalue <-1-pchisq(bp.statistic,df) 
bp.pvalue # HOMOCEDASTICO
# ARIMA(2,1,0)
y45_res_modelo2 <- residuals(y45_modelo2)
y45_df <-data.frame(cbind(y45_res_modelo2^2,lag(y45_diff,-1), lag(y45_diff,-2)))
colnames(y45_df) <- c('residuos','yt_1',"yt_2")
fitaux<-lm(y45_df$residuos~y45_df$yt_1 + y45_df$yt_2)
df<-3
bp.statistic<-(length(y45_res_modelo1)-df)*summary(fitaux)$r.squared
bp.pvalue <-1-pchisq(bp.statistic,df) 
bp.pvalue #HOMOCEDASTICO

#¿Añadimos una constante a la ARIMA(1,1,0)? No, ya que no es significativo ni mejora el BIC

######## CONCLUSION Y45 : Nos quedamos con el MODELO ARIMA(1,1,0) #########





######################################## ESTUDIO EDAD 46 ############################################

#### Idea preeliminar sobre la posible ARIMA objetivo ####
auto.arima(y_46)

#### Comprobacion estacionariedad ####

adf.test(y_46) #NO ESTACIONARIO
kpss.test(y_46) #NO ESTACIONARIO
pp.test(y_46) #NO ESTACIONARIO

# Diferimos un periodo y volvemos a comprobar la estacionariedad en media

adf.test(diff(y_46, lag = 1)) #ESTACIONARIO
kpss.test(diff(y_46, lag = 1)) # NO ESTACIONARIO
pp.test(diff(y_46, lag = 1)) #ESTACIONARIO

# La continua NO estacionariedad del KPSS puede indicar heterocedasticidad. 
# No lo resolvemos porque por un lado esta fuera del alcance de la asignatura y, por otro,
# No queremos diferenciar por segunda vez porque doblaríamos los errores y empeora la estimacion
# CONCLUSION: Hacemos el estudio con una unica diferencia

y_46_diff = diff(y_46)

#### Estudio preeliminar Box-Jenkins sobre autocorrelacion parcial y simple ####
ggAcf(y_46_diff) 
ggPacf(y_46_diff)  

# Se infieren las siguientes posibilidades: ARIMA(1,1,0) - ARIMA(2,1,0) - ARIMA(0,1,1) -ARIMA(3,1,0)

#### Estudiamos los modelo mas verosimiles segun el investigador ####

# Si estacionario 
y46_modelo1 <- Arima(y_46, c(1,1,0))
y46_modelo2 <- Arima(y_46, c(2,1,0))
y46_modelo3 <- Arima(y_46, c(0,1,1))
y46_modelo4 <- Arima(y_46, c(0,1,2))

## Siguiente FASE: mejor AIC/BIC y RMSE


y46_AIC_BIC_models = matrix(c(y46_modelo1$aic, y46_modelo2$aic,y46_modelo3$aic,y46_modelo4$aic,
                              y46_modelo1$bic,y46_modelo2$bic,y46_modelo3$bic,y46_modelo4$bic), ncol = 2, nrow = 4)

rownames(y46_AIC_BIC_models) = c("y46_modelo1","y46_modelo2","y46_modelo3","y46_modelo4")
colnames(y46_AIC_BIC_models) = c("AIC", "BIC")

y46_AIC_BIC_models = y46_AIC_BIC_models[order(y46_AIC_BIC_models[,1]),]
y46_AIC_BIC_models

## Los modelos ARIMA(1,1,0) y ARIMA(2,1,0) pasan a la siguiente fase por mejores propiedades del AIC/BIC

#### FASE SIGNIFICATIVIDAD Y DIAGNOSIS (residuos son ruido blanco) de los candidatos PREDILECTOS 

y46_modelo1_diagnosis = diagnosys_phase(y46_modelo1, y_46) 
y46_modelo1_diagnosis

#Problemas en la diagnosis: 0

y46_modelo2_diagnosis = diagnosys_phase(y46_modelo2, y_46) 
y46_modelo2_diagnosis
# Problemas en la diagnosis: 1 ->  Coeficiente NO significativo por el segundo coef. del AR(2)


## Por ultimo, nos queda contrastar la estabilidad en varianza - Homocedasticidad de los residuos

y46_diff = diff(y_46)

# ARIMA(1,1,0)
y46_res_modelo1 <- residuals(y46_modelo1)
y46_df <-data.frame(cbind(y46_res_modelo1^2,lag(y46_diff,-1)))
colnames(y46_df) <- c('residuos','yt_1')
fitaux<-lm(y46_df$residuos~y46_df$yt_1)
df<-2 
bp.statistic<-(length(y46_res_modelo1)-df)*summary(fitaux)$r.squared
bp.pvalue <-1-pchisq(bp.statistic,df) 
bp.pvalue # HOMOCEDASTICO

# ARIMA(2,1,0)
y46_res_modelo2 <- residuals(y46_modelo2)
y46_df <-data.frame(cbind(y46_res_modelo2^2,lag(y46_diff,-1), lag(y46_diff,-2)))
colnames(y46_df) <- c('residuos','yt_1',"yt_2")
fitaux<-lm(y46_df$residuos~y46_df$yt_1 + y46_df$yt_2)
df<-3
bp.statistic<-(length(y46_res_modelo1)-df)*summary(fitaux)$r.squared
bp.pvalue <-1-pchisq(bp.statistic,df) 
bp.pvalue #Homocedastico

#¿Añadimos una constante a la ARIMA(1,1,0)? No, ya que no es significativo ni mejora el BIC

######## CONCLUSION Y46 : Nos quedamos con el MODELO ARIMA(1,1,0) #########


######################################## ESTUDIO EDAD 47 ############################################

#### Idea preeliminar sobre la posible ARIMA objetivo ####
auto.arima(y_47)

#### Comprobacion estacionariedad ####

adf.test(y_47) #NO ESTACIONARIO
kpss.test(y_47) #NO ESTACIONARIO
pp.test(y_47) #NO ESTACIONARIO

# Diferimos un periodo y volvemos a comprobar la estacionariedad en media

adf.test(diff(y_47, lag = 1)) #ESTACIONARIO
kpss.test(diff(y_47, lag = 1)) # NO ESTACIONARIO
pp.test(diff(y_47, lag = 1)) #ESTACIONARIO

# La continua NO estacionariedad del KPSS puede indicar heterocedasticidad. 
# No lo resolvemos porque por un lado esta fuera del alcance de la asignatura y, por otro,
# No queremos diferenciar por segunda vez porque doblaríamos los errores y empeora la estimacion
# CONCLUSION: Hacemos el estudio con una unica diferencia

y_47_diff = diff(y_47)

#### Estudio preeliminar Box-Jenkins sobre autocorrelacion parcial y simple ####
ggAcf(y_47_diff) 
ggPacf(y_47_diff)  

# Se infieren las siguientes posibilidades: ARIMA(1,1,0) - ARIMA(2,1,0) - ARIMA(0,1,1) -ARIMA(0,1,2)

#### Estudiamos los modelo mas verosimiles segun el investigador ####

# Si estacionario 
y47_modelo1 <- Arima(y_47, c(1,1,0))
y47_modelo2 <- Arima(y_47, c(2,1,0))
y47_modelo3 <- Arima(y_47, c(0,1,1))
y47_modelo4 <- Arima(y_47, c(0,1,2))

## Siguiente FASE: mejor AIC/BIC y RMSE


y47_AIC_BIC_models = matrix(c(y47_modelo1$aic, y47_modelo2$aic,y47_modelo3$aic,y47_modelo4$aic,
                              y47_modelo1$bic,y47_modelo2$bic,y47_modelo3$bic,y47_modelo4$bic), ncol = 2, nrow = 4)

rownames(y47_AIC_BIC_models) = c("y47_modelo1","y47_modelo2","y47_modelo3","y47_modelo4")
colnames(y47_AIC_BIC_models) = c("AIC", "BIC")

y47_AIC_BIC_models = y47_AIC_BIC_models[order(y47_AIC_BIC_models[,1]),]
y47_AIC_BIC_models

## Los modelos ARIMA(1,1,0) y ARIMA(2,1,0) pasan a la siguiente fase por mejores propiedades del AIC/BIC

#### FASE SIGNIFICATIVIDAD Y DIAGNOSIS (residuos son ruido blanco) de los candidatos PREDILECTOS 

y47_modelo1_diagnosis = diagnosys_phase(y47_modelo1, y_47) 
y47_modelo1_diagnosis

#Problemas en la diagnosis: 0

y47_modelo2_diagnosis = diagnosys_phase(y47_modelo2, y_47) 
y47_modelo2_diagnosis
# Problemas en la diagnosis: 1 ->  Coeficiente NO significativo por el segundo coef. del AR(2)


## Por ultimo, nos queda contrastar la estabilidad en varianza - Homocedasticidad de los residuos

y47_diff = diff(y_47)

# ARIMA(1,1,0)
y47_res_modelo1 <- residuals(y47_modelo1)
y47_df <-data.frame(cbind(y47_res_modelo1^2,lag(y47_diff,-1)))
colnames(y47_df) <- c('residuos','yt_1')
fitaux<-lm(y47_df$residuos~y47_df$yt_1)
df<-2 
bp.statistic<-(length(y47_res_modelo1)-df)*summary(fitaux)$r.squared
bp.pvalue <-1-pchisq(bp.statistic,df) 
bp.pvalue # HOMOCEDASTICO

# ARIMA(2,1,0)
y47_res_modelo2 <- residuals(y47_modelo2)
y47_df <-data.frame(cbind(y47_res_modelo2^2,lag(y47_diff,-1), lag(y47_diff,-2)))
colnames(y47_df) <- c('residuos','yt_1',"yt_2")
fitaux<-lm(y47_df$residuos~y47_df$yt_1 + y47_df$yt_2)
df<-3
bp.statistic<-(length(y47_res_modelo1)-df)*summary(fitaux)$r.squared
bp.pvalue <-1-pchisq(bp.statistic,df) 
bp.pvalue #Heterocedastico

#¿Añadimos una constante a la ARIMA(1,1,0)? No, ya que no es significativo ni mejora el BIC

######## CONCLUSION Y47 : Nos quedamos con el MODELO ARIMA(1,1,0) #########




######################################## ESTUDIO EDAD 48 ############################################

#### Idea preeliminar sobre la posible ARIMA objetivo ####
auto.arima(y_48)

#### Comprobacion estacionariedad ####

adf.test(y_48) #NO ESTACIONARIO
kpss.test(y_48) #NO ESTACIONARIO
pp.test(y_48) #NO ESTACIONARIO

# Diferimos un periodo y volvemos a comprobar la estacionariedad en media

adf.test(diff(y_48, lag = 1)) #ESTACIONARIO
kpss.test(diff(y_48, lag = 1)) # NO ESTACIONARIO
pp.test(diff(y_48, lag = 1)) #ESTACIONARIO

# La continua NO estacionariedad del KPSS puede indicar heterocedasticidad. 
# No lo resolvemos porque por un lado esta fuera del alcance de la asignatura y, por otro,
# No queremos diferenciar por segunda vez porque doblaríamos los errores y empeora la estimacion
# CONCLUSION: Hacemos el estudio con una unica diferencia

y_48_diff = diff(y_48)

#### Estudio preeliminar Box-Jenkins sobre autocorrelacion parcial y simple ####
ggAcf(y_48_diff) 
ggPacf(y_48_diff)  

# Se infieren las siguientes posibilidades: ARIMA(1,1,0) - ARIMA(1,1,1) - ARIMA(1,1,2) -ARIMA(2,1,1)

#### Estudiamos los modelo mas verosimiles segun el investigador ####

# Si estacionario 
y48_modelo1 <- Arima(y_48, c(1,1,0))
y48_modelo2 <- Arima(y_48, c(1,1,1))
y48_modelo3 <- Arima(y_48, c(1,1,2))
y48_modelo4 <- Arima(y_48, c(2,1,1))

## Siguiente FASE: mejor AIC/BIC y RMSE


y48_AIC_BIC_models = matrix(c(y48_modelo1$aic, y48_modelo2$aic,y48_modelo3$aic,y48_modelo4$aic,
                              y48_modelo1$bic,y48_modelo2$bic,y48_modelo3$bic,y48_modelo4$bic), ncol = 2, nrow = 4)

rownames(y48_AIC_BIC_models) = c("y48_modelo1","y48_modelo2","y48_modelo3","y48_modelo4")
colnames(y48_AIC_BIC_models) = c("AIC", "BIC")

y48_AIC_BIC_models = y48_AIC_BIC_models[order(y48_AIC_BIC_models[,2]),]
y48_AIC_BIC_models

## Los modelos ARIMA(1,1,0) y ARIMA(1,1,1) pasan a la siguiente fase por mejores propiedades del AIC/BIC

#### FASE SIGNIFICATIVIDAD Y DIAGNOSIS (residuos son ruido blanco) de los candidatos PREDILECTOS 

y48_modelo1_diagnosis = diagnosys_phase(y48_modelo1, y_48) 
y48_modelo1_diagnosis

#Problemas en la diagnosis: 0

y48_modelo2_diagnosis = diagnosys_phase(y48_modelo2, y_48) 
y48_modelo2_diagnosis
# Problemas en la diagnosis: 0 


## Por ultimo, nos queda contrastar la estabilidad en varianza - Homocedasticidad de los residuos

y48_diff = diff(y_48)

# ARIMA(1,1,0)
y48_res_modelo1 <- residuals(y48_modelo1)
y48_df <-data.frame(cbind(y48_res_modelo1^2,lag(y48_diff,-1)))
colnames(y48_df) <- c('residuos','yt_1')
fitaux<-lm(y48_df$residuos~y48_df$yt_1)
df<-2 
bp.statistic<-(length(y48_res_modelo1)-df)*summary(fitaux)$r.squared
bp.pvalue <-1-pchisq(bp.statistic,df) 
bp.pvalue # HOMOCEDASTICO

# ARIMA(1,1,1)
y48_res_modelo2 <- residuals(y48_modelo2)
y48_df <-data.frame(cbind(y48_res_modelo2^2,lag(y48_diff,-1), lag(y48_diff,-2), lag(y48_res_modelo2, -1)))
colnames(y48_df) <- c('residuos','yt_1',"yt_2","et_1")
fitaux<-lm(y48_df$residuos~y48_df$yt_1 + y48_df$yt_2 + y48_df$et_1)
df<-4
bp.statistic<-(length(y48_res_modelo1)-df)*summary(fitaux)$r.squared
bp.pvalue <-1-pchisq(bp.statistic,df) 
bp.pvalue #Homocedastico

# Si similitud en el estudio de los dos modelos hace que nos quedemos con 
# ARIMA(1,1,0) por las razones de mejor BIC y por las autocorrelaciones
# simples y parciales que se asemejan mas a esta

#¿Añadimos una constante a la ARIMA(1,1,0)? No, ya que no es significativo ni mejora el BIC

######## CONCLUSION Y48 : Nos quedamos con el MODELO ARIMA(1,1,0) #########



######################################## ESTUDIO EDAD 49 ############################################

#### Idea preeliminar sobre la posible ARIMA objetivo ####
auto.arima(y_49)

#### Comprobacion estacionariedad ####

adf.test(y_49) #NO ESTACIONARIO
kpss.test(y_49) #NO ESTACIONARIO
pp.test(y_49) #NO ESTACIONARIO

# Diferimos un periodo y volvemos a comprobar la estacionariedad en media

adf.test(diff(y_49, lag = 1)) #ESTACIONARIO
kpss.test(diff(y_49, lag = 1)) # NO ESTACIONARIO
pp.test(diff(y_49, lag = 1)) #ESTACIONARIO

# La continua NO estacionariedad del KPSS puede indicar heterocedasticidad. 
# No lo resolvemos porque por un lado esta fuera del alcance de la asignatura y, por otro,
# No queremos diferenciar por segunda vez porque doblaríamos los errores y empeora la estimacion
# CONCLUSION: Hacemos el estudio con una unica diferencia

y_49_diff = diff(y_49)

#### Estudio preeliminar Box-Jenkins sobre autocorrelacion parcial y simple ####
ggAcf(y_49_diff) 
ggPacf(y_49_diff)  

# Se infieren las siguientes posibilidades: ARIMA(1,1,0) - ARIMA(1,1,1) - ARIMA(1,1,2) -ARIMA(2,1,1)

#### Estudiamos los modelo mas verosimiles segun el investigador ####

# Si estacionario 
y49_modelo1 <- Arima(y_49, c(1,1,0))
y49_modelo2 <- Arima(y_49, c(1,1,1))
y49_modelo3 <- Arima(y_49, c(1,1,2))
y49_modelo4 <- Arima(y_49, c(2,1,1))

## Siguiente FASE: mejor AIC/BIC y RMSE


y49_AIC_BIC_models = matrix(c(y49_modelo1$aic, y49_modelo2$aic,y49_modelo3$aic,y49_modelo4$aic,
                              y49_modelo1$bic,y49_modelo2$bic,y49_modelo3$bic,y49_modelo4$bic), ncol = 2, nrow = 4)

rownames(y49_AIC_BIC_models) = c("y49_modelo1","y49_modelo2","y49_modelo3","y49_modelo4")
colnames(y49_AIC_BIC_models) = c("AIC", "BIC")

y49_AIC_BIC_models = y49_AIC_BIC_models[order(y49_AIC_BIC_models[,2]),]
y49_AIC_BIC_models

## Los modelos ARIMA(1,1,0) y ARIMA(1,1,2) pasan a la siguiente fase por mejores propiedades del AIC/BIC

#### FASE SIGNIFICATIVIDAD Y DIAGNOSIS (residuos son ruido blanco) de los candidatos PREDILECTOS 

y49_modelo1_diagnosis = diagnosys_phase(y49_modelo1, y_49) 
y49_modelo1_diagnosis

#Problemas en la diagnosis: 0

y49_modelo2_diagnosis = diagnosys_phase(y49_modelo3, y_49) 
y49_modelo2_diagnosis
# Problemas en la diagnosis: 0 


## Por ultimo, nos queda contrastar la estabilidad en varianza - Homocedasticidad de los residuos

y49_diff = diff(y_49)

# ARIMA(1,1,0)
y49_res_modelo1 <- residuals(y49_modelo1)
y49_df <-data.frame(cbind(y49_res_modelo1^2,lag(y49_diff,-1)))
colnames(y49_df) <- c('residuos','yt_1')
fitaux<-lm(y49_df$residuos~y49_df$yt_1)
df<-2 
bp.statistic<-(length(y49_res_modelo1)-df)*summary(fitaux)$r.squared
bp.pvalue <-1-pchisq(bp.statistic,df) 
bp.pvalue # HOMOCEDASTICO

# ARIMA(1,1,2)
y49_res_modelo2 <- residuals(y49_modelo2)
y49_df <-data.frame(cbind(y49_res_modelo2^2,lag(y49_diff,-1), lag(y49_diff,-2), lag(y49_res_modelo2, -1), lag(y49_res_modelo2, -2)))
colnames(y49_df) <- c('residuos','yt_1',"yt_2","et_1", "et_2")
fitaux<-lm(y49_df$residuos~y49_df$yt_1 + y49_df$yt_2 + y49_df$et_1 + y49_df$et_2)
df<-4
bp.statistic<-(length(y49_res_modelo1)-df)*summary(fitaux)$r.squared
bp.pvalue <-1-pchisq(bp.statistic,df) 
bp.pvalue #HOMOCEDASTICO

#¿Añadimos una constante a la ARIMA(1,1,0)? No, ya que no es significativo ni mejora el BIC

######## CONCLUSION Y49 : Nos quedamos con el MODELO ARIMA(1,1,0) #########



######################################## ESTUDIO EDAD 50 ############################################

#### Idea preeliminar sobre la posible ARIMA objetivo ####
auto.arima(y_50)

#### Comprobacion estacionariedad ####

adf.test(y_50) #NO ESTACIONARIO
kpss.test(y_50) #NO ESTACIONARIO
pp.test(y_50) #NO ESTACIONARIO

# Diferimos un periodo y volvemos a comprobar la estacionariedad en media

adf.test(diff(y_50, lag = 1)) #ESTACIONARIO
kpss.test(diff(y_50, lag = 1)) # NO ESTACIONARIO
pp.test(diff(y_50, lag = 1)) #ESTACIONARIO

# La continua NO estacionariedad del KPSS puede indicar heterocedasticidad. 
# No lo resolvemos porque por un lado esta fuera del alcance de la asignatura y, por otro,
# No queremos diferenciar por segunda vez porque doblaríamos los errores y empeora la estimacion
# CONCLUSION: Hacemos el estudio con una unica diferencia

y_50_diff = diff(y_50)

#### Estudio preeliminar Box-Jenkins sobre autocorrelacion parcial y simple ####
ggAcf(y_50_diff) 
ggPacf(y_50_diff)  

# Se infieren las siguientes posibilidades: ARIMA(1,1,0) - ARIMA(1,1,1) - ARIMA(1,1,2) -ARIMA(2,1,1)

#### Estudiamos los modelo mas verosimiles segun el investigador ####

# Si estacionario 
y50_modelo1 <- Arima(y_50, c(1,1,0))
y50_modelo2 <- Arima(y_50, c(1,1,1))
y50_modelo3 <- Arima(y_50, c(1,1,2))
y50_modelo4 <- Arima(y_50, c(2,1,1))

## Siguiente FASE: mejor AIC/BIC y RMSE


y50_AIC_BIC_models = matrix(c(y50_modelo1$aic, y50_modelo2$aic,y50_modelo3$aic,y50_modelo4$aic,
                              y50_modelo1$bic,y50_modelo2$bic,y50_modelo3$bic,y50_modelo4$bic), ncol = 2, nrow = 4)

rownames(y50_AIC_BIC_models) = c("y50_modelo1","y50_modelo2","y50_modelo3","y50_modelo4")
colnames(y50_AIC_BIC_models) = c("AIC", "BIC")

y50_AIC_BIC_models = y50_AIC_BIC_models[order(y50_AIC_BIC_models[,2]),]
y50_AIC_BIC_models

## Los modelos ARIMA(1,1,0) y ARIMA(1,1,2) pasan a la siguiente fase por mejores propiedades del AIC/BIC

#### FASE SIGNIFICATIVIDAD Y DIAGNOSIS (residuos son ruido blanco) de los candidatos PREDILECTOS 

y50_modelo1_diagnosis = diagnosys_phase(y50_modelo1, y_50) 
y50_modelo1_diagnosis

#Problemas en la diagnosis: 0

y50_modelo2_diagnosis = diagnosys_phase(y50_modelo2, y_50) 
y50_modelo2_diagnosis
# Problemas en la diagnosis: 1


## Por ultimo, nos queda contrastar la estabilidad en varianza - Homocedasticidad de los residuos

y50_diff = diff(y_50)

# ARIMA(1,1,0)
y50_res_modelo1 <- residuals(y50_modelo1)
y50_df <-data.frame(cbind(y50_res_modelo1^2,lag(y50_diff,-1)))
colnames(y50_df) <- c('residuos','yt_1')
fitaux<-lm(y50_df$residuos~y50_df$yt_1)
df<-2 
bp.statistic<-(length(y50_res_modelo1)-df)*summary(fitaux)$r.squared
bp.pvalue <-1-pchisq(bp.statistic,df) 
bp.pvalue # HOMOCEDASTICO

# ARIMA(1,1,1)
y50_res_modelo2 <- residuals(y50_modelo2)
y50_df <-data.frame(cbind(y50_res_modelo2^2,lag(y50_diff,-1), lag(y50_diff,-2), lag(y50_res_modelo2, -1)))
colnames(y50_df) <- c('residuos','yt_1',"yt_2","et_1")
fitaux<-lm(y50_df$residuos~y50_df$yt_1 + y50_df$yt_2 + y50_df$et_1)
df<-3
bp.statistic<-(length(y50_res_modelo1)-df)*summary(fitaux)$r.squared
bp.pvalue <-1-pchisq(bp.statistic,df) 
bp.pvalue #Heterocedastico

#¿Añadimos una constante a la ARIMA(1,1,0)? No, ya que no es significativo ni mejora el BIC

######## CONCLUSION Y50 : Nos quedamos con el MODELO ARIMA(1,1,0) #########


######################################## ESTUDIO EDAD 51 ############################################

#### Idea preeliminar sobre la posible ARIMA objetivo ####
auto.arima(y_51)

#### Comprobacion estacionariedad ####

adf.test(y_51) #NO ESTACIONARIO
kpss.test(y_51) #NO ESTACIONARIO
pp.test(y_51) #NO ESTACIONARIO

# Diferimos un periodo y volvemos a comprobar la estacionariedad en media

adf.test(diff(y_51, lag = 1)) #ESTACIONARIO
kpss.test(diff(y_51, lag = 1)) # NO ESTACIONARIO
pp.test(diff(y_51, lag = 1)) #ESTACIONARIO

# La continua NO estacionariedad del KPSS puede indicar heterocedasticidad. 
# No lo resolvemos porque por un lado esta fuera del alcance de la asignatura y, por otro,
# No queremos diferenciar por segunda vez porque doblaríamos los errores y empeora la estimacion
# CONCLUSION: Hacemos el estudio con una unica diferencia

y_51_diff = diff(y_51)

#### Estudio preeliminar Box-Jenkins sobre autocorrelacion parcial y simple ####
ggAcf(y_51_diff) 
ggPacf(y_51_diff)  

# Se infieren las siguientes posibilidades: ARIMA(1,1,0) - ARIMA(1,1,1) - ARIMA(1,1,2) -ARIMA(2,1,1)

#### Estudiamos los modelo mas verosimiles segun el investigador ####

# Si estacionario 
y51_modelo1 <- Arima(y_51, c(1,1,0))
y51_modelo2 <- Arima(y_51, c(1,1,1))
y51_modelo3 <- Arima(y_51, c(1,1,2))
y51_modelo4 <- Arima(y_51, c(2,1,1))

## Siguiente FASE: mejor AIC/BIC y RMSE


y51_AIC_BIC_models = matrix(c(y51_modelo1$aic, y51_modelo2$aic,y51_modelo3$aic,y51_modelo4$aic,
                              y51_modelo1$bic,y51_modelo2$bic,y51_modelo3$bic,y51_modelo4$bic), ncol = 2, nrow = 4)

rownames(y51_AIC_BIC_models) = c("y51_modelo1","y51_modelo2","y51_modelo3","y51_modelo4")
colnames(y51_AIC_BIC_models) = c("AIC", "BIC")

y51_AIC_BIC_models = y51_AIC_BIC_models[order(y51_AIC_BIC_models[,2]),]
y51_AIC_BIC_models

## Los modelos ARIMA(1,1,0) y ARIMA(1,1,2) pasan a la siguiente fase por mejores propiedades del AIC/BIC

#### FASE SIGNIFICATIVIDAD Y DIAGNOSIS (residuos son ruido blanco) de los candidatos PREDILECTOS 

y51_modelo1_diagnosis = diagnosys_phase(y51_modelo1, y_51) 
y51_modelo1_diagnosis

#Problemas en la diagnosis: 0

y51_modelo2_diagnosis = diagnosys_phase(y51_modelo2, y_51) 
y51_modelo2_diagnosis
# Problemas en la diagnosis: 1


## Por ultimo, nos queda contrastar la estabilidad en varianza - Homocedasticidad de los residuos

y51_diff = diff(y_51)

# ARIMA(1,1,0)
y51_res_modelo1 <- residuals(y51_modelo1)
y51_df <-data.frame(cbind(y51_res_modelo1^2,lag(y51_diff,-1)))
colnames(y51_df) <- c('residuos','yt_1')
fitaux<-lm(y51_df$residuos~y51_df$yt_1)
df<-2 
bp.statistic<-(length(y51_res_modelo1)-df)*summary(fitaux)$r.squared
bp.pvalue <-1-pchisq(bp.statistic,df) 
bp.pvalue # HOMOCEDASTICO

# ARIMA(1,1,1)
y51_res_modelo2 <- residuals(y51_modelo2)
y51_df <-data.frame(cbind(y51_res_modelo2^2,lag(y51_diff,-1), lag(y51_diff,-2), lag(y51_res_modelo2, -1)))
colnames(y51_df) <- c('residuos','yt_1',"yt_2","et_1")
fitaux<-lm(y51_df$residuos~y51_df$yt_1 + y51_df$yt_2 + y51_df$et_1)
df<-3
bp.statistic<-(length(y51_res_modelo1)-df)*summary(fitaux)$r.squared
bp.pvalue <-1-pchisq(bp.statistic,df) 
bp.pvalue #Heterocedastico

#¿Añadimos una constante a la ARIMA(1,1,0)? No, ya que no es significativo ni mejora el BIC

######## CONCLUSION Y51 : Nos quedamos con el MODELO ARIMA(1,1,0) #########



######################################## ESTUDIO EDAD 52 ############################################

#### Idea preeliminar sobre la posible ARIMA objetivo ####
auto.arima(y_52)

#### Comprobacion estacionariedad ####

adf.test(y_52) #NO ESTACIONARIO
kpss.test(y_52) #NO ESTACIONARIO
pp.test(y_52) #NO ESTACIONARIO

# Diferimos un periodo y volvemos a comprobar la estacionariedad en media

adf.test(diff(y_52, lag = 1)) #ESTACIONARIO
kpss.test(diff(y_52, lag = 1)) # NO ESTACIONARIO
pp.test(diff(y_52, lag = 1)) #ESTACIONARIO

# La continua NO estacionariedad del KPSS puede indicar heterocedasticidad. 
# No lo resolvemos porque por un lado esta fuera del alcance de la asignatura y, por otro,
# No queremos diferenciar por segunda vez porque doblaríamos los errores y empeora la estimacion
# CONCLUSION: Hacemos el estudio con una unica diferencia

y_52_diff = diff(y_52)

#### Estudio preeliminar Box-Jenkins sobre autocorrelacion parcial y simple ####
ggAcf(y_52_diff) 
ggPacf(y_52_diff)  

# Se infieren las siguientes posibilidades: ARIMA(1,1,0) - ARIMA(1,1,1) - ARIMA(1,1,2) -ARIMA(2,1,1)

#### Estudiamos los modelo mas verosimiles segun el investigador ####

# Si estacionario 
y52_modelo1 <- Arima(y_52, c(1,1,0))
y52_modelo2 <- Arima(y_52, c(1,1,1))
y52_modelo3 <- Arima(y_52, c(1,1,2))
y52_modelo4 <- Arima(y_52, c(2,1,1))

## Siguiente FASE: mejor AIC/BIC y RMSE


y52_AIC_BIC_models = matrix(c(y52_modelo1$aic, y52_modelo2$aic,y52_modelo3$aic,y52_modelo4$aic,
                              y52_modelo1$bic,y52_modelo2$bic,y52_modelo3$bic,y52_modelo4$bic), ncol = 2, nrow = 4)

rownames(y52_AIC_BIC_models) = c("y52_modelo1","y52_modelo2","y52_modelo3","y52_modelo4")
colnames(y52_AIC_BIC_models) = c("AIC", "BIC")

y52_AIC_BIC_models = y52_AIC_BIC_models[order(y52_AIC_BIC_models[,2]),]
y52_AIC_BIC_models

## Los modelos ARIMA(1,1,0) y ARIMA(1,1,2) pasan a la siguiente fase por mejores propiedades del AIC/BIC

#### FASE SIGNIFICATIVIDAD Y DIAGNOSIS (residuos son ruido blanco) de los candidatos PREDILECTOS 

y52_modelo1_diagnosis = diagnosys_phase(y52_modelo1, y_52) 
y52_modelo1_diagnosis

#Problemas en la diagnosis: 0

y52_modelo2_diagnosis = diagnosys_phase(y52_modelo2, y_52) 
y52_modelo2_diagnosis
# Problemas en la diagnosis: 1


## Por ultimo, nos queda contrastar la estabilidad en varianza - Homocedasticidad de los residuos

y52_diff = diff(y_52)

# ARIMA(1,1,0)
y52_res_modelo1 <- residuals(y52_modelo1)
y52_df <-data.frame(cbind(y52_res_modelo1^2,lag(y52_diff,-1)))
colnames(y52_df) <- c('residuos','yt_1')
fitaux<-lm(y52_df$residuos~y52_df$yt_1)
df<-2 
bp.statistic<-(length(y52_res_modelo1)-df)*summary(fitaux)$r.squared
bp.pvalue <-1-pchisq(bp.statistic,df) 
bp.pvalue # HOMOCEDASTICO

# ARIMA(1,1,1)
y52_res_modelo2 <- residuals(y52_modelo2)
y52_df <-data.frame(cbind(y52_res_modelo2^2,lag(y52_diff,-1), lag(y52_diff,-2), lag(y52_res_modelo2, -1)))
colnames(y52_df) <- c('residuos','yt_1',"yt_2","et_1")
fitaux<-lm(y52_df$residuos~y52_df$yt_1 + y52_df$yt_2 + y52_df$et_1)
df<-3
bp.statistic<-(length(y52_res_modelo1)-df)*summary(fitaux)$r.squared
bp.pvalue <-1-pchisq(bp.statistic,df) 
bp.pvalue #Homocedastico

#¿Añadimos una constante a la ARIMA(1,1,0)? No, ya que no es significativo ni mejora el BIC

######## CONCLUSION Y52 : Nos quedamos con el MODELO ARIMA(1,1,0) #########


######################################## ESTUDIO EDAD 53 ############################################

#### Idea preeliminar sobre la posible ARIMA objetivo ####
auto.arima(y_53)

#### Comprobacion estacionariedad ####

adf.test(y_53) #NO ESTACIONARIO
kpss.test(y_53) #NO ESTACIONARIO
pp.test(y_53) #NO ESTACIONARIO

# Diferimos un periodo y volvemos a comprobar la estacionariedad en media

adf.test(diff(y_53, lag = 1)) #ESTACIONARIO
kpss.test(diff(y_53, lag = 1)) # NO ESTACIONARIO
pp.test(diff(y_53, lag = 1)) #ESTACIONARIO

# La continua NO estacionariedad del KPSS puede indicar heterocedasticidad. 
# No lo resolvemos porque por un lado esta fuera del alcance de la asignatura y, por otro,
# No queremos diferenciar por segunda vez porque doblaríamos los errores y empeora la estimacion
# CONCLUSION: Hacemos el estudio con una unica diferencia

y_53_diff = diff(y_53)

#### Estudio preeliminar Box-Jenkins sobre autocorrelacion parcial y simple ####
ggAcf(y_53_diff) 
ggPacf(y_53_diff)  

# Se infieren las siguientes posibilidades: ARIMA(1,1,0) - ARIMA(1,1,1) - ARIMA(1,1,2) -ARIMA(2,1,1)

#### Estudiamos los modelo mas verosimiles segun el investigador ####

# Si estacionario 
y53_modelo1 <- Arima(y_53, c(1,1,0))
y53_modelo2 <- Arima(y_53, c(1,1,1))
y53_modelo3 <- Arima(y_53, c(1,1,2))
y53_modelo4 <- Arima(y_53, c(2,1,1))

## Siguiente FASE: mejor AIC/BIC y RMSE


y53_AIC_BIC_models = matrix(c(y53_modelo1$aic, y53_modelo2$aic,y53_modelo3$aic,y53_modelo4$aic,
                              y53_modelo1$bic,y53_modelo2$bic,y53_modelo3$bic,y53_modelo4$bic), ncol = 2, nrow = 4)

rownames(y53_AIC_BIC_models) = c("y53_modelo1","y53_modelo2","y53_modelo3","y53_modelo4")
colnames(y53_AIC_BIC_models) = c("AIC", "BIC")

y53_AIC_BIC_models = y53_AIC_BIC_models[order(y53_AIC_BIC_models[,2]),]
y53_AIC_BIC_models

## Los modelos ARIMA(1,1,0) y ARIMA(1,1,2) pasan a la siguiente fase por mejores propiedades del AIC/BIC

#### FASE SIGNIFICATIVIDAD Y DIAGNOSIS (residuos son ruido blanco) de los candidatos PREDILECTOS 

y53_modelo1_diagnosis = diagnosys_phase(y53_modelo1, y_53) 
y53_modelo1_diagnosis

#Problemas en la diagnosis: 0

y53_modelo2_diagnosis = diagnosys_phase(y53_modelo1, y_53) 
y53_modelo2_diagnosis
# Problemas en la diagnosis: 1


## Por ultimo, nos queda contrastar la estabilidad en varianza - Homocedasticidad de los residuos

y53_diff = diff(y_53)

# ARIMA(1,1,0)
y53_res_modelo1 <- residuals(y53_modelo1)
y53_df <-data.frame(cbind(y53_res_modelo1^2,lag(y53_diff,-1)))
colnames(y53_df) <- c('residuos','yt_1')
fitaux<-lm(y53_df$residuos~y53_df$yt_1)
df<-2 
bp.statistic<-(length(y53_res_modelo1)-df)*summary(fitaux)$r.squared
bp.pvalue <-1-pchisq(bp.statistic,df) 
bp.pvalue # HOMOCEDASTICO

# ARIMA(1,1,1)
y53_res_modelo2 <- residuals(y53_modelo2)
y53_df <-data.frame(cbind(y53_res_modelo2^2,lag(y53_diff,-1), lag(y53_diff,-2), lag(y53_res_modelo2, -1)))
colnames(y53_df) <- c('residuos','yt_1',"yt_2","et_1")
fitaux<-lm(y53_df$residuos~y53_df$yt_1 + y53_df$yt_2 + y53_df$et_1)
df<-3
bp.statistic<-(length(y53_res_modelo1)-df)*summary(fitaux)$r.squared
bp.pvalue <-1-pchisq(bp.statistic,df) 
bp.pvalue #Homocedastico

#¿Añadimos una constante a la ARIMA(1,1,0)? No, ya que no es significativo ni mejora el BIC

######## CONCLUSION Y53 : Nos quedamos con el MODELO ARIMA(1,1,0) #########





############ Mortalidad esperada, VaR99 y TVaR99; Coste Esperado, VaR99 y TVaR99; Capital Economico ############

#### Year 45 ####
y45_fcast <- forecast(y45_modelo1, h = 1)

### Mortalidad esperada, VaR y  TVaR 
y45_esperada = y45_fcast$mean[1]
y45_Var = Time_serie_VaR_TVaR_next_point(y45_modelo1,0.99)[1]
y45_TVar = Time_serie_VaR_TVaR_next_point(y45_modelo1,0.99)[2]

### VALORES ECONOMICOS ###
## Best Estimate(Coste Esperado), VaR, TVaR ######
y45_Best_Estimate = y45_esperada[1]*personas_en_cartera[1]*suma_asegurada_individual
y45_VaR_coste = personas_en_cartera[1]*y45_Var*suma_asegurada_individual
y45_TVaR_coste = personas_en_cartera[1]*y45_TVar*suma_asegurada_individual

### CAPITAL ECONOMICO - SCR(solvency capital requirement) ###
y45_SCR = y45_VaR_coste - y45_Best_Estimate



#### Year 46####
y46_fcast <- forecast(y46_modelo1, h = 1)

### Mortalidad esperada, VaR y  TVaR 
y46_esperada = y46_fcast$mean[1]
y46_Var = Time_serie_VaR_TVaR_next_point(y46_modelo1,0.99)[1]
y46_TVar = Time_serie_VaR_TVaR_next_point(y46_modelo1,0.99)[2]

### VALORES ECONOMICOS ###
## Best Estimate(Coste Esperado), VaR, TVaR ######
y46_Best_Estimate = y46_esperada[1]*personas_en_cartera[2]*suma_asegurada_individual
y46_VaR_coste = personas_en_cartera[2]*y46_Var*suma_asegurada_individual
y46_TVaR_coste = personas_en_cartera[2]*y46_TVar*suma_asegurada_individual

### CAPITAL ECONOMICO - SCR(solvency capital requirement) ###
y46_SCR = y46_VaR_coste - y46_Best_Estimate



#### Year 47 ####
y47_fcast <- forecast(y47_modelo1, h = 1)

### Mortalidad esperada, VaR y  TVaR 
y47_esperada = y45_fcast$mean[1]
y47_Var = Time_serie_VaR_TVaR_next_point(y47_modelo1,0.99)[1]
y47_TVar = Time_serie_VaR_TVaR_next_point(y47_modelo1,0.99)[2]

### VALORES ECONOMICOS ###
## Best Estimate(Coste Esperado), VaR, TVaR ######
y47_Best_Estimate = y47_esperada[1]*personas_en_cartera[3]*suma_asegurada_individual
y47_VaR_coste = personas_en_cartera[3]*y47_Var*suma_asegurada_individual
y47_TVaR_coste = personas_en_cartera[3]*y47_TVar*suma_asegurada_individual

### CAPITAL ECONOMICO - SCR(solvency capital requirement) ###
y47_SCR = y47_VaR_coste - y47_Best_Estimate



#### Year 48 ####
y48_fcast <- forecast(y48_modelo1, h = 1)

### Mortalidad esperada, VaR y  TVaR 
y48_esperada = y48_fcast$mean[1]
y48_Var = Time_serie_VaR_TVaR_next_point(y48_modelo1,0.99)[1]
y48_TVar = Time_serie_VaR_TVaR_next_point(y48_modelo1,0.99)[2]

### VALORES ECONOMICOS ###
## Best Estimate(Coste Esperado), VaR, TVaR ######
y48_Best_Estimate = y48_esperada[1]*personas_en_cartera[4]*suma_asegurada_individual
y48_VaR_coste = personas_en_cartera[4]*y48_Var*suma_asegurada_individual
y48_TVaR_coste = personas_en_cartera[4]*y48_TVar*suma_asegurada_individual

### CAPITAL ECONOMICO - SCR(solvency capital requirement) ###
y48_SCR = y48_VaR_coste - y48_Best_Estimate



#### Year 49 ####
y49_fcast <- forecast(y49_modelo1, h = 1)

### Mortalidad esperada, VaR y  TVaR 
y49_esperada = y49_fcast$mean[1]
y49_Var = Time_serie_VaR_TVaR_next_point(y49_modelo1,0.99)[1]
y49_TVar = Time_serie_VaR_TVaR_next_point(y49_modelo1,0.99)[2]

### VALORES ECONOMICOS ###
## Best Estimate(Coste Esperado), VaR, TVaR ######
y49_Best_Estimate = y49_esperada[1]*personas_en_cartera[5]*suma_asegurada_individual
y49_VaR_coste = personas_en_cartera[5]*y49_Var*suma_asegurada_individual
y49_TVaR_coste = personas_en_cartera[5]*y49_TVar*suma_asegurada_individual

### CAPITAL ECONOMICO - SCR(solvency capital requirement) ###
y49_SCR = y49_VaR_coste - y49_Best_Estimate



#### Year 50 ####
y50_fcast <- forecast(y50_modelo1, h = 1)

### Mortalidad esperada, VaR y  TVaR 
y50_esperada = y50_fcast$mean[1]
y50_Var = Time_serie_VaR_TVaR_next_point(y50_modelo1,0.99)[1]
y50_TVar = Time_serie_VaR_TVaR_next_point(y50_modelo1,0.99)[2]

### VALORES ECONOMICOS ###
## Best Estimate(Coste Esperado), VaR, TVaR ######
y50_Best_Estimate = y50_esperada[1]*personas_en_cartera[6]*suma_asegurada_individual
y50_VaR_coste = personas_en_cartera[6]*y50_Var*suma_asegurada_individual
y50_TVaR_coste = personas_en_cartera[6]*y50_TVar*suma_asegurada_individual

### CAPITAL ECONOMICO - SCR(solvency capital requirement) ###
y50_SCR = y50_VaR_coste - y50_Best_Estimate



#### Year 51 ####
y51_fcast <- forecast(y51_modelo1, h = 1)

### Mortalidad esperada, VaR y  TVaR 
y51_esperada = y51_fcast$mean[1]
y51_Var = Time_serie_VaR_TVaR_next_point(y51_modelo1,0.99)[1]
y51_TVar = Time_serie_VaR_TVaR_next_point(y51_modelo1,0.99)[2]

### VALORES ECONOMICOS ###
## Best Estimate(Coste Esperado), VaR, TVaR ######
y51_Best_Estimate = y51_esperada[1]*personas_en_cartera[7]*suma_asegurada_individual
y51_VaR_coste = personas_en_cartera[7]*y51_Var*suma_asegurada_individual
y51_TVaR_coste = personas_en_cartera[7]*y51_TVar*suma_asegurada_individual

### CAPITAL ECONOMICO - SCR(solvency capital requirement) ###
y51_SCR = y51_VaR_coste - y51_Best_Estimate



#### Year 52 ####
y52_fcast <- forecast(y52_modelo1, h = 1)

### Mortalidad esperada, VaR y  TVaR 
y52_esperada = y52_fcast$mean[1]
y52_Var = Time_serie_VaR_TVaR_next_point(y52_modelo1,0.99)[1]
y52_TVar = Time_serie_VaR_TVaR_next_point(y52_modelo1,0.99)[2]

### VALORES ECONOMICOS ###
## Best Estimate(Coste Esperado), VaR, TVaR ######
y52_Best_Estimate = y52_esperada[1]*personas_en_cartera[8]*suma_asegurada_individual
y52_VaR_coste = personas_en_cartera[8]*y52_Var*suma_asegurada_individual
y52_TVaR_coste = personas_en_cartera[8]*y52_TVar*suma_asegurada_individual

### CAPITAL ECONOMICO - SCR(solvency capital requirement) ###
y52_SCR = y52_VaR_coste - y52_Best_Estimate



#### Year 53 ####
y53_fcast <- forecast(y53_modelo1, h = 1)

### Mortalidad esperada, VaR y  TVaR 
y53_esperada = y53_fcast$mean[1]
y53_Var = Time_serie_VaR_TVaR_next_point(y53_modelo1,0.99)[1]
y53_TVar = Time_serie_VaR_TVaR_next_point(y53_modelo1,0.99)[2]

### VALORES ECONOMICOS ###
## Best Estimate(Coste Esperado), VaR, TVaR ######
y53_Best_Estimate = y53_esperada[1]*personas_en_cartera[9]*suma_asegurada_individual
y53_VaR_coste = personas_en_cartera[9]*y53_Var*suma_asegurada_individual
y53_TVaR_coste = personas_en_cartera[9]*y53_TVar*suma_asegurada_individual

### CAPITAL ECONOMICO - SCR(solvency capital requirement) ###
y53_SCR = y53_VaR_coste - y53_Best_Estimate



# TABLA RESUMEN DE TODOS LOS years #


tabla_resumen_edades = matrix(c("ARIMA(1,1,0)","ARIMA(1,1,0)","ARIMA(1,1,0)","ARIMA(1,1,0)","ARIMA(1,1,0)","ARIMA(1,1,0)","ARIMA(1,1,0)","ARIMA(1,1,0)","ARIMA(1,1,0)",
                                as.numeric(y45_esperada),as.numeric(y46_esperada),as.numeric(y47_esperada),as.numeric(y48_esperada),as.numeric(y49_esperada),as.numeric(y50_esperada),as.numeric(y51_esperada),as.numeric(y52_esperada),as.numeric(y53_esperada),
                                as.numeric(y45_Var),as.numeric(y46_Var),as.numeric(y47_Var),as.numeric(y48_Var),as.numeric(y49_Var),as.numeric(y50_Var),as.numeric(y51_Var),as.numeric(y52_Var),as.numeric(y53_Var),
                                as.numeric(y45_TVar),as.numeric(y46_TVar),as.numeric(y47_TVar),as.numeric(y48_TVar),as.numeric(y49_TVar),as.numeric(y50_TVar),as.numeric(y51_TVar),as.numeric(y52_TVar),as.numeric(y53_TVar),
                                as.numeric(y45_Best_Estimate),as.numeric(y46_Best_Estimate),as.numeric(y47_Best_Estimate),as.numeric(y48_Best_Estimate),as.numeric(y49_Best_Estimate),as.numeric(y50_Best_Estimate),as.numeric(y51_Best_Estimate),as.numeric(y52_Best_Estimate),as.numeric(y53_Best_Estimate),
                                as.numeric(y45_VaR_coste),as.numeric(y46_VaR_coste),as.numeric(y47_VaR_coste),as.numeric(y48_VaR_coste),as.numeric(y49_VaR_coste),as.numeric(y50_VaR_coste),as.numeric(y51_VaR_coste),as.numeric(y52_VaR_coste),as.numeric(y53_VaR_coste),
                                as.numeric(y45_TVaR_coste),as.numeric(y46_TVaR_coste),as.numeric(y47_TVaR_coste),as.numeric(y48_TVaR_coste),as.numeric(y49_TVaR_coste),as.numeric(y50_TVaR_coste),as.numeric(y51_TVaR_coste),as.numeric(y52_TVaR_coste),as.numeric(y53_TVaR_coste),
                                as.numeric(y45_SCR),as.numeric(y46_SCR),as.numeric(y47_SCR),as.numeric(y48_SCR),as.numeric(y49_SCR),as.numeric(y50_SCR),as.numeric(y51_SCR),as.numeric(y52_SCR),as.numeric(y53_SCR)), ncol = 8, nrow = 9, byrow = FALSE)



rownames(tabla_resumen_edades) = c("Año 45","Año 46","Año 47","Año 48","Año 49","Año 50","Año 51","Año 52","Año 53")
colnames(tabla_resumen_edades) = c("ARIMA","Mortalidad Esperada", "VaR Mortalidad ","TVaR Mortalidad", "Coste Esperado", 
                            "VaR Coste", "TVaR Coste", "Capital Economico")



# TABLA RESUMEN MODELO COMPLETO #

Coste_Esperado_total = as.numeric(tabla_resumen_edades[,"Coste Esperado"])
VaR_Coste_total = as.numeric(tabla_resumen_edades[,"VaR Coste"])
TVar_Coste_total = as.numeric(tabla_resumen_edades[,"TVaR Coste"])
Capital_Economico_total = as.numeric(tabla_resumen_edades[,"Capital Economico"])

tabla_resumen_total = matrix(c(sum(Coste_Esperado_total),sum(VaR_Coste_total),sum(TVar_Coste_total), sum(Capital_Economico_total)), ncol = 4)
colnames(tabla_resumen_total) = c("Coste Esperado", "VaR", "TVaR", "Capital ECONOMICO")









