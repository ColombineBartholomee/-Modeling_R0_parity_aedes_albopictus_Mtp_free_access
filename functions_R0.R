# function header have been generated using ChatGPT using the following prompt:
# please generate an roxygen header for the following R function: "#code of the function#"


#' Estimate Extrinsic Incubation Period (EIP) Based on Pathogen, Temperature, and Method
#'
#' Computes the extrinsic incubation period (EIP) for a given arbovirus (DENV, CHIKV, ZIKA) 
#' based on the ambient temperature and a selected method from the literature.
#'
#' @param patho Character. Name of the pathogen. One of `"DENV"`, `"CHIKV"`, or `"ZIKA"`.
#' @param temp Numeric. Ambient temperature (in degrees Celsius).
#' @param method Character. Method or model used to compute EIP. Must match one of the published approaches 
#'   for the selected pathogen.
#'
#' @return A numeric value representing the estimated EIP at the specified temperature using the selected method.
#'
#' @details 
#' - **DENV**:
#'   - `"Caminade2016"`: Used in Caminade et al. 2016 [https://doi.org/10.1073/pnas.1614303114] and Metelman et al. 2019 [https://doi.org/10.1371/journal.pntd.0009153], based on work by Chan & Johansson, 2012 [https://doi.org/10.1371/journal.pone.0050972] (data for Aedes albopictus and Aedes aegypti from literature) and Brady et al. 2014 [https://doi.org/10.1186/1756-3305-7-338] (data for Aedes aegypti and aedes albopictus) 
#'   - `"Benkimoun2021"`: Benkimoun et al. 2021 [https://doi.org/10.1016/j.rinp.2021.104687] (Model for Aedes albopictus)
#'
#' - **CHIKV**:
#'   - `"Dumont2008"`: Dumont et al. 2008 [https://doi.org/10.1016/j.mbs.2008.02.008] (Model for Aedes albopictus)
#'   - `"Dubrulle2009"`: Dubrulle et al. 2009 [https://doi.org/10.1371/journal.pone.0005895](laboratory data for Aedes albopictus)
#'   - `"Christofferson2023"`: Based on review by Christofferson et al. 2023 [https://doi.org/10.3390/pathogens12111368] for Aedes aegypti and Aedes albopictus
#'
#' - **ZIKA**:
#'   - `"Guzetta2016"`: Expressed in Guzzetta et al. 2016 [https://doi.org/10.2807/1560-7917.ES.2016.21.15.30199] for Aedes albopictus using results from Chouin-Carneiro et al. 2016 [10.1371/journal.pntd.0004543] (laboratory data for Aedes ageypti and Aedes albopictus)
#'   - `"Caminade2016"`: Same method as for DENV.
#'   - `"Tesla2018"`: Tesla et al. 2018 [https://doi.org/10.1098/rspb.2018.0795]; data available at [https://datadryad.org/dataset/doi:10.5061/dryad.7hj6q4c] (laboratory data for Aedes aegypti)
#'   - `"Winokur2020"`: Winokur et al. 2020 [https://doi.org/10.1371/journal.pntd.0008047] (laboratory data for Aedes aegypti)
#'
#' @examples
#' compute_N("DENV", 28, "Caminade2016")
#' compute_N("CHIKV", 26, "Christofferson2023")
#' compute_N("ZIKA", 30, "Tesla2018")
#'
#' @export
compute_N <- function(patho="DENV", temp, method="Caminade2016"){
  TMN <- temp
  if (patho=="DENV"){
    if (method=="Caminade2016"){
      ## Caminade 2016 (https://doi.org/10.1073/pnas.1614303114) & Metelman 2019 (10.1371/journal.pntd.0009153)  
      N<-1.03*(4 + exp(5.15-0.123*TMN))
    }
    else if (method=="Benkimoun2021"){
      ## Benkimoun 2021 (10.1016/j.rinp.2021.104687)
      N<-0.11*TMN^2 - 7.13 * TMN+121.17
    }else {
      stop("Method non reconnue'")
    }

  }
  else if (patho == "CHIKV") {  
    #if (method=="Dumont2008"){
      ## Dumont 2008 (https://doi.org/10.1016/j.mbs.2008.02.008): Reunion studies and R0 modelling, no explanation on the value
      #n<-3
    # N<-n
    #}
    #else if (method=="Dubrulle2009"){
    # Dubrulle 2009 (https://doi.org/10.1371/journal.pone.0005895): laboratory data
    #n<-runif(n = 1, min = 2, max = 6)
    # N<-n
    #}
     if(method=="Christofferson2023") {
      ## Data from the meta analysis of  Christofferson in 2023:  10.3390/pathogens12111368 (exponential relationships)
      T_values <- c(19, 22, 28) 
      EIP_values <- c(21, 10, 4)  
      exp_model <- function(T, a, b) {
        a * exp(-b * T)
      }
      fit <- nls(EIP_values ~ exp_model(T_values, a, b), 
                 start = list(a = 20, b = 0.1), 
                 algorithm = "port")
      params <- coef(fit)
      a_opt <- params["a"]
      b_opt <- params["b"]
      N <- exp_model(TMN, a_opt, b_opt)
    } else  {
      stop("Method non reconnue'")
    }
  }
  else if (patho == "ZIKA") {  
    # if (method=="Guzetta2016"){
      ## Guzetta et al 2016 (10.2807/1560-7917.ES.2016.21.15.30199) (Data taken from other articles, where the values are not explicitly presented)
      #n<-runif(n = 1, min = 7, max = 14)
    #N<-n 
    # }
    if (method=="Caminade2016"){
      ## Caminade et al. 2016 : same as DENV (https://doi.org/10.1073/pnas.1614303114)
      N<-1.03*(4 + exp(5.15-0.123*TMN))
    }
    
    else if (method=="Winokur2020"){
      ## Data from the article of Winokur 2020:   10.1371/journal.pntd.0008047 (lab studies with Aedes aegypti)
      # Equation with temperatures
      T_values <- c(21, 26, 30) 
      EIP_values <- c(24.2, 9.6, 5.1)  
      exp_model <- function(T, a, b) {
        a * exp(-b * T)
      }
       fit <- nls(EIP_values ~ exp_model(T_values, a, b), 
        start = list(a = 20, b = 0.1), 
        algorithm = "port")
      params <- coef(fit)
      a_opt <- params["a"]
      b_opt <- params["b"]
      N <- exp_model(TMN, a_opt, b_opt)
    }
    else if (method=="Tesla2018"){
      ## Tesla et al. 2018 :For Aedes ageypti in laboratory : https://doi.org/10.1098/rspb.2018.0795, raw data  : https://datadryad.org/dataset/doi:10.5061/dryad.7hj6q4c#
      data_points <- data.frame(
        T = c(16, 20, 24, 28, 32, 34, 36, 38, 16, 20, 24, 28, 32, 34, 36, 38),
        N = c(0, 18.216,13.28399992, 7.28599998, 4.5, 4.3935, 3.8775, 1.5, 0,0, 15.86400008, 8, 6.96449998, 5.787, 3.556000002, 0)
      )
      model_N <- lm(N ~ poly(T, 2, raw = TRUE), data = data_points)
      predictions <- predict(model_N, ## Predictions avec IC de confiance de 95% a partir de la temperature souhaitee
                             newdata = data.frame(T = TMN), 
                             interval = "confidence", 
                             level = 0.95)
      N <- predictions[, "fit"] 
    } else {
      stop("Method non reconnue'")
    }

  }
  else {
    stop("Pathogene non reconnu'")
  }
  return(N)
}



#' Estimate Vector Competence Based on Pathogen, Temperature, and Method
#'
#' Computes the vector competence (probability of transmission after infection) for a given arbovirus 
#' (DENV, CHIKV, ZIKA), based on temperature and a selected method from the literature.
#'
#' @param patho Character. Name of the pathogen. Must be one of `"DENV"`, `"CHIKV"`, or `"ZIKA"`.
#' @param temp Numeric. Ambient temperature (in degrees Celsius).
#' @param method Character. Method or model used to compute vector competence. Must match one of the 
#'   published approaches for the selected pathogen.
#'
#' @return A numeric value representing the estimated vector competence (`b`) at the given temperature using the specified method.
#'
#' @details 
#' - **DENV**:
#'   - `"Verga-Rua2013"`: Verga-Ruiz et al. 2013 [https://doi.org/10.1371/journal.pone.0059716] (Laboratory data from Aedes albopictus from Montpellier) and Solimini et al. 2018 [https://doi.org/10.1038/s41598-018-34664-5] (modelling data from Aedes albopictus)
#'   - `"Metelman2019"`: Metelmann et al. 2019 [https://doi.org/10.1371/journal.pntd.0009153], based on work by Chan & Johansson, 2012 [https://doi.org/10.1371/journal.pone.0050972] (data for Aedes albopictus and Aedes aegypti from literature) and Brady et al. 2014 [https://doi.org/10.1186/1756-3305-7-338] (data for Aedes aegypti and aedes albopictus)
#'   - `"Caminade2016"`: Caminade et al. 2016 [https://doi.org/10.1073/pnas.1614303114]
#'   - `"Benkimoun2021"`: Benkimoun et al. 2021 [https://doi.org/10.1016/j.rinp.2021.104687] (Model for Aedes albopictus)
#'
#' - **CHIKV**:
#'   - `"Verga-Rua2013"`: Verga-Ruiz et al. 2013 [https://doi.org/10.1371/journal.pone.0059716] (Laboratory data from Aedes albopictus from Montpellier) and Solimini et al. 2018 [https://doi.org/10.1038/s41598-018-34664-5] (modelling data from Aedes albopictus)
#'   - `"method name"`: Based on multiple sources including Delrieu et al. 2023 [https://doi.org/10.1016/j.crpvbd.2023.100139] (systematic reveiw on Aedes aegypti and Aedes albopictus, vector competence depending on temperature), 
#'     Christofferson et al. 2023 [https://doi.org/10.3390/pathogens12111368] (systematic reveiw on Aedes aegypti and Aedes albopictus, vector competence depending on temperature), and Lühken et al. 2024 [https://doi.org/10.1186/s13071-024-06594-x] (laboratory data from German Aedes albopictus)
#'
#' - **ZIKA**:
#'   - `"DiLuca2016"`: Di Luca et al. 2016 [https://doi.org/10.2807/1560-7917.ES.2016.21.18.30223] (laboratory data from Italian Aedes albopictus) and Solimini et al. 2018 [https://doi.org/10.1038/s41598-018-34664-5]
#'   - `"Caminade2016"`: Caminade et al. 2016 (same as for DENV) [https://doi.org/10.1073/pnas.1614303114]
#'   - `"Tesla2018"`: Tesla et al. 2018 [https://doi.org/10.1098/rspb.2018.0795]; data available at [https://datadryad.org/dataset/doi:10.5061/dryad.7hj6q4c] (laboratory data for Aedes aegypti)
#'
#' @examples
#' compute_b("DENV", 28, "Caminade2016")
#' compute_b("CHIKV", 26, "Verga-Rua2013")
#' compute_b("ZIKA", 30, "Tesla2018")
#'
#' @export

compute_b <- function(patho="DENV", temp, method="Benkimoun2021"){
  TMN <- temp
  if (patho=="DENV"){
    if (method=="Verga-Rua2013"){
      ## Verga-Rua 2013 (10.1371/journal.pone.0059716) & Solimini 2018 (10.1038/s41598-018-34664-5) ## Mosquito from Montpellier
      #b<-runif(n = 1, min = 0.56, max =0.67)
      b<-0.61 # Value on the middle of the interval
    }
      # if (method=="Metelman2019"){
      ##  Metelman 2019 (10.1371/journal.pntd.0009153) Data from modelling 
      # b<-0.5
      # }
    #else if (method=="Caminade2016"){
      ## Caminade et al. 2016 : meme que DENV (https://doi.org/10.1073/pnas.1614303114) Data from modelling
      # b<-runif(n = 1, min = 0.1, max = 0.75) ## Ou 0.5
    # b<-0.5
    # }
    else if (method=="Benkimoun2021"){
      ## Benkimoun 2021 (10.1016/j.rinp.2021.104687)
      b<- 0.2593*TMN - 3.2705 -0.0043*TMN^2
      b<-ifelse(b<0,0,b)
    }
    else {
      stop("Method non reconnue'")
    }
  }
  else if (patho == "CHIKV") {
    if (method=="Verga-Rua2013"){
      ## Verga-Rua 2013 (10.1371/journal.pone.0059716) & Solimini 2018 (10.1038/s41598-018-34664-5) ## Mosquito from Montpellier
      #b<-runif(n = 1, min = 0.14, max = 0.84)
    b<-0.49
     }
   else if (method=="method name") {
      ## Varies with temperature: Delrieu 2023 (https://doi.org/10.1016/j.crpvbd.2023.100139), Christofferson (10.3390/pathogens12111368), 
      ## Probable relationship with T2 as for DENV
      ## values from Lühlen 2024 (https://doi.org/10.1186/s13071-024-06594-x), in Germany, with Fluctuating T°
      data_points <- data.frame(
        T = c(15, 18,  21,  24.5),
        b = c(0.325,0.543,  0.586, 0.325)
      )
      model_b <- lm(b ~ poly(T, 2, raw = TRUE), data = data_points)
      predictions <- predict(model_b, 
                             newdata = data.frame(T = TMN), 
                             interval = "confidence", 
                             level = 0.95)
      b <- predictions[, "fit"] 
      b<-ifelse(b<0,0,b)
    }
    else {
      stop("Method non reconnue'")
    }
    
  }
  else if (patho == "ZIKA") {  
    if (method=="DiLuca2016"){
      ## Di Luca 2016  (https://doi.org/10.2807/1560-7917.ES.2016.21.18.30223) & Solimini 2018  (10.1038/s41598-018-34664-5) ## Lab data with Aedes albopictus
      #b<-runif(n = 1, min = 0.19, max = 0.39)
      b<-0.29 
    }
    # else if (method=="Caminade2016"){
      ## Caminade et al. 2016 : meme que DENV (https://doi.org/10.1073/pnas.1614303114) # Data from modelling
    #  b<-0.5
    # }
    else if (method=="Tesla2018"){
      ## Tesla et al. 2018 : For Aedes aegypti in lab : https://doi.org/10.1098/rspb.2018.0795, raw data : https://datadryad.org/dataset/doi:10.5061/dryad.7hj6q4c#
      data_points <- data.frame(
        T = c(16, 20, 24, 28, 32, 34, 36, 38, 16, 20, 24, 28, 32, 34, 36, 38),
        b = c(0, 0.0071, 0.0429, 0.2571, 0.2, 0.1857, 0.117,0.1, 0, 0, 0.0288, 0.2, 0.1286, 0.1857, 0.0556,0)
      )
      model_b <- lm(b ~ poly(T, 2, raw = TRUE), data = data_points)
      predictions <- predict(model_b, 
                             newdata = data.frame(T = TMN), 
                             interval = "confidence", 
                             level = 0.95)
      b <- predictions[, "fit"] 
    }
    else {
      stop("Method non reconnue'")
    }

  }
  else {
    stop("Pathogene non reconnu'")
  }
  
  return(b)
}
  
#' Estimate Probability of Transmission from Human to Vector (c)
#'
#' Computes the probability of transmission from an infected human to a vector (c)
#' based on the pathogen, ambient temperature, and a selected method from the literature.
#'
#' @param patho Character. Name of the pathogen. One of `"DENV"`, `"CHIKV"`, or `"ZIKA"`.
#' @param method Character. Method or source used to compute the probability. Must match one of the available methods for the selected pathogen.
#'
#' @return A numeric value representing the estimated probability of transmission from human to vector (c).
#'
#' @details
#' - **DENV**:
#'   - `"Solimini2018"`: Random value between 0.14 and 0.39 [https://doi.org/10.1038/s41598-018-34664-5] from Lab data on Aedes albopictus of Talbalaghi et al. 2010 [https://doi.org/10.1111/j.1365-2915.2009.00853.x]
#'   - `"Metelman2019"`: Fixed value of 0.31 [https://doi.org/10.1371/journal.pntd.0009153], using lab data on Aedes albopictus from acmeroun [https://doi.org/10.1089/vbz.2009.0005]
#'
#' - **CHIKV**:
#'   - `"Solimini2018"`: Random value between 0.75 and 0.90 [https://doi.org/10.1038/s41598-018-34664-5] from Lab data on Aedes albopictus of Talbalaghi et al. 2010 [https://doi.org/10.1111/j.1365-2915.2009.00853.x]
#'
#' - **ZIKA**:
#'   - `"Solimini2018"`: Random value between 0.03 and 0.17 [https://doi.org/10.1038/s41598-018-34664-5], with Lab data From Di luca 2016 [https://doi.org/10.2807/1560-7917.ES.2016.21.18.30223] (laboratory data from Italian Aedes albopictus) 
#'   - `"Caminade2016"`: Fixed value of 0.033 [https://doi.org/10.1073/pnas.1614303114]  for Aedes albopictus using results from Chouin-Carneiro et al. 2016 [10.1371/journal.pntd.0004543] (laboratory data for Aedes ageypti and Aedes albopictus)
#'
#' @examples
#' compute_c("DENV", "Solimini2018")
#' compute_c("ZIKA", "Caminade2016")
#'
#' @export
compute_c <- function(patho="DENV", method="Metelman2019") {
  if (patho == "DENV") {
    # if (method == "Solimini2018") {
    #  c <- runif(n = 1, min = 0.14, max = 0.39)
    # } 
    if (method == "Metelman2019") { ## On the interval of Solimini 2018
      c <- 0.31
    } else {
      stop("Method non reconnue pour DENV")
    }
  } else if (patho == "CHIKV") {
    if (method == "Solimini2018") {
      # c <- runif(n = 1, min = 0.75, max = 0.90)
      c<-0.82
    } else {
      stop("Method non reconnue pour CHIKV")
    }
  } else if (patho == "ZIKA") {
    # if (method == "Solimini2018") {
      #c <- runif(n = 1, min = 0.03, max = 0.17)
      # c<-0.10
    # } 
    if (method == "Caminade2016") {  ##On the interval of Solimini 2018
      c <- 0.033
    } else {
      stop("Method non reconnue pour ZIKA")
    }
  } else {
    stop("Pathogene non reconnu")
  }
  return(c)
}

#' Estimate Daily Biting Rate (a)
#'
#' Computes the daily biting rate of mosquitoes based on ambient temperature and a selected method from the literature.
#'
#' @param temp Numeric. Ambient temperature (in degrees Celsius).
#' @param method Character. Method or literature source used to compute the biting rate.
#'               Must be one of: `"Caminade2016"`, `"Poletti2011"`, `"Fontenille2017"`, `"Solimini2018"`.
#'
#' @return A numeric value representing the estimated daily biting rate (a).
#'
#' @details
#' - **Caminade2016**: Linear function of temperature, `a = (0.0043 * T + 0.0943) / 2`. [https://doi.org/10.1073/pnas.1614303114], from lab and modelling data of Aedes aegypti of Thailand from Liu-Helmersson et al.2014[https://doi.org/10.1371/journal.pone.0089783], 
#' - **Poletti2011**: Random value between 0.05 and 0.16 [https://doi.org/10.1371/journal.pone.0018860].
#' - **Fontenille2017**: Fixed value of 0.25 [ISBN: 978-90-8686-053-1].
#' - **Solimini2018**: Random value between 0.08 and 0.10 [https://doi.org/10.1038/s41598-018-34664-5].
#'
#' @examples
#' compute_a(27, "Caminade2016")
#' compute_a(30, "Poletti2011")
#'
#' @export
compute_a <- function(temp, method="Caminade2016") { ## Prefer the a=1/gc function because gc comes from a laboratory value
  if (method == "Caminade2016") {
    a <- (0.0043 * temp + 0.0943) / 2
  } #else if (method == "Poletti2011") {
    #a <- runif(n = 1, min = 0.05, max = 0.16)
  #a<-0.11}
  #else if (method == "Fontenille2017") { # Valeur unique 
  #  a <- 0.25
  # } 
#else if (method == "Solimini2018") {
    #a <- runif(n = 1, min = 0.08, max = 0.10)
    #a<-0.09} 
  else {
    stop("Method non reconnue")
  }
  return(a)
}

#' Estimate Human Infectious Period (r)
#'
#' Returns the human infectious period (in days) based on the pathogen, using distributions or fixed values reported in the literature.
#'
#' @param patho Character. Pathogen name. Must be one of: `"DENV"`, `"CHIKV"`, `"ZIKA"`.
#'
#' @return A numeric value representing the estimated infectious period in days.
#'
#' @details
#' - **DENV**: Uniform distribution between 2 and 7 days. Sources: Solimini 2018 [https://doi.org/10.1038/s41598-018-34664-5], WHO.
#' - **CHIKV**: Uniform distribution between 2 and 6 days. Source: Solimini 2018 [https://doi.org/10.1038/s41598-018-34664-5].
#' - **ZIKA**: Uniform distribution between 4 and 7 days. Sources: Guzzetta 2016 [https://doi.org/10.2807/1560-7917.ES.2016.21.15.30199], Solimini 2018  [https://doi.org/10.1038/s41598-018-34664-5].
#'
#' @examples
#' compute_r("DENV")
#' compute_r("ZIKA")
#'
#' @export
compute_r <- function(patho="DENV") {
  if (patho == "DENV") {
    r <- 5 #runif(n = 1, min = 2, max = 7) 
  } else if (patho == "CHIKV") {
    r <- 4 #runif(n = 1, min = 2, max = 6)
  } else if (patho == "ZIKA") {
    r <- 6 #runif(n = 1, min = 4, max = 7)
  } else {
    stop("Pathogene non reconnu")
  }
  return(r)
}

#' Estimate Aedes albopictus human preference
#'
#' Returns human preference rate of Ae. albopictus, using fixed values reported in the literature from European cities.

#'
#' @return A numeric value representing the estimated infectious period in days.
#'
#' @details
#' - between 68.4% and 97.9% in Rome (Italy): Valerio et al. 2010 [10.1089/vbz.2009.0007], Valerio et al. 2008 [PMID: 18693570],
#' - 100% In Barcelona (Spain) but only in 30 mosquitoes : Munoz et al. 2011 [10.1603/ME11016]
#' @examples
#' compute_phi()
#' @export
compute_phi <- function() {
  phi<-0.842 #runif (0.684 and 1) (middle of the interval )
  return(phi)
}


#' Estimate Gonotrophic Cycle Duration
#'
#' Returns the estimated duration (in days) of the gonotrophic cycle based on temperature.
#' The model is fitted from laboratory data on Ae. albopictus, combining results from
#' Delatte et al. (2009) [10.1603/033.046.0105] and Marini et al. (2020) [10.3390/insects11110808].
#'
#' @param temp Numeric. Mean daily temperature in °C.
#' @param return_ci Logical. If TRUE, returns the full prediction with confidence interval (lower, fit, upper).
#'
#' @return Numeric. Estimated gonotrophic cycle duration in days (or full prediction if `return_ci = TRUE`).
#'
#' @examples
#' compute_gc(20)
#' compute_gc(30, return_ci = TRUE)
#'
#' @export
compute_gc <- function(temp, return_ci = FALSE) {
  data_points <- data.frame(
    T = c(20, 25, 30, 35),
    gc = c(8.1, 3.8, 3.3, 4.4)
  )
  
  model_gc <- lm(gc ~ poly(T, 2, raw = TRUE), data = data_points)
  
  predictions <- predict(model_gc, 
                         newdata = data.frame(T = temp), 
                         interval = "confidence", 
                         level = 0.95)
  
  if (return_ci) {
    return(predictions)
  } else {
    return(predictions[, "fit"])
  }
}

#' Compute simulated parity rate using GAM models
#'
#' This function estimates the parity rate (PROP_PP) via simulation from a Generalized Additive Model (GAM),
#' using either the AREA + MONTH or the ZONE + MONTH modeling structure. It returns simulated predictions
#' across multiple iterations, allowing for the computation of uncertainty (e.g., 95% confidence intervals).
#'
#' @param df A data frame containing the columns: `ID_COLLECTE`, `AREA`, `MOIS`, `ZONE`, `NB_TOT`, and `PROP_PP`,
#' which are required for modeling and simulation.
#' @param method A character string indicating the GAM structure to use. Must be one of `"area_mois"` 
#' (default) or `"zone_mois"`.
#' @param nsim An integer specifying the number of simulations to perform (default is 100).
#'
#' @return A data frame grouped by `AREA`, `MOIS`, `ZONE`, and `ID_COLLECTE`, containing:
#' - `sim_id`: Simulation index
#' - `parity_rate_sim`: Simulated parity rate values from the GAM model
#'
#' @details
#' - `"area_mois"`: fits the model `PROP_PP ~ AREA + s(MOIS, k = 5, bs = 'cr')` using binomial family with weights `NB_TOT`.
#'   Example model performance: Adjusted R² = 0.146, Deviance explained = 18.2%. See script: `'05_GAMs_For_Predictions'`.
#'
#' - `"zone_mois"`: fits the model `PROP_PP ~ ZONE + s(MOIS, k = 5, bs = 'cr') + s(MOIS, by = ZONE, k = 5)`.
#'   Example model performance: Adjusted R² = 0.14, Deviance explained = 16.6%. See script: `'05_GAMs_For_Predictions'`.
#'
#' @import dplyr tidyr mgcv
#' @export
compute_parity_rate_ic <- function(df, method = "area_mois", nsim = 100) {
  ## Prepare aggregated new data for prediction
  newdata <- df %>%
    dplyr::select(ID_COLLECTE,AREA, ZONE, MOIS, PROP_PP, NB_TOT)
    #group_by(AREA, MOIS, ZONE) %>%
    #summarise(
    #  parity_rate = mean(PROP_PP, na.rm = TRUE),
    # NB_TOT = round(mean(NB_TOT, na.rm = TRUE), 0),
    # .groups = "drop"
    # )
  
  ## Fit the appropriate model based on the chosen method
  if (method == "area_mois") {
    gam_model <- gam(PROP_PP ~ AREA + s(MOIS, k = 5, bs = 'tp'),
                     data = df,
                     family = binomial,
                     weights = NB_TOT,
                     # method = "REML"
                     gamma=1)
  } else if (method == "zone_mois") {
    gam_model <- gam(PROP_PP ~ ZONE + s(MOIS, k=5, bs='tp')+s(MOIS, by=ZONE, k=5),
                     data = df,
                     family = binomial,
                     weights = NB_TOT,
                     # method = "REML"
                     gamma=1)
  }
  else {
    stop("Invalid method. Please use 'area_mois' or 'zone_mois'.")
  }
  ## Simulate nsim draws from the fitted model using new data
  sim <- simulate(gam_model, data = newdata, nsim = nsim, weights = newdata$NB_TOT)
  colnames(sim) <- paste0("sim_", 1:nsim)
  sim_data <- cbind(newdata, sim)
  
  ## Convert simulations to long format
  sim_long <- sim_data %>%
    pivot_longer(
      cols = starts_with("sim_"),
      names_to = "sim_id",
      values_to = "parity_rate_sim")%>%
    mutate(parity_rate_sim = ifelse(is.nan(parity_rate_sim), 0, parity_rate_sim)) 
  
  ## Compute summary statistics: mean, lower and upper bounds of 95% CI
  #ic_data <- sim_long %>%
  # group_by(ID_COLLECTE,AREA, ZONE, MOIS,  PROP_PP, NB_TOT)%>%
  #group_by( MOIS, ZONE) %>%
  # summarise(
  #  parity_rate_sim = mean(sim_value),
  # parity_rate_low = quantile(sim_value, 0.025),
  # parity_rate_high = quantile(sim_value, 0.975),
  # .groups = "drop"
  #)

  return(sim_long)
}



#' Simulate count predictions for Aedes albopictus with 95% confidence intervals using a GAM
#'
#' This function fits a negative binomial GAM to predict the number of female Aedes albopictus 
#' mosquitoes based on spatial (AREA or ZONE) and temporal (MOIS) covariates, and includes random effects. 
#' Simulations are used to estimate predicted values.
#'
#' @param df A data frame containing at least the columns: AREA, ZONE, MOIS, NB_ALBO_F, and num_session.
#' @param method A character string: either `"zone_mois"` or `"area_mois"` to choose the model structure.
#' @param nsim An integer for the number of simulations to perform. Default is 100.
#'
#' @return  data frame grouped by `AREA`, `MOIS`, `ZONE`, and `ID_COLLECTE`, containing:
#' - `sim_id`: Simulation index
#' - `ma_sim`: Simulated abundance from the GAM model

#'#' @details
#' - `"area_mois"`: fits the model ` NB_ALBO_F ~  AREA +s(MOIS, k=5, bs='tp')+s(num_session, bs='re')`.  
#'   Example: R²(adj) =  0.367   Deviance explained = 50.3%. See script `05_GAMs_For_Predictions`.
#' - `"zone_mois"`: fits the model `NB_ALBO_F ~ ZONE + s(MOIS, k = 5, bs = 'cr') + s(num_session, bs = 're') + s(AREA, bs = 're')`. 
#'   Example:  R²(adj) =  0.354   Deviance explained =   49%. See script `05_GAMs_For_Predictions`.
#' @import dplyr tidyr mgcv MASS
#' @export
compute_ma_counts <- function(df, method = "area_mois", nsim = 100) {
  ## Estimate theta (dispersion parameter) using a simple GLM to provide a starting point
  theta_hat <- MASS::glm.nb(NB_ALBO_F_DIS ~ 1, data = df)$theta
  
  ## Prepare newdata for simulation: average per combination
  newdata <- df %>%
    dplyr::select(ID_COLLECTE,AREA, MOIS, ZONE, num_session,NB_ALBO_F)
  #group_by(AREA, ZONE, MOIS, num_session) %>%
  #summarise(
  #  NB_ALBO_F = round(mean(NB_ALBO_F, na.rm = TRUE), 0),
  #  .groups = "drop"
  # )
  
  ## Fit the appropriate GAM model with random effects and estimated theta
  if (method == "zone_mois") {
    model <- gam(NB_ALBO_F_DIS ~ ZONE + s(MOIS, k = 7, bs = 'cr') + 
                   s(num_session, bs = 're') + s(AREA, bs = 're'),
                 data = df,
                 family = negbin(theta = theta_hat),
                 # method = "REML"
                 gamma=1)
  } else if (method == "area_mois") {
    model <- gam(NB_ALBO_F_DIS ~ AREA + s(MOIS, k = 7, bs = 'tp') + 
                   s(num_session, bs = 're'),
                 data = df,
                 family = negbin(theta = theta_hat),
                 # method = "REML"
                 gamma=1)
  } else {
    stop("Method must be either 'zone_mois' or 'area_mois'")
  }
  ## Perform manual simulations based on predicted means
  mu <- predict(model, newdata = newdata, type = "response")
  theta <- model$family$getTheta()  
  
  sim_matrix <- replicate(nsim, rnbinom(n = length(mu), mu = mu, size = theta))
  sim_df <- as.data.frame(sim_matrix)
  colnames(sim_df) <- paste0("sim_", 1:nsim)
  
  ## Combine with identifiers
  sim_data <- cbind(newdata, sim_df)
  
  ## Reshape simulations to long format
  sim_long <- sim_data %>%
    pivot_longer(
      cols = starts_with("sim_"),
      names_to = "sim_id",
      values_to = "ma_sim"
    )
  return(sim_long)
}

#' Calculate Vectorial Capacity (Garrett-Jones 1964)
#'
#' Computes vectorial capacity based on biting rate, mosquito-to-human ratio,
#' survival probability, and extrinsic incubation period.
#'
#' @param ma Numeric. Mosquito-to-human ratio (m × a).
#' @param a Numeric. Daily biting rate.
#' @param p Numeric. Daily survival probability of the mosquito (between 0 and 1).
#' @param expo Numeric. Human exposure to Ae. albopictus in a day (between 0 and 1).(from Richards et al. 2006 : [10.1093/jmedent/43.3.543])
#' @param N Numeric. Duration of the extrinsic incubation period (in days).
#'
#' @return Numeric. Vectorial capacity.
#'
#' @examples
#' compute_capacity(ma = 2.5, a = 0.3, p = 0.9, N = 10)
#'
#' @export
compute_capacity <- function(ma, a, p,expo, N) {
  if (any(is.na(c(ma, a, p, N)))) {
    return(NA_real_)}
  else if (any(p < 0 | p > 1)) {
    stop("Survival probability 'p' must be in (0, 1].")
  }
  else if (p==0){
    vc<-0 ## to prevent _ log (0)
  }
  else if (p==1){
    #vc<-(expo*ma * a * ((0.99)^N)) / (1-(0.99)) ## Taylor development −log(p)≈1−p(when p --> 1) and to prevent / 0. Approximation of p~0.99
    vc<-(expo*ma * a * ((0.99)^N)) / (1-(0.99)) ## Taylor development −log(p)≈1−p(when p --> 1) and to prevent / 0. Approximation of p~0.99
  }
  else {
    vc <- (expo*ma * a * ((p)^N)) / (-log(p))
  }
  return(vc)
}

#' Calculate Basic Reproduction Number R₀
#'
#' Computes R₀ using vectorial capacity and parameters for transmission probabilities
#' and duration of human infectiousness.
#'
#' @param capacity_vect Numeric. Vectorial capacity (Garrett-Jones formula).
#' @param r Numeric. Duration of human infectiousness (in days).
#' @param b Numeric. Vector competence (probability of transmission from mosquito to human).
#' @param c Numeric. Probability of transmission from human to mosquito.
#'
#' @return Numeric. Basic reproduction number R₀.
#'
#' @examples
#' compute_R0(capacity_vect = 2.5, r = 5, b = 0.3, c = 0.2)
#'
#' @export
compute_R0 <- function(capacity_vect, r, b, c) {
  r0 <- capacity_vect * r * b * c
  return(r0)
}

#' Compute Basic Reproduction Number R₀ (Poletti 2011)
#'
#' Calculates R₀ based on the formulation by Poletti et al. (2011), 
#' separating transmission from vector to human (R0VH) and human to vector (R0HV).
#'
#' @param a Numeric. Daily biting rate.
#' @param b Numeric. Vector competence (probability of transmission from mosquito to human).
#' @param c Numeric. Probability of transmission from human to mosquito.
#' @param r Numeric. Duration of infectiousness in humans.
#' @param N Numeric. Duration of extrinsic incubation period in days.
#' @param p Numeric. Daily mosquito survival probability.
#' @param ma Numeric. Vector-to-host ratio (mosquitoes per human).
#' @param expo Numeric. Human exposure to Ae. albopictus in a day (between 0 and 1). (from Richards et al. 2006 : [10.1093/jmedent/43.3.543])
#' @return Numeric. Basic reproduction number R₀.
#'
#' @examples
#' compute_R0_poletti(a = 0.3, b = 0.5, c = 0.5, r = 5, N = 10, p = 0.9, ma = 3)
#'
#' @export
compute_R0_poletti <- function(a, b, c, r, N, p,expo, ma) {
  # Avoid log(0)
  to <- ifelse(p == 0, -log(p + 0.01), -log(p))
  
  # R0 vector to human
  r0VH <- ifelse(to == 0, a * b / (to + 0.01), a * b / to)
  
  # R0 human to vector
  r0HV <- expo * ma * c * r * ((1 / N) / ((1 / N) + to)) 
  
  # Total R0
  r0 <- r0VH * r0HV
  
  return(list(r0 = r0, r0VH = r0VH, r0HV = r0HV))
}

# Daily survival probability p
compute_p <- function(parity_rate, gc) {
  parity_rate^(1 / gc)
}


#' Compute EIP using Caminade's formula for Sensibility analysis
#' @title Get N (Survival Probability) based on temperature using Caminade's formula
#' @param T Temperature in degrees Celsius
#' @return Numeric value of N (survival probability)
get_N_from_T_DENV_Caminade <- function(T) {
  return(1.03 * (4 + exp(5.15 - 0.123 * T)))  # Caminade 2016 formula
}

#' Compute EIP using Benkimoun's formula for Sensibility analysis
#' @title Get N (Survival Probability) based on temperature using Benkimoun's formula
#' @param T Temperature in degrees Celsius
#' @return Numeric value of N (survival probability)
get_N_from_T_DENV_Benkimoun <- function(T) {
  return(0.11 * T^2 - 7.13 * T + 121.17)  # Benkimoun 2021 formula
}

#' Compute vector competence using Benkimoun's formula for Sensibility analysis
#' @title Get Vector Competence 'b' based on temperature using Benkimoun's formula
#' @param T Temperature in degrees Celsius
#' @return Numeric value of b (vector competence)
get_b_from_T_DENV_Benkimoun <- function(T) {
  return(0.2593 * T - 3.2705 - 0.0043 * T^2)
}

#' Compute vector capacity for sobol analysis indices 
#' @title Compute Vector Competence
#' @param X A list of parameters (including temperature, parity, etc.)
#' @param method_N The method to compute N ("Caminade2016" or "Benkimoun2021")
#' @return Numeric value of vector competence (vc)
compute_capacity_sobol <- function(X, method_N) {
  with(as.list(X), {
    # Compute N based on the selected method
    N <- switch(method_N,
                "Caminade2016" = get_N_from_T_DENV_Caminade(T),
                "Benkimoun2021" = get_N_from_T_DENV_Benkimoun(T),
                stop("Unrecognized method for N"))
    
    gc <- compute_gc(T)  # Ensure compute_gc() is properly defined elsewhere
    p <- parity^(1 / gc)  # Probability of survival
    a <- phi / gc  
    
    # Validate probability p (survival probability)
    if (any(p < 0 | p > 1)) {
      stop("Survival probability 'p' must be between 0 and 1.")
    }
    
    # Compute vector competence (vc) based on survival probability p
    if (any(p == 0)) {
      vc <- 0  # No competence if survival probability is zero
    } else if (any(p == 1)) {
      vc <- (expo*ma * a * (p^N)) / (1 - 0.99)  # Check the formula here
    } else {
      vc <- (expo*ma * a * (p^N)) / (-log(p))  # Classical calculation
    }
    return(vc)
  })
}

#' Compute R0 for sobol analysis indices 
#' @title Compute R0 based on vector competence
#' @param X A list of parameters (including temperature, parity, etc.)
#' @param method_N The method to compute N ("Caminade2016" or "Benkimoun2021")
#' @param method_b The method to compute b ("Benkimoun2021" or others)
#' @return Numeric value of R0 (basic reproduction number)
compute_R0_sobol <- function(X, method_N, method_b) {
  with(as.list(X), {
    capacity_vect <- compute_capacity_sobol(X, method_N)  # Get vector competence
    if (method_b %in% c("Verga-Rua2013", "Metelman2019", "Caminade2016")) { ## if the method used for b is one of these, b is constant
      b <- b  # b is in the 6th column of X
    } else {
      b <- switch(method_b,
                  "Benkimoun2021" = max(0, get_b_from_T_DENV_Benkimoun(T)), ## b depends of temperature
                  stop("Unrecognized method for b"))
    }
    r0 <- capacity_vect * r * b * c
    return(r0)
  })
}

#' Compute R0 for sobol analysis indices 
#' @title Compute R0 using Poletti's method
#' @param X A list of parameters (including temperature, parity, etc.)
#' @param method_N The method to compute N ("Caminade2016" or "Benkimoun2021")
#' @param method_b The method to compute b ("Benkimoun2021" or others)
#' @return Numeric value of R0 using Poletti's method
compute_R0_poletti_sobol <- function(X, method_N, method_b) {
  with(as.list(X), {
    # Handle 'b' value based on selected method
    if (method_b %in% c("Verga-Rua2013", "Metelman2019", "Caminade2016")) {## if the method used for b is one of these, b is constant
      b <- b  # b is in the 6th column of X
    } else {
      b <- switch(method_b,
                  "Benkimoun2021" = max(0, get_b_from_T_DENV_Benkimoun(T)),
                  stop("Unrecognized method for b"))
    }
    
    # Compute N based on the selected method
    N <- switch(method_N,
                "Caminade2016" = get_N_from_T_DENV_Caminade(T),
                "Benkimoun2021" = get_N_from_T_DENV_Benkimoun(T),
                stop("Unrecognized method for N"))
    
    # Compute gonotrophic cycle
    gc <- compute_gc(T)
    p <- parity^(1 / gc)
    a <- phi / gc
    to <- ifelse(p == 0, -log(p + 1e-10), -log(p))
    
    # Compute R0 using Poletti's method
    r0VH <- ifelse(to == 0, a * b / (to + 1e-3), a * b / to)
    r0HV <- expo*ma * c * r * ((1 / N) / ((1 / N) + to))
    r0 <- r0VH * r0HV
    return(r0)
  })
}

#' Calculate Sobol indices for the parameters
#' @title Run Sobol Sensitivity Analysis for different combinations
#' @param method_N The method to compute N ("Caminade2016" or "Benkimoun2021")
#' @param method_b The method to compute b ("Benkimoun2021" or others)
#' @param method_R0 The method to compute R0 ("R0_basic" or "R0_poletti")
#' @param n Number of simulations (default is 1000)
#' @return A tibble summarizing the first-order and total-order sensitivity indices
run_sensitivity_for_combo_sobol <- function(method_N, method_b, method_R0, n = 1000) {
  # Define parameter ranges
  param_ranges <- base_param_ranges
  
  # Generate Sobol sensitivity matrices
  X1 <- data.frame(lapply(param_ranges, function(range) runif(n, min = range[1], max = range[2])))
  X2 <- data.frame(lapply(param_ranges, function(range) runif(n, min = range[1], max = range[2])))
  
  # Run Sobol-Jansen sensitivity analysis
  sobol_res <- sensitivity::soboljansen(
    model = function(X) {
      if (method_R0 == "R0_poletti") {
        return(compute_R0_poletti_sobol(X, method_N, method_b))  # Using Poletti's method
      } else if (method_R0 == "R0_basic") {
        return(compute_R0_sobol(X, method_N, method_b))  # Using basic method
      } else {
        stop("Unrecognized method for R0")
      }
    },
    X1 = X1, X2 = X2,
    nboot = 100  # Bootstrapping for precision
  )
  
  # Summarize Sobol results
  sobol_df <- tibble(
    parameter = names(param_ranges),
    first_order = sobol_res$S$original,
    total_order = sobol_res$T$original
  )
  
  return(sobol_df)
}

#' Compute R0 for sensibility analysis indices using the method fast99
#' @title Compute R0 based on vector competence
#' @param X A list of parameters (including temperature, parity, etc.)
#' @param method_N The method to compute N ("Caminade2016" or "Benkimoun2021")
#' @param method_b The method to compute b ("Benkimoun2021" or others)
#' @return Numeric value of R0 (basic reproduction number)

compute_R0_fast99 <- function(X, method_N, method_b) {
  with(as.list(X), {
    gc <- compute_gc(T)
    a <- phi / gc
    p <- compute_p(parity, gc)
    
    N <- switch(method_N,
                "Caminade2016" = get_N_from_T_DENV_Caminade(T),
                "Benkimoun2021" = get_N_from_T_DENV_Benkimoun(T),
                stop("Unknown N method"))
    capacity_vect <- compute_capacity(ma, a, p, expo,N)
    
    b <- if (method_b %in% c("Verga-Rua2013", "Metelman2019", "Caminade2016")) {## if the method used for b is one of these, b is constant
      b
    } else {
      switch(method_b,
             "Benkimoun2021" = max(0, get_b_from_T_DENV_Benkimoun(T)),
             stop("Unrecognized method for b"))
    }
    
    r0 <- capacity_vect * r * b * c
    return(r0)
  })
}

#' Compute R0 of POletti for sensibility analysis indices using the method fast99
#' @title Compute R0 using Poletti's method
#' @param X A list of parameters (including temperature, parity, etc.)
#' @param method_N The method to compute N ("Caminade2016" or "Benkimoun2021")
#' @param method_b The method to compute b ("Benkimoun2021" or others)
#' @return Numeric value of R0 using Poletti's method
compute_R0_poletti_fast99 <- function(X, method_N, method_b) {
  with(as.list(X), {
    # Handle 'b' value based on selected method
    if (method_b %in% c("Verga-Rua2013", "Metelman2019", "Caminade2016")) {
      b <- b  # b is in the 8th column of X
    } else {
      b <- switch(method_b,
                  "Benkimoun2021" = max(0, get_b_from_T_DENV_Benkimoun(T)),
                  stop("Unrecognized method for b"))
    }
    
    # Compute N based on the selected method
    N <- switch(method_N,
                "Caminade2016" = get_N_from_T_DENV_Caminade(T),
                "Benkimoun2021" = get_N_from_T_DENV_Benkimoun(T),
                stop("Unrecognized method for N"))
    
    # Compute gonotrophic cycle
    gc <- compute_gc(T)
    p <- parity^(1 / gc)
    a <- phi / gc
    to <- ifelse(p == 0, -log(p + 1e-10), -log(p))
    # Compute R0 using Poletti's method
    r0VH <- ifelse(to == 0, a * b / (to + 1e-3), a * b / to)
    r0HV <-expo* ma * c * r * ((1 / N) / ((1 / N) + to))
    r0 <- r0VH * r0HV
    return(r0)
  })
}


#' Compute R0 for sensitivity analysis (FAST99-compatible wrapper)
#'
#' @title Adaptive wrapper to compute R0 (basic or Poletti) for DENV models
#' @description This function prepares parameter sets for the computation of R0 using either a basic method or Poletti's formulation.
#'              It handles different methods for estimating vector competence (b) and extrinsic incubation period (N).
#'              Intended for use with sensitivity analysis (e.g. FAST99).
#'
#' @param X A matrix of input parameters (each row is one parameter set; columns must be in order: T, c, r, parity, ma, and possibly b)
#' @param method_N Character. Method used to compute the extrinsic incubation period N. Options: "Caminade2016", "Benkimoun2021".
#' @param method_b Character. Method used to compute vector competence b. Options: "Caminade2016", "Benkimoun2021", "Verga-Rua2013", "Metelman2019".
#' @param method_R0 Character. Method used to compute R0. Options: "R0_basic", "R0_poletti".
#'
#' @return A numeric vector of R0 values computed for each parameter set.
#' @examples
#' # Used internally in sensitivity functions such as run_sensitivity_fast99()
#' # Not intended for standalone use unless inputs are well-formed.
wrapper_R0_DENV_adapt <- function(X, method_N, method_b, method_R0) {
  print(paste("method_N:", method_N, "method_b:", method_b, "method_R0:", method_R0)) 
  apply(X, 1, function(x) {
    # Assemble input parameter list
    x_list <- list(
      T = x[1],                                # Temperature
      c = x[2],                                # Transmission probability host→vector
      r = x[3],                                # Human recovery rate
      parity = x[4],                           # Parity rate
      ma = x[5],                                # Vector-to-host ratio
      expo=x[6],
      phi=x[7],
      b = if (method_b %in% c("Verga-Rua2013", "Metelman2019", "Caminade2016")) x[8] else NA  # Optional: b
    )
    
    # Compute 'b' dynamically if needed
    if (method_b == "Benkimoun2021") {
      x_list$b <- max(0, get_b_from_T_DENV_Benkimoun(x_list$T))
    }
    
    # Choose R0 calculation method
    if (method_R0 == "R0_poletti") {
      return(compute_R0_poletti_fast99(x_list, method_N, method_b))
    } else if (method_R0 == "R0_basic") {
      return(compute_R0_fast99(x_list, method_N, method_b))
    } else {
      stop("Unrecognized method for R0")
    }
  })
}

#' Compute FAST99 Sensitivity Analysis on R0 Model
#' @title Run FAST99 Sensitivity Analysis on R0 Model
#' @description Executes the FAST99 global sensitivity analysis method using specified
#' combinations of transmission model components (methods for N, b, and R0).
#'
#' @param method_N Character. Method used to compute N (e.g., "Caminade2016" or "Benkimoun2021").
#' @param method_b Character. Method used to compute b (e.g., "Caminade2016", "Verga-Rua2013", or "Benkimoun2021").
#' @param method_R0 Character. R0 model used (e.g., "R0_poletti" or "R0_basic").
#' @param n Integer. Number of samples per input factor (default is 1000).
#'
#' @return A list with two elements:
#' \describe{
#'   \item{fast_res}{A list object returned by `fast99()` containing raw sensitivity analysis results:
#'     \itemize{
#'       \item \code{D1}: First-order (main effect) variances.
#'       \item \code{Dt}: Total-effect variances.
#'       \item \code{V}: Total output variance.
#'     }
#'   }
#'   \item{summary}{A tibble summarizing normalized sensitivity indices:
#'     \itemize{
#'       \item \code{parameter}: The name of the input factor.
#'       \item \code{S_main}: Normalized first-order index (\code{D1 / V}).
#'       \item \code{S_total}: Normalized total-effect index (\code{Dt / V}).
#'     }
#'   }
#' }
#'
#' @details Adds 'b' to the parameter ranges if the selected `method_b` treats it as a free input.
#' Sampling is done via a uniform distribution using `q = "qunif"`.
#'
#' @examples
#' run_sensitivity_fast99("Caminade2016", "Benkimoun2021", "R0_poletti", n = 500)
run_sensitivity_fast99 <- function(method_N, method_b, method_R0, n = 1000) {
  param_ranges <- base_param_ranges
  if (method_b %in% c("Verga-Rua2013", "Metelman2019", "Caminade2016")) {
    param_ranges <- rbind(param_ranges, b = c(0.1, 0.75))
  }
  param_names <- rownames(param_ranges)
  q.arg <- as.list(param_ranges %>% pmap(~ list(min = ..1, max = ..2)))
  set.seed(123)
  fast_res <- fast99(
    model = function(X) wrapper_R0_DENV_adapt(X, method_N, method_b, method_R0),
    factors = param_names,
    n = n,
    q = "qunif",
    q.arg = q.arg
  )
  summary_res <- tibble(
    parameter = names(fast_res$D1),
    S_main = fast_res$D1 / fast_res$V,
    S_total = fast_res$Dt / fast_res$V
  ) %>% arrange(desc(S_main))
  list(fast_res = fast_res, summary = summary_res)
}

