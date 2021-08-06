#' Computes and tabulates day-time and night-time sleep statistics
#'
#' @description
#' This function allows users to estimate day-time and night-time average sleep bout duration, number and latency. Sleep bout latency is defined as the time taken (in minutes) for the occurrence of the first sleep bout since respective transitions. The input for this function must be the output from the trimData() function. The output of this function is a matrix which contains fly-wise (each row) data. Sleep statistics for several sleep data can be examined by users; (i) data = "Sleep" will provide sleep statistics for standard sleep data wherein any contiguous bout of inactivity that is longer than or equal to 5-minutes is considered sleep, (ii) data = "DeepSleep" will provide sleep statistics for only all the deep sleep bouts, (iii) data = "LightSleep" will provide sleep statistics for only all the light sleep bouts, (iv) data = "LongDeepSleep" and (v) data = "ShortDeepSleep" will provide sleep statistics for the long bouts of deep sleep and short bouts of deep sleep, respectively, (vi) data = "ShortLightSleep", (vii) data = "InterLightSleep", and (viii) data = "LongLightSleep" will provide sleep statistics for short, intermediate and long bouts of light sleep, respectively. See ??deepSleepData() and ??lightSleepData() for information on deep sleep and light sleep categories and definitions.
#'
#' @param input Input data file. The input for this function must be the output of the function trimData(). See ??trimData().
#' @param sleep.def Definition of sleep. Traditionally, a single bout of sleep is defined as any duration of inactivity that is equal to or greater than 5-minutes. However, sometimes it may be of interest to examine longer bouts of sleep or specific bout durations; sleep.def allows users to change the definition of sleep. The default input is a single value vector of value 5. If users wish to analyse sleep only between 5 to 20 mins, the input must be c(5,20).
#' @param t.cycle Define the period of the environmental cycle or a single day in hours. This defaults to 24.
#' @param photoperiod Duration (in hours) of what can be considered day-phase. This defaults to 12.
#'
#' @importFrom grDevices rgb
#' @importFrom stats aggregate fitted lm na.omit sd
#' 
#' @return A \code{matrix} \code{array} matrix with 32 rows (one for each fly) and 7 columns:
#' \describe{
#' \item{Channel}{Fly identity.}
#' \item{Day.BoutNumber}{Number of sleep bouts in the user defined day time.}
#' \item{Day.BoutDuration}{Mean sleep duration in the user defined day time.}
#' \item{Day.Latency}{Time taken for the first sleep bout to occur in the user defined day time.}
#' \item{Night.BoutNumber}{Number of sleep bouts in the user defined night time.}
#' \item{Night.BoutDuration}{Mean sleep duration in the user defined night time.}
#' \item{Night.Latency}{Time taken for the first sleep bout to occur in the user defined night time.}
#' }
#'
#' @export sleepStat
#'
#' @examples
#' td <- trimData(data = df, start.date = "19 Dec 20", start.time = "21:00",
#' n.days = 10, bin = 1, t.cycle = 24)
#' slp.stat <- sleepStat(input = td)

sleepStat <- function(input, sleep.def = c(5), t.cycle = 24, photoperiod = 12) {
  raw <- input[,-c(1:10)]
  
  s_per_day <- 60*t.cycle
  n.days <- length(raw[,1])/s_per_day
  
  raw[,1:length(raw[1,])][raw[,1:length(raw[1,])] >= 1] <- -1
  raw[,1:length(raw[1,])][raw[,1:length(raw[1,])] == 0] <- 1
  raw[,1:length(raw[1,])][raw[,1:length(raw[1,])] == -1] <- 0
  
  cyc.wise.split <- list()
  cyc.wise.index <- seq(1, length(raw[,1]), by = (60*t.cycle))
  for (i in 1:n.days) {
    cyc.wise.split[[i]] <- raw[cyc.wise.index[i]:(cyc.wise.index[i]+((60*t.cycle)-1)),]
  }
  
  output <- matrix(NA, nrow = 32, ncol = 7)
  colnames(output) <- c("Channel", "Day.BoutNumber", "Day.BoutDuration", "Day.Latency",
                        "Night.BoutNumber", "Night.BoutDuration", "Night.Latency")
  for.row.nam <- seq(1, 32, by = 1)
  
  output[,"Channel"] <- for.row.nam
  
  day.cyc.boutnum <- matrix(NA, nrow = 32, ncol = n.days)
  day.cyc.boutdur <- matrix(NA, nrow = 32, ncol = n.days)
  day.cyc.boutlatency <- matrix(NA, nrow = 32, ncol = n.days)
  
  night.cyc.boutnum <- matrix(NA, nrow = 32, ncol = n.days)
  night.cyc.boutdur <- matrix(NA, nrow = 32, ncol = n.days)
  night.cyc.boutlatency <- matrix(NA, nrow = 32, ncol = n.days)
  
  if (length(sleep.def) == 1) {
    for (i in 1:length(cyc.wise.split)) {
      day <- cyc.wise.split[[i]][(1:(photoperiod*60)),]
      night <- cyc.wise.split[[i]][(((photoperiod*60) + 1):(60*t.cycle)),]
      
      for (j in 1:length(day[1,])) {
        d <- day[,j]
        d.rle <- rle(d)
        df.d <- as.data.frame((unclass(d.rle)))
        df.d$end <- cumsum(df.d$lengths)
        df.d$start <- df.d$end - df.d$lengths + 1
        ddf.d <- subset(df.d, df.d$values == 1 & df.d$lengths >= sleep.def[1])
        
        day.cyc.boutnum[j,i] <- length(ddf.d[,1])
        
        temp.boutdur.day <- mean(ddf.d$lengths)
        if (is.nan(temp.boutdur.day) | is.na(temp.boutdur.day)) {
          day.cyc.boutdur[j,i] <- NA
        } else {
          day.cyc.boutdur[j,i] <- temp.boutdur.day
        }
        
        temp.bout.lat.day <- ddf.d[1,"start"]
        if (is.nan(temp.bout.lat.day) | is.na(temp.bout.lat.day)) {
          day.cyc.boutlatency[j,i] <- NA
        } else {
          day.cyc.boutlatency[j,i] <- temp.bout.lat.day
        }
        
        
        n <- night[,j]
        n.rle <- rle(n)
        df.n <- as.data.frame((unclass(n.rle)))
        df.n$end <- cumsum(df.n$lengths)
        df.n$start <- df.n$end - df.n$lengths + 1
        ddf.n <- subset(df.n, df.n$values == 1 & df.n$lengths >= sleep.def[1])
        
        night.cyc.boutnum[j,i] <- length(ddf.n[,1])
        
        temp.boutdur.night <- mean(ddf.n$lengths)
        if (is.nan(temp.boutdur.night) | is.na(temp.boutdur.night)) {
          night.cyc.boutdur[j,i] <- NA
        } else {
          night.cyc.boutdur[j,i] <- temp.boutdur.night
        }
        
        temp.bout.lat.night <- ddf.n[1,"start"]
        if (is.nan(temp.bout.lat.night) | is.na(temp.bout.lat.night)) {
          night.cyc.boutlatency[j,i] <- NA
        } else {
          night.cyc.boutlatency[j,i] <- temp.bout.lat.night
        }
        
        
      }
    }
    
    for (i in 1:length(raw[1,])) {
      output[i,"Day.BoutNumber"] = mean(day.cyc.boutnum[i,], na.rm = T)
      output[i,"Night.BoutNumber"] = mean(night.cyc.boutnum[i,], na.rm = T)
      
      output[i,"Day.BoutDuration"] = mean(day.cyc.boutdur[i,], na.rm = T)
      output[i,"Night.BoutDuration"] = mean(night.cyc.boutdur[i,], na.rm = T)
      
      output[i,"Day.Latency"] = mean(day.cyc.boutlatency[i,], na.rm = T)
      output[i,"Night.Latency"] = mean(night.cyc.boutlatency[i,], na.rm = T)
    }
  } else if (length(sleep.def) == 2) {
    for (i in 1:length(cyc.wise.split)) {
      day <- cyc.wise.split[[i]][(1:(photoperiod*60)),]
      night <- cyc.wise.split[[i]][(((photoperiod*60) + 1):(60*t.cycle)),]
      
      for (j in 1:length(day[1,])) {
        d <- day[,j]
        d.rle <- rle(d)
        df.d <- as.data.frame((unclass(d.rle)))
        df.d$end <- cumsum(df.d$lengths)
        df.d$start <- df.d$end - df.d$lengths + 1
        ddf.d <- subset(df.d, df.d$values == 1 & df.d$lengths >= sleep.def[1] & df.d$lengths <= sleep.def[2])
        
        day.cyc.boutnum[j,i] <- length(ddf.d[,1])
        
        temp.boutdur.day <- mean(ddf.d$lengths)
        if (is.nan(temp.boutdur.day) | is.na(temp.boutdur.day)) {
          day.cyc.boutdur[j,i] <- NA
        } else {
          day.cyc.boutdur[j,i] <- temp.boutdur.day
        }
        
        temp.bout.lat.day <- ddf.d[1,"start"]
        if (is.nan(temp.bout.lat.day) | is.na(temp.bout.lat.day)) {
          day.cyc.boutlatency[j,i] <- NA
        } else {
          day.cyc.boutlatency[j,i] <- temp.bout.lat.day
        }
        
        n <- night[,j]
        n.rle <- rle(n)
        df.n <- as.data.frame((unclass(n.rle)))
        df.n$end <- cumsum(df.n$lengths)
        df.n$start <- df.n$end - df.n$lengths + 1
        ddf.n <- subset(df.n, df.n$values == 1 & df.n$lengths >= sleep.def[1] & df.n$lengths <= sleep.def[2])
        
        night.cyc.boutnum[j,i] <- length(ddf.n[,1])
        
        temp.boutdur.night <- mean(ddf.n$lengths)
        if (is.nan(temp.boutdur.night) | is.na(temp.boutdur.night)) {
          night.cyc.boutdur[j,i] <- NA
        } else {
          night.cyc.boutdur[j,i] <- temp.boutdur.night
        }
        
        temp.bout.lat.night <- ddf.n[1,"start"]
        if (is.nan(temp.bout.lat.night) | is.na(temp.bout.lat.night)) {
          night.cyc.boutlatency[j,i] <- NA
        } else {
          night.cyc.boutlatency[j,i] <- temp.bout.lat.night
        }
        
      }
    }
    
    for (i in 1:length(raw[1,])) {
      output[i,"Day.BoutNumber"] = mean(day.cyc.boutnum[i,], na.rm = T)
      output[i,"Night.BoutNumber"] = mean(night.cyc.boutnum[i,], na.rm = T)
      
      output[i,"Day.BoutDuration"] = mean(day.cyc.boutdur[i,], na.rm = T)
      output[i,"Night.BoutDuration"] = mean(night.cyc.boutdur[i,], na.rm = T)
      
      output[i,"Day.Latency"] = mean(day.cyc.boutlatency[i,], na.rm = T)
      output[i,"Night.Latency"] = mean(night.cyc.boutlatency[i,], na.rm = T)
    }
  }
  
  return(output)
}
