#' Computes and tabulates day-time and night-time sleep statistics
#'
#' @description
#' This function allows users to estimate day-time and night-time average sleep bout duration, number and latency. Sleep bout latency is defined as the time taken (in minutes) for the occurrence of the first sleep bout since respective transitions. The input for this function must be the output from the trimData() function. The output of this function is a matrix which contains fly-wise (each row) data. Sleep statistics for several sleep data can be examined by users; (i) data = "Sleep" will provide sleep statistics for standard sleep data wherein any contiguous bout of inactivity that is longer than or equal to 5-minutes is considered sleep, (ii) data = "DeepSleep" will provide sleep statistics for only all the deep sleep bouts, (iii) data = "LightSleep" will provide sleep statistics for only all the light sleep bouts, (iv) data = "LongDeepSleep" and (v) data = "ShortDeepSleep" will provide sleep statistics for the long bouts of deep sleep and short bouts of deep sleep, respectively, (vi) data = "ShortLightSleep", (vii) data = "InterLightSleep", and (viii) data = "LongLightSleep" will provide sleep statistics for short, intermediate and long bouts of light sleep, respectively. See ??deepSleepData() and ??lightSleepData() for information on deep sleep and light sleep categories and definitions.
#'
#' @param input Input data file. The input for this function must be the output of the function trimData(). See ??trimData().
#' @param sleep.def Definition of sleep. Traditionally, a single bout of sleep is defined as any duration of inactivity that is equal to or greater than 5-minutes. However, sometimes it may be of interest to examine longer bouts of sleep; sleep.def allows users to change the definition of sleep. This defaults to 5.
#' @param t.cycle Define the period of the environmental cycle or a single day in hours. This defaults to 24.
#' @param photoperiod Duration (in hours) of what can be considered day-phase. This defaults to 12.
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
  cyc.wise.index <- seq(1, length(raw[,1]), by = 1440)
  for (i in 1:n.days) {
    cyc.wise.split[[i]] <- raw[cyc.wise.index[i]:(cyc.wise.index[i]+(1440-1)),]
  }

  output <- matrix(NA, nrow = 32, ncol = 7)
  # row.nam <- c()
  # for (i in 1:32) {
  #   row.nam[i] <- paste("Ch-", i, sep = "")
  # }
  # rownames(output) <- row.nam
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
      night <- cyc.wise.split[[i]][(((photoperiod*60) + 1):1440),]

      for (j in 1:length(day[1,])) {
        d <- day[,j]
        d.rle <- rle(d)
        df.d <- as.data.frame((unclass(d.rle)))
        df.d$end <- cumsum(df.d$lengths)
        df.d$start <- df.d$end - df.d$lengths + 1
        ddf.d <- subset(df.d, values == 1 & lengths >= sleep.def[1])

        day.cyc.boutnum[j,i] <- length(ddf.d[,1])
        day.cyc.boutdur[j,i] <- mean(ddf.d$lengths)
        day.cyc.boutlatency[j,i] <- ddf.d[1,"start"]

        n <- night[,j]
        n.rle <- rle(n)
        df.n <- as.data.frame((unclass(n.rle)))
        df.n$end <- cumsum(df.n$lengths)
        df.n$start <- df.n$end - df.n$lengths + 1
        ddf.n <- subset(df.n, values == 1 & lengths >= sleep.def[1])

        night.cyc.boutnum[j,i] <- length(ddf.n[,1])
        night.cyc.boutdur[j,i] <- mean(ddf.n$lengths)
        night.cyc.boutlatency[j,i] <- ddf.n[1,"start"]

      }
    }

    for (i in 1:length(raw[1,])) {
      output[i,"Day.BoutNumber"] = mean(day.cyc.boutnum[i,])
      output[i,"Night.BoutNumber"] = mean(night.cyc.boutnum[i,])

      output[i,"Day.BoutDuration"] = mean(day.cyc.boutdur[i,])
      output[i,"Night.BoutDuration"] = mean(night.cyc.boutdur[i,])

      output[i,"Day.Latency"] = mean(day.cyc.boutlatency[i,])
      output[i,"Night.Latency"] = mean(night.cyc.boutlatency[i,])
    }
  } else if (length(sleep.def) == 2) {
    for (i in 1:length(cyc.wise.split)) {
      day <- cyc.wise.split[[i]][(1:(photoperiod*60)),]
      night <- cyc.wise.split[[i]][(((photoperiod*60) + 1):1440),]

      for (j in 1:length(day[1,])) {
        d <- day[,j]
        d.rle <- rle(d)
        df.d <- as.data.frame((unclass(d.rle)))
        df.d$end <- cumsum(df.d$lengths)
        df.d$start <- df.d$end - df.d$lengths + 1
        ddf.d <- subset(df.d, values == 1 & lengths >= sleep.def[1] & lengths < sleep.def[2])

        day.cyc.boutnum[j,i] <- length(ddf.d[,1])
        day.cyc.boutdur[j,i] <- mean(ddf.d$lengths)
        day.cyc.boutlatency[j,i] <- ddf.d[1,"start"]

        n <- night[,j]
        n.rle <- rle(n)
        df.n <- as.data.frame((unclass(n.rle)))
        df.n$end <- cumsum(df.n$lengths)
        df.n$start <- df.n$end - df.n$lengths + 1
        ddf.n <- subset(df.n, values == 1 & lengths >= sleep.def[1] & lengths < sleep.def[2])

        night.cyc.boutnum[j,i] <- length(ddf.n[,1])
        night.cyc.boutdur[j,i] <- mean(ddf.n$lengths)
        night.cyc.boutlatency[j,i] <- ddf.n[1,"start"]

      }
    }

    for (i in 1:length(raw[1,])) {
      output[i,"Day.BoutNumber"] = mean(day.cyc.boutnum[i,])
      output[i,"Night.BoutNumber"] = mean(night.cyc.boutnum[i,])

      output[i,"Day.BoutDuration"] = mean(day.cyc.boutdur[i,])
      output[i,"Night.BoutDuration"] = mean(night.cyc.boutdur[i,])

      output[i,"Day.Latency"] = mean(day.cyc.boutlatency[i,])
      output[i,"Night.Latency"] = mean(night.cyc.boutlatency[i,])
    }
  }

  return(output)
}
