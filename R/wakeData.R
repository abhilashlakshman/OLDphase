#' Equivalent of binData() but for "wake" data
#'
#' @description
#' Allows users to bin data sets into intervals of time different from the data collection interval. The input for this function must be the output of the trimData() function. The output of this function is a data frame. The first column of which stores Zeitgeber Time values (assuming that the start.time in the trimData() function was set at Zeitgeber Time 00). All subsequent columns have binned "wake" data for each fly.
#' In a particular bin, sleep is calculated as the total minutes of inactivity equal to or greater than the defined threshold (sleep.def; typically, 5-minutes). Wake is defined as the total time spent by the fly not sleeping. See also ??wakeData().
#'
#' @param data Input data file. The input for this function must be the output of the function trimData(). See ??trimData().
#' @param sleep.def Definition of sleep. Traditionally, a single bout of sleep is defined as any duration of inactivity that is equal to or greater than 5-minutes. However, sometimes it may be of interest to examine longer bouts of sleep; sleep.def allows users to change the definition of sleep. This defaults to 5.
#' @param bin Intervals in which data are saved (in minutes). This defaults to 30. The value of bin cannot be lower than that of sleep.def.
#' @param t.cycle Define the period of the environmental cycle or a single day in hours. This defaults to 24.
#'
#' @importFrom grDevices rgb
#' @importFrom stats aggregate fitted lm na.omit sd
#'
#' @export wakeData
#'
#' @examples
#' td <- trimData(data = df, start.date = "19 Dec 20", start.time = "21:00",
#' n.days = 4, bin = 1, t.cycle = 24)
#' wd <- wakeData(data = td)
#' wd <- wakeData(data = td, sleep.def = 20, bin = 60, t.cycle = 24)

wakeData <- function(data, sleep.def = c(5), bin = 30, t.cycle = 24) {
  raw <- data[,-c(1:10)]
  s_per_day <- (60/bin)*t.cycle

  raw[,1:length(raw[1,])][raw[,1:length(raw[1,])] >= 1] <- -1
  raw[,1:length(raw[1,])][raw[,1:length(raw[1,])] == 0] <- 1
  raw[,1:length(raw[1,])][raw[,1:length(raw[1,])] == -1] <- 0

  # binned_full_run.sleep <- (length(raw[,1])/1440)*s_per_day
  binned_full_run.wake <- (length(raw[,1])/1440)*s_per_day

  # sleep <- matrix(NA, nrow = binned_full_run.sleep, ncol = 32)
  wake <- matrix(NA, nrow = binned_full_run.wake, ncol = 32)

  index.wake <- seq(1, length(raw[,1]), by = bin)

  if (length(sleep.def) == 1) {
    for (i in 1:length(index.wake)) {
      for (j in 1:length(raw[1,])) {
        x <- raw[index.wake[i]:(index.wake[i]+bin-1),j]
        y <- rle(x)
        d_y <- as.data.frame(unclass(y))
        # dd_y_sleep <- subset(d_y, values == 1 & lengths >= sleep.def)
        dd_y_wake <- subset(d_y, (d_y$values == 1 & d_y$lengths < sleep.def[1]) | (d_y$values == 0))
        # sleep[i,j] <- sum(dd_y_sleep$lengths)
        wake[i,j] <- sum(dd_y_wake$lengths)
      }
    }
  } else if (length(sleep.def) == 2) {
    for (i in 1:length(index.wake)) {
      for (j in 1:length(raw[1,])) {
        x <- raw[index.wake[i]:(index.wake[i]+bin-1),j]
        y <- rle(x)
        d_y <- as.data.frame(unclass(y))
        # dd_y_sleep <- subset(d_y, values == 1 & lengths >= sleep.def)
        dd_y_wake <- subset(d_y, (d_y$values == 1 & d_y$lengths < sleep.def[1]) | (d_y$values == 0) | d_y$lengths > sleep.def[2])
        # sleep[i,j] <- sum(dd_y_sleep$lengths)
        wake[i,j] <- sum(dd_y_wake$lengths)
      }
    }
  }

  column.names <- c()
  for (i in 1:length(wake[1,])) {
    column.names[i] <- paste("I",i, sep = "")
  }
  # colnames(sleep) <- column.names
  colnames(wake) <- column.names

  t <- seq((bin/60), t.cycle, by = (bin/60))
  zt <- as.data.frame(rep(t, length(wake[,1])/(t.cycle*(60/bin))))
  colnames(zt) <- c("ZT")

  # data.plot.sleep <- cbind(zt,sleep)
  data.plot.wake <- cbind(zt,wake)

  output <- list(
    # "Sleep" = data.plot.sleep,
    "Wake" = data.plot.wake
  )
  return(data.plot.wake)
}
