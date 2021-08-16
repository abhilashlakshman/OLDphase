#' Bin activity data to desired intervals
#'
#' @description
#' Allows users to bin data sets into intervals of time different from the data collection interval. The input for this function must be the output of the trimData() function. The output of this function is a data frame. The first column of which stores Zeitgeber Time values (assuming that the start.time in the trimData() function was set at Zeitgeber Time 00). All subsequent columns have binned activity data for each fly.
#'
#' @param data Input data file. The input for this function must be the output of the function trimData(). See ??trimData().
#' @param input.bin Define the bin size (in minutes) in the input file. This defaults to 1.
#' @param output.bin Define the desired output bin size (in minutes). This defaults to 30.
#' @param t.cycle Define the period of the environmental cycle or a single day in hours. This defaults to 24.
#'
#' @importFrom grDevices rgb
#' @importFrom stats aggregate fitted lm na.omit sd
#' 
#' @return A \code{matrix} \code{array} with 33 columns (number of rows depends on number of days, and the input parameters of this function):
#' \describe{
#' \item{ZT/CT/Time}{ZT/CT values starting at ZT/CT00 (time at which light turn ON or onset of subjective day). In case start time is not ZT/CT00, Time refers to the first window since the user defined start window, which will repeat every user-defined T-cycle.}
#' \item{I1:I32}{Columns of binned locomotor activity data (each column represents a single fly).}
#' }
#'
#' @export binData
#'
#' @examples
#' td <- trimData(data = df, start.date = "19 Dec 20", start.time = "21:00",
#' n.days = 10, bin = 1, t.cycle = 24)
#' bd <- binData(data = td)
#' bd <- binData(data = td, input.bin = 1, output.bin = 15, t.cycle = 24)

binData <- function(data, input.bin = 1, output.bin = 30, t.cycle = 24) {

  raw <- data[,-c(1:10)]
  s_per_day <- (60/output.bin)*t.cycle
  indices <- seq(1, length(raw[,1]), by = output.bin/input.bin)
  binned <- matrix(NA, nrow = (s_per_day*length(raw[,1])/((60/input.bin)*t.cycle)), ncol = length(raw[1,]))

  for (i in 1:length(indices)){
    for (j in 1:length(raw[1,])){
      binned[i,j] <- sum(raw[indices[i]:(indices[i] + output.bin - 1),j])
    }
  }

  column.names <- c()
  for (i in 1:length(binned[1,])) {
    column.names[i] <- paste("I",i, sep = "")
  }
  colnames(binned) <- column.names

  time <- seq(output.bin/60, t.cycle, by = output.bin/60)
  t <- as.matrix(rep(time, length(binned[,1])/s_per_day))
  colnames(t) <- "ZT/CT/Time"
  output <- cbind(t,binned)

  return(output)
}
