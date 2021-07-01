#' Trim data for downstream analyses
#'
#' @description
#' trimData() allows users to define start date, start time and number of days for analysing data sets.
#'
#' @param data Input data file. This must be DAM scanned monitor files, saved in 1-min intervals starting at time 00:00.
#' @param start.date Define starting date for analysis. Date inputs must be in this form: "DD Mon YY"; August 6, 2020 must be input as "6 Aug 20". Any other form will result in computation errors.
#' @param start.time Define start time for analysis. Times must be input in the 24-h format and in the following form: "HH:MM"; 10 AM must be "10:00" and 10 PM must be "22:00".
#' @param n.days Number of days to analyse.
#' @param bin Intervals in which data are sampled (in minutes). This defaults to 1.
#' @param t.cycle Define the period of the environmental cycle or a single day in hours. This defaults to 24.
#'
#' @importFrom lubridate hour minute
#' @importFrom grDevices rgb
#' @importFrom stats aggregate fitted lm na.omit sd
#' 
#' @return A \code{data.frame} containing trimmed data.
#'
#' @export trimData
#'
#' @examples
#' td <- trimData(data = df, start.date = "19 Dec 20", start.time = "21:00",
#' n.days = 10, bin = 1, t.cycle = 24)

trimData <- function(data, start.date, start.time, n.days, bin = 1, t.cycle = 24) {
  
  requireNamespace("lubridate")

  # library(lubridate)

  startdate <- as.Date(start.date, format = '%d %b %y')
  colnames(data) <- c("x","date","time",seq(1, (length(data[1,])-3), by = 1))
  data[,"date"] <- as.Date(data[,"date"], format = '%d %b %y')

  df <- subset(data, date > (startdate-1))

  starttime <- as.POSIXct(paste(startdate, start.time, sep = " "))
  s_per_day <- (60/bin)*t.cycle
  h <- lubridate::hour(starttime)
  m <- lubridate::minute(starttime)
  h_in_min <- h*60
  min <- h_in_min + m
  ind.start <- (min/bin) + 1 + 1

  df <- df[ind.start:((ind.start - 1) + (s_per_day*n.days)),]
  return(df)
}
