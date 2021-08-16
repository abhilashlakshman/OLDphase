#' Computes, tabulates and plots day-time and night-time onset of sleep bouts and their respective lengths
#'
#' @description
#' This function allows users to compute day-time and night-time onset of sleep bouts and their respective bout lengths. Also generates a plot to visualise and compare within and across day and night-time windows. The input for this function must be the output from the trimData() function. Number of days to analyse must be at least 2 days. If the number of days is more than 2 days, the function will compute statistics for the first day only. The output of this function is a list with three elements, i.e., "Daytime.Data" is a data frame which has the onset of sleep bout (minutes since the start of the day window) and the length of that bout (in minutes) for each fly, "Nighttime.Data" which is the same as "Daytime.Data" but for the night window, and "Plots" which allows the visualisation of the two data sets.
#'
#' @param data Input data file. The input for this function must be the output of the function trimData(). See ??trimData().
#' @param sleep.def Definition of sleep. Traditionally, a single bout of sleep is defined as any duration of inactivity that is equal to or greater than 5-minutes. However, sometimes it may be of interest to examine longer bouts of sleep or specific bout durations; sleep.def allows users to change the definition of sleep. The default input is a single value vector of value 5. If users wish to analyse sleep only between 5 to 20 mins, the input must be c(5,20).
#' @param t.cycle Define the period of the environmental cycle or a single day in hours. This defaults to 24.
#' @param photoperiod Duration (in hours) of what can be considered day-phase. This defaults to 12.
#'
#' @importFrom grDevices rgb
#' @importFrom stats aggregate fitted lm na.omit sd
#' 
#' @return A \code{list} with three items:
#' \describe{
#' \item{Daytime.Data}{
#' \describe{
#' \item{Channel}{Fly identity.}
#' \item{Start}{Onset of sleep bout (in minutes since the start of day window)}
#' \item{BoutLength}{Length of the sleep bout (in minutes)}
#' }
#' }
#' 
#' \item{Nighttim.Data}{
#' \describe{
#' \item{Channel}{Fly identity.}
#' \item{Start}{Onset of sleep bout (in minutes since the start of night window)}
#' \item{BoutLength}{Length of the sleep bout (in minutes)}
#' }
#' }
#' }
#'
#' @export sleepOnsetBoutLength
#'
#' @examples
#' td <- trimData(data = df, start.date = "28 Dec 20", start.time = "21:00",
#' n.days = 2, bin = 1, t.cycle = 24)
#' bout.onset.vs.length <- sleepOnsetBoutLength(data = td)

sleepOnsetBoutLength <- function(data, sleep.def = c(5), t.cycle = 24, photoperiod = 12) {
  requireNamespace("plotly")
  if (requireNamespace("plotly", quietly = T)) {
    
    raw <- data[,-c(1:10)]
    
    raw[,1:length(raw[1,])][raw[,1:length(raw[1,])] >= 1] <- -1
    raw[,1:length(raw[1,])][raw[,1:length(raw[1,])] == 0] <- 1
    raw[,1:length(raw[1,])][raw[,1:length(raw[1,])] == -1] <- 0
    
    day.raw <- raw[1:1440,]
    night.raw <- raw[(((photoperiod*60) + 1):((photoperiod*60)+(60*t.cycle))),]
    
    day.flies.dat <- list()
    night.flies.dat <- list()
    
    if (length(sleep.def) == 1) {
      for (i in 1:length(day.raw[1,])) {
        x <- day.raw[,i]
        y <- as.data.frame(unclass(rle(x)))
        y$end <- cumsum(y$lengths)
        y$start <- y$end - y$lengths + 1
        z <- subset(y, y$values == 1 & y$lengths >= sleep.def[1])
        z$Channel <- rep(i, length(z[,1]))
        z$Start <- z$start
        z$BoutLength <- z$lengths
        day.flies.dat[[i]] <- subset(z, z$start <= 720, select = c("Channel", "Start", "BoutLength"))
      }
      
      for (i in 1:length(night.raw[1,])) {
        x <- night.raw[,i]
        y <- as.data.frame(unclass(rle(x)))
        y$end <- cumsum(y$lengths)
        y$start <- y$end - y$lengths + 1
        z <- subset(y, y$values == 1 & y$lengths >= sleep.def[1])
        z$Channel <- rep(i, length(z[,1]))
        z$Start <- z$start
        z$BoutLength <- z$lengths
        night.flies.dat[[i]] <- subset(z, z$start <= 720, select = c("Channel", "Start", "BoutLength"))
      }
    } else if (length(sleep.def) == 2) {
      for (i in 1:length(day.raw[1,])) {
        x <- day.raw[,i]
        y <- as.data.frame(unclass(rle(x)))
        y$end <- cumsum(y$lengths)
        y$start <- y$end - y$lengths + 1
        z <- subset(y, y$values == 1 & y$lengths >= sleep.def[1] & y$lengths <= sleep.def[2])
        z$Channel <- rep(i, length(z[,1]))
        z$Start <- z$start
        z$BoutLength <- z$lengths
        day.flies.dat[[i]] <- subset(z, z$start <= 720, select = c("Channel", "Start", "BoutLength"))
      }
      
      for (i in 1:length(night.raw[1,])) {
        x <- night.raw[,i]
        y <- as.data.frame(unclass(rle(x)))
        y$end <- cumsum(y$lengths)
        y$start <- y$end - y$lengths + 1
        z <- subset(y, y$values == 1 & y$lengths >= sleep.def[1] & y$lengths <= sleep.def[2])
        z$Channel <- rep(i, length(z[,1]))
        z$Start <- z$start
        z$BoutLength <- z$lengths
        night.flies.dat[[i]] <- subset(z, z$start <= 720, select = c("Channel", "Start", "BoutLength"))
      }
    }
    
    day.flies.compiled <- do.call(rbind, day.flies.dat)
    night.flies.compiled <- do.call(rbind, night.flies.dat)
    
    f1 <- list(
      family = "Arial, sans-serif",
      size = 20,
      color = "black"
    )
    f2 <- list(
      family = "Arial, sans-serif",
      size = 14,
      color = "black"
    )
    f3 <- list(
      family = "Arial, sans-serif",
      size = 14,
      color = rgb(0,0,0,0)
    )
    
    day.plot <- plot_ly(
      x = day.flies.compiled$Start,
      y = day.flies.compiled$BoutLength,
      type = "scatter",
      mode = "markers",
      marker = list(
        color = rgb(0,0,1,0.6)
      ),
      name = "Day-time"
    )
    
    night.plot <- plot_ly(
      x = night.flies.compiled$Start,
      y = night.flies.compiled$BoutLength,
      type = "scatter",
      mode = "markers",
      marker = list(
        color = rgb(1,0,0,0.6)
      ),
      name = "Night-time"
    )
    
    output.plot <- subplot(
      day.plot, night.plot,
      nrows = 1,
      shareX = F, shareY = F
    )%>%
      layout(
        showlegend = T,
        legend = list(
          font = f2,
          bgcolor = rgb(0,0,0,0)
        ),
        xaxis = list(
          showgrid = F,
          showline = T,
          titlefont = f1,
          tickfont = f2,
          title = "Onset of sleep bout (min)",
          linecolor = "black",
          # linewidth = 4,
          mirror = TRUE,
          autotick = TRUE,
          ticks = "inside",
          tickcolor = "black",
          # tickwidth = 4,
          range = c(0, 720+0.5)
        ),
        xaxis2 = list(
          showgrid = F,
          showline = T,
          titlefont = f1,
          tickfont = f2,
          title = "Onset of sleep bout (min)",
          linecolor = "black",
          # linewidth = 4,
          mirror = TRUE,
          autotick = TRUE,
          ticks = "inside",
          tickcolor = "black",
          # tickwidth = 4,
          range = c(0, 720+0.5)
        ),
        yaxis = list(
          showgrid = F,
          showline = T,
          titlefont = f1,
          tickfont = f2,
          title = "Length of sleep bout (min)",
          linecolor = "black",
          # linewidth = 4,
          mirror = TRUE,
          autotick = TRUE,
          ticks = "inside",
          tickcolor = "black",
          range = c(0, 720+5)
          # tickwidth = 4,
        ),
        yaxis2 = list(
          showgrid = F,
          showline = T,
          titlefont = f1,
          tickfont = f3,
          title = "",
          linecolor = "black",
          # linewidth = 4,
          mirror = TRUE,
          autotick = TRUE,
          ticks = "inside",
          tickcolor = "black",
          range = c(0, 720+5)
          # tickwidth = 4,
        )
      )
    
    final.output <- list(
      "Daytime.Data" = day.flies.compiled,
      "Nighttime.Data" = night.flies.compiled,
      "Plots" = output.plot
    )
    
  }
  return(final.output)
}