#' Generate actogram for "wake" data for individual fly
#'
#' @description
#' This function allows users to generate individual actograms for "wake" data. Input for this function must be an output from the trimData() function. The output of this function is a plotly object.
#' In a particular bin, sleep is calculated as the total minutes of inactivity equal to or greater than the defined threshold (sleep.def; typically, 5-minutes). Wake is defined as the total time spent by the fly not sleeping. See also wakeData().
#'
#' @param data Input data file. The input for this function must be the output of the function trimData(). See ??trimData().
#' @param sleep.def Definition of sleep. Traditionally, a single bout of sleep is defined as any duration of inactivity that is equal to or greater than 5-minutes. However, sometimes it may be of interest to examine longer bouts of sleep; sleep.def allows users to change the definition of sleep. This defaults to 5.
#' @param bin Intervals in which data are saved (in minutes). This defaults to 30. The value of bin cannot be lower than that of sleep.def.
#' @param t.cycle Define the period of the environmental cycle or a single day in hours. This defaults to 24.
#' @param ind The channel number (or individual) whose actogram must be plotted.
#'
#'
#' @export indActogramWake
#'
#' @examples
#' td <- trimData(data = df, start.date = "19 Dec 20", start.time = "21:00",
#' n.days = 10, bin = 1, t.cycle = 24)
#' ind.actogram.wake <- indActogramWake(data = td, ind = 7)
#' ind.actogram.wake <- indActogramWake(data = td, sleep.def = 20, bin = 60, t.cycle = 24, ind = 14)

indActogramWake <- function(data, sleep.def = 5, bin = 30, t.cycle = 24, ind = 1, key = "k.iaw") {

  library(plotly)
  library(zoo)

  n.plot = 2
  s_per_day <- (60/bin)*t.cycle
  dummy <- matrix(0, nrow = s_per_day*(n.plot-1), ncol = 1)

  pre.raw <- data[,-c(1:10)]
  raw <- as.matrix(pre.raw[,c(ind)])

  raw[,1:length(raw[1,])][raw[,1:length(raw[1,])] >= 1] <- -1
  raw[,1:length(raw[1,])][raw[,1:length(raw[1,])] == 0] <- 1
  raw[,1:length(raw[1,])][raw[,1:length(raw[1,])] == -1] <- 0

  binned_full_run.wake <- (length(raw[,1])/1440)*s_per_day
  wake <- matrix(NA, nrow = binned_full_run.wake, ncol = 1)
  index.wake <- seq(1, length(raw[,1]), by = bin)

  for (i in 1:length(index.wake)) {
    x <- raw[index.wake[i]:(index.wake[i]+bin-1),1]
    y <- rle(x)
    d_y <- as.data.frame(unclass(y))
    # dd_y <- subset(d_y, values == 1 & lengths >= sleep.def)
    dd_y <- subset(d_y, (values == 1 & lengths < sleep.def) | (values == 0))
    wake[i,1] <- sum(dd_y$lengths)
  }
  data <- rbind(dummy, wake, dummy)

  p <- list()

  f1 <- list(
    family = "Arial, sans-serif",
    size = 24,
    color = "black"
  )
  f2 <- list(
    family = "Arial, sans-serif",
    size = 20,
    color = "black"
  )

  a <- t(as.matrix(zoo::rollapply(data[,1],
                             width = s_per_day*n.plot,
                             by = s_per_day, as.numeric)))
  for (j in 1:length(a[1,])){
    p[[j]] <- plot_ly(
      # x = seq(0, ((length(a[,j])*(bin/60))-(bin/60)), by = bin/60),
      y = a[,j]/max(a[,j]),
      type = "bar",
      marker = list(
        color = "darkred",
        line = list(
          color = "darkred"
        )
      ),
      source = "actogram.phases",
      key = key
    )%>%
      layout(
        barmode = 'overlay',
        bargap = 0,
        yaxis = list(
          showticklabels = F,
          showline = T,
          showgrid = F,
          linecolor = "black"
        ),
        xaxis  = list(
          showgrid = F,
          showline = T,
          titlefont = f1,
          tickfont = f2,
          title = "",
          linecolor = "black",
          autotick = FALSE,
          ticks = "outside",
          tick0 = 0,
          dtick = length(a[,1])/6,
          ticklen = 7,
          tickcolor = toRGB("black"),
          range = c(0, (length(a[,1])+1))
        )
      )
  }
  plot.ind.wakogram <- subplot(
    p,
    nrows = length(a[1,]),
    shareX = T,
    margin = 0.0
  )%>%
    layout(
      showlegend = F,
      yaxis = list(
        showticklabels = F,
        showline = T,
        showgrid = F,
        linecolor = "black"
        # range = c(y.min, y.max)
      ),
      xaxis = list(
        showline = T,
        showgrid = F,
        titlefont = f1,
        tickfont = f2,
        title = "Time index",
        linecolor = "black",
        # linewidth = 4,
        # mirror = TRUE,
        autotick = FALSE,
        ticks = "outside",
        tick0 = 0,
        dtick = length(a[,1])/6,
        ticklen = 7,
        tickcolor = toRGB("black"),
        # tickwidth = 4,
        range = c(0, (length(a[,1])+1))
      )
    )

  return(plot.ind.wakogram)
}
