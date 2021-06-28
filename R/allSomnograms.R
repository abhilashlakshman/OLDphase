#' Generate actograms for sleep data (Somnograms)
#'
#' @description
#' This function generates a composite figure with actograms (referred to as somnograms, here) for all flies in a DAM scanned monitor file. Input for this function must be an output from the trimData() function. The output of this function is a large plotly object. This function requires the packages "plotly" and "zoo".
#' In a particular bin, sleep is calculated as the total minutes of inactivity equal to or greater than the defined threshold (sleep.def; typically, 5-minutes).
#'
#' @param data Input data file. The input for this function must be the output of the function trimData(). See ??trimData().
#' @param sleep.def Definition of sleep. Traditionally, a single bout of sleep is defined as any duration of inactivity that is equal to or greater than 5-minutes. However, sometimes it may be of interest to examine longer bouts of sleep; sleep.def allows users to change the definition of sleep. This defaults to 5.
#' @param bin Intervals in which data are saved (in minutes). This defaults to 30. The value of bin cannot be lower than that of sleep.def.
#' @param t.cycle Define the period of the environmental cycle or a single day in hours. This defaults to 24.
#'
#' @importFrom zoo rollapply
#' @importFrom plotly plot_ly add_trace layout %>% subplot
#' @importFrom grDevices rgb
#' @importFrom stats aggregate fitted lm na.omit sd
#'
#' @export allSomnograms
#'
#' @examples
#' td <- trimData(data = df, start.date = "19 Dec 20", start.time = "21:00",
#' n.days = 3, bin = 1, t.cycle = 24)
#' somnograms <- allSomnograms(data = td[,1:15])

allSomnograms <- function(data, sleep.def = 5, bin = 30, t.cycle = 24) {
  requireNamespace("plotly")
  requireNamespace("zoo")
  
  if (requireNamespace("plotly", quietly = T)) {
    # library(plotly)
    # library(zoo)
    
    n.plot = 2
    s_per_day <- (60/bin)*t.cycle
    dummy <- matrix(0, nrow = s_per_day*(n.plot-1), ncol = 32)
    
    raw <- data[,-c(1:10)]
    
    raw[,1:length(raw[1,])][raw[,1:length(raw[1,])] >= 1] <- -1
    raw[,1:length(raw[1,])][raw[,1:length(raw[1,])] == 0] <- 1
    raw[,1:length(raw[1,])][raw[,1:length(raw[1,])] == -1] <- 0
    
    binned_full_run.sleep <- (length(raw[,1])/1440)*s_per_day
    sleep <- matrix(NA, nrow = binned_full_run.sleep, ncol = 32)
    index.sleep <- seq(1, length(raw[,1]), by = bin)
    
    for (i in 1:length(index.sleep)) {
      for (j in 1:length(raw[1,])) {
        x <- raw[index.sleep[i]:(index.sleep[i]+bin-1),j]
        y <- rle(x)
        d_y <- as.data.frame(unclass(y))
        dd_y <- subset(d_y, d_y$values == 1 & d_y$lengths >= sleep.def)
        sleep[i,j] <- sum(dd_y$lengths)
      }
    }
    
    data <- rbind(dummy, sleep, dummy)
    
    plots <- list()
    p <- list()
    
    f1 <- list(
      family = "Arial, sans-serif",
      size = 24,
      color = "black"
    )
    f2 <- list(
      family = "Arial, sans-serif",
      # size = 18,
      color = "black"
    )
    
    for (i in 1:length(data[1,])){
      a <- t(as.matrix(zoo::rollapply(data[,i],
                                      width = s_per_day*n.plot,
                                      by = s_per_day, as.numeric)))
      for (j in 1:length(a[1,])){
        p[[j]] <- plot_ly(
          # x = seq(0, ((length(a[,j])*(bin/60))-(bin/60)), by = bin/60),
          y = a[,j]/max(a[,j]),
          type = "bar",
          marker = list(
            color = "darkblue",
            line = list(
              color = "darkblue"
            )
          )
        )%>%
          layout(
            barmode = 'overlay',
            bargap = 0,
            yaxis = list(
              showticklabels = F,
              showline = T,
              showgrid = F,
              linecolor = "black"
              # range = c(y.min, y.max)
            ),
            xaxis  = list(
              showgrid = F,
              showline = T,
              titlefont = f1,
              tickfont = f2,
              title = "",
              linecolor = "black",
              # linewidth = 4,
              # mirror = TRUE,
              autotick = FALSE,
              ticks = "outside",
              tick0 = 0,
              dtick = length(a[,1])/6,
              ticklen = 7,
              tickcolor = "black",
              # tickwidth = 4,
              range = c(0, (length(a[,1])+1))
            )
          )
      }
      plots[[i]] <- subplot(
        p,
        nrows = length(a[1,]),
        shareX = T,
        margin = 0.0
      )%>%
        layout(
          yaxis = list(
            showticklabels = F,
            showline = T,
            showgrid = F,
            linecolor = "black",
            range = c(0, 1.1)
          ),
          xaxis = list(
            showline = T,
            showgrid = F,
            titlefont = f1,
            tickfont = f2,
            title = "",
            linecolor = "black",
            # linewidth = 4,
            # mirror = TRUE,
            autotick = FALSE,
            ticks = "outside",
            tick0 = 0,
            dtick = length(a[,1])/6,
            ticklen = 7,
            tickcolor = "black",
            # tickwidth = 4,
            range = c(0, (length(a[,1])+1))
          )
        )
    }
    
    ann1 <- list(x = -0.4, y = 0.65,
                 text = "[1..8]", showarrow = F,
                 textangle = -90,
                 xref='paper', yref='paper',
                 font = f1
    )
    ann2 <- list(x = -0.4, y = 0.55,
                 text = "[9..16]", showarrow = F,
                 textangle = -90,
                 xref='paper', yref='paper',
                 font = f1
    )
    ann3 <- list(x = -0.4, y = 0.5,
                 text = "[17..24]", showarrow = F,
                 textangle = -90,
                 xref='paper', yref='paper',
                 font = f1
    )
    ann4 <- list(x = -0.4, y = 0.3,
                 text = "[25..32]", showarrow = F,
                 textangle = -90,
                 xref='paper', yref='paper',
                 font = f1
    )
    
    plots[[1]] <- plots[[1]]%>%
      layout(
        annotations = ann1
      )
    plots[[9]] <- plots[[9]]%>%
      layout(
        annotations = ann2
      )
    plots[[17]] <- plots[[17]]%>%
      layout(
        annotations = ann3
      )
    plots[[25]] <- plots[[25]]%>%
      layout(
        annotations = ann4
      )
    
    final <- subplot(
      plots,
      nrows = 4
    )%>%
      layout(showlegend = F,
             # autosize = F,
             # height = 900,
             # width = 1300,
             yaxis = list(
               showticklabels = F,
               showgrid = F,
               showline = T,
               linecolor = "black"
             ),
             xaxis = list(
               showgrid = F,
               showline = T,
               titlefont = f1,
               tickfont = f2,
               title = "",
               linecolor = "black",
               # linewidth = 4,
               # mirror = TRUE,
               autotick = FALSE,
               ticks = "outside",
               tick0 = 0,
               dtick = length(a[,1])/6,
               ticklen = 7,
               tickcolor = "black",
               # tickwidth = 4,
               range = c(0, (length(a[,1])+1))
             )
      )
    return(final)
  }
}
