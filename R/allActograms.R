#' Generate actograms
#'
#' @description
#' This function generates a composite figure with actograms for all flies in a DAM scanned monitor file. Input for this function must be an output from the binData() function. The output of this function is a large plotly object. This function requires the packages "plotly" and "zoo".
#'
#' @param data Input data file. The input for this function must be the output of the function binData(). See ??binData().
#' @param bin Intervals in which data are saved (in minutes). This defaults to 30.
#' @param t.cycle Define the period of the environmental cycle or a single day in hours. This defaults to 24.
#' 
#' @importFrom zoo rollapply
#' @importFrom plotly plot_ly add_trace layout %>% subplot
#'
#'
#' @export allActograms
#'
#' @examples
#' td <- trimData(data = df, start.date = "19 Dec 20", start.time = "21:00",
#' n.days = 10, bin = 1, t.cycle = 24)
#' bd <- binData(data = td)
#' actograms <- allActograms(data = bd)


allActograms <- function(data, bin = 30, t.cycle = 24) {
  library(plotly)
  library(zoo)

  n.plot = 2
  s_per_day <- (60/bin)*t.cycle
  dummy <- matrix(0, nrow = s_per_day*(n.plot-1), ncol = 32)

  raw <- data[,-c(1)]
  data <- rbind(dummy, raw, dummy)

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
          color = "green",
          line = list(
            color = "green"
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
            tickcolor = toRGB("black"),
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
          linecolor = "black"
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
          tickcolor = toRGB("black"),
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
             tickcolor = toRGB("black"),
             # tickwidth = 4,
             range = c(0, (length(a[,1])+1))
           )
    )
  return(final)
}
