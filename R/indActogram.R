#' Plot actogram of individual fly
#'
#' @description
#' Allows users to generate individual actograms. The input for this function must be the output of the binData() function. The output of this function is a plotly object.
#'
#' @param data Input data file. The input for this function must be the output of the function binData(). See ??binData().
#' @param bin Define the bin size (in minutes) in which input data is saved.
#' @param t.cycle Define the period of the environmental cycle or a single day in hours. This defaults to 24.
#' @param ind The channel number (or individual) whose actogram must be plotted.
#' @param key.acto Key for reactive input tables in the shiny app.
#'
#' @importFrom zoo rollapply
#' @importFrom plotly plot_ly add_trace layout %>% subplot
#' @importFrom grDevices rgb
#' @importFrom stats aggregate fitted lm na.omit sd
#' 
#' @return A \code{plotly} \code{htmlwidget} with the actogram of a user defined fly.
#'
#' @export indActogram
#'
#' @examples
#' td <- trimData(data = df, start.date = "19 Dec 20", start.time = "21:00",
#' n.days = 10, bin = 1, t.cycle = 24)
#' bd <- binData(td)
#' ind.actogram <- indActogram(data = bd, bin = 30, ind = 2)

indActogram <- function(data, bin = 30, t.cycle = 24, ind = 1, key.acto = 1) {

  requireNamespace("plotly")
  requireNamespace("zoo")
  
  if (requireNamespace("plotly", quietly = T)) {
    # library(plotly)
    # library(zoo)
    
    n.plot = 2
    s_per_day <- (60/bin)*t.cycle
    dummy <- matrix(0, nrow = s_per_day*(n.plot-1), ncol = 1)
    
    raw <- as.matrix(data[,c(1+ind)])
    data <- rbind(dummy, raw, dummy)
    
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
          color = "green",
          line = list(
            color = "green"
          )
        ),
        source = "actogram.phases",
        key = key.acto
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
            title = "Time index",
            linecolor = "black",
            autotick = FALSE,
            ticks = "outside",
            tick0 = 0,
            dtick = length(a[,1])/6,
            ticklen = 7,
            tickcolor = "black",
            range = c(0, (length(a[,1]+1))
            )
          )
        )
    }
    
    n.row = length(a[1,])
    
    plot.ind.actogram <- subplot(
      p,
      nrows = n.row,
      shareX = T,
      margin = 0.0, heights = rep(1/n.row, n.row)
    )%>%
      layout(
        showlegend = F,
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
          title = "Time index",
          linecolor = "black",
          autotick = FALSE,
          ticks = "outside",
          tick0 = 0,
          dtick = length(a[,1])/6,
          ticklen = 7,
          tickcolor = "black",
          range = c(0, (length(a[,1]+1))
          )
        )
      )
    
    return(plot.ind.actogram) 
  }
}
