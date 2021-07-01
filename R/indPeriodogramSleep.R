#' Periodogram plot for sleep data of individual flies
#'
#' @description
#' This function generates a periodogram for the sleep data of a single fly. The output of this function is a plotly object.
#'
#' @param data Input data file. The input for this function must be an output from one of either sleepData(), deepSleepData() or lightSleepData() functions. See ??sleepData(), ??deepSleepData() and/or ??lightSleepData().
#' @param bin Intervals in which input data is saved (in minutes). This defaults to 30.
#' @param method Choose the method for performing time-series analysis. Currently, three methods are implemented for analysis - "ChiSquare", "Autocorrelation", and "LombScargle". This defaults to "ChiSquare".
#' @param low.per Choose the lowest period (in hours) for analysis. This defaults to 16.
#' @param high.per Choose the highest period (in hours) for analysis. This defaults to 32.
#' @param alpha Choose the significance level for periodogram analysis. This defaults to 0.05.
#' @param time.res Resolution of periods (in minutes) to analyse while using the ChiSquare periodogram. For instance, if users wish to scan periods from low.per to high.per in the following manner: 16, 16.5, 17, 17.5, and so on, then time.res must be 30. This defaults to 20.
#' @param ind The channel number (or individual) whose periodogram must be plotted.
#'
#' @importFrom zeitgebr chi_sq_periodogram ac_periodogram ls_periodogram
#' @importFrom behavr hours mins
#' @importFrom plotly plot_ly add_trace layout %>% subplot
#' @importFrom grDevices rgb
#' @importFrom stats aggregate fitted lm na.omit sd
#' 
#' @return A \code{plotly} \code{htmlwidget} with the individual periodogram of a user defined fly.
#'
#' @export indPeriodogramSleep
#'
#' @examples
#' td <- trimData(data = df, start.date = "19 Dec 20", start.time = "21:00",
#' n.days = 3, bin = 1, t.cycle = 24)
#' sd <- sleepData(td)
#' ind.periodogram.sleep <- indPeriodogramSleep(data = sd, ind = 10)

indPeriodogramSleep <- function(data, bin = 30, method = "ChiSquare", low.per = 16, high.per = 32, alpha = 0.05, time.res = 20, ind = 1) {
  
  requireNamespace("zeitgebr")
  requireNamespace("behavr")
  requireNamespace("plotly")
  
  if (requireNamespace("plotly", quietly = T)) {
    # library(zeitgebr)
    # library(behavr)
    # library(plotly)
    
    raw <- data[,(1+ind)]
    sr <- 1/behavr::mins(bin)
    
    
    if ("ChiSquare" %in% method){
      x <- as.matrix(zeitgebr::chi_sq_periodogram(as.vector(raw),
                                                  period_range = as.integer(c(behavr::hours(low.per),
                                                                              behavr::hours(high.per))),
                                                  sampling_rate = sr, alpha = alpha,
                                                  time_resolution = behavr::mins(time.res)))
      
      adj.power <- matrix(x[,"power"] - x[,"signif_threshold"], ncol = 1)
      y <- cbind(x, adj.power)
      colnames(y)[5] <- c("adj.power")
      
      ind_data <- y
      
    } else if ("Autocorrelation" %in% method) {
      x <- as.matrix(zeitgebr::ac_periodogram(as.vector(raw),
                                              period_range = as.integer(c(behavr::hours(low.per),
                                                                          behavr::hours(high.per))),
                                              sampling_rate = sr, alpha = alpha))
      
      adj.power <- matrix(x[,"power"] - x[,"signif_threshold"], ncol = 1)
      y <- cbind(x, adj.power)
      colnames(y)[5] <- c("adj.power")
      
      ind_data <- y
      
    } else if ("LombScargle" %in% method) {
      x <- as.matrix(zeitgebr::ls_periodogram(as.vector(raw),
                                              period_range = as.integer(c(behavr::hours(low.per),
                                                                          behavr::hours(high.per))),
                                              sampling_rate = sr, alpha = alpha))
      
      adj.power <- matrix(x[,"power"] - x[,"signif_threshold"], ncol = 1)
      y <- cbind(x, adj.power)
      colnames(y)[5] <- c("adj.power")
      
      ind_data <- y
      
    }
    
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
    ax <- list(
      showgrid = F,
      showline = T,
      titlefont = f1,
      tickfont = f2,
      title = "Period (h)",
      linecolor = "black",
      linewidth = 4,
      mirror = TRUE,
      autotick = FALSE,
      ticks = "inside",
      tick0 = low.per,
      dtick = 4,
      ticklen = 7,
      tickcolor = "black",
      tickwidth = 4,
      range = c(low.per-1, high.per+1)
    )
    ay <- list(
      showgrid = F,
      showline = T,
      titlefont = f1,
      tickfont = f2,
      title = "Power",
      linecolor = "black",
      linewidth = 4,
      mirror = TRUE,
      autotick = TRUE,
      ticks = "inside",
      ticklen = 7,
      tickcolor = "black",
      tickwidth = 4
      # tick0 = 0,
      # dtick = max(ind_data$power)/6
    )
    
    ind.plot <- plot_ly(x = ind_data[,"period"]/3600,
                        y = ind_data[,"adj.power"],
                        name = paste("Channel: ", ind, sep = ""),
                        type = "scatter",
                        mode = "lines",
                        line = list(color = "black", width = 2),
                        source = "period_power_ind_act"
    )%>%
      add_trace(
        y = rep(0, length(ind_data[,1])),
        line = list(
          color = "red",
          width = 2
        )
      )%>%
      layout(
        showlegend = FALSE,
        xaxis = ax,
        yaxis = ay
      )
    
    return(ind.plot)
  }
}
