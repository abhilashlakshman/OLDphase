#' Periodogram analysis for sleep data
#'
#' @description
#' This function generates a composite figure with periodogram plots for all flies in a DAM scanned monitor file. Input for this function must be an output from the sleepData() function. The output of this function is a list with two components - (a) large plotly object with periodogram plots for all flies, and (b) a table which has channel wise information of significant period and adjusted power values, from the chosen time-series analysis method. This function requires the packages "plotly" and "zeitgebr".
#'
#' @param data Input data file. The input for this function must be an output from one of either sleepData(). See ??sleepData().
#' @param bin Intervals in which input data is saved (in minutes). This defaults to 30.
#' @param method Choose the method for performing time-series analysis. Currently, three methods are implemented for analysis - "ChiSquare", "Autocorrelation", and "LombScargle". This defaults to "ChiSquare".
#' @param low.per Choose the lowest period (in hours) for analysis. This defaults to 16.
#' @param high.per Choose the highest period (in hours) for analysis. This defaults to 32.
#' @param alpha Choose the significance level for periodogram analysis. This defaults to 0.05.
#' @param time.res Resolution of periods (in minutes) to analyse while using the ChiSquare periodogram. For instance, if users wish to scan periods from low.per to high.per in the following manner: 16, 16.5, 17, 17.5, and so on, then time.res must be 30. This defaults to 20.
#' 
#' @importFrom zeitgebr chi_sq_periodogram ac_periodogram ls_periodogram
#' @importFrom behavr hours mins
#' @importFrom plotly plot_ly add_trace layout %>% subplot
#' @importFrom grDevices rgb
#' @importFrom stats aggregate fitted lm na.omit sd
#' 
#' @return A \code{list} with two items:
#' \describe{
#' \item{Plots}{A \code{plotly} \code{htmlwidget} with all periodograms in a 4-by-8 array.}
#' \item{Data}{A \code{matrix} \code{array} with 32 rows (one for each fly) and 2 columns (Period and Adjusted Power).}
#' }
#'
#' @export allPeriodogramsSleep
#'
#' @examples
#' td <- trimData(data = df, start.date = "19 Dec 20", start.time = "21:00",
#' n.days = 5, bin = 1, t.cycle = 24)
#' sd <- sleepData(td)
#' all.periodograms.sleep <- allPeriodogramsSleep(data = sd[,1:6])

allPeriodogramsSleep <- function(data, bin = 30, method = "ChiSquare", low.per = 16, high.per = 32, alpha = 0.05, time.res = 20) {
  
  requireNamespace("zeitgebr")
  requireNamespace("plotly")

  if (requireNamespace("plotly", quietly = T)) {
    # library(zeitgebr)
    # library(plotly)
    
    raw <- data[,-1]
    
    sr <- 1/behavr::mins(bin)
    
    ind_data <- list()
    table.period.power <- matrix(NA, nrow = length(raw[1,]), ncol = 2)
    colnames(table.period.power) <- c("Period", "Adjusted.Power")
    row_names <- c()
    for (i in 1:length(raw[1,])) {
      row_names[i] <- paste("Ch-", i, sep = "")
    }
    rownames(table.period.power) <- row_names
    ind_plot <- list()
    
    for (i in 1:length(raw[1,])){
      if ("ChiSquare" %in% method){
        x <- as.matrix(zeitgebr::chi_sq_periodogram(raw[,i],
                                                    period_range = as.integer(c(behavr::hours(low.per),
                                                                                behavr::hours(high.per))),
                                                    sampling_rate = sr, alpha = alpha,
                                                    time_resolution = behavr::mins(time.res)))
        adj.power <- matrix(x[,"power"] - x[,"signif_threshold"], ncol = 1)
        y <- cbind(x, adj.power)
        colnames(y)[5] <- c("adj.power")
        
        yy <- subset(y, adj.power > 0)
        
        max_power <- which(yy[,"adj.power"] == max(yy[,"adj.power"]))
        
        if (is.null(max_power) | length(max_power) == 0) {
          table.period.power[i,"Period"] <- NA
          table.period.power[i,"Adjusted.Power"] <- NA
        } else if (length(max_power) > 1) {
          temp.period <- list()
          temp.adj.power <- list()
          for (j in 1:length(max_power)) {
            temp.period[[j]] <- yy[max_power[j],"period"]
            temp.adj.power[[j]] <- yy[max_power[j],"adj.power"]
          }
          table.period.power[i,"Period"] <- (mean(unlist(temp.period)))/(60*60)
          table.period.power[i,"Adjusted.Power"] <- mean(unlist(temp.adj.power))
        } else {
          table.period.power[i,"Period"] <- yy[max_power,"period"]/(60*60)
          table.period.power[i,"Adjusted.Power"] <- yy[max_power,"adj.power"]
        }
        
        ind_data[[i]] <- y
        
      } else if ("Autocorrelation" %in% method) {
        x <- as.matrix(zeitgebr::ac_periodogram(as.vector(raw[,i]),
                                                period_range = as.integer(c(behavr::hours(low.per),
                                                                            behavr::hours(high.per))),
                                                sampling_rate = sr, alpha = alpha))
        
        adj.power <- matrix(x[,"power"] - x[,"signif_threshold"], ncol = 1)
        y <- cbind(x, adj.power)
        colnames(y)[5] <- c("adj.power")
        
        yy <- subset(y, adj.power > 0)
        
        max_power <- which(yy[,"adj.power"] == max(yy[,"adj.power"]))
        
        if (is.null(max_power) | length(max_power) == 0) {
          table.period.power[i,"Period"] <- NA
          table.period.power[i,"Adjusted.Power"] <- NA
        } else if (length(max_power) > 1) {
          temp.period <- list()
          temp.adj.power <- list()
          for (j in 1:length(max_power)) {
            temp.period[[j]] <- yy[max_power[j],"period"]
            temp.adj.power[[j]] <- yy[max_power[j],"adj.power"]
          }
          table.period.power[i,"Period"] <- (mean(unlist(temp.period)))/(60*60)
          table.period.power[i,"Adjusted.Power"] <- mean(unlist(temp.adj.power))
        } else {
          table.period.power[i,"Period"] <- yy[max_power,"period"]/(60*60)
          table.period.power[i,"Adjusted.Power"] <- yy[max_power,"adj.power"]
        }
        
        ind_data[[i]] <- y
        
      } else if ("LombScargle" %in% method) {
        x <- as.matrix(zeitgebr::ls_periodogram(as.vector(raw[,i]),
                                                period_range = as.integer(c(behavr::hours(low.per),
                                                                            behavr::hours(high.per))),
                                                sampling_rate = sr, alpha = alpha))
        
        adj.power <- matrix(x[,"power"] - x[,"signif_threshold"], ncol = 1)
        y <- cbind(x, adj.power)
        colnames(y)[5] <- c("adj.power")
        
        yy <- subset(y, adj.power > 0)
        
        max_power <- which(yy[,"adj.power"] == max(yy[,"adj.power"]))
        
        if (is.null(max_power) | length(max_power) == 0) {
          table.period.power[i,"Period"] <- NA
          table.period.power[i,"Adjusted.Power"] <- NA
        } else if (length(max_power) > 1) {
          temp.period <- list()
          temp.adj.power <- list()
          for (j in 1:length(max_power)) {
            temp.period[[j]] <- yy[max_power[j],"period"]
            temp.adj.power[[j]] <- yy[max_power[j],"adj.power"]
          }
          table.period.power[i,"Period"] <- (mean(unlist(temp.period)))/(60*60)
          table.period.power[i,"Adjusted.Power"] <- mean(unlist(temp.adj.power))
        } else {
          table.period.power[i,"Period"] <- yy[max_power,"period"]/(60*60)
          table.period.power[i,"Adjusted.Power"] <- yy[max_power,"adj.power"]
        }
        
        ind_data[[i]] <- y
        
      }
    }
    
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
    ax <- list(
      showgrid = F,
      showline = T,
      titlefont = f1,
      tickfont = f2,
      title = "Period (h)",
      linecolor = "black",
      # linewidth = 4,
      # mirror = TRUE,
      autotick = FALSE,
      ticks = "inside",
      tick0 = low.per,
      dtick = 4,
      ticklen = 7,
      tickcolor = "black",
      # tickwidth = 4,
      range = c(low.per-1, high.per+1)
    )
    ay <- list(
      showgrid = F,
      showline = T,
      titlefont = f1,
      tickfont = f2,
      title = "Power",
      linecolor = "black",
      # linewidth = 4,
      # mirror = TRUE,
      autotick = TRUE,
      ticks = "inside",
      tick0 = 0,
      # dtick = max(table.period.power[,"Power"], na.rm = T)/6,
      ticklen = 7,
      tickcolor = "black"
      # range = c(0, max(table.period.power[,"Power"]))
      # tickwidth = 4,
    )
    
    for (i in 1:length(ind_data)) {
      
      a <- as.matrix(ind_data[[i]])
      
      ind_plot[[i]] <- plot_ly(x = a[,"period"]/3600,
                               y = a[,"adj.power"],
                               name = i,
                               type = "scatter",
                               mode = "lines",
                               line = list(color = "black", width = 2),
                               source = "period_power_act"
      )%>%
        add_trace(
          y = rep(0, length(a[,1])),
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
    }
    
    plot.periodograms <- subplot(ind_plot, nrows = 4, shareX = TRUE, shareY = TRUE)%>%
      layout(
        showlegend=FALSE
      )
    
    output <- list(
      "Plots" = plot.periodograms,
      "Data" = table.period.power
    )
    
    return(output)
  }
}
