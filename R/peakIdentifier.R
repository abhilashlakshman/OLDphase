#' Phase identifier for activity data
#'
#' @description
#' This function generates a list of outputs. The first element is a plot of raw activity along with the smoothed data and objectively estimated phases of peak of activity data. Smoothing of activity is done using a Savitzky-Golay filter. The input of this function must be the output of the trimData function. This function requires the packages "plotly", "pracma" and "signal".
#'
#' @param data Input data file. The input for this function must be the output of the function trimData(). See ??trimData(). This assumes that the output from trimData has data binned in 1-min intervals.
#' @param filt.order The filter order. This defaults to 3.
#' @param filt.length The length of filter in indices. This defaults to 51.
#' @param min.peak.dist The minimum distance between peaks to be picked (in indices). This defaults to 100.
#' @param peak.ht.scal A scaling factor to identify peaks. This defaults to 0.5. A value of 0.5 will only find peaks that are at least half as tall as that of the tallest peak.
#' @param windows A list of vectors defining the start and end of the windows (in ZT) for morning and evening peaks, respectively around which the algorithm should look for peaks. This defaults to list(c(18,6), c(6,18)). The first window is treated as the window for morning peak, and the second as the window for evening peak.
#' @param rm.channels All the channels that users want to remove from their averaging. This must be a vector, i.e., channels must be separated by commas. For instance, if users choose to remove channels 1 to 5, 25 and 32, then the input should be either c(1,2,3,4,5,25,32) or c(1:5,25,32). This defaults to an empty vector, meaning no individuals are removed from analysis.
#'
#' @importFrom plotly plot_ly add_trace subplot %>% layout
#' @importFrom pracma findpeaks isempty
#' @importFrom signal sgolayfilt
#' @importFrom grDevices rgb
#' @importFrom stats aggregate fitted lm na.omit sd
#' 
#' @return A \code{list} with two items:
#' \describe{
#' \item{Plots}{A \code{plotly} \code{htmlwidget} with all the averaged activity overlayed with the smoothed data and markers to point out identified peaks in a 4-by-8 array.}
#' \item{Data}{A \code{matrix} \code{array} with 32 rows (one for each fly) and 3 columns (Channel/Fly identity, Morning peak phase and Evening peak phase (measured in ZT)).}
#' }
#' 
#'
#' @export peakIdentifier
#'
#' @examples
#' td <- trimData(data = df, start.date = "19 Dec 20", start.time = "21:00",
#' n.days = 3, bin = 1, t.cycle = 24)
#' pks <- peakIdentifier(data = td)

peakIdentifier <- function(data, filt.order = 3, filt.length = 51, min.peak.dist = 100, peak.ht.scal = 0.5, windows = list(c(18,6), c(18,6)), rm.channels = c()) {
  
  requireNamespace("plotly")
  requireNamespace("pracma")
  requireNamespace("signal")
  
  if (requireNamespace("plotly", quietly = T)) {
    # library(plotly)
    # library(pracma)
    # library(signal)
    
    bd <- binData(data = data, input.bin = 1, output.bin = 5, t.cycle = 24)
    
    pro <- profilesAct(data = bd, bin = 5, t.cycle = 24, average.type = "Days", rm.channels = rm.channels)
    
    pre.dat <- pro$Profiles
    dat <- rbind(subset(pre.dat, pre.dat$ZT > 18), subset(pre.dat, pre.dat$ZT < 18.01))
    
    p <- list()
    
    for (i in 1:32) {
      ind = i
      if (requireNamespace("signal", quietly = T)) {
        xx = signal::sgolayfilt(x = na.omit(dat[,1+ind]), p = filt.order, n = filt.length)
      }
      pks <- pracma::findpeaks(xx, minpeakdistance = min.peak.dist, minpeakheight = max(xx)*peak.ht.scal)
      
      
      p[[i]] <- plot_ly(
      )%>%
        add_trace(
          x = 1:length(dat[,1]),
          y = dat[,1+ind],
          type = "scatter",
          mode = "lines",
          line = list(
            color = "black",
            dash = "dash",
            width = 1
          )
        )%>%
        add_trace(
          x = 1:length(xx),
          y = xx,
          type = "scatter",
          mode = "lines",
          line = list(
            color = "red",
            dash = "solid",
            width = 2
          )
        )%>%
        add_trace(
          x = pks[,2],
          y = pks[,1]+3,
          type = "scatter",
          mode = "markers",
          marker = list(
            color = "blue",
            symbol = "triangle-down",
            size = 20
          )
        )
    }
    
    sp <- subplot(p, nrows = 4, shareX = T, shareY = T, margin = 0.01)%>%
      layout(
        showlegend = F
      )
    
    
    phase <- matrix(NA, nrow = 32, ncol = 3)
    colnames(phase) <- c("Channel", "Morning Peak", "Evening Peak")
    
    phase[1:32,"Channel"] <- 1:32
    
    for (i in 1:32) {
      ind = i
      if (requireNamespace("signal", quietly = T)) {
        xx = signal::sgolayfilt(x = na.omit(dat[,1+ind]), p = filt.order, n = filt.length)
      }
      pks <- pracma::findpeaks(xx, minpeakdistance = min.peak.dist, minpeakheight = max(xx)*peak.ht.scal)
      
      pks.zt <- dat[pks[,2],1]
      morn.peak.times <- windows[[1]]
      eve.peak.times <- windows[[2]]
      
      pks.morn <- pks.zt[pks.zt < morn.peak.times[2] | pks.zt > morn.peak.times[1]]
      pks.eve <- pks.zt[pks.zt < eve.peak.times[1] & pks.zt > eve.peak.times[2]]
      
      if (pracma::isempty(pks.eve)) {
        phase[i, "Evening Peak"] = NA
      } else {
        phase[i, "Evening Peak"] = pks.eve
      }
      
      if (pracma::isempty(pks.morn)) {
        phase[i, "Morning Peak"] = NA
      } else {
        phase[i, "Morning Peak"] = pks.morn
      }
    }
    
    output <- list(
      "Plots" = sp,
      "Data" = phase
    )
    
    return(output)
  }

}
