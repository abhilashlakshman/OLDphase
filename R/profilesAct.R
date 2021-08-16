#' Compute and plot activity profiles
#'
#' @description
#' Users can compute average profiles and visualise the same. Averages can be performed either over Flies, Days, Both or None. Except in the case of "None" the output of this function will be a list with two elements. One is the generated plot and the other is a table with the activity values and corresponding standard errors of the mean.
#'
#' @param data Input data file. The input for this function must be the output of the function binData(). See ??binData().
#' @param bin Intervals in which data are saved (in minutes). This defaults to 30. This value must be the same as that for binData().
#' @param t.cycle Define the period of the environmental cycle or a single day in hours. This defaults to 24. This value must be the same as that for binData().
#' @param average.type Define how the averaging must be done for computing profiles. There are 4 choices; i) "Flies": This will average over all the flies for the number of days of data that there are, and produce an averaged time-series, ii) "Days": This will average activity over all the days but for each fly and provide 32 averaged plots, iii) "Both": This will average over both flies and days and provide one composite average profile, iv) "None": This will not average and produce any plots; the output will be the same as the input file. This defaults to "Both".
#' @param rm.channels All the channels that users want to remove from their averaging. This must be a vector, i.e., channels must be separated by commas. For instance, if users choose to remove channels 1 to 5, 25 and 32, then the input should be either c(1,2,3,4,5,25,32) or c(1:5,25,32). This defaults to an empty vector, meaning no individuals are removed from analysis.
#'
#' @importFrom plotly plot_ly add_trace layout %>% subplot
#' @importFrom grDevices rgb
#' @importFrom stats aggregate fitted lm na.omit sd
#' 
#' @return Except when average.type = "None", a \code{list} with two items. When average.type = "None", input file is returned.
#' \describe{
#' If average.type = "Days":
#' \item{Profiles}{
#' \describe{
#' \item{ZT}{Column with ZT values.}
#' \item{I1:I32}{Data averaged over days for each of 32 flies.}
#' \item{ZT}{Column with ZT values.}
#' \item{I1:I32}{SEM (across days) for each of 32 flies.}
#' }
#' }
#' \item{Plot}{A \code{plotly} \code{htmlwidget} with the activity profiles in a 4-by-8 array.}
#' 
#' If average.type = "Flies":
#' \item{Profiles}{
#' \describe{
#' \item{ZT}{Column with ZT values.}
#' \item{Mean}{Data averaged over all 32 flies for the entire duration of chosen days.}
#' \item{SEM}{SEM (across flies).}
#' }
#' }
#' \item{Plot}{A \code{plotly} \code{htmlwidget} with the activity time-series.}
#' 
#' If average.type = "Both":
#' \item{Profiles}{
#' \describe{
#' \item{ZT}{Column with ZT values.}
#' \item{Mean}{Data averaged over all days and all 32 flies.}
#' \item{SEM}{SEM (across flies).}
#' }
#' }
#' \item{Plot}{A \code{plotly} \code{htmlwidget} with the activity profile.}
#' }
#'
#' @export profilesAct
#'
#' @examples
#' td <- trimData(data = df, start.date = "19 Dec 20", start.time = "21:00",
#' n.days = 10, bin = 1, t.cycle = 24)
#' bd <- binData(td)
#' pro <- profilesAct(data = bd)
#' pro <- profilesAct(data = bd, bin = 60, t.cycle = 32,
#' average.type = "Days", rm.channels = c(1:5,25,32))

profilesAct <- function(data, bin = 30, t.cycle = 24, average.type = "Both", rm.channels = c()) {

  requireNamespace("plotly")
  
  if (requireNamespace("plotly", quietly = T)) {
    # library(plotly)
    
    if ("Flies" %in% average.type) {
      df <- data[,-c(1)]
      if (is.null(rm.channels)) {
        df <- df
      } else {
        df <- df[,-c(rm.channels)]
      }
      averaged.time.series <- as.matrix(rowMeans(df))
      sd.across.flies <- as.matrix(apply(df, 1, sd))
      sem <- sd.across.flies/sqrt(length(df[1,]))
      
      output <- cbind(data[,1], averaged.time.series, sem)
      colnames(output) <- c("ZT", "Mean", "SEM")
      
      t = table(data[,"ZT"])
      n.cyc <- as.numeric(t[[1]])
      p <- plot_ly(
        x = seq(1, length(output[,1])),
        y = output[,"Mean"],
        type = "scatter",
        mode = "lines",
        line = list(
          width = 2,
          dash = "solid",
          color = rgb(0,0,0,1)
        ),
        name = "Mean Time Series",
        error_y = list(
          array = output[,"SEM"],
          color = rgb(0,0,0,0.6),
          width = 0
        )
      )%>%
        layout(
          showlegend = F,
          xaxis = list(
            showgrid = F,
            showline = T,
            titlefont = list(
              family = "Arial, sans-serif",
              size = 20,
              color = "black"
            ),
            tickfont = list(
              family = "Arial, sans-serif",
              size = 14,
              color = "black"
            ),
            title = "Time (h)",
            linecolor = "black",
            mirror = TRUE,
            autotick = FALSE,
            ticks = "inside",
            tick0 = 0,
            dtick = 12*60/bin,
            ticklen = 7,
            tickcolor = "black",
            ticktext = as.list(
              c("0/24", rep(c("12", "0/24"), n.cyc))
            ),
            tickvals = as.list(
              c(0, seq(12*60/bin, length(output[,1]), by = 12*60/bin))
            ),
            tickmode = "array",
            range = c(0, length(output[,1])+1)
          ),
          yaxis = list(
            showgrid = F,
            showline = T,
            titlefont = list(
              family = "Arial, sans-serif",
              size = 20,
              color = "black"
            ),
            tickfont = list(
              family = "Arial, sans-serif",
              size = 14,
              color = "black"
            ),
            title = paste("Activity counts/", bin, "-min", sep = ""),
            linecolor = "black",
            mirror = TRUE,
            autotick = TRUE,
            ticks = "inside",
            tick0 = 0,
            dtick = max(output[,"Mean"])/6,
            ticklen = 7,
            tickcolor = "black"
            # range = c(0, max(output[,"Mean"]+output[,"SEM"])+5)
          )
        )
      out <- list(
        "Profiles" = output,
        "Plot" = p
      )
      return(out)
    } else if ("Days" %in% average.type) {
      if (is.null(rm.channels)) {
        df <- data
      } else {
        data[,c(1+rm.channels)] = NA
        df <- data
      }
      pre.averaged.vals <- aggregate(df, by = list(df[,"ZT"]), FUN = mean)
      averaged.vals <- pre.averaged.vals[,-c(1)]
      
      pre.sd.vals <- aggregate(df, by = list(df[,"ZT"]), FUN = sd)
      sd.vals <- pre.sd.vals[,-c(1,2)]
      t = table(data[,"ZT"])
      sem <- sd.vals/sqrt(as.numeric(t[[1]]))
      
      output <- cbind(averaged.vals, averaged.vals[,"ZT"], sem)
      colnames(output) <- c(
        colnames(averaged.vals),
        colnames(averaged.vals)
      )
      
      max.val <- max(output, na.rm = T)
      
      p <- list()
      for (i in 2:length(averaged.vals[1,])) {
        p[[i-1]] <- plot_ly(
          x = seq(1, length(output[,1])),
          y = output[,i],
          type = "scatter",
          mode = "lines",
          line = list(
            width = 2,
            dash = "solid",
            color = rgb(0,0,0,1)
          ),
          name = paste("Ch-", i-1, sep = ""),
          error_y = list(
            array = output[,(length(averaged.vals[1,])+i)],
            color = rgb(0,0,0,0.6),
            width = 0
          )
        )%>%
          layout(
            showlegend = F,
            xaxis = list(
              showgrid = F,
              showline = T,
              titlefont = list(
                family = "Arial, sans-serif",
                size = 20,
                color = "black"
              ),
              tickfont = list(
                family = "Arial, sans-serif",
                size = 14,
                color = "black"
              ),
              title = "Time (h)",
              linecolor = "black",
              mirror = F,
              autotick = FALSE,
              ticks = "inside",
              tick0 = 0,
              dtick = 12*60/bin,
              ticklen = 7,
              tickcolor = "black",
              ticktext = as.list(
                c("00", "06", "12", "18", "24")
              ),
              tickvals = as.list(
                c(0, seq(6*60/bin, length(output[,1]), by = 6*60/bin))
              ),
              tickmode = "array",
              range = c(0, length(output[,1])+1)
            ),
            yaxis = list(
              showgrid = F,
              showline = T,
              titlefont = list(
                family = "Arial, sans-serif",
                size = 20,
                color = "black"
              ),
              tickfont = list(
                family = "Arial, sans-serif",
                size = 14,
                color = "black"
              ),
              # title = paste("Activity counts/", bin, "-min", sep = ""),
              title = "Activity",
              linecolor = "black",
              mirror = F,
              autotick = TRUE,
              ticks = "inside",
              tick0 = 0,
              dtick = max.val/5,
              ticklen = 7,
              tickcolor = "black"
              # range = c(0, max.val+5)
            )
          )
      }
      pp <- subplot(p, nrows = 4, shareX = T, shareY = T, margin = 0.01)
      out <- list(
        "Profiles" = output,
        "Plot" = pp
      )
      return(out)
    } else if ("Both" %in% average.type) {
      s_per_day <- (60/bin)*t.cycle
      if (is.null(rm.channels)) {
        df <- data
      } else {
        df <- data[,-c(1+rm.channels)]
      }
      pre.averaged.vals <- aggregate(df, by = list(df[,"ZT"]), FUN = mean)
      averaged.vals <- pre.averaged.vals[,-c(1,2)]
      
      mean.profile <- as.matrix(rowMeans(averaged.vals))
      sd.across.flies <- as.matrix(apply(averaged.vals, 1, sd))
      sem <- sd.across.flies/sqrt(length(averaged.vals[1,]))
      
      output <- cbind(data[1:s_per_day,"ZT"], mean.profile, sem)
      colnames(output) <- c("ZT", "Mean", "SEM")
      
      p <- plot_ly(
        x = seq(1, length(output[,1])),
        y = output[,"Mean"],
        type = "scatter",
        mode = "lines",
        line = list(
          width = 2,
          dash = "solid",
          color = rgb(0,0,0,1)
        ),
        name = "Mean Profile",
        error_y = list(
          array = output[,"SEM"],
          color = rgb(0,0,0,0.6),
          width = 0
        )
      )%>%
        layout(
          showlegend = F,
          xaxis = list(
            showgrid = F,
            showline = T,
            titlefont = list(
              family = "Arial, sans-serif",
              size = 20,
              color = "black"
            ),
            tickfont = list(
              family = "Arial, sans-serif",
              size = 14,
              color = "black"
            ),
            title = "Time (h)",
            linecolor = "black",
            mirror = TRUE,
            autotick = FALSE,
            ticks = "inside",
            tick0 = 0,
            dtick = 6*60/bin,
            ticklen = 7,
            tickcolor = "black",
            ticktext = as.list(
              c("00", "06", "12", "18", "24")
            ),
            tickvals = as.list(
              c(0, seq(6*60/bin, length(output[,1]), by = 6*60/bin))
            ),
            tickmode = "array",
            range = c(0, length(output[,1])+1)
          ),
          yaxis = list(
            showgrid = F,
            showline = T,
            titlefont = list(
              family = "Arial, sans-serif",
              size = 20,
              color = "black"
            ),
            tickfont = list(
              family = "Arial, sans-serif",
              size = 14,
              color = "black"
            ),
            title = paste("Activity counts/", bin, "-min", sep = ""),
            linecolor = "black",
            mirror = TRUE,
            autotick = TRUE,
            ticks = "inside",
            tick0 = 0,
            dtick = max(output[,"Mean"])/6,
            ticklen = 7,
            tickcolor = "black"
            # range = c(0, max(output[,"Mean"]+output[,"SEM"])+5)
          )
        )
      
      out <- list(
        "Profiles" = output,
        "Plot" = p
      )
      return(out)
    } else if ("None" %in% average.type) {
      return(data)
    }
  }
}
