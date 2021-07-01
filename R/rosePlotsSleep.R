#' Generate rose plots for averaged sleep data
#'
#' @description
#' Users can generate rose plots for the averaged sleep data. The input for this function must be output from sleepData(). The output of this function is a plotly object. By default the averaging of profiles is done over all flies and days.
#'
#' @param data Input data file. The input for this function must be the output of the function sleepData(). See ??sleepData().
#' @param bin Intervals in which data are saved (in minutes). This defaults to 30. This value must be the same as that for sleepData().
#' @param t.cycle Define the period of the environmental cycle or a single day in hours. This defaults to 24. This value must be the same as that for sleepData().
#' @param rm.channels All the channels that users want to remove from their averaging. This must be a vector, i.e., channels must be separated by commas. For instance, if users choose to remove channels 1 to 5, 25 and 32, then the input should be either c(1,2,3,4,5,25,32) or c(1:5,25,32). This defaults to an empty vector, meaning no individuals are removed from analysis.
#'
#' @importFrom plotly plot_ly add_trace layout %>% subplot
#' @importFrom grDevices rgb
#' @importFrom stats aggregate fitted lm na.omit sd
#' 
#' @return #' @return A \code{plotly} {htmlwidget} with rose plots for sleep data.
#'
#' @export rosePlotsSleep
#'
#' @examples
#' td <- trimData(data = df, start.date = "19 Dec 20", start.time = "21:00",
#' n.days = 1, bin = 1, t.cycle = 24)
#' sd <- sleepData(td)
#' r.plot <- rosePlotsSleep(data = sd)

rosePlotsSleep <- function (data, bin = 30, t.cycle = 24, rm.channels = c()) {

  requireNamespace("plotly")
  if (requireNamespace("plotly", quietly = T)) {
    # library(plotly)
    
    pre.df <- profilesSleep(data = data, bin = bin, t.cycle = t.cycle, average.type = "Both", rm.channels = rm.channels)
    
    df <- pre.df$Profiles
    
    theta <- as.matrix(df[,"ZT"]*360/t.cycle)
    
    p <- plot_ly()%>%
      layout(showlegend = FALSE,
             margin = list(
               t = 50,
               b = 50,
               l = 100,
               r = 100,
               pad = 2
             ),
             polar = list(
               angularaxis = list(
                 direction = 'clockwise',
                 showgrid = T,
                 gridcolor = rgb(0,0,0,1),
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
                 title = "",
                 linecolor = "black",
                 ticks = "outside",
                 ticklen = 7,
                 tickcolor = "black",
                 ticktext = as.list(
                   c("00", "03", "06", "09", "12", "15", "18", "21")
                 ),
                 tickvals = as.list(
                   seq(0, 320, by = 45)
                 ),
                 tickmode = "array"
               ),
               radialaxis = list(
                 showgrid = T,
                 gridcolor = rgb(0,0,0,0.7),
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
                 title = "",
                 linecolor = "black",
                 ticks = "outside",
                 ticklen = 7,
                 tickcolor = "black"
               )
             )
      )
    for (j in 1:length(df[,1])){
      x = df[j,"Mean"]
      y = theta[j,1]
      
      p <- add_trace(p,
                     r = c(0, x, x, 0),
                     theta = c(0, y-((bin/60)*360/t.cycle), y, 0),
                     type = 'scatterpolar',
                     mode = "lines",
                     line = list(color="black"),
                     fill = 'toself',
                     fillcolor = rgb(0,0,0.5,0.5)
      )
    }
    
    return(p)
  }
}
