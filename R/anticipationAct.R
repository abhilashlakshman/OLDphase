#' Compute anticipation for activity data
#'
#' @description
#' This function calculates anticipation of activity in user-defined windows. The input for this function must be the output of the trimData() function. There are three choices for the estimation of anticipation. The "Slope" method simply finds best-fitting lines to activity data in the defined windows (by using linear regressions). The slope of this line, for each fly, is used as a measure of anticipation. The "Stoleru" method is based on the estimates used in the paper by Stoleru et al., 2004 (https://doi.org/10.1038/nature02926). The authors estimate progressive increase in activity before environmental transitions, scaled by the startle response (immediately post-transition). The "Harrisingh" method is based on the method used by Harrisingh et al., 2007 (10.1523/JNEUROSCI.3680-07.2007). In this method, the authors simply examine activity in a 3-hour interval before environmental transitions, scaled by activity in the 6-hour interval before the transitions. Only if method = "Slope", the output of this function is a list with three elements, i.e., plots of the regression line along with the activity in the morning window, plots of the regression line along with the activity in the evening window, and a table with the morning and evening anticipation values. In the other two cases, the output of this function is only a table with the morning and evening anticipation values. No plots will be generated.
#'
#' @param data Input data file. The input for this function must be the output of the function trimData(). See ??trimData().
#' @param method Choose the method for estimating anticipation. The three possible choices are, (i) "Slope", (ii) "Stoleru" and (iii) "Harrisingh". This defaults to "Slope". See Description for more details on these methods.
#' @param t.cycle Define the period of the environmental cycle or a single day in hours. This defaults to 24.
#' @param morn.win.start Define the start of morning window in ZT. ZT00 is the onset of zeitgeber; in case of light/dark cycles, ZT00 is the time at which lights turn ON. The end of the morning window is considered as ZT00.
#' @param eve.win.start Define the start of evening window in ZT.
#' @param eve.win.end Define the end of evening window.
#' @param rm.channels A vector of channels from a DAM monitor that must be discarded from analysis. If channels 1 to 5 must be removed, type in c(1:5). If channels 1 to 5 and 10 to 13 and 15 and 17 must be removed, type in c(1:5,10:13,15,17). Default is to include all individuals.
#' @param max.y.morn Set the upper limit of the y-axis in the plot for morning window.
#' @param max.y.eve Set the upper limit of the y-axis in the plot for evening window.
#'
#' @importFrom plotly plot_ly add_trace layout %>% subplot add_lines
#' @importFrom grDevices rgb
#' @importFrom stats aggregate fitted lm na.omit sd
#' 
#' @return If method = "Slope", a \code{list} with two items, else a \code{matrix} \code{array} with 32 rows (one for each fly) and 3 columns (Channel/Fly identity, Morning anticipation index and Evening anticipation index).
#' \describe{
#' If method = "Slope":
#' \item{Plot.morn}{A \code{plotly} \code{htmlwidget} with the anticipation estimates for the morning window.}
#' \item{Plot.eve}{A \code{plotly} \code{htmlwidget} with the anticipation estimates for the evening window.}
#' \item{Data}{A \code{matrix} \code{array} with 32 rows (one for each fly) and 3 columns (Channel/Fly identity, Morning anticipation index and Evening anticipation index).}
#' }
#'
#' @export anticipationAct
#'
#'
#' @examples
#' td <- trimData(data = df, start.date = "19 Dec 20", start.time = "21:00",
#' n.days = 10, bin = 1, t.cycle = 24)
#' anticip <- anticipationAct(data = td, method = "Stoleru", t.cycle = 24,
#'                            morn.win.start = 21,
#'                            eve.win.start = 9, eve.win.end = 12,
#'                            rm.channels = c())

anticipationAct <- function(data, method = "Slope", t.cycle = 24, morn.win.start, eve.win.start, eve.win.end, rm.channels = c(), max.y.morn = "auto", max.y.eve = "auto") {

  requireNamespace("plotly")
  # library(plotly)
  
  if (requireNamespace("plotly", quietly = T)) {
    if ("Slope" %in% method) {
      bd <- binData(data = data, input.bin = 1, output.bin = 15, t.cycle = t.cycle)
    } else if ("Stoleru" %in% method) {
      bd <- binData(data = data, input.bin = 1, output.bin = 60, t.cycle = t.cycle)
    } else if("Harrisingh" %in% method) {
      bd <- binData(data = data, input.bin = 1, output.bin = 60, t.cycle = t.cycle)
    }
    
    if (is.null(rm.channels)) {
      df <- bd
    } else {
      df <- bd[,-c(1+rm.channels)]
    }
    pre.averaged.vals <- aggregate(df, by = list(df[,"ZT"]), FUN = mean)
    averaged.vals <- pre.averaged.vals[,-c(1)]
    
    if ("Slope" %in% method) {
      
      df.ant.morn <- subset(averaged.vals, averaged.vals$ZT > morn.win.start - (15/60))
      df.ant.eve <- subset(averaged.vals, averaged.vals$ZT > eve.win.start - (15/60) & averaged.vals$ZT < eve.win.end + (15/60))
      
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
      
      if (max.y.morn == "auto") {
        p.morn <- plot_ly(
          x = df.ant.morn[,1],
          y = rowMeans(df.ant.morn[,-1]),
          type = "scatter",
          mode = "lines",
          line = list(
            color = "red",
            width = 2
          ),
          name = "Mean Activity"
        )%>%
          layout(
            showlegend = F,
            xaxis = list(
              showgrid = F,
              showline = T,
              titlefont = f1,
              tickfont = f2,
              title = "ZT (h)",
              linecolor = "black",
              # linewidth = 4,
              mirror = TRUE,
              autotick = TRUE,
              ticks = "inside",
              tickcolor = "black",
              # tickwidth = 4,
              range = c(df.ant.morn[1,1]-0.5, df.ant.morn[length(df.ant.morn[,1]),1]+0.5)
            ),
            yaxis = list(
              showgrid = F,
              showline = T,
              titlefont = f1,
              tickfont = f2,
              title = "Activity",
              linecolor = "black",
              # linewidth = 4,
              mirror = TRUE,
              autotick = TRUE,
              ticks = "inside",
              tickcolor = "black",
              range = c(0, max(c(rowMeans(df.ant.morn[,-1])))+5)
              # tickwidth = 4,
            )
          )
      } else {
        p.morn <- plot_ly(
          x = df.ant.morn[,1],
          y = rowMeans(df.ant.morn[,-1]),
          type = "scatter",
          mode = "lines",
          line = list(
            color = "red",
            width = 2
          ),
          name = "Mean Activity"
        )%>%
          layout(
            showlegend = F,
            xaxis = list(
              showgrid = F,
              showline = T,
              titlefont = f1,
              tickfont = f2,
              title = "ZT (h)",
              linecolor = "black",
              # linewidth = 4,
              mirror = TRUE,
              autotick = TRUE,
              ticks = "inside",
              tickcolor = "black",
              # tickwidth = 4,
              range = c(df.ant.morn[1,1]-0.5, df.ant.morn[length(df.ant.morn[,1]),1]+0.5)
            ),
            yaxis = list(
              showgrid = F,
              showline = T,
              titlefont = f1,
              tickfont = f2,
              title = "Activity",
              linecolor = "black",
              # linewidth = 4,
              mirror = TRUE,
              autotick = TRUE,
              ticks = "inside",
              tickcolor = "black",
              range = c(0, max.y.morn)
              # tickwidth = 4,
            )
          )
      }
      
      if (max.y.eve == "auto") {
        p.eve <- plot_ly(
          x = df.ant.eve[,1],
          y = rowMeans(df.ant.eve[,-1]),
          type = "scatter",
          mode = "lines",
          line = list(
            color = "red",
            width = 2
          ),
          name = "Mean Activity"
        )%>%
          layout(
            showlegend = F,
            xaxis = list(
              showgrid = F,
              showline = T,
              titlefont = f1,
              tickfont = f2,
              title = "ZT (h)",
              linecolor = "black",
              # linewidth = 4,
              mirror = TRUE,
              autotick = TRUE,
              ticks = "inside",
              tickcolor = "black",
              # tickwidth = 4,
              range = c(df.ant.eve[1,1]-0.5, df.ant.eve[length(df.ant.eve[,1]),1]+0.5)
            ),
            yaxis = list(
              showgrid = F,
              showline = T,
              titlefont = f1,
              tickfont = f2,
              title = "Activity",
              linecolor = "black",
              # linewidth = 4,
              mirror = TRUE,
              autotick = TRUE,
              ticks = "inside",
              tickcolor = "black",
              range = c(0, max(c(rowMeans(df.ant.eve[,-1])))+5)
              # tickwidth = 4,
            )
          )
      } else {
        p.eve <- plot_ly(
          x = df.ant.eve[,1],
          y = rowMeans(df.ant.eve[,-1]),
          type = "scatter",
          mode = "lines",
          line = list(
            color = "red",
            width = 2
          ),
          name = "Mean Activity"
        )%>%
          layout(
            showlegend = F,
            xaxis = list(
              showgrid = F,
              showline = T,
              titlefont = f1,
              tickfont = f2,
              title = "ZT (h)",
              linecolor = "black",
              # linewidth = 4,
              mirror = TRUE,
              autotick = TRUE,
              ticks = "inside",
              tickcolor = "black",
              # tickwidth = 4,
              range = c(df.ant.eve[1,1]-0.5, df.ant.eve[length(df.ant.eve[,1]),1]+0.5)
            ),
            yaxis = list(
              showgrid = F,
              showline = T,
              titlefont = f1,
              tickfont = f2,
              title = "Activity",
              linecolor = "black",
              # linewidth = 4,
              mirror = TRUE,
              autotick = TRUE,
              ticks = "inside",
              tickcolor = "black",
              range = c(0, max.y.eve)
              # tickwidth = 4,
            )
          )
      }
      
      
      
      anticipation <- matrix(NA, nrow = (length(averaged.vals[1,])-1), ncol = 3)
      colnames(anticipation) <- c("Channel", "Morning", "Evening")
      for.row.nam <- averaged.vals[,-1]
      nam <- colnames(for.row.nam)
      row.nam <- as.numeric(unlist(strsplit(nam, "I")))
      row.nam = row.nam[!is.na(row.nam)]
      
      for (i in 1:length(anticipation[,1])) {
        anticipation[i,"Channel"] = row.nam[i]
      }
      
      for (i in 2:length(averaged.vals[1,])) {
        morn.lin.mod <- lm(df.ant.morn[,i]~df.ant.morn[,"ZT"])
        eve.lin.mod <- lm(df.ant.eve[,i]~df.ant.eve[,"ZT"])
        
        p.morn <- p.morn%>%
          add_lines(
            x = df.ant.morn[,1],
            y = fitted(morn.lin.mod),
            line = list(
              color = rgb(0,0,0,0.5),
              width = 1
            ),
            name = colnames(df.ant.morn)[i]
          )
        
        p.eve <- p.eve%>%
          add_lines(
            x = df.ant.eve[,1],
            y = fitted(eve.lin.mod),
            line = list(
              color = rgb(0,0,0,0.5),
              width = 1
            ),
            name = colnames(df.ant.eve)[i]
          )
        
        anticipation[i-1,"Morning"] <- as.numeric(morn.lin.mod$coefficients[2])
        anticipation[i-1,"Evening"] <- as.numeric(eve.lin.mod$coefficients[2])
      }
      
      
      output <- list(
        "Plot.morn" = p.morn,
        "Plot.eve" = p.eve,
        "Data" = anticipation
      )
      
      return(output)
      
    } else if ("Stoleru" %in% method) {
      df.ant.morn <- rbind(
        subset(averaged.vals, averaged.vals$ZT > t.cycle - 3),
        averaged.vals[1,]
      )
      df.ant.eve <- subset(averaged.vals, averaged.vals$ZT > eve.win.end - 3 & averaged.vals$ZT < eve.win.end + 2)
      
      anticipation <- matrix(NA, nrow = (length(averaged.vals[1,])-1), ncol = 3)
      colnames(anticipation) <- c("Channel", "Morning", "Evening")
      for.row.nam <- averaged.vals[,-1]
      nam <- colnames(for.row.nam)
      row.nam <- as.numeric(unlist(strsplit(nam, "I")))
      row.nam = row.nam[!is.na(row.nam)]
      
      for (i in 1:length(anticipation[,1])) {
        anticipation[i,"Channel"] = row.nam[i]
      }
      
      for (i in 2:length(averaged.vals[1,])) {
        morn.ant <- (df.ant.morn[3,i] * ((df.ant.morn[3,i] - df.ant.morn[2,i]) * (df.ant.morn[2,i] - df.ant.morn[1,i])))/df.ant.morn[4,i]
        anticipation[i-1,"Morning"] <- as.numeric(morn.ant)
        
        eve.ant <- (df.ant.eve[3,i] * ((df.ant.eve[3,i] - df.ant.eve[2,i]) * (df.ant.eve[2,i] - df.ant.eve[1,i])))/df.ant.eve[4,i]
        anticipation[i-1,"Evening"] <- as.numeric(eve.ant)
      }
      
      return(anticipation)
      
    } else if ("Harrisingh" %in% method) {
      df.ant.morn <- subset(averaged.vals, averaged.vals$ZT > t.cycle - 6)
      df.ant.eve <- subset(averaged.vals, averaged.vals$ZT > eve.win.end - 6 & averaged.vals$ZT < eve.win.end + 1)
      
      anticipation <- matrix(NA, nrow = (length(averaged.vals[1,])-1), ncol = 3)
      colnames(anticipation) <- c("Channel", "Morning", "Evening")
      for.row.nam <- averaged.vals[,-1]
      nam <- colnames(for.row.nam)
      row.nam <- as.numeric(unlist(strsplit(nam, "I")))
      row.nam = row.nam[!is.na(row.nam)]
      
      for (i in 1:length(anticipation[,1])) {
        anticipation[i,"Channel"] = row.nam[i]
      }
      
      for (i in 2:length(averaged.vals[1,])) {
        morn.ant <- sum(df.ant.morn[4:6,i])/sum(df.ant.morn[,i])
        anticipation[i-1,"Morning"] <- as.numeric(morn.ant)
        
        eve.ant <- sum(df.ant.eve[4:6,i])/sum(df.ant.eve[,i])
        anticipation[i-1,"Evening"] <- as.numeric(eve.ant)
      }
      
      return(anticipation)
    }
  }

}
