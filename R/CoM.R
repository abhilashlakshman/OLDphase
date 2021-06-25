#' Objectively quantify and visualise phases and calculate consolidation in pre-defined time-windows
#'
#' @description
#' This function calculates phase of center of mass of either activity or sleep data, and generates polar plots. The function also computes consolidation of activity or sleep in the defined windows of time. The output of this function is a list, the elements of which are a polar plot with the phases in each window for each fly, depicted, and a table which has the window specific times in ZT units and the consolidation values. For this function to work as expected, users must start the analysis at ZT00.
#'
#' @param input Input data file. The input for this function must be the output of the function binData(), sleepData() or ??wakeData(). See ??binData(), ??sleepData() and ??wakeData().
#' @param data If the input for this function is the output of binData(), sleepData() or wakeData(), then this must be "Activity", "Sleep" or "Wakefulness", respectively.
#' @param bin Intervals in which input data are saved (in minutes). This defaults to 30.
#' @param t.cycle Define the period of the environmental cycle or a single day in hours. This defaults to 24.
#' @param window A list of vectors, each defining the start and end ZT values of the window of interest. For instance, if users want to define the morning window as ZT21 to ZT03 and evening window as ZT09 to ZT15, then this should be list(c(21,3), c(9,15)). This defaults to list(c(18,6), c(6,18)). For the purposes of the polar plot, data from only the first two windows are used. But, for the data table, users can define as many windows as they like, and the function should work as expected.
#' @param rm.channels A vector of channels from a DAM monitor that must be discarded from analysis. If channels 1 to 5 must be removed, type in c(1:5). If channels 1 to 5 and 10 to 13 and 15 and 17 must be removed, type in c(1:5,10:13,15,17). Default is to include all individuals.
#'
#' @importFrom circular rho.circular circular
#' @importFrom plotly plot_ly add_trace layout %>% subplot
#' @importFrom grDevices rgb
#' @importFrom stats aggregate fitted lm na.omit sd
#'
#' @export CoM
#'
#'
#' @examples
#' td <- trimData(data = df, start.date = "19 Dec 20", start.time = "21:00",
#' n.days = 10, bin = 1, t.cycle = 24)
#' bd <- binData(td)
#' phase <- CoM(input = bd)

CoM <- function (input, data = "Activity", bin = 30, t.cycle = 24, window = list(c(21, 3), c(9, 15)), rm.channels = c()) {
  requireNamespace("plotly")
  requireNamespace("circular")
  
  if (requireNamespace("plotly", quietly = T)) {
    # library(plotly)
    # library(circular)
    
    if ("Activity" %in% data) {
      
      out <- list()
      
      for (i in 1:length(window)) {
        win.start <- window[[i]][1]
        win.end <- window[[i]][2]
        
        pre.pre.df <- profilesAct(data = input, bin = bin, t.cycle = t.cycle, average.type = "Days", rm.channels = rm.channels)
        pre.df <- pre.pre.df$Profiles[,c(1:33)]
        
        if (win.start > win.end) {
          df1 <- subset(pre.df, pre.df$ZT > win.start)
          df2 <- subset(pre.df, pre.df$ZT <= win.end)
          df <- rbind(df1, df2)
        } else {
          df <- subset(pre.df, pre.df$ZT > win.start & pre.df$ZT <= win.end)
        }
        
        win.size <- length(df[,1] * (bin/60))
        
        theta <- as.matrix(seq(0, length(df[,1]), by = (bin/60))*360/win.size)
        
        sin_theta <- as.matrix(sin(theta*(pi/180)))
        cos_theta <- as.matrix(cos(theta*(pi/180)))
        
        fsintheta <- matrix(0, nrow = length(df[,1]), ncol = (length(df[1,])-1))
        fcostheta <- matrix(0, nrow = length(df[,1]), ncol = (length(df[1,])-1))
        
        for (k in 1:length(df[,1])){
          for (j in 2:length(df[1,])){
            fsintheta[k,j-1] <- sin_theta[k,1]*df[k,j]
            fcostheta[k,j-1] <- cos_theta[k,1]*df[k,j]
          }
        }
        sum_fsintheta <- matrix(colSums(fsintheta), nrow = 1)
        sum_fcostheta <- matrix(colSums(fcostheta), nrow = 1)
        x = df[,-1]
        n_samples <- matrix(colSums(x), nrow = 1)
        X <- matrix(0, nrow = 1, ncol = length(x[1,]))
        Y <- matrix(0, nrow = 1, ncol = length(x[1,]))
        mean_theta <- matrix(0, nrow = 1, ncol = length(x[1,]))
        mean_r <- matrix(0, nrow = 1, ncol = length(x[1,]))
        
        for (l in 1:length(n_samples[1,])){
          X[1,l] <- sum_fcostheta[1,l]/n_samples[1,l]
          Y[1,l] <- sum_fsintheta[1,l]/n_samples[1,l]
          if (is.na(X[1,l])) {
            mean_theta[1,l] <- NA
            mean_r[1,l] <- NA
          } else if (X[1,l] < 0){
            mean_theta[1,l] <- (pi + atan(Y[1,l]/X[1,l])) * 180/pi
            mean_r[1,l] <- sqrt((X[1,l]^2) + (Y[1,l]^2))
          } else {
            mean_theta[1,l] <- (atan(Y[1,l]/X[1,l])) * 180/pi
            mean_r[1,l] <- sqrt((X[1,l]^2) + (Y[1,l]^2))
          }
        }
        
        for.out <- matrix(NA, nrow = 32, ncol = 2)
        
        for (a in 1:length(for.out[,1])) {
          if (is.na(mean_theta[1,a])) {
            for.out[a,1] <- NA
            for.out[a,2] <- NA
          } else if (mean_theta[1,a] > 0) {
            for.out[a,1] <- (mean_theta[1,a] * win.size/360) + (win.start + (bin/60))
            for.out[a,2] <- mean_r[1,a]
          } else if (mean_theta[1,a] < 0) {
            for.out[a,1] <- (mean_theta[1,a] * win.size/360) + win.end
            for.out[a,2] <- mean_r[1,a]
          }
        }
        
        out[[i]] <- for.out
        
      }
      output <- matrix(unlist(out), ncol = length(window) * 2, byrow = F)
      colnames(output) <- rep("nam", length(output[1,]))
      
      phase.index <- seq(1, length(output[1,]), by = 2)
      consol.index <- seq(2, length(output[1,]), by = 2)
      
      for (ii in 1:length(phase.index)) {
        colnames(output)[phase.index[ii]] <- paste("Phase_", "Window.", ii, sep = "")
      }
      for (ii in 1:length(consol.index)) {
        colnames(output)[consol.index[ii]] <- paste("Consolidation_", "Window.", ii, sep = "")
      }
      
      chan.num.col <- data.frame(
        "Channel" = seq(1, 32, by = 1)
      )
      
      new.out <- cbind(chan.num.col, output)
      
      p <- plot_ly(
        type = "scatterpolar",
        mode = "lines+markers"
      )
      p <- p%>%
        add_trace(
          r = rep(1, length(output[,1])),
          theta = output[,1]*360/t.cycle,
          mode = "markers",
          marker = list(
            color = "black"
          ),
          name = "Window-1"
        )%>%
        add_trace(
          r = rep(1, length(output[,1])),
          theta = output[,3]*360/t.cycle,
          mode = "markers",
          marker = list(
            color = "red"
          ),
          name = "Window-2"
        )%>%
        add_trace(
          r = c(0,
                circular::rho.circular(circular::circular(((output[,1]*360/t.cycle)*pi/180)), na.rm = T)),
          theta = c((mean(output[,1], na.rm = T))*360/t.cycle, (mean(output[,1], na.rm = T))*360/t.cycle),
          mode = "lines+markers",
          marker = list(
            color = rgb(0,0,0,0)
          ),
          line = list(
            color = "black",
            width = 2
          ),
          showlegend = F
        )%>%
        add_trace(
          r = c(0,
                circular::rho.circular(circular::circular(((output[,3]*360/t.cycle)*pi/180)), na.rm = T)),
          theta = c((mean(output[,3], na.rm = T))*360/t.cycle, (mean(output[,3], na.rm = T))*360/t.cycle),
          mode = "lines+markers",
          marker = list(
            color = rgb(0,0,0,0)
          ),
          line = list(
            color = "red",
            width = 2
          ),
          showlegend = F
        )%>%
        layout(
          polar = list(
            angularaxis = list(
              direction = 'clockwise',
              showgrid = T,
              gridcolor = rgb(0,0,0,0.4),
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
              gridcolor = rgb(0,0,0,0.4),
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
              range = c(0, 1.1)
            )
          )
        )
      
      final <- list(
        "Plot" = p,
        "Data" = new.out
      )
      
      return(final)
      
    } else if ("Sleep" %in% data) {
      out <- list()
      
      for (i in 1:length(window)) {
        win.start <- window[[i]][1]
        win.end <- window[[i]][2]
        
        pre.pre.df <- profilesSleep(data = input, bin = bin, t.cycle = t.cycle, average.type = "Days", rm.channels = rm.channels)
        pre.df <- pre.pre.df$Profiles[,c(1:33)]
        
        if (win.start > win.end) {
          df1 <- subset(pre.df, pre.df$ZT > win.start)
          df2 <- subset(pre.df, pre.df$ZT <= win.end)
          df <- rbind(df1, df2)
        } else {
          df <- subset(pre.df, pre.df$ZT > win.start & pre.df$ZT <= win.end)
        }
        
        win.size <- length(df[,1] * (bin/60))
        
        theta <- as.matrix(seq(0, length(df[,1]), by = (bin/60))*360/win.size)
        
        sin_theta <- as.matrix(sin(theta*(pi/180)))
        cos_theta <- as.matrix(cos(theta*(pi/180)))
        
        fsintheta <- matrix(0, nrow = length(df[,1]), ncol = (length(df[1,])-1))
        fcostheta <- matrix(0, nrow = length(df[,1]), ncol = (length(df[1,])-1))
        
        for (k in 1:length(df[,1])){
          for (j in 2:length(df[1,])){
            fsintheta[k,j-1] <- sin_theta[k,1]*df[k,j]
            fcostheta[k,j-1] <- cos_theta[k,1]*df[k,j]
          }
        }
        sum_fsintheta <- matrix(colSums(fsintheta), nrow = 1)
        sum_fcostheta <- matrix(colSums(fcostheta), nrow = 1)
        x = df[,-1]
        n_samples <- matrix(colSums(x), nrow = 1)
        X <- matrix(0, nrow = 1, ncol = length(x[1,]))
        Y <- matrix(0, nrow = 1, ncol = length(x[1,]))
        mean_theta <- matrix(0, nrow = 1, ncol = length(x[1,]))
        mean_r <- matrix(0, nrow = 1, ncol = length(x[1,]))
        
        for (l in 1:length(n_samples[1,])){
          X[1,l] <- sum_fcostheta[1,l]/n_samples[1,l]
          Y[1,l] <- sum_fsintheta[1,l]/n_samples[1,l]
          if (is.na(X[1,l])) {
            mean_theta[1,l] <- NA
            mean_r[1,l] <- NA
          } else if (X[1,l] < 0){
            mean_theta[1,l] <- (pi + atan(Y[1,l]/X[1,l])) * 180/pi
            mean_r[1,l] <- sqrt((X[1,l]^2) + (Y[1,l]^2))
          } else {
            mean_theta[1,l] <- (atan(Y[1,l]/X[1,l])) * 180/pi
            mean_r[1,l] <- sqrt((X[1,l]^2) + (Y[1,l]^2))
          }
        }
        
        for.out <- matrix(NA, nrow = 32, ncol = 2)
        
        for (a in 1:length(for.out[,1])) {
          if (is.na(mean_theta[1,a])) {
            for.out[a,1] <- NA
            for.out[a,2] <- NA
          } else if (mean_theta[1,a] > 0) {
            for.out[a,1] <- (mean_theta[1,a] * win.size/360) + (win.start + (bin/60))
            for.out[a,2] <- mean_r[1,a]
          } else if (mean_theta[1,a] < 0) {
            for.out[a,1] <- (mean_theta[1,a] * win.size/360) + win.end
            for.out[a,2] <- mean_r[1,a]
          }
        }
        
        out[[i]] <- for.out
        
      }
      output <- matrix(unlist(out), ncol = length(window) * 2, byrow = F)
      colnames(output) <- rep("nam", length(output[1,]))
      
      phase.index <- seq(1, length(output[1,]), by = 2)
      consol.index <- seq(2, length(output[1,]), by = 2)
      
      for (ii in 1:length(phase.index)) {
        colnames(output)[phase.index[ii]] <- paste("Phase_", "Window.", ii, sep = "")
      }
      for (ii in 1:length(consol.index)) {
        colnames(output)[consol.index[ii]] <- paste("Consolidation_", "Window.", ii, sep = "")
      }
      
      chan.num.col <- data.frame(
        "Channel" = seq(1, 32, by = 1)
      )
      
      new.out <- cbind(chan.num.col, output)
      
      p <- plot_ly(
        type = "scatterpolar",
        mode = "lines+markers"
      )
      p <- p%>%
        add_trace(
          r = rep(1, length(output[,1])),
          theta = output[,1]*360/t.cycle,
          mode = "markers",
          marker = list(
            color = "black"
          ),
          name = "Window-1"
        )%>%
        add_trace(
          r = rep(1, length(output[,1])),
          theta = output[,3]*360/t.cycle,
          mode = "markers",
          marker = list(
            color = "red"
          ),
          name = "Window-2"
        )%>%
        add_trace(
          r = c(0,
                circular::rho.circular(circular::circular(((output[,1]*360/t.cycle)*pi/180)), na.rm = T)),
          theta = c((mean(output[,1], na.rm = T))*360/t.cycle, (mean(output[,1], na.rm = T))*360/t.cycle),
          mode = "lines+markers",
          marker = list(
            color = rgb(0,0,0,0)
          ),
          line = list(
            color = "black",
            width = 2
          ),
          showlegend = F
        )%>%
        add_trace(
          r = c(0,
                circular::rho.circular(circular::circular(((output[,3]*360/t.cycle)*pi/180)), na.rm = T)),
          theta = c((mean(output[,3], na.rm = T))*360/t.cycle, (mean(output[,3], na.rm = T))*360/t.cycle),
          mode = "lines+markers",
          marker = list(
            color = rgb(0,0,0,0)
          ),
          line = list(
            color = "red",
            width = 2
          ),
          showlegend = F
        )%>%
        layout(
          polar = list(
            angularaxis = list(
              direction = 'clockwise',
              showgrid = T,
              gridcolor = rgb(0,0,0,0.4),
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
              gridcolor = rgb(0,0,0,0.4),
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
              range = c(0, 1.1)
            )
          )
        )
      
      final <- list(
        "Plot" = p,
        "Data" = new.out
      )
      
      return(final)
      
    } else if ("Wakefulness" %in% data) {
      out <- list()
      
      for (i in 1:length(window)) {
        win.start <- window[[i]][1]
        win.end <- window[[i]][2]
        
        pre.pre.df <- profilesWake(data = input, bin = bin, t.cycle = t.cycle, average.type = "Days", rm.channels = rm.channels)
        pre.df <- pre.pre.df$Profiles[,c(1:33)]
        
        if (win.start > win.end) {
          df1 <- subset(pre.df, pre.df$ZT > win.start)
          df2 <- subset(pre.df, pre.df$ZT <= win.end)
          df <- rbind(df1, df2)
        } else {
          df <- subset(pre.df, pre.df$ZT > win.start & pre.df$ZT <= win.end)
        }
        
        win.size <- length(df[,1] * (bin/60))
        
        theta <- as.matrix(seq(0, length(df[,1]), by = (bin/60))*360/win.size)
        
        sin_theta <- as.matrix(sin(theta*(pi/180)))
        cos_theta <- as.matrix(cos(theta*(pi/180)))
        
        fsintheta <- matrix(0, nrow = length(df[,1]), ncol = (length(df[1,])-1))
        fcostheta <- matrix(0, nrow = length(df[,1]), ncol = (length(df[1,])-1))
        
        for (k in 1:length(df[,1])){
          for (j in 2:length(df[1,])){
            fsintheta[k,j-1] <- sin_theta[k,1]*df[k,j]
            fcostheta[k,j-1] <- cos_theta[k,1]*df[k,j]
          }
        }
        sum_fsintheta <- matrix(colSums(fsintheta), nrow = 1)
        sum_fcostheta <- matrix(colSums(fcostheta), nrow = 1)
        x = df[,-1]
        n_samples <- matrix(colSums(x), nrow = 1)
        X <- matrix(0, nrow = 1, ncol = length(x[1,]))
        Y <- matrix(0, nrow = 1, ncol = length(x[1,]))
        mean_theta <- matrix(0, nrow = 1, ncol = length(x[1,]))
        mean_r <- matrix(0, nrow = 1, ncol = length(x[1,]))
        
        for (l in 1:length(n_samples[1,])){
          X[1,l] <- sum_fcostheta[1,l]/n_samples[1,l]
          Y[1,l] <- sum_fsintheta[1,l]/n_samples[1,l]
          if (is.na(X[1,l])) {
            mean_theta[1,l] <- NA
            mean_r[1,l] <- NA
          } else if (X[1,l] < 0){
            mean_theta[1,l] <- (pi + atan(Y[1,l]/X[1,l])) * 180/pi
            mean_r[1,l] <- sqrt((X[1,l]^2) + (Y[1,l]^2))
          } else {
            mean_theta[1,l] <- (atan(Y[1,l]/X[1,l])) * 180/pi
            mean_r[1,l] <- sqrt((X[1,l]^2) + (Y[1,l]^2))
          }
        }
        
        for.out <- matrix(NA, nrow = 32, ncol = 2)
        
        for (a in 1:length(for.out[,1])) {
          if (is.na(mean_theta[1,a])) {
            for.out[a,1] <- NA
            for.out[a,2] <- NA
          } else if (mean_theta[1,a] > 0) {
            for.out[a,1] <- (mean_theta[1,a] * win.size/360) + (win.start + (bin/60))
            for.out[a,2] <- mean_r[1,a]
          } else if (mean_theta[1,a] < 0) {
            for.out[a,1] <- (mean_theta[1,a] * win.size/360) + win.end
            for.out[a,2] <- mean_r[1,a]
          }
        }
        
        out[[i]] <- for.out
        
      }
      output <- matrix(unlist(out), ncol = length(window) * 2, byrow = F)
      colnames(output) <- rep("nam", length(output[1,]))
      
      phase.index <- seq(1, length(output[1,]), by = 2)
      consol.index <- seq(2, length(output[1,]), by = 2)
      
      for (ii in 1:length(phase.index)) {
        colnames(output)[phase.index[ii]] <- paste("Phase_", "Window.", ii, sep = "")
      }
      for (ii in 1:length(consol.index)) {
        colnames(output)[consol.index[ii]] <- paste("Consolidation_", "Window.", ii, sep = "")
      }
      
      chan.num.col <- data.frame(
        "Channel" = seq(1, 32, by = 1)
      )
      
      new.out <- cbind(chan.num.col, output)
      
      p <- plot_ly(
        type = "scatterpolar",
        mode = "lines+markers"
      )
      p <- p%>%
        add_trace(
          r = rep(1, length(output[,1])),
          theta = output[,1]*360/t.cycle,
          mode = "markers",
          marker = list(
            color = "black"
          ),
          name = "Window-1"
        )%>%
        add_trace(
          r = rep(1, length(output[,1])),
          theta = output[,3]*360/t.cycle,
          mode = "markers",
          marker = list(
            color = "red"
          ),
          name = "Window-2"
        )%>%
        add_trace(
          r = c(0,
                circular::rho.circular(circular::circular(((output[,1]*360/t.cycle)*pi/180)), na.rm = T)),
          theta = c((mean(output[,1], na.rm = T))*360/t.cycle, (mean(output[,1], na.rm = T))*360/t.cycle),
          mode = "lines+markers",
          marker = list(
            color = rgb(0,0,0,0)
          ),
          line = list(
            color = "black",
            width = 2
          ),
          showlegend = F
        )%>%
        add_trace(
          r = c(0,
                circular::rho.circular(circular::circular(((output[,3]*360/t.cycle)*pi/180)), na.rm = T)),
          theta = c((mean(output[,3], na.rm = T))*360/t.cycle, (mean(output[,3], na.rm = T))*360/t.cycle),
          mode = "lines+markers",
          marker = list(
            color = rgb(0,0,0,0)
          ),
          line = list(
            color = "red",
            width = 2
          ),
          showlegend = F
        )%>%
        layout(
          polar = list(
            angularaxis = list(
              direction = 'clockwise',
              showgrid = T,
              gridcolor = rgb(0,0,0,0.4),
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
              gridcolor = rgb(0,0,0,0.4),
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
              range = c(0, 1.1)
            )
          )
        )
      
      final <- list(
        "Plot" = p,
        "Data" = new.out
      )
      
      return(final)
      
    }
  }
}
