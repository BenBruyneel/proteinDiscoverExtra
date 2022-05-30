library(dplyr)
library(lubridate)
library(pool)
library(grid)
library(proteinDiscover)
library(BBPersonalR)

#' helper function to retrieve the accession code (uniprot) from a short
#'  protein annotation (eg "His" or "Met6")
#'  
#' @param short protein short annotation
#' @param proteinInfo data.frame containing protein info like the one coming
#'  from the function \code{\link{knockoutProteins}}
#'  
#' @returns character vector
#' @export
shortToAccession <- function(short = "His4", proteinInfo = knockOutProteins()){
  return(proteinInfo$Accession[proteinInfo$short == short])
}

#' helper function to retrieve the short annotation from a (uniprot) protein
#'  code (eg P00815 for His4)
#'  
#' @param Accession protein Accession code
#' @param proteinInfo data.frame containing protein info like the one coming
#'  from the function \code{\link{knockoutProteins}}
#'  
#' @returns character vector
#' @export
accessionToShort <- function(Accession = "P00815", proteinInfo = knockOutProteins()){
  return(proteinInfo$short[proteinInfo$Accession == Accession])
}


#' function that calculates the average reporter signal to noise (ARSN) from
#'  the data in a psm table in a pdResult file from a TMT plex data
#'
#' @param db database access 'handle'
#' @param data default = NA, can be a data.frame or a numeric vector. If present
#'  then this data will be taken to calculate the arsn in stead of what's in the
#'  database
#' @param cutOff the level which is used to calculate the arsn, ie 10 means that
#'  we get the percentage of values that are above or equal to 10
#' @param removeNA if TRUE then prior to any calculation, the NA's are removed
#'  from the data
#'  
#' @returns a numeric value
#' @export
aboveARSN <- function(db = NA, data = NA, cutOff = 10, removeNA = TRUE){
  if (identical(db,NA) & identical(data,NA)){
    return(NA)
  }
  if (identical(data, NA)){
    data <- dbGetPsmTable(db = db, columnNames = "AverageReporterSN")
  }
  if (removeNA){
    data <- data[!is.na(data)]
  }
  return((sum(data >= cutOff)/length(data))*100)
}

#' function that plots the Average Reporter Signal to Noise (arsn) vs the
#'  retention time (rt) of all rows in a psms table in a pdResult file. The
#'  result will be a combined plot that shows a scatterplot together with a
#'  vertical density graph to show the diestribution of the arsn. Please note
#'  that the y-axis of the scatterplot and the x-axis of the density graph are
#'  logartihmic (log10)
#'
#' @param db database access 'handle'
#' @param cutOff the level which is used to draw the 'cutoff level', usually 10
#' @param pointColor defines the color of the border of the data points
#' @param pointFill defines the color of the data points themselves
#' @param pointAlpha alpha ('see through' value) of the data points
#' @param pointShape shape of the data points
#' @param pointSize size of the data points
#' @param title specifies the title
#' @param widths horizontal: two number (integer) vector specifying the amount
#'  of the plot to be used for the plots (horizontally)
#' @param drawPlot logical, if FALSE the plot will not be draw. Note that
#'  clearPlots() should be used before executing this function
#' @param returnData logical, if TRUE then a list with the following data wil be
#'  returned: dateTime, mean arsn, sd arsn, median arsn, mad arsn, above cutoff
#'  percentage, cutoff, graph (use grid.draw to display the graph)
#' @param dateTimeFormat default is NA, can be used to set the format of the
#'  dateTime item of the returned list (returnData = TRUE). If not NA, then
#'  should be a character vector specifying the format of the resulting
#'  POSIXct/POSIXt object. See \code{\link[base]{strptime}} for more info
#' @param returnRaw if TRUE (and if returnData = TRUE) then a data.frame will
#'  be added to the returned list. This data.frame contains the data for the
#'  scatterplot (retention time & arsn)
#'  
#' @returns a grid.draw object (gtable) or a list of information
#' @epxort
arsnPlot <- function(db, cutOff = 10,
                     pointAlpha = 0.75, pointColor = "black",
                     pointFill = "red", pointShape = 21,
                     pointSize = 2,
                     title = paste0("Average Reporter S/N vs Retention Time - ",
                                    as.character(
                                      ifelseProper(identical(dateTimeFormat,NA),
                                                   getAcquistionDateTime(db = db),
                                                   getAcquistionDateTime(db = db,
                                                                         format = dateTimeFormat)))),
                     widths = c(10,90),
                     drawPlot = TRUE,
                     returnData = FALSE,
                     dateTimeFormat = NA,
                     returnRaw = FALSE){
  tpsms <- dbGetPsmTable(db =db, columnNames = c("RetentionTime", "AverageReporterSN"))
  thePlot <- plotPlusMatrix(
    sPlot = scatterPlot(tpsms,
                        pointAlpha = pointAlpha, pointColor = pointColor,
                        pointFill = pointFill, pointShape = pointShape,
                        pointSize = pointSize,
                        xColumn = "RetentionTime", yColumn = "AverageReporterSN",
                        yLog = TRUE, xLabel = "Retention Time (min)",
                        yLabel = "Average Reporter S/N",
                        title = title) %>%
      lineMarks(hlines = cutOff,
                hlinesAttributes = lineAttributes(color = "blue",
                                                  linetype = "dotted", size = 0.4)),
    yPlot = statDensity(tpsms, column = "AverageReporterSN", xLog = TRUE,
                        vertical = TRUE, yExpand = c(0,0,0.025,0),
                        fillColor = "lightgray", outlineColor = "darkgray",
                        xAxis = FALSE, yAxis = FALSE, title = NULL,
                        gridLinesX = FALSE) %>%
      lineMarks(vlines = cutOff,
                vlinesAttributes = lineAttributes(color = "blue",
                                                  linetype = "dotted", size = 0.4)),
    widths = widths, heights = c(99,1))
    if (drawPlot){
      grid.draw(thePlot) # doing first clearPlot() doesn't work...
    }
  if (returnData){
    if (returnRaw){
      return(list(
        data = tpsms,
        dateTime = ifelseProper(identical(dateTimeFormat,NA),
                                getAcquistionDateTime(db = db),
                                getAcquistionDateTime(db = db,
                                                      format = dateTimeFormat)),
        meanARSN = mean(tpsms$AverageReporterSN, na.rm = TRUE),
        medianARSN = median(tpsms$AverageReporterSN, na.rm = TRUE),
        aboveCutoffPerc = aboveARSN(data = tpsms$AverageReporterSN,
                                    cutOff = cutOff),
        cutOff = cutOff,
        graph = thePlot
      ))
    } else {
      return(list(
        dateTime = ifelseProper(identical(dateTimeFormat,NA),
                                getAcquistionDateTime(db = db),
                                getAcquistionDateTime(db = db,
                                                      format = dateTimeFormat)),
        meanARSN = mean(tpsms$AverageReporterSN, na.rm = TRUE),
        sdARSN = sd(tpsms$AverageReporterSN, na.rm = TRUE),
        medianARSN = median(tpsms$AverageReporterSN, na.rm = TRUE),
        madARSN = mad(tpsms$AverageReporterSN, na.rm = TRUE),
        aboveARSNperc = aboveARSN(data = tpsms$AverageReporterSN,
                                  cutOff = cutOff),
        cutOff = cutOff,
        graph = thePlot
      ))
    }
  }
}

#' Essentially a wrapper around the plotPlusMatrix function to quickly generate
#'  plots from data in a pdResult file. These will be scatterplots with or
#'  without density plots for the x-axis and y-axis
#'  
#' @param db database access 'handle'
#' @param tableName table name of the data in the database
#' @param tData if NA, then data will re retrieved from the database, otherwise
#'  tData should be a data.frame (or similair) to use as a source for the data
#'  to be plotted
#' @param xColumn column name (or number of the column) for the x-axis
#' @param yColumn column name (or number of the column) for the y-axis
#' @param xLimits 2 element numeric vector specifying the minimum/maximum of the
#'  range of the x-axis
#' @param xLog if TRUE, then x-axis will be set to logarithmic
#' @param yLimits 2 element numeric vector specifying the minimum/maximum of the
#'  range of the y-axis
#' @param yLog if TRUE, then y-axis will be set to logarithmic
#' @param xLabel label for the x-axis
#' @param yLabel label for the y-axis
#' @param title specifies the title
#' @param xDensity if TRUE then a densityplot for the x-axis will be added
#' @param yDensity if TRUE then a densityplot for the y-axis will be added
#' @param widths horizontal: two number (integer) vector specifying the amount
#'  of the plot to be used for the plots (horizontally)
#' @param heights vertical: two number (integer) vector specifying the amount of
#'  the plot to be used for the plots (vertically)
#' @param xMarkers numeric vector specifying the vertical marker lines to be
#'  drawn
#' @param yMarkers numeric vector specifying the vertical marker lines to be
#'  drawn
#' @param pointColor defines the color of the border of the data points
#' @param pointFill defines the color of the data points themselves
#' @param pointAlpha alpha ('see through' value) of the data points
#' @param pointShape shape of the data points
#' @param pointSize size of the data points
#' @param widths horizontal: two number (integer) vector specifying the amount
#'  of the plot to be used for the plots (horizontally)
#' @param drawPlot logical, if FALSE the plot will not be draw. Note that
#'  clearPlots() should be used before executing this function
#' @param returnData logical, if TRUE then a list with the following data wil be
#'  returned: dateTime(s) of the data files & graph (use grid.draw to display
#'  the graph)
#' @param returnRaw if TRUE (and if returnData = TRUE) then a data.frame will
#'  be added to the returned list. This data.frame contains the data for the
#'  plot
#'  
#' @returns a grid.draw object (gtable) or a list of information
#' @epxort
pPlot <- function(db,
                  tableName = c("psms","peptides","proteins")[1],
                  tData = NA,
                  xColumn = "RetentionTime",
                  yColumn = "DeltaMassInPPM",
                  xLimits = NA, xLog = FALSE,
                  yLimits = NA, yLog = FALSE,
                  xLabel = NA, yLabel = NA,
                  title = NA,
                  xDensity = FALSE, yDensity =TRUE,
                  xWidths = c(10,90), yWidths = c(90,10),
                  xMarkers = NA, yMarkers = c(-5,0,5),
                  pointAlpha = 0.75, pointColor = "black",
                  pointFill = "red", pointShape = 21,
                  pointSize = 2,
                  drawPlot = TRUE,
                  returnData = FALSE,
                  returnRaw = FALSE){
  if (!identical(tData,NA)){
    tData <- tData %>% dplyr::select(all_of(c(xColumn, yColumn)))
  } else {
    if (tableName == "psms"){
      tData <- dbGetPsmTable(db = db, columnNames = c(xColumn, yColumn))
    } else {
      if (tableName == "peptides"){
        tData <- dbGetPeptideTable(db = db, columnNames = c(xColumn, yColumn))
      } else {
        if (tableName == "proteins"){
          tData <- dbGetProteinTable(db = db, columnNames = c(xColumn, yColumn))
        } else {
          stop("Invalid table was specified")
        }
      }
    }
  }
  if (!identical(xLimits,NA)){
    tData <- tData[tData[,xColumn] >= xLimits[1] & tData[,xColumn] <= xLimits[2], ]
  }
  if (!identical(yLimits,NA)){
    tData <- tData[tData[,yColumn] >= yLimits[1] & tData[,yColumn] <= yLimits[2], ]
  }
  sPlot <- scatterPlot(tData, xColumn = xColumn, yColumn = yColumn,
                       xLabel = ifelse(is.na(xLabel), xColumn, xLabel),
                       yLabel = ifelse(is.na(yLabel), yColumn, yLabel),
                       xDefault = ifelse(identical(xLimits,NA), TRUE, FALSE),
                       xLimits = xLimits, xLog = xLog,
                       yDefault = ifelse(identical(yLimits,NA), TRUE, FALSE),
                       yLimits = yLimits, yLog = yLog,
                       pointAlpha = pointAlpha, pointColor = pointColor,
                       pointFill = pointFill, pointShape = pointShape,
                       pointSize = pointSize,
                       title = ifelse(is.na(title),
                                      paste(c(yColumn," vs ",xColumn), collapse = ""),
                                      title))
  if (!identical(yMarkers,NA)){
    sPlot <- sPlot %>%
      lineMarks(hlines = yMarkers,
                hlinesAttributes = lineAttributes(color = "blue",
                                                  linetype = "dotted", size = 0.4))
  }
  if (!identical(xMarkers, NA)){
    sPlot <- sPlot %>%
      lineMarks(vlines = xMarkers,
                vlinesAttributes = lineAttributes(color = "blue",
                                                  linetype = "dotted", size = 0.4))
  }
  if (xDensity){
    xPlot <- statDensity(tData, column = xColumn,
                         vertical = FALSE, yExpand = c(0,0,0.025,0),
                         fillColor = "lightgray", outlineColor = "darkgray",
                         xAxis = FALSE, yAxis = FALSE, title = NULL,
                         xDefault = ifelse(identical(xLimits,NA), TRUE, FALSE),
                         xLimits = xLimits, xLog = xLog,
                         gridLinesX = FALSE)
    if (!identical(xMarkers, NA)){
      xPlot <- xPlot %>%
        lineMarks(vlines = xMarkers,
                  vlinesAttributes = lineAttributes(color = "blue",
                                                    linetype = "dotted", size = 0.4))
    }
  } else {
    xPlot <- clearPlot()
  }
  if (yDensity){
    yPlot <- statDensity(tData, column = yColumn,
                         vertical = TRUE, yExpand = c(0,0,0.025,0),
                         fillColor = "lightgray", outlineColor = "darkgray",
                         xAxis = FALSE, yAxis = FALSE, title = NULL,
                         xDefault = ifelse(identical(yLimits,NA), TRUE, FALSE),
                         xLimits = yLimits, xLog = yLog,
                         gridLinesX = FALSE)
    if (!identical(yMarkers,NA)){
      yPlot <- yPlot %>%
        lineMarks(vlines = yMarkers,
                  vlinesAttributes = lineAttributes(color = "blue",
                                                    linetype = "dotted", size = 0.4))
    }
  } else {
    yPlot <- clearPlot()
  }
  if (!xDensity){
    yWidths <- c(99,1)
  }
  if (!yDensity){
    xWidths <- c(1,99)
  }
  fullPlot <- plotPlusMatrix(sPlot = sPlot,
                            xPlot = xPlot, yPlot = yPlot,
                            widths = xWidths, heights = yWidths)
  if (drawPlot){
    grid.draw(fullPlot) # doing first clearPlot() doesn't work...
  }
  if (returnData){
    if (returnRaw){
      return(list(
        data = tData,
        date = getAcquistionDate(db = db),
        graph = fullPlot
      ))
    } else {
      return(list(
        date = getAcquistionDate(db = db),
        graph = fullPlot
      ))
    }
  }
}

