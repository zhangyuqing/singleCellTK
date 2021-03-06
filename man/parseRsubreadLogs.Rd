% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/alignSingleCellData.R
\name{parseRsubreadLogs}
\alias{parseRsubreadLogs}
\title{Parse Rsubread Logs for Mapping and Feature Count Statistics}
\usage{
parseRsubreadLogs(alignLog = NULL, featurecountLog = NULL,
  sampleName = NULL)
}
\arguments{
\item{alignLog}{Path to a log file created by the Rsubread align function}

\item{featurecountLog}{Path to a log file created by the Rsubread feature
count function}

\item{sampleName}{Sample name corresponding to the two log files}
}
\value{
A single line of a data frame with alignment and feature count
information
}
\description{
Parse Rsubread Logs for Mapping and Feature Count Statistics
}
