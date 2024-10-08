% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/weights_scatterplot.R
\name{weights_scatterplot}
\alias{weights_scatterplot}
\title{Draw a comparative scatterplot for original weights and robust weights.}
\usage{
weights_scatterplot(
  robust_estimate,
  save_pdf = FALSE,
  file = "nak_weights",
  height = 9,
  width = 9
)
}
\arguments{
\item{robust_estimate}{The object to use for building the plot.}

\item{save_pdf}{If True, the scatterplot is saved and not displayed.}

\item{file}{If save_pdf is true, string giving the file path (including the file name).}

\item{height}{The height of the plot.}

\item{width}{The width of the plot.}
}
\description{
Draw a comparative scatterplot for original weights and robust weights.
}
