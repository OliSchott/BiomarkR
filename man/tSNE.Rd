% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/R.R
\name{tSNE}
\alias{tSNE}
\title{tSNE}
\usage{
tSNE(dataset, plotname = "", show_ellipse = F, ColPalette = "custom_vibrant")
}
\arguments{
\item{dataset}{The dataset to be plotted}

\item{plotname}{The name to be displayed on created plots}

\item{show_ellipse}{Logical value indicating if the 95\% confidence ellipse should be plotted}
}
\value{
A list object containing the results of the tSNE calculations and the tSNE plot
}
\description{
t-distributed Stochastic Neighbor Embedding for the specified dataset.
}
