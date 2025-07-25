% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/R.R
\name{assign_colors}
\alias{assign_colors}
\title{assign_colors}
\usage{
assign_colors(labels, palette = "custom_vibrant")
}
\arguments{
\item{palette}{Color pellets c("custom_vibrant", "viridis", "plasma", "magma", "inferno")}

\item{unique_entries}{A vector of unique entries (like unique(dataset$Status))}
}
\value{
A color mapping
}
\description{
This function assigns colors from different pellets to a list of labels
}
