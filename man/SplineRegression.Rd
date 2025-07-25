% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/R.R
\name{SplineRegression}
\alias{SplineRegression}
\title{SplineRegression}
\usage{
SplineRegression(
  dataset,
  PoIs,
  Timecol,
  alpha = 0.05,
  split_status = FALSE,
  plotname = "",
  B = 100,
  col_pellet = "custom_vibrant",
  show_all_proteins = T
)
}
\arguments{
\item{dataset}{The dataset to be tested}

\item{PoIs}{A vector containing the Proteins of interest. Example: c(""Q8TF72_SHROOM3" , "Q9ULZ3_PYCARD") or (unique(dataset$Protein))}

\item{Timecol}{The column containing the time points. Example "GestationalAge"}

\item{alpha}{The confidence level for the confidence intervals}

\item{split_status}{A boolean indicating if the model should be split by status}

\item{plotname}{The name to be displayed on created plots}

\item{B}{The number of bootstrap samples to use for estimating confidence intervals}

\item{show_all_proteins}{A boolean indicating if all proteins should be shown on the plot}
}
\value{
A list object containing the results of the spline regression, the confidence intervals and the plot
}
\description{
Spline regression for the specified dataset.
}
