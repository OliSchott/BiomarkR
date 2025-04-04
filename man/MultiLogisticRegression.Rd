% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/R.R
\name{MultiLogisticRegression}
\alias{MultiLogisticRegression}
\title{MultiLogisticRegression}
\usage{
MultiLogisticRegression(dataset, PoIs, nIterations = 10)
}
\arguments{
\item{dataset}{The dataset to be tested}

\item{PoIs}{A vector containing the Proteins of interest. Example: c("Q8TF72_SHROOM3" , "Q9ULZ3_PYCARD") or (unique(dataset$Protein))}

\item{nIterations}{The number of iterations for the cross validation}
}
\value{
A list object containing the results of the logistic regression calculations and the ROC plot
}
\description{
Calculates the logistic regression for multiple features in the specified dataset.
}
