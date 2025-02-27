% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/R.R
\name{GBM}
\alias{GBM}
\title{GBM}
\usage{
GBM(dataset, PoIs, plotname = "")
}
\arguments{
\item{dataset}{The dataset to be tested}

\item{PoIs}{A vector containing the Proteins of interest. Example: c(""Q8TF72_SHROOM3" , "Q9ULZ3_PYCARD") or (unique(dataset$Protein))}

\item{plotname}{The name to be displayed on created plots}
}
\value{
A list object containing the results of the GBM model, the confusion matrix and the ROC plot and the AUC if the dataset has 2 classes
}
\description{
Gradient boosting machine for the specified dataset.
}
