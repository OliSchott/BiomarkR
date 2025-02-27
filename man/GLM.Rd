% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/R.R
\name{GLM}
\alias{GLM}
\title{GLM}
\usage{
GLM(dataset, PoIs, crossvalidation = F, plotname = "")
}
\arguments{
\item{dataset}{The dataset to be tested}

\item{PoIs}{A vector containing the Proteins of interest. Example: c(""Q8TF72_SHROOM3" , "Q9ULZ3_PYCARD") or (unique(dataset$Protein))}

\item{crossvalidation}{A boolean value indicating whether to use cross-validation> I recommend false for desriptive analysese and true for predictife tasks}

\item{plotname}{The name to be displayed on created plots}
}
\value{
A list object containing the results of the GLM model, the confusion matrix and the ROC plot and the AUC if the dataset has 2 classes
}
\description{
Generalized linear model for the specified dataset.
}
