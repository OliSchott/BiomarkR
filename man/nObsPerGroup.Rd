% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/R.R
\name{nObsPerGroup}
\alias{nObsPerGroup}
\title{nObsPerGroup}
\usage{
nObsPerGroup(dataset, groupVar, n = 10)
}
\arguments{
\item{dataset}{The dataset to be filtered}

\item{groupVar}{The grouping variable (Input without quotation makrs)}

\item{n}{The number of observations}
}
\value{
The filtered dataset
}
\description{
This pipeline friendly function filters out every Protein with less than n observations in a specified grouping variable
}
