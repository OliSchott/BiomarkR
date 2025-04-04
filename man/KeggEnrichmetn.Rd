% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/R.R
\name{KeggEnrichmetn}
\alias{KeggEnrichmetn}
\title{KEGGEnrichment}
\usage{
KeggEnrichmetn(dataset, PoIs, folder = NA, plotname = "")
}
\arguments{
\item{dataset}{The dataset to be tested}

\item{PoIs}{A vector containing the Proteins of interest. Example: c("Q8TF72_SHROOM3" , "Q9ULZ3_PYCARD")}

\item{folder}{The folder to save the pathway plots}

\item{plotname}{The name to be displayed on created plots}
}
\value{
A list object containing the results of the KEGG enrichment analysis, the KEGG enrichment plot. The pictures of the pathways are saved in the specified folder
}
\description{
KEGG enrichment analysis for the specified dataset.
}
