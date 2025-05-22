##Function for basic bottom up Biomarker discovery

## The point of these functions is to make the workflow easier and to provide
## sensible default values for quick (first pass) analysis of biomarker discovery experiments

## Loading dependencies
#' @import stats
#' @import tidyverse
#' @import readr
#' @import readxl
#' @import rstatix
#' @import ComplexHeatmap
#' @import InteractiveComplexHeatmap
#' @import pROC
#' @import caret
#' @import rstatix
#' @import stats
#' @import circlize
#' @import umap
#' @import stringi
#' @import gtools
#' @import seqinr
#' @import sva
#' @import dendextend
#' @import STRINGdb
#' @import MLmetrics
#' @import plotly
#' @import ggrepel
#' @import magrittr
#' @import igraph
#' @import Rtsne
#' @import grid
#' @import MEGENA
#' @import clusterProfiler
#' @import org.Hs.eg.db
#' @import pathview

## Data Import and Management
## add roxygen comments
#' @title SampleRandomization
#' @description This function designes the layout of a 96 well plate given a dataset and grouping variables whos distribution needs to be stable across plates
#' @param df The dataset to be randomized
#' @param group_vars The vector containing the grouping variables
#' @param IDColumn The column containing the sample IDs
#' @return A list containing the plate layout, summary statistics, and plots of the group_vars
#' @export
SampleRandomization <- function(df, group_vars, IDColumn) {

  df <- df

  nPlates <- ceiling(nrow(df) / 96)
  nSamplesPerPlate <- nrow(df) / nPlates

  # Create a list to store the results
  Results <- list()

  for (iteration in 1:1000) {

    ## randomly assign samples into groups no larger than 96
    random_split <- df %>%
      ## select necessary columns
      dplyr::select(group_vars) %>%
      ## make every column except ID column a factor
      dplyr::mutate(dplyr::across(where(is.character), as.factor)) %>%
      ## make the factors numeric
      dplyr::mutate(dplyr::across(where(is.factor), as.numeric)) %>%
      ## put back ID column
      dplyr::mutate(Sample = df[[IDColumn]]) %>%
      ## create random ID column
      dplyr::mutate(RandomID = sample(1:nrow(df), nrow(df), replace = FALSE)) %>%
      ## shuffle
      dplyr::arrange(-RandomID) %>%
      dplyr::select(-RandomID) %>%
      ## create a "Plate" column that increases by 1 every 96 rows
      dplyr::mutate(Plate = rep(1:ceiling(n() / nSamplesPerPlate), each = nSamplesPerPlate, length.out = n()))

    ## check significance between the grouping variables
    Significance <- data.frame(matrix(ncol = 2, nrow = 0))

    for (j in 1:length(group_vars)) {  # Use a different variable 'j' for the inner loop
      ANOVA_Results <- stats::anova(stats::lm(random_split[[group_vars[j]]] ~ Plate, data = random_split))
      p_value <- ANOVA_Results$`Pr(>F)`[1]

      ## add group variable and p-value into the Significance data.frame
      Significance <- rbind(Significance, data.frame(Variable = group_vars[j], p_value = p_value))
    }

    ## check if all p_values in Significance are < 0.025
    if (all(Significance$p_value > 0.5)) {
      print("Solution found!")
      ## assign solution
      final_splits <- random_split

      ## put final split into results
      Results$PlateLookup <- final_splits %>% dplyr::select(Sample, Plate)

      ## split into sub dataframes based on Plate
      final_splits <- split(final_splits, final_splits$Plate)
      ## stop loop
      break
    }

    ## stop execution if no solution found after 10000 iterations
    if (iteration == 1000) {
      stop("No solution found after 10k iterations")
    }
  }

  PlateLookup <- data.frame()

  # Create 96 well plate layout for each group
  for (i in 1:length(final_splits)) {

    ## Get the Study IDs of both lists
    Plate <- final_splits[[i]]$Sample %>%
      ## Shuffle the Study IDs
      sample() %>%
      data.frame(SampleID = .) %>%
      ## if less than 96 observations fill in NA until we reach a multiple of 96
      dplyr::add_row(SampleID = rep(NA_character_, 96 - nrow(.))) %>%
      ## Add row and column number to ensure A1, A2, A3,... order
      dplyr::mutate(
        row = rep(LETTERS[1:8], each = 12, length.out = n()),  # Row repeats every 12 samples
        column = rep(1:12, times = 8, length.out = n())        # Column cycles every 12 samples
      ) %>%
      ## add plate number
      dplyr::mutate(Plate = i)

    colnames(Plate)[1] <- "SampleID"

    ## append PlateLookup with plate
    PlateLookup <- rbind(PlateLookup, Plate)

    ## add PlateLookup to results
    Results$PlateLookup <- PlateLookup

    ## plot layout on plate
    Plot <- ggplot2::ggplot(Plate, ggplot2::aes(x = column, y = factor(row, levels = rev(unique(row))))) +  # Reverse row order
      ggplot2::geom_label(label = Plate$SampleID, size = 3) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 1)) +
      ggplot2::labs(title = paste("Plate", i), x = "Column", y = "Row") +
      ## make x axis display every value and put x-axis on top
      ggplot2::scale_x_continuous(breaks = 1:12, position = "top")

    ## store plot in results
    Results$Plots[paste("Plate", i)] <- list(Plot)

  }

  ## prepare data for summary statistics
  SummaryData <- dplyr::bind_rows(final_splits, .id = "Plate") %>%
    ### Check for non-numeric columns in final_splits and make them factors
    dplyr::mutate(dplyr::across(where(is.character), as.factor)) %>%
    ## make the factor columns numeric
    dplyr::mutate(dplyr::across(where(is.factor), as.numeric))

  ## create summary statistics for each plate
  Summary <- SummaryData %>%
    dplyr::group_by(Plate) %>%
    ## summarise the non-NA values for each column
    dplyr::summarise(mean = dplyr::across(dplyr::everything(), ~ mean(., na.rm = TRUE)),
                     sd = dplyr::across(dplyr::everything(), ~ sd(., na.rm = TRUE))) %>%
    ## data transformation
    t() %>%
    data.frame() %>%
    ## filter out rows that are entirely NA
    dplyr::filter(!is.na(X1))

  ## add summary to results
  Results$Summary <- Summary

  ## plot distributions of the group_vars
  for (i in 1:length(group_vars)) {

    ## Calculate significance metric
    Significance <- stats::anova(stats::lm(SummaryData[[group_vars[i]]] ~ Plate, data = SummaryData))

    p_value <- Significance$`Pr(>F)`[1]

    # Use sym to convert the string to a symbol, and !! to evaluate it in aes
    Histogram <- ggplot2::ggplot(SummaryData, ggplot2::aes(x = !!rlang::sym(group_vars[i]))) +
      ggplot2::facet_wrap(~Plate) +
      ggplot2::geom_histogram(binwidth = 1, ggplot2::aes(fill = factor(Plate))) +
      ## add minimal theme
      ggplot2::theme_minimal() +
      ## add boxes around the facets
      ggplot2::theme(panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = 1)) +
      ## add p-value to title
      ggplot2::ggtitle(paste("Distribution of", group_vars[i]), paste("p-value =", round(p_value, 3)))

    ## Add Histogram to Results under "Histograms"
    Results$Histograms[[group_vars[i]]] <- Histogram

    ## add significance to results
    Results$Significance[[group_vars[i]]] <- Significance
  }

  return(Results)
}


## add roxygen comments
#' @title GenerateSampleIDsFromFilePath
#' @description This function generates random sample IDs for a given dataset
#' @param filepath The path to the dataset
#' @return A dataset with random sample IDs
#' @export
GenerateSampleIDsFromFilePath <- function(filepath){
  ## reading excel sheet
  dataset <- readxl::read_excel(filepath)

  ## extract document name from filepath
  docname <- stringr::str_split(filepath, pattern = "\\\\") %>% data.frame()
  ## select last entry in docname
  docname <- docname[nrow(docname),1]
  ## delete file extension
  docname <- stringr::str_split_i(docname, pattern = "\\.", 1)

  ## Create an empty dataframe with the same structure as the input dataset
  SampleIDs <- data.frame(matrix(NA, nrow = nrow(dataset), ncol = 1))
  ## Genedate random String
  for (i in 1:nrow(dataset)) {
    ## Write random string and put two underscores in front and after it
    string <- paste0("__",stringi::stri_rand_strings(1, length = 5, pattern = "[A-Z]"),"__")
    ## Check if string has been used before
    if(!string %in% dataset){
      ## if not -> put in the data frame
      SampleIDs[i, 1] <- string
      ## if so -> create new string
    } else {i <- i-1}
  }

  ## Populate Sample frame
  SampleFrame <- cbind(SampleIDs,dataset)
  ## Rename column to "Sample ID"
  colnames(SampleFrame)[1] <- "Sample ID"
  ## write excel readable file into working directory
  write_excel_csv(SampleFrame, file = paste(docname,"with Sample","IDs"))
  ## confirm that everything went ok
  print(paste("Sample ID table safet to:", getwd()))
}
## Example: __AOGEL__

## Data Input and management

## ImportMSData
## Importing MS Data from DIA-NN output
## add roxygen comments
#' @title ImportMSData
#' @description This function imports MS data from DIA-NN output
#' @param filepath The path to the dataset
#' @param programm The program used to generate the data (for now only DIA-NN is supported)
#' @param SampleID A boolean value indicating whether the data has sample IDs generated by the function GenerateSampleIDsFromFilePath
#' @param feature The feature to be imported (Protein or Peptide)
#' @return A dataset with the imported data
#' @export
ImportMSData <- function(filepath, SampleID = FALSE, feature = "Protein", programm = "DIA-NN") {
  ## Filepath is the directory where the Output is stored
  filepath <- base::file.path(filepath)

  ## Reading DIA-NN output file
  if (programm == "DIA-NN" | programm == "dia-nn" | programm == "diann") {

    ## Importing Protein Data
    if (feature == "protein" | feature == "Protein") {
      ## Checking if SampleID was TRUE or FALSE
      if(SampleID == FALSE) {

        ## Import DIA-NN Datafile
        dataset <- readr::read_tsv(filepath)

        ## Cleaning data
        dataset <- dataset %>%
          dplyr::group_by(Protein.Group) %>%
          dplyr::mutate(Genes = ifelse(base::is.na(Genes), paste0("XYX", stringi::stri_rand_strings(1, length = 3, pattern = "[A-Z]")), Genes)) %>%
          dplyr::ungroup() %>%
          dplyr::mutate(Protein = base::paste0(Protein.Group, "_", Genes)) %>%
          dplyr::select(c("Protein", base::colnames(dataset)[6:base::ncol(dataset)]))

        ## Make the Proteins the Row Names
        dataset <- tibble::column_to_rownames(dataset, "Protein")

        ## Transpose, remove zeros, log2
        dataset <- base::data.frame(base::t(dataset))
        dataset[dataset == 0] <- NA
        dataset <- base::log2(dataset)

        ## Long format data
        dataset <- tibble::rownames_to_column(dataset, "Sample")
        dataset <- tidyr::gather(dataset, key = "Protein", value = "Intensity", base::colnames(dataset)[2:base::ncol(dataset)])

      }
      if (SampleID == TRUE) {

        ## Import DIA-NN Datafile
        dataset <- readr::read_tsv(filepath)

        ## Cleaning data
        dataset <- dataset %>%
          dplyr::group_by(Protein.Group) %>%
          dplyr::mutate(Genes = ifelse(base::is.na(Genes), paste0("XYX", stringi::stri_rand_strings(1, length = 3, pattern = "[A-Z]")), Genes)) %>%
          dplyr::ungroup() %>%
          dplyr::mutate(Protein = base::paste0(Protein.Group, "_", Genes)) %>%
          dplyr::select(c("Protein", base::colnames(dataset)[6:base::ncol(dataset)]))

        ## Make the Proteins the Row Names
        dataset <- tibble::column_to_rownames(dataset, "Protein")

        ## Transpose, remove zeros, log2
        dataset <- base::data.frame(base::t(dataset))
        dataset[dataset == 0] <- NA
        dataset <- base::log2(dataset)

        ## Long format data
        dataset <- tibble::rownames_to_column(dataset, "Sample")
        dataset <- tidyr::gather(dataset, key = "Protein", value = "Intensity", base::colnames(dataset)[2:base::ncol(dataset)])

        ## Parsing sample names
        dataset <- dataset %>%
          dplyr::mutate(Sample = stringr::str_split_i(Sample, pattern = "__", 2)) %>%
          dplyr::select("Sample", "Protein", "Intensity")

      }

      ## Checking if import was successful
      if (base::ncol(dataset) == 3 &
          base::sum(base::colnames(dataset) == c("Sample", "Protein", "Intensity")) == 3) {
        base::print("Import successful")
      } else {
        base::print("Something went wrong during data Import")
      }

    }

    ## Importing Peptide Data
    if(feature == "peptide" | feature == "Peptide") {
      ## Checking if SampleID was FALSE
      if(SampleID == FALSE) {

        ## Loading In Diann_Output
        dataset <- readr::read_tsv(filepath)

        ## Creating output dataframe
        dataset <- dataset %>%
          dplyr::group_by(Modified.Sequence) %>%
          ## if there is no gene make one up
          dplyr::mutate(Genes = ifelse(base::is.na(Genes), paste0("XYX", stringi::stri_rand_strings(1, length = 3, pattern = "[A-Z]")), Genes)) %>%
          dplyr::ungroup() %>%
          dplyr::mutate(Peptide = base::paste0(Protein.Group,"_", Modified.Sequence, "_", Genes)) %>%
          dplyr::select(-c("Protein.Names", "First.Protein.Description", "Proteotypic", "Stripped.Sequence", "Precursor.Id", "Protein.Group", "Protein.Ids")) %>%
          dplyr::select(-c(Modified.Sequence, Genes)) %>%
          tidyr::pivot_longer(names_to = "Sample", values_to = "Intensity", cols = -c("Precursor.Charge", "Peptide")) %>%
          dplyr::mutate(Peptide = base::paste0(Peptide, "_", Precursor.Charge)) %>%
          dplyr::select("Sample", "Peptide", "Intensity") %>%
          dplyr::mutate(Intensity = ifelse(Intensity == 0, NA, Intensity)) %>%
          dplyr::mutate(Intensity = base::log2(Intensity))
      }

      ## Checking if import was successful
      if (base::ncol(dataset) == 3 &
          base::sum(base::colnames(dataset) == c("Sample", "Peptide", "Intensity")) == 3) {
        base::print("Import successful")
      } else {
        base::print("Something went wrong during data Import")
      }

    }

    return(base::invisible(dataset))

  }

  ## Checking if import was successful
  if (base::ncol(dataset) == 3 &
      base::sum(base::colnames(dataset) == c("Sample", "Protein", "Intensity")) == 3) {
    base::print("Import successful")
  } else {
    base::print("Something went wrong during data Import")
  }

  return(base::invisible(dataset))

}


## only DIA-NN Data can be imported ATM

## CORE FUNCTIONS
## Ranked Intensities
## Optional Function that plots the ranked mean intensities for each Protein
## add roxygen comments
#' @title RankedIntensities
#' @description This function plots the ranked mean intensities for each Protein
#' @param dataset The dataset to be plotted
#' @param plotname The name to be displayed on created plots
#' @return A list containing the log2 and normal plot and a table with the ranked intensities
#' @export
RankedIntensities <- function(dataset, plotname = ""){
  if("Protein" %in% colnames(dataset)){
    ## Calculating mean intensities for each protein
    rankedIntensities <- dataset %>%
      group_by(Protein) %>%
      summarise(meanInt = mean(Intensity, na.rm = TRUE)) %>%
      arrange(-meanInt)

      ## Create plot with normal Intensity scale
      IntensityPlot <- rankedIntensities %>%
        ggplot(mapping = aes(x = reorder(Protein,- meanInt), y = (2^meanInt),colour = Protein))+
        geom_point() +
        xlab(paste("Protein Rank")) +
        ylab("mean Intensity") +
        theme(legend.position = "none")+
        ggtitle(paste("Ranked Protein intensities",plotname)) +
        theme(axis.text.x = element_blank())  # Hide x-axis text labels

      ## Creating plot with logarithmic Intensity scale
      Log2IntensityPlot <- rankedIntensities %>%
        ggplot(mapping = aes(x = reorder(Protein,- meanInt), y = meanInt,colour = Protein))+
        geom_point() +
        xlab(paste("Peptide Rank")) +
        ylab("log2 of mean Intensity") +
        theme(legend.position = "none")+
        ggtitle(paste("Ranked Protein intensities", plotname)) +
        theme(axis.text.x = element_blank()) # Hide x-axis text labels

      ## Print the number of unique proteins in the dataset
      nProteins <- dataset %>% distinct(Protein) %>% nrow()
      print(paste(nProteins, "unique Proteins Identified"))

  }
  if("Peptide" %in% colnames(dataset)){

    ## Calculating mean intensities for each protein
    rankedIntensities <<- dataset %>%
      group_by(Peptide) %>%
      summarise(meanInt = mean(Intensity, na.rm = TRUE)) %>%
      arrange(-meanInt)



      ## Create plot with normal Intensity scale
      Log2IntensityPlot <- rankedIntensities %>%
        ggplot(mapping = aes(x = reorder(Peptide,- meanInt), y = (2^meanInt),colour = Peptide))+
        geom_point() +
        xlab(paste(feature, "Rank")) +
        ylab("mean Intensity") +
        theme(legend.position = "none")+
        ggtitle(paste("Ranked Peptide intensities", plotname)) +
        theme(axis.text.x = element_blank())  # Hide x-axis text labels

      plot(IntensityPlot)

      ## Creating plot with logarithmic Intensity scale
      IntensityPlot <- rankedIntensities %>%
        ggplot(mapping = aes(x = reorder(Peptide,- meanInt), y = meanInt,colour = Peptide))+
        geom_point() +
        xlab(paste(feature, "Rank")) +
        ylab("log2 of mean Intensity") +
        theme(legend.position = "none")+
        ggtitle(paste("Ranked Peptide intensities",plotname)) +
        theme(axis.text.x = element_blank()) # Hide x-axis text labels

      plot(IntensityPlot)

      ## Print the number of unique proteins in the dataset
      nProteins <- dataset %>% distinct(Peptide) %>% nrow()
      print(paste(nProteins, "unique Proteins Identified"))
  }

  ## Generating Output object

  Output <- list()
  Output$log2Plot <- Log2IntensityPlot
  Output$Intensityplot <- IntensityPlot
  Output$Table <- rankedIntensities


  ## Return Output
  return(Output)

}

## NA CutOff
## Get rid of every protein that has more than x % missing values (default = 70)
## add roxygen comments
#' @title NaCutoff
#' @description This pipeline friendly function gets rid of every protein that has more than x % missing values (default = 70)
#' @param dataset The dataset to be filtered
#' @param cutoffvalue The cutoff value for the percentage of missing values
#' @return A dataset with the filtered proteins
#' @export
NaCutoff <- function(dataset, cutoffvalue = 70) {


  if ("Protein" %in% colnames(dataset)) {

    dataset <- dataset %>%
      tidyr::pivot_wider(names_from = "Protein", values_from = "Intensity") %>%
      tidyr::pivot_longer(cols = contains("_"), names_to = "Protein", values_to = "Intensity")

    Results <- dataset %>%
      dplyr::group_by(Protein) %>%
      dplyr::mutate(ISNA = is.na(Intensity)) %>%
      dplyr::mutate(SUMNA = sum(ISNA)) %>%
      dplyr::mutate(ProtTotal = sum(ISNA) + sum(!ISNA)) %>%
      dplyr::mutate(PercentMissing = (SUMNA / ProtTotal) * 100) %>%
      dplyr::filter((100 - PercentMissing) >= cutoffvalue) %>%
      dplyr::select(-c(ISNA, SUMNA, ProtTotal, PercentMissing))

  }
  if ("Peptide" %in% colnames(dataset)) {

    dataset <- dataset %>%
      tidyr::pivot_wider(names_from = "Peptide", values_from = "Intensity") %>%
      tidyr::pivot_longer(cols = contains("_"), names_to = "Peptide", values_to = "Intensity")


    Results <- dataset %>%
      dplyr::group_by(Peptide) %>%
      dplyr::mutate(ISNA = is.na(Intensity)) %>%
      dplyr::mutate(SUMNA = sum(ISNA)) %>%
      dplyr::mutate(ProtTotal = sum(ISNA) + sum(!ISNA)) %>%
      dplyr::mutate(PercentMissing = (SUMNA / ProtTotal) * 100) %>%
      dplyr::filter((100 - PercentMissing) >= cutoffvalue) %>%
      dplyr::select(-c(ISNA, SUMNA, ProtTotal, PercentMissing))
    }
  return(Results)
}


## Pipeline functions
## The pipeline functions are a collection of functions that are for easy data manipulation

## filter out every Protein with less than n observations total
## add roxygen comments
#' @title nObsPerFeature
#' @description This pipeline friendly function filters out every Protein with less than n observations total
#' @param dataset The dataset to be filtered
#' @param n The number of observations
#' @return The filtered dataset
#' @export
nObsPerFeature <- function(dataset, n = 10){
  if("Protein" %in% colnames(dataset)){

    results <- dataset %>%
      dplyr::group_by(Protein) %>%
      ## column with the number of observations of each Protein
      dplyr::mutate(nObs = n()) %>%
      ## nunmber of NA Observations of each Protein
      dplyr::mutate(NAObs = sum(is.na(Intensity))) %>%
      ## nunber of valid entires
      dplyr::mutate(ValidObs = nObs - NAObs) %>%
      ## filtering Proteins that have less than 5 Observation
      dplyr::filter(ValidObs > n) %>%
      ## Removing unnecessary columns
      dplyr::select(-c("nObs", "NAObs", "ValidObs"))
  }
  if("Peptide" %in% colnames(dataset)){

    results <- dataset %>%
      dplyr::group_by(Peptide) %>%
      ## column with the number of observations of each Protein
      dplyr::mutate(nObs = n()) %>%
      ## nunmber of NA Observations of each Protein
      dplyr::mutate(NAObs = sum(is.na(Intensity))) %>%
      ## nunber of valid entires
      dplyr::mutate(ValidObs = nObs - NAObs) %>%
      ## filtering Proteins that have less than 5 Observation
      dplyr::filter(ValidObs > n) %>%
      ## Removing unnecessary columns
      dplyr::select(-c("nObs", "NAObs", "ValidObs"))
  }
  return(invisible(results))
}

## filter out every Protein with less than n observations in a specified grouping variable
## add roxygen comments
#' @title nObsPerGroup
#' @description This pipeline friendly function filters out every Protein with less than n observations in a specified grouping variable
#' @param dataset The dataset to be filtered
#' @param groupVar The grouping variable (Input without quotation makrs)
#' @param n The number of observations
#' @return The filtered dataset
#' @export
nObsPerGroup <- function(dataset, groupVar, n = 10) {

  ## Error if groupVar is not a column name of dataset
  if (!deparse(substitute(groupVar)) %in% colnames(dataset)) {
    stop("groupVar must be a column name of dataset")
  }

  if("Protein" %in% colnames(dataset)){

    Filter <- dataset %>%
      dplyr::filter(!is.na(Intensity)) %>%
      dplyr::group_by({{groupVar}}, Protein) %>%
      dplyr::summarise(nObs = n()) %>%
      dplyr::ungroup() %>%
      tidyr::pivot_wider(names_from = {{groupVar}}, values_from = nObs) %>%
      ## convert so True if > n
      dplyr::mutate(across(-Protein, ~ . >= n)) %>%
      ## calculate row sum if True == 1
      dplyr::mutate(Total = rowSums(select(., -Protein))) %>%
      ## filter for rows where Total == ncolums - 2
      dplyr::filter(Total == ncol(.) - 2) %>%
      ## get Protein names
      dplyr::pull(Protein)

    FilteredData <- dataset %>%
      dplyr::filter(Protein %in% Filter)

  }


  if("Peptide" %in% colnames(dataset)){

    Filter <- dataset %>%
      dplyr::filter(!is.na(Intensity)) %>%
      dplyr::group_by({{groupVar}}, Peptide) %>%
      dplyr::summarise(nObs = n()) %>%
      dplyr::ungroup() %>%
      tidyr::pivot_wider(names_from = {{groupVar}}, values_from = nObs) %>%
      ## convert so True if > n
      dplyr::mutate(across(-Peptide, ~ . >=n)) %>%
      ## calculate row sum if True == 1
      dplyr::mutate(Total = rowSums(select(., -Peptide))) %>%
      ## filter for rows where Total == ncolums - 2
      dplyr::filter(Total == ncol(.) - 2) %>%
      ## get Peptide names
      dplyr::pull(Peptide)

    FilteredData <- dataset %>%
      dplyr::filter(Peptide %in% Filter)

  }

  return(FilteredData)
}

## assign_colors
## a function to assign colors from a Brewer’s palette to entries in a vector
## add roxygen comments
#' @title assign_colors
#' @description This function assigns colors from different pellets to a list of labels
#' @param unique_entries A vector of unique entries (like unique(dataset$Status))
#' @param palette Color pellets c("custom_vibrant", "viridis", "plasma", "magma", "inferno")
#' @return A color mapping
#' @export
assign_colors <- function(labels, palette = "custom_vibrant") {
  labels <- unique(labels)
  n <- length(labels)

  # Generate colors
  colors <- switch(
    palette,
    "custom_vibrant" = rep(c("#FFAA00",  # golden orange
                        "#9932CC",  # purple (avoid bluish ones)
                        "#00C19F",  # teal / sea green
                        "#ADFF2F",  # lime green (replaces yellow)
                        "#8B8000",  # olive
                        "#FF69B4",  # hot pink
                        "#A0522D",  # sienna brown
                        "#40E0D0"   # turquoise (on the greenish side)
    ), length.out = n),
    "viridis"   = viridisLite::viridis(n),
    "plasma"    = viridisLite::plasma(n),
    "magma"     = viridisLite::magma(n),
    "inferno"   = viridisLite::inferno(n),
    stop("Unknown palette")
  )

  # Return named vector
  setNames(colors, labels)
}



## Imputation of missing value (default method = "mean")
## add roxygen comments
#' @title ImputeFeatureIntensity
#' @description This pipeline friendly function imputes missing values in the dataset
#' @param dataset The dataset to be imputed
#' @param method The method to be used for imputation (mean, zero, half_min, median)
#' @return The imputed dataset
#' @export
ImputeFeatureIntensity <- function(dataset, method = "half_min"){


  MissingValue <- base::sum(base::is.na(dataset$Intensity))

  if("Protein" %in% base::colnames(dataset)){

    dataset <- dataset %>%
      tidyr::pivot_wider(names_from = "Protein", values_from = "Intensity") %>%
      tidyr::pivot_longer(cols = contains("_"), names_to = "Protein", values_to = "Intensity")


    if (method == "mean" | method == "Mean"){
      ## Creating temporary columns
      dataset <- dataset %>%
        dplyr::group_by(Protein) %>%
        dplyr::mutate(meanInt = base::mean(Intensity, na.rm = TRUE)) %>%
        dplyr::ungroup() %>%
        ## Check for missing values and impute them
        dplyr::mutate(Intensity = ifelse(base::is.na(Intensity), meanInt, Intensity))
      ## Selecting necessary columns
      dataset <- dataset %>%
        dplyr::select(-c("meanInt"))

      base::print(base::paste("Imputed", MissingValue, "values"))
      return(dataset)
    }
    if (method == "Zero" | method == "zero" | method == "0"){
      dataset <- dataset %>% dplyr::mutate(Intensity = ifelse(base::is.na(Intensity), 0, Intensity))

      base::print(base::paste("Imputed", MissingValue, "values"))
      return(dataset)
    }
    if (method == "half_min" | method == "halfmin" | method == "HalfMin") {
      dataset <- dataset %>%
        dplyr::group_by(Protein) %>%
        dplyr::mutate(half_min = (base::min(Intensity, na.rm = TRUE) / 2)) %>%
        dplyr::mutate(Intensity = ifelse(base::is.na(Intensity), half_min, Intensity)) %>%
        dplyr::select(-c("half_min"))

    }
    ## Median imputation
    if (method == "median" | method == "Median") {
      dataset <- dataset %>%
        dplyr::group_by(Protein) %>%
        dplyr::mutate(medianInt = base::median(Intensity, na.rm = TRUE)) %>%
        dplyr::mutate(Intensity = ifelse(base::is.na(Intensity), medianInt, Intensity)) %>%
        dplyr::select(-c("medianInt"))

    }
  }

  if("Peptide" %in% base::colnames(dataset)){

    dataset <- dataset %>%
      tidyr::pivot_wider(names_from = "Peptide", values_from = "Intensity") %>%
      tidyr::pivot_longer(cols = contains("_"), names_to = "Peptide", values_to = "Intensity")


    if (method == "mean" | method == "Mean"){
      ## Creating temporary columns
      dataset <- dataset %>%
        dplyr::group_by(Peptide) %>%
        dplyr::mutate(meanInt = base::mean(Intensity, na.rm = TRUE)) %>%
        dplyr::ungroup() %>%
        ## Check for missing values and impute them
        dplyr::mutate(Intensity = ifelse(base::is.na(Intensity), meanInt, Intensity))
      ## Selecting necessary columns
      dataset <- dataset %>%
        dplyr::select(-c("meanInt"))

      base::print(base::paste("Imputed", MissingValue, "values"))
      return(dataset)
    }
    if (method == "Zero" | method == "zero" | method == "0"){
      dataset <- dataset %>% dplyr::mutate(Intensity = ifelse(base::is.na(Intensity), 0, Intensity))

      base::print(base::paste("Imputed", MissingValue, "values"))
      return(dataset)
    }
    if (method == "half_min" | method == "halfmin" | method == "HalfMin") {
      dataset <- dataset %>%
        dplyr::group_by(Peptide) %>%
        dplyr::mutate(half_min = (base::min(Intensity, na.rm = TRUE) / 2)) %>%
        dplyr::mutate(Intensity = ifelse(base::is.na(Intensity), half_min, Intensity)) %>%
        dplyr::select(-c("half_min"))

    }

    if (method == "median" | method == "Median") {
      dataset <- dataset %>%
        dplyr::group_by(Peptide) %>%
        dplyr::mutate(medianInt = base::median(Intensity, na.rm = TRUE)) %>%
        dplyr::mutate(Intensity = ifelse(base::is.na(Intensity), medianInt, Intensity)) %>%
        dplyr::select(-c("medianInt"))

    }
  }

  return(base::invisible(dataset))
}

## method c("mean", "zero", "half_min")

## Intensity normalization on Intensities grouped by Protein (default method ="median')
## add roxygen comments
#' @title normalizeIntensityOnFeature
#' @description This function normalizes MS Intensities on a feature level
#' @param dataset The dataset to be normalized
#' @param method The method to be used for normalization (mean, median, z-score, minmax)
#' @param plot A boolean value indicating whether a plot should be generated
#' @return The normalized dataset
#' @export
normalizeIntensityOnFeature <- function(dataset , method = "minmax", plot = FALSE){

  print(paste("Normalizing MS Intensities using", method))

  if("Protein" %in% colnames(dataset)){

    ## Mean normalization
    if(method == "mean" | method == "Mean"){
      dataset_norm <- dataset %>%
        group_by(Protein) %>%
        mutate(meanInt = mean(Intensity, na.rm = TRUE)) %>%
        group_by(Sample) %>%
        mutate(Intensity = Intensity - meanInt) %>%
        select(-c("medianInt"))
    }
    ## Median normalization
    if(method == "median" | method == "Median"){
      dataset_norm <- dataset %>%
        group_by(Protein) %>%
        mutate(medianInt = median(Intensity, na.rm = TRUE)) %>%
        group_by(Sample) %>%
        mutate(Intensity = Intensity - medianInt) %>%
        select(-c("medianInt"))
    }
    ## z-score normalization
    if(method == "z-score" | method == "zscore" |method == "Z-score"){
      dataset_norm <- dataset %>%
        group_by(Protein) %>%
        mutate(meanInt = mean(Intensity, na.rm = TRUE)) %>%
        mutate(stdev = sd(Intensity, na.rm = TRUE)) %>%
        group_by(Sample) %>%
        mutate(Intensity = ((Intensity-meanInt)/stdev)) %>%
        select(-c("meanInt","stdev"))

    }
    ## MinMax normalization
    if (method == "MinMax" | method == "minmax") {
      dataset_norm <- dataset %>%
        group_by(Protein) %>%
        ## calculating mean intensiy
        mutate(meanInt = mean(Intensity, na.rm = TRUE)) %>%
        mutate(Intensity = Intensity - meanInt) %>%
        ## Rescaling so values are between -1 and 1
        mutate(intMax = max(Intensity, na.rm = TRUE)) %>%
        mutate(intMin = min(Intensity, na.rm = TRUE)) %>%
        mutate(Intensity = (Intensity - intMin) / (intMax - intMin)) %>%
        ## Selecting columns
        select(-c("meanInt","intMax", "intMin"))
    }

    ## Plotting Data before and after normalization
    if (plot == TRUE){
      beforedata <- dataset
      afterdata <- dataset_norm

      ## Plot data before normalization
      beforeplot <- beforedata %>% ggplot(mapping = aes(x = reorder(Protein, -Intensity), y = Intensity)) +
        stat_boxplot() +
        ylab("log 2 Intensity") +
        xlab("Protein") +
        ggtitle("Intensities before normalization on protein intensity") +
        theme(axis.text.x = element_blank())
      ## Plot data after normalization
      afterplot <- afterdata %>% ggplot(mapping = aes(x = reorder(Protein, -Intensity), y = Intensity)) +
        stat_boxplot() +
        ylab("Normalized intensity") +
        xlab("Protein") +
        ggtitle(paste("Intensities after normalization on protein intensity using", method, "normalization")) +
        theme(axis.text.x = element_blank())

      plot(beforeplot)
      plot(afterplot)

    }

  }
  if("Peptide" %in% colnames(dataset)){

    ## Mean normalization
    if(method == "mean" | method == "Mean"){
      dataset_norm <- dataset %>%
        group_by(Peptide) %>%
        mutate(meanInt = mean(Intensity, na.rm = TRUE)) %>%
        group_by(Sample) %>%
        mutate(Intensity = Intensity - meanInt) %>%
        select(-c("medianInt"))
    }
    ## Median normalization
    if(method == "median" | method == "Median"){
      dataset_norm <- dataset %>%
        group_by(Peptide) %>%
        mutate(medianInt = median(Intensity, na.rm = TRUE)) %>%
        group_by(Sample) %>%
        mutate(Intensity = Intensity - medianInt) %>%
        select(-c("medianInt"))
    }
    ## z-score normalization
    if(method == "z-score" | method == "zscore" |method == "Z-score"){
      dataset_norm <- dataset %>%
        group_by(Peptide) %>%
        mutate(meanInt = mean(Intensity, na.rm = TRUE)) %>%
        mutate(stdev = sd(Intensity, na.rm = TRUE)) %>%
        group_by(Sample) %>%
        mutate(Intensity = ((Intensity-meanInt)/stdev)) %>%
        select(-c("meanInt","stdev"))

    }
    ## MinMax normalization
    if (method == "MinMax" | method == "minmax") {
      dataset_norm <- dataset %>%
        group_by(Peptide) %>%
        ## calculating mean intensiy
        mutate(meanInt = mean(Intensity, na.rm = TRUE)) %>%
        mutate(Intensity = Intensity - meanInt) %>%
        ## Rescaling so values are between -1 and 1
        mutate(intMax = max(Intensity, na.rm = TRUE)) %>%
        mutate(intMin = min(Intensity, na.rm = TRUE)) %>%
        mutate(Intensity = (Intensity - intMin) / (intMax - intMin)) %>%
        ## Selecting columns
        select(-c("meanInt","intMax", "intMin"))
    }

    ## Plotting Data before and after normalization
    if (plot == TRUE){
      beforedata <- dataset
      afterdata <- dataset_norm

      ## Plot data before normalization
      beforeplot <- beforedata %>% ggplot(mapping = aes(x = reorder(Peptide, -Intensity), y = Intensity)) +
        stat_boxplot() +
        ylab("log 2 Intensity") +
        xlab("Peptide") +
        ggtitle("Intensities before normalization on protein intensity") +
        theme(axis.text.x = element_blank())
      ## Plot data after normalization
      afterplot <- afterdata %>% ggplot(mapping = aes(x = reorder(Peptide, -Intensity), y = Intensity)) +
        stat_boxplot() +
        ylab("Normalized intensity") +
        xlab("Peptide") +
        ggtitle(paste("Intensities after normalization on protein intensity using", method, "normalization")) +
        theme(axis.text.x = element_blank())

      plot(beforeplot)
      plot(afterplot)

    }

  }





  ## Returning normalized dataset
  return(invisible(dataset_norm))
}
## method <- c("median", "mean", "z-score", MinMax)

## Intensity normalization on Intensities grouped by Samples
## add roxygen comments
#' @title normalizeIntensityOnSample
#' @description This pipeline friendly function normalizes MS Intensities on a sample level
#' @param dataset The dataset to be normalized
#' @param method The method to be used for normalization (mean, median, z-score, minmax)
#' @param plot A boolean value indicating whether a plot should be generated
#' @return The normalized dataset
#' @export
normalizeIntensityOnSample <- function(dataset, method = "minmax", plot = FALSE){

  ## z-score normalization
  if (method == "z-score" | method == "zscore"){
    dataNorm <- dataset %>%
      group_by(Sample) %>%
      mutate(meanInt = mean(Intensity, na.rm = TRUE)) %>%
      mutate(IntNorm = (Intensity - meanInt)/sd(Intensity, na.rm = TRUE)) %>%
      mutate(Intensity = IntNorm) %>%
      select(-c("meanInt","IntNorm"))

    ## Median normalization
  }
  if(method == "median" | method == "Median"){
    dataNorm <- dataset %>%
      group_by(Sample) %>%
      mutate(medianInt = median(Intensity, na.rm = TRUE)) %>%
      mutate(IntNorm = (Intensity - medianInt)) %>%
      mutate(Intensity = IntNorm) %>%
      select(-c("medianInt", "IntNorm"))
  }

  ## Mean normalization
  if(method == "mean" | method == "Mean"){
    dataNorm <- dataset %>%
      group_by(Sample) %>%
      mutate(meanInt = mean(Intensity, na.rm = TRUE)) %>%
      mutate(IntNorm = (Intensity - meanInt)) %>%
      mutate(Intensity = IntNorm) %>%
      select(-c("meanInt", "IntNorm"))
  }
  ## MinMax normalization
  if (method == "MinMax" | method == "minmax") {
    dataNorm <- dataset %>%
      ## calculating difference from the mean intensity
      mutate(meanInt = mean(Intensity, na.rm = TRUE)) %>%
      mutate(Intensity = Intensity - meanInt) %>%
      ## Rescaling so values are between -1 and 1
      mutate(intMax = max(Intensity, na.rm = TRUE)) %>%
      mutate(intMin = min(Intensity, na.rm = TRUE)) %>%
      mutate(Intensity = (Intensity - intMin) / (intMax - intMin)) %>%
      ## Selecting columns
      select(-c("meanInt", "intMax", "intMin"))
  }

  if (plot == TRUE){

    ## plotting Results
    ## before normalization
    plot(ggplot(dataset, aes(x = Sample, y = Intensity))+
           stat_boxplot() +
           ggtitle("Intensities before normalization on Samples") +
           theme(axis.text.x = element_blank()))

    ## after normalization
    plot(ggplot(dataNorm, aes(x = Sample, y = Intensity))+
           stat_boxplot() +
           ggtitle(paste("Intensities after normalization on Samples using", method, "normalization")) +
           theme(axis.text.x = element_blank()))

  }

  ## returning normalized dataset
  return(invisible(dataNorm))

}
## method <- c("median", "mean", "z-score", MinMax)

## CompletentessAnlysis
## add roxygen comments
#' @title CompletenessAnalysis
#' @description This function calculates the number of unique proteins at different completeness levels
#' @param dataset The dataset to be analyzed
#' @return A list containing the results and a plot
#' @export
CompletenessAnalysis <- function(dataset){

  ## define Results dataframe
  Results <- base::data.frame(Cutoff = base::numeric(),
                              UniqueProteins = base::numeric())

  ## define sequence from 0 to 100 in increments of 10
  seq <- base::seq(0, 100, by = 10)

  ## loop through sequence
  for (i in seq){

    CutoffData <- dataset %>% NaCutoff(i)

    ## calculate unique numbers of Proteins
    uniqueProteins <- base::length(base::unique(CutoffData$Protein))

    ## add to Results dataframe
    Results <- base::rbind(Results, base::data.frame(Cutoff = i, UniqueProteins = uniqueProteins))
  }

  ## plot Results using ggplot2
  Plot <- ggplot2::ggplot(Results, ggplot2::aes(x = Cutoff, y = UniqueProteins)) +
    ggplot2::geom_point(size = 5) +
    ggplot2::geom_line(linewidth = 1, color = "blue") +
    ## remove ribbon
    ggplot2::geom_ribbon(aes(ymin = UniqueProteins, ymax = UniqueProteins), alpha = 0) +
    ggplot2::labs(title = "Completeness Analysis",
                  x = "Percent completeness",
                  y = "Unique Proteins") +
    ggplot2::theme_minimal() +
    ## let y axis start at 0
    ggplot2::scale_y_continuous(limits = c(0, base::max(Results$UniqueProteins) + 10))

  ## define Output List
  Output <- base::list(Results = Results,
                       Plot = Plot)

  return(Output)
}


## Non-pipeline functions
## Outlier Removal (from Patrick`s markdown code)
## add roxygen comments
#' @title RemoveOutliers
#' @description This function removes outliers from the dataset
#' @param dataset The dataset to be filtered
#' @param Stdev The standard deviation for the outlier removal
#' @param plotname The name to be displayed on created plots
#' @param plot A boolean value indicating whether a plot should be generated
#' @return A list containing the outlier plot, a table with the outliers and the dataset without the outliers
#' @export
RemoveOutliers <- function(dataset, Stdev = 2, plotname = "Outlier plot", plot = TRUE){
  if("Protein" %in% colnames(dataset)){

    #Determine correlations between sample and theoretical
    SampleCorrelations <- dataset%>%
      group_by(Protein) %>%
      dplyr::mutate(medianProtein = median(Intensity, na.rm = T)) %>%
      ungroup() %>%
      group_by(Sample) %>%
      dplyr::mutate(SampleCor = stats::cor(Intensity, medianProtein, method = "pearson", use="pairwise.complete.obs")) %>%
      ungroup() %>%
      distinct(Sample, .keep_all = T)

    ## Calculate correlation between the Samples and the mean Sample
    SampleCorrelations <- SampleCorrelations %>%
      dplyr::mutate(Outlier = SampleCor < (1 - Stdev*sd(unique(SampleCorrelations$SampleCor), na.rm = T)))

    #Make the DF that we will return (i.e. without outliers)
    returnDataset <- dataset %>%
      dplyr::filter(Sample %in% dplyr::filter(SampleCorrelations, Outlier == F)$Sample)

    #Purely for visualization: Also attached PCA Plots of the removed samples

    #Half minimum value (per protein) imputation
    PCA_DF <- dataset %>%
      group_by(Protein) %>%
      dplyr::mutate(IntensImputed = replace_na(Intensity, mean(Intensity, na.rm = T)/2)) %>%
      dplyr::select(Protein, Sample, IntensImputed) %>%
      spread(key = "Protein", value = "IntensImputed") %>%
      column_to_rownames("Sample")

    #scale and PCA
    PCA_DF_Results <- prcomp(scale(as.matrix(PCA_DF)))
    eigs <- PCA_DF_Results$sdev^2
    variance_percentage <- (eigs / sum(eigs))*100
    pc1var <- round(variance_percentage[1],digits=0)
    pc2var <- round(variance_percentage[2],digits=0)

    #Plot
    OutlierPCAPlot <- left_join(SampleCorrelations, rownames_to_column(data.frame(PCA_DF_Results$x)), by = c("Sample" = "rowname")) %>%
      ggplot(aes(x = PC1, y = PC2, colour = Outlier, shape = Outlier)) +
      geom_point(size = 2.5, alpha = 0.6) +
      theme_bw() +
      xlab(paste('PC1',' (',pc1var,'% variance)',sep='')) +
      ylab(paste('PC2',' (',pc2var,'% variance)',sep='')) +
      scale_color_manual(values = c("black", "red")) +
      ggtitle(plotname)

  }

  if("Peptide" %in% colnames(dataset)){

    #Determine correlations between sample and theoretical
    SampleCorrelations <- dataset%>%
      group_by(Peptide) %>%
      dplyr::mutate(medianPeptide = median(Intensity, na.rm = T)) %>%
      ungroup() %>%
      group_by(Sample) %>%
      dplyr::mutate(SampleCor = stats::cor(Intensity, medianPeptide, method = "pearson", use="pairwise.complete.obs")) %>%
      ungroup() %>%
      distinct(Sample, .keep_all = T)

    ## Calculate correlation between the Samples and the mean Sample
    SampleCorrelations <- SampleCorrelations %>%
      dplyr::mutate(Outlier = SampleCor < (1 - Stdev*sd(unique(SampleCorrelations$SampleCor), na.rm = T)))

    #Make the DF that we will return (i.e. without outliers)
    returnDataset <- dataset %>%
      dplyr::filter(Sample %in% dplyr::filter(SampleCorrelations, Outlier == F)$Sample) %>%
      dplyr::select(Sample, Peptide, Intensity)

    #Purely for visualization: Also attached PCA Plots of the removed samples

    #Half minimum value (per Peptide) imputation
    PCA_DF <- dataset %>%
      group_by(Peptide) %>%
      dplyr::mutate(IntensImputed = replace_na(Intensity, mean(Intensity, na.rm = T)/2)) %>%
      dplyr::select(Peptide, Sample, IntensImputed) %>%
      spread(key = "Peptide", value = "IntensImputed") %>%
      column_to_rownames("Sample")

    #scale and PCA
    PCA_DF_Results <- prcomp(scale(as.matrix(PCA_DF)))
    eigs <- PCA_DF_Results$sdev^2
    variance_percentage <- (eigs / sum(eigs))*100
    pc1var <- round(variance_percentage[1],digits=0)
    pc2var <- round(variance_percentage[2],digits=0)

    #Plot
    OutlierPCAPlot <- left_join(SampleCorrelations, rownames_to_column(data.frame(PCA_DF_Results$x)), by = c("Sample" = "rowname")) %>%
      ggplot(aes(x = PC1, y = PC2, colour = Outlier, shape = Outlier)) +
      geom_point(size = 2.5, alpha = 0.6) +
      theme_bw() +
      xlab(paste('PC1',' (',pc1var,'% variance)',sep='')) +
      ylab(paste('PC2',' (',pc2var,'% variance)',sep='')) +
      scale_color_manual(values = c("black", "red")) +
      ggtitle(plotname)


  }

  ## creating output object
  Output <- list()
  Output$OutlierPlot <- OutlierPCAPlot
  Output$Table <- SampleCorrelations
  Output$Dataset <- returnDataset

  #Print and Return
  print(paste("The following were removed from:", toString(dplyr::filter(SampleCorrelations, Outlier == T)$Sample)))
  return(Output)

}

## Batch correction using ComBat
## add roxygen comments
#' @title ComBat
#' @description This function corrects batch effects using ComBat
#' @param dataset The dataset to be corrected
#' @return The corrected dataset
#' @export
ComBat <- function(dataset){

  ## Preparing Quant Data
  ComBatDataQaunt <- dataset %>%
    ## Impute missing values
    NaCutoff(50) %>%
    ImputeFeatureIntensity() %>%
    ## Pivot wider
    pivot_wider(names_from = "Protein", values_from = "Intensity") %>%
    ## Select Proteins
    select(contains("_")) %>%
    t() %>% as.matrix()

  ComBatDataClin <- dataset %>%
    pivot_wider(names_from = Protein, values_from = Intensity) %>%
    select(!contains("_"))

  colnames(ComBatDataQaunt) <- ComBatDataClin$Sample

  PlateVector <- ComBatDataClin$Plate

  CorrectedData <- sva::ComBat(dat = ComBatDataQaunt, batch = PlateVector, prior.plots = TRUE) %>%
    t() %>%
    data.frame() %>%
    rownames_to_column(var = "Sample") %>%
    pivot_longer(cols = -"Sample", names_to = "Protein", values_to = "Intensity")

  CorrectedData <- merge(ComBatDataClin, CorrectedData , by = "Sample")

  Output <- list()
  Output$CorrectedData <- CorrectedData

  return(Output)
}

## Significance Tests
## Finding NA Cutoff
## add roxygen comments
#' @title FindNACutoff
#' @description This function finds the optimal NA cutoff value to maximize the number of significant features by t-test
#' @param dataset The dataset to be filtered
#' @param plotname The name to be displayed on created plots
#' @param p.adj.method The method to be used for p-value adjustment ("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
#' @param Test The test to be used for the analysis ("Wilcox", "Student")
#' @export
FindNACutoff <- function(dataset, plotname = "", p.adj.method = "BH", Test = "Wilcox"){

  Results <- data.frame()

  ## define range of steps (1 to 100 in steps of 10)
  Steps <- seq(0, 100, by = 10)

  for(i in Steps){

    TestData <- dataset %>% NaCutoff(i)


    if(Test == "Wilcox"){

      TestResults <- BiomarkR::WTest(dataset = TestData, plotname = plotname, p.adj.method = p.adj.method)

    }
    if(Test == "Student"){

      TestResults <- BiomarkR::TTest(dataset = TestData, plotname = plotname, p.adj.method = p.adj.method)
    }

    nprot <- length(unique(TestData$Protein))
    SigProt <- nrow(TestResults$Significant)

    ## put results into data frame
    Results <- rbind(Results, data.frame(NACutoff = i, nprot = nprot, SigProt = SigProt))

  }

  ## plot Results
  Plot <- Results %>%
    ggplot(aes(x = NACutoff, y = SigProt)) +
    geom_point() +
    geom_line() +
    theme_minimal() +
    ggtitle("NA Cutoff Optimization") +
    xlab("NA Cutoff") + ylab("Number of significant Proteins") +
    ## fix min(y) at zero
    scale_y_continuous(limits = c(0, max(Results$SigProt)))

  return(list(Results = Results, Plot = Plot))
}


## T-Test, volcano plot and Heat Map
## add roxygen comments
#' @title TTest
#' @description This function performs differential expression analysis using a T-Test.
#' @param dataset The dataset to be tested
#' @param plotname The name to be displayed on created plots
#' @param method The method to be used for the Heatmap (unsupervised, supervised)
#' @param clustDist The distance metric to be used for clustering in the Heatmap ("euclidean", "maximum", "man-hattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall")
#' @param p.adj.method The method to be used for p-value adjustment ("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
#' @export
TTest <- function(dataset, plotname = "", method = "unsupervised", clustDist = "euclidean", p.adj.method = "BH") {
  datasetT <- dataset %>% dplyr::arrange(Status)

  Status1 <- unique(datasetT$Status)[1] %>% as.character()
  Status2 <- unique(datasetT$Status)[2] %>% as.character()

  ## error if there are more than 2 groups
  if(length(unique(datasetT$Status)) > 2) {
    stop("Only two groups (Status) are allowed for Wilcox-Test")
  }

  ## Running the W-Test
  if("Protein" %in% colnames(dataset)) {
    ## Making sure we have at least 2 Observations per group
    filter <- datasetT %>%
      dplyr::select(c("Sample", "Status", "Protein", "Intensity")) %>%
      tidyr::pivot_wider(names_from = "Status", values_from = "Intensity") %>%
      dplyr::group_by(Protein) %>%
      dplyr::summarise(
        ObsInStatus1 = sum(!is.na(.data[[Status1]])),
        ObsInStatus2 = sum(!is.na(.data[[Status2]]))
      ) %>%
      dplyr::mutate(possible = ifelse(ObsInStatus1 < 2 | ObsInStatus2 < 2, FALSE, TRUE)) %>%
      dplyr::filter(possible)

    TResults <- datasetT %>%
      dplyr::filter(Protein %in% filter$Protein) %>%
      dplyr::group_by(Protein) %>%
      rstatix::t_test(Intensity ~ Status, detailed = TRUE) %>%
      rstatix::adjust_pvalue(method = p.adj.method) %>%
      dplyr::mutate(
        UniprotID = stringr::str_split_i(Protein, pattern = "_", 1),
        Gene = stringr::str_split_i(Protein, pattern = "_", 2),
        log10adjustP = -1 * log10(p.adj)
      )

    ## Create dataframe of significant Proteins in global environment
    TSignificantFeatures <- TResults %>%
      dplyr::filter(p.adj < 0.05) %>%
      dplyr::arrange(p.adj) %>%
      as.data.frame() %>%
      dplyr::mutate(Protein = paste0(UniprotID, "_", Gene))

    print(paste(nrow(TSignificantFeatures), "Significant proteins have been identified"))

    FoldChangeData <- dataset %>%
      dplyr::group_by(Status, Protein) %>%
      dplyr::summarise(meanInt = mean(Intensity, na.rm = TRUE)) %>%
      tidyr::pivot_wider(names_from = "Status", values_from = "meanInt") %>%
      dplyr::mutate(FC = (.[[Status1]] - .[[Status2]]))

    VulconaoPlotData <- merge(TResults, FoldChangeData, by = "Protein") %>%
      dplyr::mutate(Direction = ifelse(p.adj > 0.05, "NotSignificant", ifelse(FC < 0, "Down", "Up")))

    ## Volcano plot of results
    vulcanoPlot <- ggplot2::ggplot(data = VulconaoPlotData) +
      ggplot2::geom_point(
        size = 3.5, shape = 21,
        data = subset(VulconaoPlotData, Direction == "NotSignificant"),
        ggplot2::aes(x = FC, y = log10adjustP, fill = "NotSignificant")
      ) +
      ggplot2::geom_point(
        size = 3.5, shape = 21,
        data = subset(VulconaoPlotData, Direction == "Up"),
        ggplot2::aes(x = FC, y = log10adjustP, fill = "Up")
      ) +
      ggplot2::geom_point(
        size = 3.5, shape = 21,
        data = subset(VulconaoPlotData, Direction == "Down"),
        ggplot2::aes(x = FC, y = log10adjustP, fill = "Down")
      ) +
      ggplot2::scale_fill_manual(values = c("Up" = "red", "Down" = "blue", "NotSignificant" = "grey")) +
      ggplot2::geom_hline(yintercept = -log10(0.05), alpha = 0.7, linetype = 2) +
      ggplot2::geom_hline(yintercept = -log10(0.01), alpha = 0.7, linetype = 2, col = "red") +
      ## Add Protein names using ggrepel
      ggrepel::geom_text_repel(
        data = subset(VulconaoPlotData, log10adjustP > 1.3),
        ggplot2::aes(label = Gene, x = FC, y = log10adjustP),
        box.padding = 0.3,
        point.padding = 0.3,
        segment.color = "grey50",
        segment.size = 0.1,
        segment.alpha = 0.1,
        size = 3,
        angle = 30
      ) +
      ggplot2::ggtitle(plotname) +
      ggplot2::xlab(paste("Fold change in", unique(TResults$group1))) +
      ggplot2::theme_light(base_size = 13) +
      ## rename y axis
      ggplot2::ylab(ifelse(p.adj.method != "none", "log10 (p.adj)", "log10 (p)"))
  }
  if("Peptide" %in% colnames(dataset)) {
    ## Making sure we have at least 2 Observations per group
    filter <- datasetT %>%
      dplyr::select(c("Sample", "Status", "Peptide", "Intensity")) %>%
      tidyr::pivot_wider(names_from = "Status", values_from = "Intensity") %>%
      dplyr::group_by(Peptide) %>%
      dplyr::summarise(
        ObsInStatus1 = sum(!is.na(.data[[Status1]])),
        ObsInStatus2 = sum(!is.na(.data[[Status2]]))
      ) %>%
      dplyr::mutate(possible = ifelse(ObsInStatus1 < 2 | ObsInStatus2 < 2, FALSE, TRUE)) %>%
      dplyr::filter(possible)

    TResults <- datasetT %>%
      dplyr::filter(Peptide %in% filter$Peptide) %>%
      dplyr::group_by(Peptide) %>%
      rstatix::t_test(Intensity ~ Status, detailed = TRUE) %>%
      ## Adjusting p-values for multiple testing
      dplyr::mutate(Gene = stringr::str_split_i(Peptide, pattern = "_", 2)) %>%
      dplyr::mutate(p.adj = p * length(unique(Gene))) %>%
      dplyr::mutate(
        UniprotID = stringr::str_split_i(Peptide, pattern = "_", 1),
        Gene = stringr::str_split_i(Peptide, pattern = "_", 3),
        log10adjustP = -1 * log10(p.adj)
      )

    ## Create dataframe of significant Peptides in global environment
    TSignificantFeatures <- TResults %>%
      dplyr::filter(p.adj < 0.05) %>%
      dplyr::arrange(p.adj) %>%
      as.data.frame()

    print(paste(nrow(TSignificantFeatures), "Significant Peptides have been identified"))

    FoldChangeData <- dataset %>%
      dplyr::group_by(Status, Peptide) %>%
      dplyr::summarise(meanInt = mean(Intensity, na.rm = TRUE)) %>%
      tidyr::pivot_wider(names_from = "Status", values_from = "meanInt") %>%
      dplyr::mutate(FC = (.[[Status1]] - .[[Status2]]))

    VulconaoPlotData <- merge(TResults, FoldChangeData, by = "Peptide") %>%
      dplyr::mutate(Direction = ifelse(p.adj > 0.05, "NotSignificant", ifelse(FC < 0, "Down", "Up")))

    ## Volcano plot of results
    vulcanoPlot <- ggplot2::ggplot(data = VulconaoPlotData) +
      ggplot2::geom_point(
        size = 3.5, shape = 21,
        data = subset(VulconaoPlotData, Direction == "NotSignificant"),
        ggplot2::aes(x = FC, y = log10adjustP, fill = "NotSignificant")
      ) +
      ggplot2::geom_point(
        size = 3.5, shape = 21,
        data = subset(VulconaoPlotData, Direction == "Up"),
        ggplot2::aes(x = FC, y = log10adjustP, fill = "Up")
      ) +
      ggplot2::geom_point(
        size = 3.5, shape = 21,
        data = subset(VulconaoPlotData, Direction == "Down"),
        ggplot2::aes(x = FC, y = log10adjustP, fill = "Down")
      ) +
      ggplot2::scale_fill_manual(values = c("Up" = "red", "Down" = "blue", "NotSignificant" = "grey")) +
      ggplot2::geom_hline(yintercept = -log10(0.05), alpha = 0.7, linetype = 2) +
      ggplot2::geom_hline(yintercept = -log10(0.01), alpha = 0.7, linetype = 2, col = "red") +
      ggplot2::geom_text(
        data = subset(VulconaoPlotData, log10adjustP > 1.3),
        ggplot2::aes(label = Gene, x = FC, y = log10adjustP),
        vjust = 0.5, hjust = -0.2, size = 3, angle = 30
      ) +
      ggplot2::ggtitle(plotname) +
      ggplot2::xlab(paste("Fold change in", unique(TResults$group1))) +
      ggplot2::theme_light(base_size = 13) +
      ggplot2::ylab(ifelse(p.adj.method != "none", "log10 (p.adj)", "log10 (p)"))
  }

  ## plotting Heatmap using sigificant Proteins or peptides
  ## NOTE clustering distance c("euclidean", "maximum", "man-hattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall")

  if("Protein" %in% colnames(dataset)){

    Heatmap <- BiomarkR::HeatMap(datasetT, PoIs = TSignificantFeatures$Protein, method = method, clustDist = clustDist, show_column_names = F, show_row_names = F, plotname = plotname)

  }

  if("Peptide" %in% colnames(dataset)){

    Heatmap <- BiomarkR::HeatMap(datasetT, PoIs = TSignificantFeatures$Peptide, method = method, clustDist = clustDist, show_column_names = F, show_row_names = F,plotname = plotname)


  }


  ## Preparing Output object
  Output <- list()
  Output$raw <- TResults
  Output$Significant <- TSignificantFeatures
  Output$Vulcanoplot <- vulcanoPlot
  Output$Heatmap <- Heatmap

  return(Output)
}

## Wilcox test for significance
## add roxygen comments
#' @title WTest
#' @description This function performs differential expression analysis using a Wilcox test.
#' @param dataset The dataset to be tested
#' @param plotname The name to be displayed on created plots
#' @param method The method to be used for the Heatmap (unsupervised, supervised)
#' @param clustDist The distance metric to be used for clustering in the Heatmap ("euclidean", "maximum", "man-hattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall")
#' @param p.adj.method The method to be used for p-value adjustment ("BH", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
#' @return A list object containing the results of the Wilcox test, the significant features and a volcano plot
#' @export
WTest <- function(dataset, plotname = "", method = "unsupervised", clustDist = "euclidean", p.adj.method = "BH", ColPalette = "custom_vibrant") {
  datasetW <- dataset %>% dplyr::arrange(Status)

  Status1 <- unique(datasetW$Status)[1] %>% as.character()
  Status2 <- unique(datasetW$Status)[2] %>% as.character()

  ## error if there are more than 2 groups
  if(length(unique(datasetW$Status)) > 2) {
    stop("Only two groups (Status) are allowed for Wilcox-Test")
  }

  ## Running the W-Test
  if("Protein" %in% colnames(dataset)) {
    ## Making sure we have at least 2 Observations per group
    filter <- datasetW %>%
      dplyr::select(c("Sample", "Status", "Protein", "Intensity")) %>%
      tidyr::pivot_wider(names_from = "Status", values_from = "Intensity") %>%
      dplyr::group_by(Protein) %>%
      dplyr::summarise(
        ObsInStatus1 = sum(!is.na(.data[[Status1]])),
        ObsInStatus2 = sum(!is.na(.data[[Status2]]))
      ) %>%
      dplyr::mutate(possible = ifelse(ObsInStatus1 < 2 | ObsInStatus2 < 2, FALSE, TRUE)) %>%
      dplyr::filter(possible)

    WResults <- datasetW %>%
      dplyr::filter(Protein %in% filter$Protein) %>%
      dplyr::group_by(Protein) %>%
      rstatix::wilcox_test(Intensity ~ Status, detailed = TRUE) %>%
      rstatix::adjust_pvalue(method = p.adj.method) %>%
      dplyr::mutate(
        UniprotID = stringr::str_split_i(Protein, pattern = "_", 1),
        Gene = stringr::str_split_i(Protein, pattern = "_", 2),
        log10adjustP = -1 * log10(p.adj)
      )

    ## Create dataframe of significant Proteins in global environment
    WilcoxSignificantFeatures <- WResults %>%
      dplyr::filter(p.adj < 0.05) %>%
      dplyr::arrange(p.adj) %>%
      as.data.frame() %>%
      dplyr::mutate(Protein = paste0(UniprotID, "_", Gene))

    print(paste(nrow(WilcoxSignificantFeatures), "Significant proteins have been identified"))

    FoldChangeData <- dataset %>%
      dplyr::group_by(Status, Protein) %>%
      dplyr::summarise(meanInt = mean(Intensity, na.rm = TRUE)) %>%
      tidyr::pivot_wider(names_from = "Status", values_from = "meanInt") %>%
      dplyr::mutate(FC = (.[[Status1]] - .[[Status2]]))

    VulconaoPlotData <- merge(WResults, FoldChangeData, by = "Protein") %>%
      dplyr::mutate(Direction = ifelse(p.adj > 0.05, "NotSignificant", ifelse(FC < 0, "Down", "Up")))

    ## Volcano plot of results
    vulcanoPlot <- ggplot2::ggplot(data = VulconaoPlotData) +
      ggplot2::geom_point(
        size = 3.5, shape = 21,
        data = subset(VulconaoPlotData, Direction == "NotSignificant"),
        ggplot2::aes(x = FC, y = log10adjustP, fill = "NotSignificant")
      ) +
      ggplot2::geom_point(
        size = 3.5, shape = 21,
        data = subset(VulconaoPlotData, Direction == "Up"),
        ggplot2::aes(x = FC, y = log10adjustP, fill = "Up")
      ) +
      ggplot2::geom_point(
        size = 3.5, shape = 21,
        data = subset(VulconaoPlotData, Direction == "Down"),
        ggplot2::aes(x = FC, y = log10adjustP, fill = "Down")
      ) +
      ggplot2::scale_fill_manual(values = c("Up" = "red", "Down" = "blue", "NotSignificant" = "grey")) +
      ggplot2::geom_hline(yintercept = -log10(0.05), alpha = 0.7, linetype = 2) +
      ggplot2::geom_hline(yintercept = -log10(0.01), alpha = 0.7, linetype = 2, col = "red") +
      ## Add Protein names using ggrepel
      ggrepel::geom_text_repel(
        data = subset(VulconaoPlotData, log10adjustP > 1.3),
        ggplot2::aes(label = Gene, x = FC, y = log10adjustP),
        box.padding = 0.3,
        point.padding = 0.3,
        segment.color = "grey50",
        segment.size = 0.1,
        segment.alpha = 0.1,
        size = 3,
        angle = 30
      ) +
      ggplot2::ggtitle(plotname) +
      ggplot2::xlab(paste("Fold change in", unique(WResults$group1))) +
      ggplot2::theme_light(base_size = 13) +
      ## rename y axis
      ggplot2::ylab(ifelse(p.adj.method != "none", "log10 (p.adj)", "log10 (p)"))
  }
  if("Peptide" %in% colnames(dataset)) {
    ## Making sure we have at least 2 Observations per group
    filter <- datasetW %>%
      dplyr::select(c("Sample", "Status", "Peptide", "Intensity")) %>%
      tidyr::pivot_wider(names_from = "Status", values_from = "Intensity") %>%
      dplyr::group_by(Peptide) %>%
      dplyr::summarise(
        ObsInStatus1 = sum(!is.na(.data[[Status1]])),
        ObsInStatus2 = sum(!is.na(.data[[Status2]]))
      ) %>%
      dplyr::mutate(possible = ifelse(ObsInStatus1 < 2 | ObsInStatus2 < 2, FALSE, TRUE)) %>%
      dplyr::filter(possible)

    WResults <- datasetW %>%
      dplyr::filter(Peptide %in% filter$Peptide) %>%
      dplyr::group_by(Peptide) %>%
      rstatix::wilcox_test(Intensity ~ Status, detailed = TRUE) %>%
      ## Adjusting p-values for multiple testing
      dplyr::mutate(Gene = stringr::str_split_i(Peptide, pattern = "_", 2)) %>%
      rstatix::adjust_pvalue(method = p.adj.method) %>%
      dplyr::mutate(
        UniprotID = stringr::str_split_i(Peptide, pattern = "_", 1),
        Gene = stringr::str_split_i(Peptide, pattern = "_", 3),
        log10adjustP = -1 * log10(p.adj)
      )

    ## Create dataframe of significant Peptides in global environment
    WilcoxSignificantFeatures <- WResults %>%
      dplyr::filter(p.adj < 0.05) %>%
      dplyr::arrange(p.adj) %>%
      as.data.frame()

    print(paste(nrow(WilcoxSignificantFeatures), "Significant Peptides have been identified"))

    FoldChangeData <- dataset %>%
      dplyr::group_by(Status, Peptide) %>%
      dplyr::summarise(meanInt = mean(Intensity, na.rm = TRUE)) %>%
      tidyr::pivot_wider(names_from = "Status", values_from = "meanInt") %>%
      dplyr::mutate(FC = (.[[Status1]] - .[[Status2]]))

    VulconaoPlotData <- merge(WResults, FoldChangeData, by = "Peptide") %>%
      dplyr::mutate(Direction = ifelse(p.adj > 0.05, "NotSignificant", ifelse(FC < 0, "Down", "Up")))

    ## Volcano plot of results
    vulcanoPlot <- ggplot2::ggplot(data = VulconaoPlotData) +
      ggplot2::geom_point(
        size = 3.5, shape = 21,
        data = subset(VulconaoPlotData, Direction == "NotSignificant"),
        ggplot2::aes(x = FC, y = log10adjustP, fill = "NotSignificant")
      ) +
      ggplot2::geom_point(
        size = 3.5, shape = 21,
        data = subset(VulconaoPlotData, Direction == "Up"),
        ggplot2::aes(x = FC, y = log10adjustP, fill = "Up")
      ) +
      ggplot2::geom_point(
        size = 3.5, shape = 21,
        data = subset(VulconaoPlotData, Direction == "Down"),
        ggplot2::aes(x = FC, y = log10adjustP, fill = "Down")
      ) +
      ggplot2::scale_fill_manual(values = c("Up" = "red", "Down" = "blue", "NotSignificant" = "grey")) +
      ggplot2::geom_hline(yintercept = -log10(0.05), alpha = 0.7, linetype = 2) +
      ggplot2::geom_hline(yintercept = -log10(0.01), alpha = 0.7, linetype = 2, col = "red") +
      ggplot2::geom_text(
        data = subset(VulconaoPlotData, log10adjustP > 1.3),
        ggplot2::aes(label = Gene, x = FC, y = log10adjustP),
        vjust = 0.5, hjust = -0.2, size = 3, angle = 30
      ) +
      ggplot2::ggtitle(plotname) +
      ggplot2::xlab(paste("Fold change in", unique(WResults$group1))) +
      ggplot2::theme_light(base_size = 13) +
      ggplot2::ylab(ifelse(p.adj.method != "none", "log10 (p.adj)", "log10 (p)"))
  }

  ## plotting Heatmap using sigificant Proteins or peptides
  ## NOTE clustering distance c("euclidean", "maximum", "man-hattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall")

  if("Protein" %in% colnames(dataset)){

    Heatmap <- BiomarkR::HeatMap(datasetW %>% select(c(Sample,Protein,Intensity,Status)), PoIs = WilcoxSignificantFeatures$Protein, method = method, clustDist = clustDist, show_column_names = F, show_row_names = F, plotname = plotname, ColPalette = ColPalette)

  }

  if("Peptide" %in% colnames(dataset)){

    Heatmap <- BiomarkR::HeatMap(datasetW %>% select(c(Sample,Peptide,Intensity,Status)), PoIs = WilcoxSignificantFeatures$Peptide, method = method, clustDist = clustDist, show_column_names = F, show_row_names = F,plotname = plotname, ColPalette = ColPalette)


  }


  ## Preparing Output object
  Output <- list()
  Output$raw <- WResults
  Output$Significant <- WilcoxSignificantFeatures
  Output$Vulcanoplot <- vulcanoPlot
  Output$Heatmap <- Heatmap

  return(Output)
}

## Fisher Test
## add roxygen comments
#' @title FisherTest
#' @description This function performs a Fisher test for differential expression analysis.
#' @param dataset The dataset to be tested
#' @return A list object containing the results of the Fisher test
#' @export
FisherTest <- function(dataset, p.adjust.method = "BH"){

  ## Prepare Data

  if("Protein" %in% colnames(dataset)){

    PoIs <- unique(dataset$Protein)

    DataForFisher <- dataset %>%
      dplyr::filter(Protein %in% PoIs) %>%
      dplyr::select(Sample, Protein, Intensity, Status) %>%
      tidyr::pivot_wider(names_from = Protein, values_from = Intensity) %>%
      tidyr::pivot_longer(cols = -c(Sample, Status), names_to = "Protein", values_to = "Intensity")

    ## Create contingency table
    ContingencyTable <- DataForFisher %>%
      dplyr::group_by(Protein, Status) %>%
      dplyr::summarise(Total = dplyr::n(),
                       Missing = sum(is.na(Intensity)), .groups = 'drop',
                       Count = Total - Missing) %>%
      split(.$Protein)

    ## Run Fisher Test
    Results <- list()

    for (i in 1:length(ContingencyTable)) {
      ConfusionMatrix <- as.data.frame(ContingencyTable[[i]]) %>%
        dplyr::select(-c("Protein", "Total")) %>%
        tibble::column_to_rownames("Status") %>%
        as.matrix()

      FisherResults <- rstatix::fisher_test(ConfusionMatrix, detailed = TRUE,)

      Protein <- ContingencyTable[[i]]$Protein[1]
      p <- FisherResults$p
      n <- FisherResults$n

      Results[[i]] <- data.frame(Protein, p, n)

    }
  }

  if("Peptide" %in% colnames(dataset)){

    PoIs <- unique(dataset$Peptide)

    DataForFisher <- dataset %>%
      dplyr::filter(Peptide %in% PoIs) %>%
      dplyr::select(Sample, Peptide, Intensity, Status) %>%
      tidyr::pivot_wider(names_from = Peptide, values_from = Intensity) %>%
      tidyr::pivot_longer(cols = -c(Sample, Status), names_to = "Peptide", values_to = "Intensity")

    ## Create contingency table
    ContingencyTable <- DataForFisher %>%
      dplyr::group_by(Peptide, Status) %>%
      dplyr::summarise(Total = dplyr::n(),
                       Missing = sum(is.na(Intensity)), .groups = 'drop',
                       Count = Total - Missing) %>%
      split(.$Peptide)

    ## Run Fisher Test
    Results <- list()

    for (i in 1:length(ContingencyTable)) {
      ConfusionMatrix <- as.data.frame(ContingencyTable[[i]]) %>%
        dplyr::select(-c("Peptide", "Total")) %>%
        tibble::column_to_rownames("Status") %>%
        as.matrix() %>% t()

      FisherResults <- rstatix::fisher_test(ConfusionMatrix, detailed = TRUE)

      Peptide <- ContingencyTable[[i]]$Peptide[1]
      p <- FisherResults$p
      n <- FisherResults$n

      Results[[i]] <- data.frame(Peptide, p, n)

  }


  }

  ## Create output

  if("Protein" %in% colnames(dataset)){

    ## calculate frequency of observations per group
    Frequency <- DataForFisher %>%
      dplyr::filter(!is.na(Intensity)) %>%
      dplyr::count(Protein, Status) %>%
      tidyr::pivot_wider(names_from = Status, values_from = n)

    Results <- dplyr::bind_rows(Results) %>%
      ## Adjusting p-values for multiple testing
      dplyr::mutate(Gene = stringr::str_split_i(Protein, pattern = "_", 2)) %>%
      rstatix::adjust_pvalue(method = p.adjust.method) %>%
      dplyr::arrange(p.adj) %>%
      dplyr::mutate(log10adjustP = -log10(p.adj)) %>%
      dplyr::mutate(Gene = stringr::str_split_i(Protein, pattern = "_", 2))

  }

  if("Peptide" %in% colnames(dataset)){

    ## calculate frequency of observations per group
    Frequency <- DataForFisher %>%
      dplyr::filter(!is.na(Intensity)) %>%
      dplyr::count(Peptide, Status) %>%
      tidyr::pivot_wider(names_from = Status, values_from = n)

    Results <- dplyr::bind_rows(Results) %>%
      ## Adjusting p-values for multiple testing
      dplyr::mutate(Gene = stringr::str_split_i(Peptide, pattern = "_", 2)) %>%
      rstatix::adjust_pvalue(method = p.adjust.method) %>%
      dplyr::arrange(p.adj) %>%
      dplyr::mutate(log10adjustP = -log10(p.adj)) %>%
      dplyr::mutate(Gene = stringr::str_split_i(Peptide, pattern = "_", 2))

  }

  ## make dataframe with significant entries
  Significant <- Results %>% dplyr::filter(p.adj < 0.05) %>% dplyr::arrange(p.adj)

  ## return output
  return(list(Results = Results,
              Frequency = Frequency,
              Significant = Significant))



}


## ANOVA
## add roxygen comments
#' @title ANOVA
#' @description This function performs differential expression analysis using an ANOVA.
#' @param dataset The dataset to be tested
#' @param plotname The name to be displayed on created plots
#' @param clustDist The distance metric to be used for clustering in the Heatmap ("euclidean", "maximum", "man-hattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall")
#' @param method The method to be used for the Heatmap (unsupervised, supervised)
#' @return A list object containing the results of the ANOVA, the significant features and a heatmap
#' @export
ANOVA <- function(dataset, plotname= "", clustDist = "euclidean", method = "unsupervised"){

  ## Initializing Dataframe
  datasetANOVA <- dataset

  ## Making the "Status" column into a factor
  datasetANOVA$Status <- datasetANOVA$Status %>% as.factor()

  ## Cheking if we are operating on Peptides or Proteins
  if("Protein" %in% colnames(dataset)){

    ## Running the ANOVA
    ANOVAResults <<- datasetANOVA %>%
      ## Preparing Dataframe
      drop_na(Intensity) %>%
      group_by(Protein) %>%
      anova_test(Intensity ~ Status, detailed = T) %>%
      magrittr::set_class(c("anova_test", "rstatix_test", "data.frame")) %>%
      adjust_pvalue(method = "BH") %>%
      separate(Protein, into = c("Uniprot", "Gene"), sep = "_", remove = F)

    ANOVASignificantFeatures <<- ANOVAResults %>%
      filter(p.adj < 0.05) %>%
      arrange(p.adj)


  }
  if("Peptide" %in% colnames(dataset)){

    ## Running the ANOVA
    ANOVAResults <<- datasetANOVA %>%
      ## Preparing Dataframe
      drop_na(Intensity) %>%
      group_by(Peptide) %>%
      anova_test(Intensity ~ Status, detailed = T) %>%
      magrittr::set_class(c("anova_test", "rstatix_test", "data.frame")) %>%
      adjust_pvalue(method = "BH") %>%
      separate(Protein, into = c("Uniprot", "Gene"), sep = "_", remove = F)

    ANOVASignificantFeatures <- ANOVAResults %>%
      filter(p.adj < 0.05) %>%
      arrange(p.adj)

  }
  ## Running the ANOVA
  ANOVAResults <- datasetANOVA %>%
    ## Preparing Dataframe
    drop_na(Intensity) %>%
    group_by(Protein) %>%
    anova_test(Intensity ~ Status, detailed = T) %>%
    magrittr::set_class(c("anova_test", "rstatix_test", "data.frame")) %>%
    adjust_pvalue(method = "BH") %>%
    separate(Protein, into = c("Uniprot", "Gene"), sep = "_", remove = F)

  ANOVASignificantFeatures <- ANOVAResults %>%
    filter(p.adj < 0.05) %>%
    arrange(- p.adj)

    ## creating heat map Data
    if("Protein" %in% colnames(dataset)){

      HeatMapData <- dataset %>%
        filter(Protein %in% ANOVASignificantFeatures$Protein) %>%
        group_by(Protein) %>%
        mutate(Intensity = scale(Intensity, center = TRUE, scale = TRUE))

      ## Quantitative heat map data
      HeatMapDataQuant <- HeatMapData %>%
        pivot_wider(names_from = "Protein", values_from = "Intensity") %>%
        select(contains("_")) %>%
        t() %>% as.matrix()

      ## clinical heat map data
      HeatMapDataClin <- HeatMapData %>%
        pivot_wider(names_from = "Protein", values_from = "Intensity") %>%
        select(!contains("_"))
    }
    if("Peptide" %in% colnames(dataset)){

      HeatMapData <- dataset %>%
        filter(Peptide %in% ANOVASignificantFeatures$Peptide) %>%
        group_by(Peptide) %>%
        mutate(Intensity = scale(Intensity, center = TRUE, scale = TRUE))

      ## Quantitative heat map data
      HeatMapDataQuant <- HeatMapData %>%
        pivot_wider(names_from = "Peptide", values_from = "Intensity") %>%
        select(contains("_")) %>%
        t() %>% as.matrix()

      ## clinical heat map data
      HeatMapDataClin <- HeatMapData %>%
        pivot_wider(names_from = "Peptide", values_from = "Intensity") %>%
        select(!contains("_"))

    }

    ## Annotations
    ## for now only status is annotated
    ## maybe i can find a general way to annotate all clinical variables
    "colnames(HeatMapDataClin[1]) = HeatMapDataClin[1] works; need to find a way to generalize"

    Annotation <- HeatmapAnnotation(
      Status = HeatMapDataClin$Status
    )

    if(method == "supervised" | method == "Supervised"){

      ANOVAHeatMap <- ComplexHeatmap::Heatmap(HeatMapDataQuant,

                                              ## Annotation stuff
                                              top_annotation = Annotation,

                                              ## clustering specifics
                                              ## Clustering columns
                                              cluster_columns = TRUE,
                                              clustering_distance_columns = clustDist,

                                              ## clustering Rows
                                              cluster_rows = TRUE,
                                              clustering_distance_rows = clustDist,
                                              show_row_names = FALSE,

                                              ## specify pacient status as main cluster
                                              column_split = HeatMapDataClin$Status,

                                              ## Changing Legend title
                                              name = "z-score Int",

                                              ## naming Plot
                                              column_title = paste(plotname,"Supervised Heat map clustered by", clustDist)
      )


    }
    if(method == "unsupervised" | method == "Unsupervised"){

      ANOVAHeatMap <- ComplexHeatmap::Heatmap(HeatMapDataQuant,

                                              ## Annotations Stuff
                                              top_annotation = Annotation,

                                              ## clustering specifics
                                              ## Clustering columns
                                              cluster_columns = TRUE,
                                              clustering_distance_columns = "euclidean",


                                              ## clustering Rows
                                              cluster_rows = TRUE,
                                              clustering_distance_rows = clustDist,
                                              show_row_names = FALSE,

                                              ## Change legend title
                                              name = "z-score Int",

                                              ## naming Plot
                                              column_title = paste("Unsupervised Heat map", plotname, "clustered by", clustDist)

      )

  }
  ## Preparing Output object
  Output <- list()

  Output$raw <- ANOVAResults
  Output$Significant <- ANOVASignificantFeatures
  Output$Heatmap <- ANOVAHeatMap

  return(Output)
}

## Kruskal Test
## add roxygen comments
#' @title KruskalTest
#' @description This function performs differential expression analysis using a Kruskal test.
#' @param dataset The dataset to be tested
#' @param plotname The name to be displayed on created plots
#' @param clustDist The distance metric to be used for clustering in the Heatmap ("euclidean", "maximum", "man-hattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall")
#' @param method The method to be used for the Heatmap (unsupervised, supervised)
#' #' @param p.adj.method The method to be used for p-value adjustment ("BH", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
#' @return A list object containing the results of the Kruskal test, the significant features and a heatmap
#' @export
KruskalTest <- function(dataset, plotname= "", clustDist = "euclidean", method = "unsupervised", p.adj.method = "BH"){


  ## Initializing Dataframe
  datasetkruskal <- dataset

  ## Making the "Status" column into a factor
  datasetkruskal$Status <- datasetkruskal$Status %>% as.factor()


  ## Cheking if we are operating on Peptides or Proteins
  if("Protein" %in% colnames(dataset)){

    ## Running the kruskal
    kruskalResults <- datasetkruskal %>%
      ## Preparing Dataframe
      drop_na(Intensity) %>%
      dplyr::group_by(Protein) %>%
      rstatix::kruskal_test(Intensity ~ Status) %>%
      rstatix::adjust_pvalue(method = p.adj.method) %>%
      tidyr::separate(Protein, into = c("Uniprot", "Gene"), sep = "_", remove = F)

    KruskalSignificantFeatures <- kruskalResults %>%
      dplyr::filter(p.adj < 0.05) %>%
      dplyr::arrange(p.adj)


  }
  if("Peptide" %in% colnames(dataset)){

    ## Running the kruskal
    kruskalResults <- datasetkruskal %>%
      ## Preparing Dataframe
      drop_na(Intensity) %>%
      dplyr::group_by(Peptide) %>%
      rstatix::kruskal_test(Intensity ~ Status) %>%
      rstatix::adjust_pvalue(method = p.adj.method) %>%
      tidyr::separate(Peptide, into = c("Uniprot", "Gene"), sep = "_", remove = F)

    KruskalSignificantFeatures <- kruskalResults %>%
      dplyr::filter(p.adj < 0.05) %>%
      dplyr::arrange(p.adj)

    KruskalSignificantFeatures <- kruskalResults %>% dplyr::filter(p.adj < 0.05)

  }
  ## Making Heatmap of significant features
  if("Protein" %in% colnames(dataset)){

    ## supervised Heatmap
    if(method %in% c("supervised", "Supervised")){

      KurskalHeatMap <- HeatMap(dataset, PoIs = KruskalSignificantFeatures$Protein, plotname = plotname, method = "supervised", clustDist = clustDist, show_row_names = F)

    }
    ## unsupervised Heatmap
    if(method %in% c("unsupervised", "Unsupervised")){

      KurskalHeatMap <- HeatMap(dataset, PoIs = KruskalSignificantFeatures$Protein, plotname = plotname, method = "unsupervised", clustDist = clustDist, show_row_names = F)

    }

  }
  if("Peptide" %in% colnames(dataset)){

    ## supervised Heatmap
    if(method %in% c("supervised", "Supervised")){

      KurskalHeatMap <- HeatMap(dataset, PoIs = KruskalsignificantFeatures$Peptide, plotname = plotname, method = "supervised", clustDist = clustDist, show_row_names = F)

    }
    ## unsupervised Heatmap
    if(method %in% c("unsupervised", "Unsupervised")){

      KurskalHeatMap <- HeatMap(dataset, PoIs = KruskalsignificantFeatures$Peptide, plotname = plotname, method = "unsupervised", clustDist = clustDist, show_row_names = F)

    }

  }


  ##  Preparing Output object
  Output <- list()
  Output$raw <- kruskalResults
  Output$Significant <- KruskalSignificantFeatures
  Output$Heatmap <- KurskalHeatMap

  return(Output)
}

## Uni variate feature analysis

## Plots logistic regression of one feature
## add roxygen comments
#' @title LogisticRegressionSingleFeature
#' @description Plots the logistic regression of a single feature.
#' @param dataset The dataset to be tested
#' @param PoI The feature of interest
#' @return A plot object
#' @export
LogisticRegressionSingleFeature <- function(dataset, PoI){

  if("Protein" %in% colnames(dataset)){

    Status1 <- dataset$Status %>% unique()
    Status1 <- Status1[1]

    Status2 <- dataset$Status %>% unique()
    Status2 <- Status2[2]

    plotData <- dataset %>%
      filter(Protein %in% PoI) %>%
      ## setting up dummy variable
      mutate(DStatus = ifelse(Status == Status1 , 1, 0))}

  if("Peptide" %in% colnames(dataset)){

    Status1 <- dataset$Status %>% unique()
    Status1 <- Status1[1]

    Status2 <- dataset$Status %>% unique()
    Status2 <- Status2[2]

    plotData <- dataset %>%
      filter(Peptide %in% PoI) %>%
      ## setting up dummy variable
      mutate(DStatus = ifelse(Status == Status1 , 1, 0))}

  # Plot the data points and logistic regression line

  plot <- ggplot(plotData, aes(x = Intensity, y = DStatus)) +
    geom_point(aes(col = Status)) +  # Raw data points
    stat_smooth(method="glm", color="green", se=FALSE,
                method.args = list(family=binomial)) +
    labs(x = "Log2 Intensity", y = "Status") +
    ggtitle(paste("Logistic Regression for", PoI))+
    theme_light(base_size = 13)


  return(plot)

}

## Calculating the AUCs for a list of PoIs
## add roxygen comments
#' @title AUCs
#' @description Calculates the AUCs for a list of PoIs in the provided dataset.
#' @param dataset The dataset to be tested
#' @param PoIs A vector containing the Proteins of interest. Example: c(""Q8TF72_SHROOM3" , "Q9ULZ3_PYCARD") or (unique(dataset$Protein))
#' @param plotname The name to be displayed on created plots
#' @return A list object containing the results of the AUCs and a plot
#' @export
AUCs <- function(dataset, PoIs, plotname = "") {

  dataset <- dataset %>% dplyr::arrange(Status)
  ## Calculate difference betwen mean exrpession in status
  ### Check if dataset has more than 2 unique entries in Status
  if(length(unique(dataset$Status)) > 2){
    stop("More than 2 unique entries in Status")
  } else {

    Status1 <- unique(dataset$Status)[1]
    Status2 <- unique(dataset$Status)[2]

  }

  # Create an empty data frame to store results
  AUCResults <- data.frame(PoI = character(0), AUC = numeric(0))

  # Loop through each PoI (Protein or Peptide of Interest)
  for (i in 1:length(PoIs)) {
    PoI <- PoIs[i]

    # Filter dataset based on whether it's Protein or Peptide
    if ("Protein" %in% colnames(dataset)) {
      ROCData <- dplyr::filter(dataset, Protein %in% PoI)
    }

    if ("Peptide" %in% colnames(dataset)) {
      ROCData <- dplyr::filter(dataset, Peptide %in% PoI)
    }

    # Ensure Status is treated as a factor
    ROCData$Status <- ROCData$Status %>% as.factor()

    # Calculate ROC curve and AUC
    ROC <- pROC::roc(data = ROCData, response = "Status", predictor = "Intensity")

    # Extract AUC value
    AUC <- ROC$auc

    ## count number of observation in Stauts[1]
    n1 <- sum(ROCData$Status == Status1)
    ## count number of observation in Stauts[2]
    n2 <- sum(ROCData$Status == Status2)

    # Populate AUCResults data frame
    AUCResults[i, 1] <- PoI
    AUCResults[i, 2] <- AUC
    AUCResults[i, 3] <- n1
    AUCResults[i, 4] <- n2

    # Arrange AUCResults by decreasing AUC values
    AUCResults <- dplyr::arrange(AUCResults, desc(AUC))

  }

  colnames(AUCResults) <- c("PoI", "AUC", Status1, Status2)

  # Create AUC plot using ggplot2
  AUCPlot <- ggplot2::ggplot(data = AUCResults) +
    ggplot2::geom_col(mapping = ggplot2::aes(x = reorder(PoI, -AUC), y = AUC)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::ggtitle(label = plotname) +
    ggplot2::xlab(label = "PoI") +
    ggplot2::geom_hline(yintercept = 0.5)

  if("Protein" %in% colnames(dataset)){

    Diff <- dataset %>%
      group_by(Protein) %>%
      summarise(Diff = mean(Intensity[Status == Status1], na.rm = T) - mean(Intensity[Status == Status2], na.rm = T))

  }

  if("Peptide" %in% colnames(dataset)){

    Diff <- dataset %>%
      group_by(Peptide) %>%
      summarise(Diff = mean(Intensity[Status == Status1], na.rm = T) - mean(Intensity[Status == Status2], na.rm = T))

  }

  colnames(AUCResults)[1] <- ifelse("Protein" %in% colnames(dataset), "Protein", "Peptide")

  if("Protein" %in% colnames(dataset)){
    ## Combine the two dataframes
    Vulcanoplotdata <- merge(Diff, AUCResults, by = ifelse("Protein" %in% colnames(dataset),"Protein", "Peptide")) %>%
      ## Add AUC to 0.5 if value is less than 0.5
      dplyr::mutate(AUC = ifelse(AUC < 0.5, 0.5 + 0.5- AUC, AUC)) %>%
      ## make gene column for later plot annotation
      dplyr::mutate(Gene = str_split_i(Protein, "_", 2)) %>%
      ## make direction column for coloring
      dplyr::mutate(Direction = ifelse(Diff > 0, "Up", "Down")) %>%
      dplyr::arrange(desc(AUC))

  }

  if("Peptide" %in% colnames(dataset)){
    ## Combine the two dataframes
    Vulcanoplotdata <- merge(Diff, AUCResults, by = ifelse("Peptide" %in% colnames(dataset),"Peptide", "Peptide")) %>%
      ## Add AUC to 0.5 if value is less than 0.5
      dplyr::mutate(AUC = ifelse(AUC < 0.5, 0.5 + 0.5- AUC, AUC)) %>%
      ## make gene column for later plot annotation
      dplyr::mutate(Gene = str_split_i(Peptide, "_", 2)) %>%
      ## make direction column for coloring
      dplyr::mutate(Direction = ifelse(Diff > 0, "Up", "Down")) %>%
      dplyr::arrange(desc(AUC))

  }

  VulcanoPlot <- ## volcano plot of Results
    ggplot2::ggplot(Vulcanoplotdata, aes(x = Diff, y = AUC)) +
    ggplot2::geom_point(aes(col = ifelse(Direction == "Up", "blue", "red"))) +
    ggrepel::geom_text_repel(aes(label = Gene), box.padding = 0.5) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ## fix y axis from 0 to 1
    ggplot2::scale_y_continuous(limits = c(0.5, 1)) +
    ## rename axis
    ggplot2::xlab(paste("Fold Change in ", Status1)) +
    ggplot2::ylab("AUC") +
    ## add title
    ggplot2::ggtitle(plotname) +
    ## remove legend
    ggplot2::theme(legend.position = "none")

  ## create Histogram plot of AUCs
  Histogram <- ggplot2::ggplot(AUCResults, ggplot2::aes(x = AUC)) +
    ggplot2::geom_histogram(fill = "blue", alpha = 0.5, binwidth = 0.01) +
    ggplot2::theme_minimal() +
    ggplot2::ggtitle(paste("Histogram of AUCs for", plotname)) +
    ggplot2::xlab("AUC") +
    ggplot2::ylab("Count")

  # Create output list containing results and plot
  Output <- list()
  Output$results <- AUCResults
  Output$plot <- AUCPlot
  Output$VulcanoPlot <- VulcanoPlot
  Output$Histogram <- Histogram
  return(Output)
}

## Plots the ROC for a single feature (adressed by name) in the specified dataset
## Uses n-fold cross validation
## add roxygen comments
#' @title ROC
#' @description Plots the ROC for a single feature in the specified dataset.
#' @param dataset The dataset to be tested
#' @param PoI The feature of interest
#' @param nCross The number of iterations for the cross validation
#' @return A list object containing the results of the ROC calculations and a plot
#' @export
ROC <- function(dataset, PoI, plotname = ""){

  ## error of more than 2 classes
  if(length(unique(dataset$Status)) > 2){
    stop("ROC curve only works for two classes")
  }

  ## error if more than one PoI
  if(length(PoI) > 1){
    stop("ROC curve only works for one Protein of Interest")
  }

  ## filter data and prepare it for analysis
  MLData <- dataset %>%
    dplyr::filter(Protein %in% PoI) %>%
    dplyr::mutate(Status = as.numeric(as.factor(Status)))

  ## calculate glm model using
  glm_model <- stats::glm(Status ~ Intensity, data = MLData)

  ## calculate ROC curve using pROC
  roc_curve <- pROC::roc(predictor = stats::predict(glm_model, type = "response"),
                         response = MLData$Status)

  ## plot ROC curve using ggplot2
  AUC <- round(pROC::auc(roc_curve), 2)

  # Customizing aesthetics within ggroc
  roc_plot <- pROC::ggroc(roc_curve, size = 1.2, color = "black") + # Set line size and color
    theme_minimal() +                         # Use a minimal theme
    labs(
      title = "ROC Curve",                    # Add a title
      x = "1 - Specificity",
      y = "Sensitivity"
    ) +
    ggplot2::geom_abline(intercept = 1, slope = 1, linetype = "dashed", colour = "red", linewidth = 1.5) + # Add a diagonal line
    ggplot2::ggtitle(plotname, paste("AUC:", AUC))

  ## create output list
  output <- list(ROC_Plot = roc_plot,
                 AUC = AUC,
                 PoI = PoI)

  return(output)
}


## Creating Boxplots of Any Protein from a vector containing Protein names
## add roxygen comments
#' @title BoxPlotsFeatures
#' @description Creates boxplots of any protein from a vector containing protein or peptide names.
#' @param dataset The dataset to be tested
#' @param PoIs The vector containing the protein names
#' @param plotname The name to be displayed on created plots
#' @return A plot object
#' @export
BoxPlotsFeatures <- function(dataset, PoIs, plotname = "", pellet = "custom_vibrant") {

  dataset <- dataset %>% arrange(Status)

  ## generate colors
  colors <- assign_colors(unique(dataset$Status), palette = pellet)

  if ("Protein" %in% colnames(dataset)) {
    BoxPlotData <- dataset %>%
      dplyr::filter(Protein %in% PoIs)



    Boxplot <- ggplot2::ggplot(data = BoxPlotData, ggplot2::aes(x = Status, y = Intensity, col = Status)) +
      ## use colors
      ggplot2::scale_color_manual(values = colors) +
      ggplot2::geom_violin() +
      ggplot2::geom_jitter() +
      ggplot2::facet_wrap(~Protein) +
      ggplot2::ggtitle(paste("Expressions of Proteins", plotname)) +
      ggplot2::theme(legend.position = "none") +
      ggplot2::ylab("log 2 Intensity") +
      ggplot2::xlab("") +
      ggplot2::theme_light(base_size = 13)
  }

  if ("Peptide" %in% colnames(dataset)) {
    BoxPlotData <- dataset %>%
      dplyr::filter(Peptide %in% PoIs)

    Boxplot <- ggplot2::ggplot(data = BoxPlotData, ggplot2::aes(x = Status, y = Intensity, col = Status)) +
      ## use colors
      ggplot2::scale_color_manual(values = colors) +
      ggplot2::geom_violin() +
      ggplot2::geom_jitter() +
      ggplot2::facet_wrap(~Peptide) +
      ggplot2::ggtitle(paste("Expressions of Proteins", plotname)) +
      ggplot2::theme(legend.position = "none") +
      ggplot2::ylab("log 2 Intensity") +
      ggplot2::xlab("") +
      ggplot2::theme_light(base_size = 13)
  }

  return(Boxplot)
}

## Creating a Histogram of the Intensity of a given Protein
## add roxygen comments
#' @title PlotHistogram
#' @description Creates a histogram of the intensity of a given protein.
#' @param dataset The dataset to be tested
#' @param PoIs The protein or peptide of interest
#' @param plotname The name to be displayed on created plots
#' @return A plot object
#' @export
PlotHistogram <- function(dataset, PoIs, plotname = ""){

  for(k in 1:length(PoIs)){

    if("Protein" %in% colnames(dataset)){
      HistData <- dataset %>%
        filter(Protein %in% PoIs[k])

      Histogram <- ggplot(data = HistData, aes(x = Intensity, fill = Status)) +
        geom_histogram(alpha = 0.5, bins =  25) +
        ggtitle(paste("Histogram of Protein Intensity", PoIs[k], plotname)) +
        xlab("log2 Intensity")

      plot(Histogram)

    }
    if("Peptide" %in% colnames(dataset)){
      HistData <- dataset %>%
        filter(Peptide %in% PoIs[k])

      Histogram <- ggplot(data = HistData, aes(x = Intensity, fill = Status)) +
        geom_histogram(alpha = 0.5, bins =  25) +
        ggtitle(paste("Histogram of peptide Intensity", PoIs[k], plotname)) +
        xlab("log2 Intensity")

      plot(Histogram)
    }

  }

  return(Histogram)

}

## Multivariate Analysis

## Heatmap
## plots a Heatmap of the specified dataset
## add roxygen comments
#' @title HeatMap
#' @description veratile Heatmap function that can be used for unsupervised and supervised clustering of a list of proteins or peptides
#' @param dataset The dataset to be plottet
#' @param PoIs The list of proteins or peptides of interest
#' @param method The method to be used for the Heatmap (unsupervised, supervised)
#' @param clustDist The distance metric to be used for clustering in the Heatmap ("euclidean", "maximum", "man-hattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall")
#' @param plotname The name to be displayed on created plots
#' @param show_column_names Logical value indicating if column names should be shown
#' @param show_row_names Logical value indicating if row names should be shown
#' @param column_split A numeric indicating into how many clusters the columns should be split
#' @param row_split A numeric indicating into how many clusters the rows should be split
#' @param cluster_columns Logical value indicating if columns should be clustered at all
#' @param Annotations comumn names of the dataset that should be annotated
#' @param ColPalette The color palette to be used for the annotations (Specified using the names of a color palette from the RColorBrewer package): https://www.aptech.com/releases/gauss18/graphics-updates/color-brewer-palettes/
#' @param contColors The colors to be used for continuous annotations of numeric columns (Input is a vector of colors, example c("green", "red"))
#' @param fontsize The fontsize to be used for the annotations (numeric value)
#' @return A Heatmap object
#' @export
HeatMap <- function(dataset, PoIs, method = "unsupervised", clustDist = "euclidean", plotname = "", show_column_names = FALSE, show_row_names = FALSE, column_split = NULL, row_split = NULL, cluster_columns = TRUE, Annotations = NULL, ColPalette = "custom_vibrant", contColors = c("yellow", "purple"), fontsize = 8) {

  dataset <- dataset %>% dplyr::arrange(Status)

  ## creating heat map Data
  if ("Protein" %in% colnames(dataset)) {
    HeatMapData <- dataset %>%
      dplyr::filter(Protein %in% PoIs) %>%
      dplyr::group_by(Protein) %>%
      dplyr::mutate(Intensity = scale(Intensity, center = TRUE, scale = TRUE))

    ## Quantitative heat map data
    HeatMapDataQuant <- HeatMapData %>%
      tidyr::pivot_wider(names_from = "Protein", values_from = "Intensity") %>%
      column_to_rownames(var = "Sample") %>%
      dplyr::select(dplyr::contains("_")) %>%
      t() %>%
      as.matrix()

    ## clinical heat map data
    HeatMapDataClin <- HeatMapData %>%
      tidyr::pivot_wider(names_from = "Protein", values_from = "Intensity") %>%
      dplyr::select(!dplyr::contains("_"))

  } else if ("Peptide" %in% colnames(dataset)) {
    HeatMapData <- dataset %>%
      dplyr::filter(Peptide %in% PoIs) %>%
      dplyr::group_by(Peptide) %>%
      dplyr::mutate(Intensity = scale(Intensity, center = TRUE, scale = TRUE))

    ## Quantitative heat map data
    HeatMapDataQuant <- HeatMapData %>%
      tidyr::pivot_wider(names_from = "Peptide", values_from = "Intensity") %>%
      column_to_rownames(var = "Sample") %>%
      dplyr::select(dplyr::contains("_")) %>%
      t() %>%
      as.matrix()

    ## clinical heat map data
    HeatMapDataClin <- HeatMapData %>%
      tidyr::pivot_wider(names_from = "Peptide", values_from = "Intensity") %>%
      dplyr::select(!dplyr::contains("_"))
  } else {
    stop("Dataset must contain either 'Protein' or 'Peptide' column.")
  }


  # Function to generate a color palette for annotations
  generate_annotation_colors <- function(annotation_levels) {
    # Check if the annotation levels are numeric
    if (is.numeric(annotation_levels)) {
      # Generate a continuous color gradient for numeric data
      min_val <- min(annotation_levels, na.rm = TRUE)
      max_val <- max(annotation_levels, na.rm = TRUE)
      colors <- colorRamp2(c(min_val, max_val), c(contColors[1], contColors[2]))
    } else {
      # For categorical data, use qualitative color palette
      colors <- assign_colors(annotation_levels, ColPalette)
    }

    return(colors)
  }


  # Generate annotation colors and create annotation list dynamically
  annotation_list <- list()
  annotation_colors <- list()

  annotation_columns <- c("Status", Annotations)

  ## Generate palette for every annotation column
  for (col_name in annotation_columns) {
    annotation_levels <- unique(HeatMapDataClin[[col_name]])
    colors <- generate_annotation_colors(annotation_levels)
    annotation_list[[col_name]] <- HeatMapDataClin[[col_name]]
    annotation_colors[[col_name]] <- colors
  }

  # Create the annotation object
  Annotation <- ComplexHeatmap::HeatmapAnnotation(
    df = annotation_list,
    col = annotation_colors
  )

  if (tolower(method) == "supervised") {
    column_split <- HeatMapDataClin$Status

    HeatMapPlot <- ComplexHeatmap::Heatmap(
      HeatMapDataQuant,
      ## Annotation stuff
      top_annotation = Annotation,
      ## clustering specifics
      ## Clustering columns
      cluster_columns = cluster_columns,
      clustering_distance_columns = clustDist,
      show_column_names = show_column_names,
      ## clustering Rows
      cluster_rows = TRUE,
      clustering_distance_rows = clustDist,
      show_row_names = show_row_names,
      row_split = row_split,
      ## specify patient status as main cluster
      column_split = HeatMapDataClin$Status,
      ## Changing Legend title
      name = "z-score Int",
      ## naming Plot
      column_title = paste(plotname),
      ## adjust font size for row names
      row_names_gp = grid::gpar(fontsize = fontsize),
      column_names_gp = grid::gpar(fontsize = fontsize)
    )
  } else if (tolower(method) == "unsupervised") {

    HeatMapPlot <- ComplexHeatmap::Heatmap(
      HeatMapDataQuant,
      ## Annotations Stuff
      top_annotation = Annotation,
      ## clustering specifics
      ## Clustering columns
      cluster_columns = cluster_columns,
      clustering_distance_columns = "euclidean",
      show_column_names = show_column_names,
      column_split = column_split,
      ## clustering Rows
      cluster_rows = TRUE,
      clustering_distance_rows = clustDist,
      show_row_names = show_row_names,
      row_split = row_split,
      ## Change legend title
      name = "z-score Int",
      ## naming Plot
      column_title = paste(plotname),
      ## adjust font size for row names
      row_names_gp = grid::gpar(fontsize = fontsize),
      column_names_gp = grid::gpar(fontsize = fontsize)
    )
  } else {
    stop("Invalid method. Choose either 'supervised' or 'unsupervised'.")
  }

  return(HeatMapPlot)
}

## logistic regression for multiple features
## add roxygen comments
#' @title MultiLogisticRegression
#' @description Calculates the logistic regression for multiple features in the specified dataset.
#' @param dataset The dataset to be tested
#' @param PoIs A vector containing the Proteins of interest. Example: c("Q8TF72_SHROOM3" , "Q9ULZ3_PYCARD") or (unique(dataset$Protein))
#' @param nIterations The number of iterations for the cross validation
#' @return A list object containing the results of the logistic regression calculations and the ROC plot
#' @export
MultiLogisticRegression <- function(dataset, PoIs, nIterations = 10){
  ## Preparing Data
  Status1 <- unique(dataset$Status)[1]
  Status2 <- unique(dataset$Status)[2]

  Sensitivities <- data.frame()
  Specificities <- data.frame()
  ModelAUCs <- data.frame()

  if("Protein" %in% colnames(dataset)){

    MLData <- dataset %>%
      ## filtering for PoIs
      filter(Protein %in% PoIs) %>%
      ## creating dummy status variable
      mutate(DStatus = ifelse(Status == Status1, 1,0)) %>%
      drop_na()



    for(k in 1 : nIterations){
      ## Select Training and validation data
      ## Creating Training Data
      ## Select 80 % of the Samples the Entries in both statuses
      Data1 <- MLData %>%
        filter(Status == Status1) %>%
        group_by(Sample) %>%
        summarise() %>%
        sample_n(size = round(nrow(.)*0.8))

      Data2 <- MLData %>%
        filter(Status == Status2) %>%
        group_by(Sample) %>%
        summarise() %>%
        sample_n(size = round(nrow(.)*0.8))

      TrainingData <- MLData %>% filter(Sample %in% Data1$Sample | Sample %in% Data2$Sample)


      ValData <- MLData %>%
        anti_join(TrainingData, by = NULL)


      ## calculating multiple logistic regression model
      Model <- glm(DStatus ~ Protein:Intensity, data = TrainingData)

      ## Predicting Outcome

      Predictions <- predict(Model, newdata = ValData) %>%
        data.frame() %>%
        cbind(ValData)

      colnames(Predictions)[1] <- "Pred"

      ## defining ROC thersholds

      ROC <- roc(data = Predictions, response = DStatus, predictor = Pred)

      AUC <- ROC$auc

      sensitivity <- ROC$sensitivities %>% data.frame()
      specificity <- ROC$specificities %>% data.frame()

      ## First iteration
      if (k == 1){
        Sensitivities <- sensitivity$. %>% data.frame()
        Specificities <- specificity$. %>% data.frame()
      }
      ## Other iterations
      if (k!= 1){
        if(nrow(sensitivity) == nrow(Sensitivities)) {
          Sensitivities <- Sensitivities %>% cbind(sensitivity)
        }
        if (nrow(sensitivity) != nrow(Sensitivities)){
          k <- k-1
          next
        }
        if(nrow(specificity) == nrow(Specificities)){
          Specificities <- Specificities %>% cbind(specificity)
        }

      }

      ModelAUCs[k,1] <- AUC
    }



  }
  if("Peptide" %in% colnames(dataset)){

    for(k in 1 : nIterations){
      ## Select Training and validation data

      MLData <- dataset %>%
        ## filtering for PoIs
        filter(Peptide %in% PoIs) %>%
        ## creating dummy status variable
        mutate(DStatus = ifelse(Status == Status1, 1,0))

      ## Creating Training Data
      ## Here we select 80 % on the Entires in both stati
      Data1 <- MLData %>%
        filter(Status == Status1) %>%
        group_by(Sample) %>%
        summarise() %>%
        sample_n(size = round(nrow(.)*0.8))

      Data2 <- MLData %>%
        filter(Status == Status2) %>%
        group_by(Sample) %>%
        summarise() %>%
        sample_n(size = round(nrow(.)*0.8))

      TrainingData <- MLData %>% filter(Sample %in% Data1$Sample | Sample %in% Data2$Sample)


      ValData <- MLData %>%
        anti_join(TrainingData, by = NULL)


      ## calculating multiple logistic regression model
      Model <- glm(DStatus ~ Peptide:Intensity, data = TrainingData)

      ## Predicting Outcome

      Predictions <- predict(Model, newdata = ValData) %>%
        data.frame() %>%
        cbind(ValData)

      colnames(Predictions)[1] <- "Pred"


      ROC <- roc(data = Predictions, response = Status, predictor = Pred)

      AUC <- ROC$auc

      sensitivity <- ROC$sensitivities %>% data.frame()
      specificity <- ROC$specificities %>% data.frame()

      if (k == 1){

        Sensitivities <- sensitivity$. %>% data.frame()
        Specificities <- specificity$. %>% data.frame()
      }
      if (k!= 1){
        if(nrow(sensitivity) == nrow(Sensitivities)) {
          Sensitivities <- Sensitivities %>% cbind(sensitivity)
        }
        if (nrow(sensitivity) != nrow(Sensitivities)){
          k <- k-1
          next
        }
        if(nrow(specificity) == nrow(Specificities)){
          Specificities <- Specificities %>% cbind(specificity)
        }

      }

      ModelAUCs[k,1] <- AUC
    }
  }

  ## Plotting a given model Results
  ## Summarizing results of logistical regressions
  ## Calculate mean AUC of the model
  MeanAUC = colMeans(ModelAUCs, na.rm = TRUE) %>% round(digits = 3)


  PloTPoIs <- PoIs %>% str_split_i(pattern = "_", 2)
  Sensummary <- Sensitivities %>% t() %>% data.frame()
  MeanSen <- colMeans(Sensummary, na.rm = TRUE)
  minSen <- apply(Sensummary,2, min)
  maxSen <- apply(Sensummary,2, max)
  sdSen <- apply(Sensummary, 2, sd)
  TValue <- qt(p = 0.05, df = nIterations-1)
  CISen <- TValue*(sdSen/sqrt(nIterations))
  Spesummary <- Specificities %>% t() %>% data.frame()
  MeanSpe <- colMeans(Spesummary, na.rm = TRUE)

  PlotData <- rbind(MeanSen,minSen,maxSen, MeanSpe, CISen) %>% data.frame() %>% t() %>% data.frame()


  plot <- ggplot(data = PlotData)+
    geom_ribbon(aes(x =1-MeanSpe, ymin = MeanSen - CISen, ymax = MeanSen + CISen), alpha = 0.5) +
    geom_point(aes(x = 1-MeanSpe, y = MeanSen))+
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed")+
    ggtitle(paste("Multiple logistic regression; mean AUC = ", MeanAUC), paste("PoIs = ", paste(PloTPoIs, " ", collapse = ","), "   |   ", nIterations, "Iterations")) +
    xlab("False positive rate") +
    ylab("Sensitivity") +
    xlim(0,1) +
    ylim(0,1)+
    theme_light(base_size = 13)


  ## create output object
  Output <- list()
  Output$Results <- ModelAUCs
  Output$ROCPlot <- plot

}


## PCA
## The Function selects Proteins in the combined clinical and quantitative dataframe based on
##  colnames contains("_"), Make sure non of the clinical names have an Underscore in their name
## add roxygen comments
#' @title PCA
#' @description Principal compenent analysis for the specified dataset.
#' @param dataset The dataset to be tested
#' @param nPcs The number of principal components to be calculated
#' @param plotname The name to be displayed on created plots
#' @param PoIs The protein or peptide of interest
#' @param plotTopNLoading Logical value indicating if the top 3 loadings should be plotted
#' @param topNLoading The number of top loadings to be plotted
#' @param show_ellipse Logical value indicating if the 95% confidence ellipse should be plotted
#' @return A list object containing the results of the PCA calculations and the PCA plot
#' @export
PCA <- function(dataset, nPcs = 3, plotname = "PCA", PoIs = "", plotTopNLoading = TRUE, topNLoading = 3, show_ellipse = F, ColPalette = "custom_vibrant"){

  if("Protein" %in% colnames(dataset)){

    PCAData <- dataset %>%
      ## making sure every column has entries
      BiomarkR::NaCutoff(0.1) %>%
      ## imputing missing values
      BiomarkR::ImputeFeatureIntensity(method = "halfmin") %>%
      ## Normalize on sample
      BiomarkR::normalizeIntensityOnSample(plot = F) %>%
      ## Selecting necessary columns
      dplyr::select("Protein", "Intensity", "Status", "Sample") %>%
      ## pivoting wider
      tidyr::pivot_wider(names_from = "Protein", values_from = "Intensity")

    ## Calculating PCA
    PCA <- pcaMethods::pca(select(PCAData, contains("_")), nPcs = 3)


  }

  if("Peptide" %in% colnames(dataset)){

    PCAData <- dataset %>%
      ## making sure every column has entries
      BiomarkR::NaCutoff(0.1) %>%
      ## imputing missing values
      BiomarkR::ImputeFeatureIntensity(method = "halfmin") %>%
      ## Normalize on sample
      BiomarkR::normalizeIntensityOnSample(plot = F) %>%
      ## Selecting necessary columns
      dplyr::select("Peptide", "Intensity", "Status", "Sample") %>%
      ## pivoting wider
      tidyr::pivot_wider(names_from = "Peptide", values_from = "Intensity")

    ## Calculating PCA
    PCA <- pcaMethods::pca(select(PCAData, contains("_")), nPcs = 3)

  }

  ## Plotting Results
  ## Faceted Score Plot
  PCAPlotData <- merge(PCA@scores, PCAData , by = 0)
  ## splitting into three dataframes for plotting reasons
  ## Only PC1 and PC 2
  PCAPlotData12 <- PCAPlotData %>%
    dplyr::select(-c("PC3")) %>%
    dplyr::mutate(facet = "PC1 (x) on PC2 (y)")
  colnames(PCAPlotData12)[2:3] <- c("PlotPC1", "PlotPC2")

  ## Only PC1 and PC 3
  PCAPlotData13 <- PCAPlotData %>%
    dplyr::select(-c("PC2")) %>%
    dplyr::mutate(facet = "PC1 (x) on PC3 (y)")
  colnames(PCAPlotData13)[2:3] <- c("PlotPC1", "PlotPC2")

  ## Only PC2 and PC 3
  PCAPlotData23 <- PCAPlotData %>%
    dplyr::select(-c("PC1")) %>%
    dplyr::mutate(facet = "PC2 (x) on PC3 (y)")
  colnames(PCAPlotData23)[2:3] <- c("PlotPC1", "PlotPC2")

  ## Recombining Dataframes in long format
  PCAPlotData <- rbind(PCAPlotData12, PCAPlotData13, PCAPlotData23)

  ## get variance explained
  VarianceExplained <- PCA@R2

  ## assign color
  if(! is.numeric(PCAPlotData$Status)){
    colors <- BiomarkR::assign_colors(unique(PCAPlotData$Status), palette = ColPalette)
    }




  ## making Score plot
  scorePlot <- ggplot(PCAPlotData, aes(x = PlotPC1, y = PlotPC2, colour = Status)) +
    geom_jitter() +
    ## use colors
    facet_wrap(~ PCAPlotData$facet, axis.labels = "all") +
    ggtitle(plotname, paste0(
      "PC1: ",round(VarianceExplained[1],2)*100, " %, ",
      "PC2: ",round(VarianceExplained[2],2)*100, " %, ",
      "PC3: ",round(VarianceExplained[3],2)*100, " % ")
    ) +
    ylab("") +
    ## name x axis properly
    xlab("") +
    theme_light(base_size = 13) +
    ## make facet title background white and text black
    theme(strip.background = element_rect(fill = "white"),
          strip.text = element_text(colour = "black")) +
    if(show_ellipse == T){
      ggplot2::stat_ellipse(aes(x = PlotPC1, y = PlotPC2, color = Status), level = 0.95, linetype = 1)
    }

  ## change the colors if Status is numeric
  if(is.numeric(PCAPlotData$Status)){
    scorePlot <- scorePlot + scale_color_viridis_c()
  } else {
    scorePlot <- scorePlot + scale_color_manual(values = colors)
  }

  ScorePlot_12 <- ggplot(PCAPlotData12, aes(x = PlotPC1, y = PlotPC2, colour = Status)) +
    geom_jitter() +
    ggtitle(plotname) +
    ylab(paste0("PC2 (", round(VarianceExplained[2],2)*100, " %)")) +
    ## name x axis properly
    xlab(paste0("PC1 (", round(VarianceExplained[1],2)*100, " %)")) +
    theme_light(base_size = 13) +
    ## make facet title background white and text black
    theme(strip.background = element_rect(fill = "white"),
          strip.text = element_text(colour = "black")) +
    if(show_ellipse == T){
      geom_ellipse(aes(x = PlotPC1, y = PlotPC2, color = Status), level = 0.95, linetype = 2)
    }

  ## change the colors if Status is numeric
  if(is.numeric(PCAPlotData$Status)){
    ScorePlot_12 <- ScorePlot_12 + scale_color_viridis_c()
  } else {
    ScorePlot_12 <- ScorePlot_12 + scale_color_manual(values = colors)
  }


  ## Making 3D score plot
  ## extracting scores for the first three principal components
  # Create the 3D Score Plot
  ScorePlot3DData <- PCA@scores %>%
    ## Combine with sample information
    cbind(PCAData)

  ## Make 3D score plot using plotly
  # Check if Status is a factor
  ScorePlot3D <- if (!is.numeric(ScorePlot3DData$Status)) {
    # Use categorical color mapping
    plotly::plot_ly(
      data = ScorePlot3DData,
      x = ~PC1,
      y = ~PC2,
      z = ~PC3,
      color = ~Status,  # Categorical color mapping
      colors = assign_colors(unique(ScorePlot3DData$Status)), # Define colors for levels
      type = "scatter3d",
      mode = "markers",
      marker = list(size = 5)
    )
  } else {
    # Use continuous color mapping
    plotly::plot_ly(
      data = ScorePlot3DData,
      x = ~PC1,
      y = ~PC2,
      z = ~PC3,
      marker = list(
        size = 5,
        color = ~Status,   # Continuous color scale
        colorscale = "Viridis", # Apply Viridis color scale
        showscale = TRUE       # Show the color legend
      ),
      type = "scatter3d",
      mode = "markers"
    )
  }
    ## Add a plot title
  ScorePlot3D <- ScorePlot3D %>% plotly::layout(title = paste("3D Score plot", plotname))


  ## loading plot
  ## How many Proteins should be on the loading plot

  PCALoadings <- PCA@loadings
  Top5Data <- PCALoadings %>% data.frame() %>% abs()

  Top5PC1 <- Top5Data %>%
    dplyr::arrange(-PC1) %>%
    utils::head(topNLoading) %>% row.names()

  Top5PC2 <- Top5Data %>%
    dplyr::arrange(-PC2) %>%
    utils::head(topNLoading) %>% row.names()

  Top5PC3 <- Top5Data %>%
    dplyr::arrange(-PC3) %>%
    utils::head(topNLoading) %>% row.names()

  ## get loadings for Top5Proteins
  Top5Proteins <- c(Top5PC1, Top5PC2, Top5PC3)

  Top5Loadings <- PCALoadings %>% data.frame() %>% rownames_to_column(var = "Protein") %>%
    dplyr::filter(Protein %in% Top5Proteins) %>% mutate(var = "Top5")

  PoILoadings <- PCALoadings %>% data.frame() %>% rownames_to_column(var = "Protein") %>%
    dplyr::filter(Protein %in% PoIs) %>% mutate(var = "PoIs")

  ifelse(plotTopNLoading == T,
         LoadingPlotData <- rbind(Top5Loadings,PoILoadings),
         LoadingPlotData <- PoILoadings)

  ## Only PC1 and PC 2
  LoadingPlotData12 <- LoadingPlotData %>%
    dplyr::select(-c("PC3")) %>%
    dplyr::mutate(facet = "PC1 (x) on PC2 (y)")
  colnames(LoadingPlotData12)[2:3] <- c("PlotPC1", "PlotPC2")

  ## Only PC1 and PC 3
  LoadingPlotData13 <- LoadingPlotData %>%
    dplyr::select(-c("PC2")) %>%
    dplyr::mutate(facet = "PC1 (x) on PC3 (y)")
  colnames(LoadingPlotData13)[2:3] <- c("PlotPC1", "PlotPC2")

  ## Only PC2 and PC 3
  LoadingPlotData23 <- LoadingPlotData %>%
    dplyr::select(-c("PC1")) %>%
    dplyr::mutate(facet = "PC2 (x) on PC3 (y)")
  colnames(LoadingPlotData23)[2:3] <- c("PlotPC1", "PlotPC2")

  ## Recombining Dataframes in long format
  LoadingPlotData <- rbind(LoadingPlotData12, LoadingPlotData13, LoadingPlotData23) %>%
    dplyr::mutate(Gene = str_split_i(Protein, "_", 2))

  loadingPlot <- ggplot2::ggplot(LoadingPlotData, aes(x = PlotPC1, y = PlotPC2, label = Protein, colour = var)) +
    ggplot2::facet_wrap(~ facet) +
    ggplot2::geom_segment(aes(x = 0, y = 0, xend = PlotPC1, yend = PlotPC2), arrow = arrow(type = "closed", length = unit(0.05, "inches"))) +  # Lines from origin
    ggrepel::geom_text_repel(aes(label = Gene), size = 4) +  # Optional: Add text labels
    ggplot2::ggtitle(plotname, paste0(
      "PC1: ",round(VarianceExplained[1],2)*100, " %, ",
      "PC2: ",round(VarianceExplained[2],2)*100, " %, ",
      "PC3: ",round(VarianceExplained[3],2)*100, " % ")
    ) +
    ggplot2::ylab("") +
    ## name x axis properly
    ggplot2::xlab("") +
    ggplot2::theme_light(base_size = 13) +
    ## make facet title background white and text black
    ggplot2::theme(strip.background = element_rect(fill = "white"),
                   strip.text = element_text(colour = "black"))


  ## Creating Output Object
  Output <- list()
  Output$ScorePlot2D <- scorePlot
  Output$Loadingplot <- loadingPlot
  Output$PCA <- PCA
  Output$ScorePlot_12 <- ScorePlot_12
  Output$ScorePlot3D <- ScorePlot3D


  return(Output)
}

## Umap for visualization
## add roxygen comments
#' @title UMAP
#' @description Uniform Manifold Approximation and Projection for the specified dataset.
#' @param dataset The dataset to be plotted
#' @param plotname The name to be displayed on created plots
#' @param show_ellipse Logical value indicating if the 95% confidence ellipse should be plotted
#' @return A list object containing the results of the UMAP calculations and the UMAP plot
#' @export
UMAP <- function(dataset, plotname = "", show_ellipse = F,ColPalette = "custom_vibrant") {

  ## error of missing values in Intensity column
  if (base::any(base::is.na(dataset$Intensity))) {
    base::stop("Missing values in dataset, please impute")
  }

  if ("Protein" %in% base::colnames(dataset)) {

    # Reshape the dataset and normalize intensity
    dataset <- dataset %>%
      tidyr::pivot_wider(names_from = Protein, values_from = Intensity) %>%
      tidyr::pivot_longer(cols = dplyr::contains("_"), names_to = "Protein", values_to = "Intensity") %>%
      normalizeIntensityOnSample(plot = FALSE)

    # Prepare clinical data for UMAP
    UmapDataClin <- dataset %>%
      normalizeIntensityOnSample(plot = FALSE) %>%
      tidyr::pivot_wider(names_from = Protein, values_from = Intensity) %>%
      dplyr::select(!dplyr::contains("_"))

    # Prepare quantitative data for UMAP
    UmapDataQuant <- dataset %>%
      normalizeIntensityOnSample(plot = FALSE) %>%
      tidyr::pivot_wider(names_from = Protein, values_from = Intensity) %>%
      tibble::column_to_rownames("Sample") %>%
      dplyr::select(dplyr::contains("_"))
  }

  if ("Peptide" %in% base::colnames(dataset)) {

    # Reshape the dataset and normalize intensity
    dataset <- dataset %>%
      tidyr::pivot_wider(names_from = Peptide, values_from = Intensity) %>%
      tidyr::pivot_longer(cols = dplyr::contains("_"), names_to = "Peptide", values_to = "Intensity") %>%
      normalizeIntensityOnSample(plot = FALSE)

    # Prepare clinical data for UMAP
    UmapDataClin <- dataset %>%
      tidyr::pivot_wider(names_from = Peptide, values_from = Intensity) %>%
      dplyr::select(!dplyr::contains("_"))

    # Prepare quantitative data for UMAP
    UmapDataQuant <- dataset %>%
      normalizeIntensityOnSample(plot = FALSE) %>%
      tidyr::pivot_wider(names_from = Peptide, values_from = Intensity) %>%
      tibble::column_to_rownames("Sample") %>%
      dplyr::select(dplyr::contains("_"))
  }

  set.seed(2105)  # Use any fixed number for reproducability

  # Perform UMAP with three components
  UMAP <- umap::umap(UmapDataQuant, n_components = 3)

  UMAPPlotData <- base::data.frame(Sample = base::row.names(UMAP$layout),
                                   UMAP1 = UMAP$layout[, 1],
                                   UMAP2 = UMAP$layout[, 2],
                                   UMAP3 = UMAP$layout[, 3])

  ## Only UMAP1 and UMAP2
  UMAPPlotData12 <- UMAPPlotData %>%
    dplyr::select(-c("UMAP3")) %>%
    dplyr::mutate(facet = "UMAP1 (x) on UMAP2 (y)")
  base::colnames(UMAPPlotData12)[2:3] <- c("PlotUMAP1", "PlotUMAP2")

  ## Only UMAP1 and UMAP3
  UMAPPlotData13 <- UMAPPlotData %>%
    dplyr::select(-c("UMAP2")) %>%
    dplyr::mutate(facet = "UMAP1 (x) on UMAP3 (y)")
  base::colnames(UMAPPlotData13)[2:3] <- c("PlotUMAP1", "PlotUMAP2")

  ## Only UMAP2 and UMAP3
  UMAPPlotData23 <- UMAPPlotData %>%
    dplyr::select(-c("UMAP1")) %>%
    dplyr::mutate(facet = "UMAP2 (x) on UMAP3 (y)")
  base::colnames(UMAPPlotData23)[2:3] <- c("PlotUMAP1", "PlotUMAP2")

  ## Recombining Dataframes in long format
  UMAPPlotData2D <- base::rbind(UMAPPlotData12, UMAPPlotData13, UMAPPlotData23) %>%
    base::merge(UmapDataClin, by = "Sample")

  ## assign colors to the status variable
  colors <- assign_colors(unique(dataset$Status),palette = ColPalette)

  ## Making Score plot
  UMAPPlot2D <- ggplot2::ggplot(UMAPPlotData2D, ggplot2::aes(x = PlotUMAP1, y = PlotUMAP2, colour = Status)) +
    ggplot2::geom_jitter() +
    ## use colors
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::facet_wrap(~facet, axis.labels = "all") +
    ggplot2::ylab("") +
    ggplot2::xlab("") +
    ggplot2::theme_light(base_size = 13) +
    ggplot2::theme(strip.background = ggplot2::element_rect(fill = "white"),
                   strip.text = ggplot2::element_text(colour = "black")) +
    ggtitle(plotname)+
    ## add ellipses
    if(show_ellipse == T){
      ggplot2::geom_ellipse(ggplot2::aes(x = PlotUMAP1, y = PlotUMAP2, color = Status), level = 0.95, linetype = 2)
    }


  UMAPPlotData12 <- UMAPPlotData2D %>% dplyr::filter(facet == "UMAP1 (x) on UMAP2 (y)")

  UMAPPlot_12 <- ggplot(UMAPPlotData12, aes(x =PlotUMAP1 , y = PlotUMAP2, colour = Status)) +
    ggplot2::scale_color_manual(values = colors) +
    geom_jitter() +
    ggtitle(plotname) +
    ylab("UMAP2") +
    ## name x axis properly
    xlab("UMAP1") +
    theme_light(base_size = 13) +
    ## make facet title background white and text black
    theme(strip.background = element_rect(fill = "white"),
          strip.text = element_text(colour = "black")) +
    if(show_ellipse == T){
      geom_ellipse(aes(x = PlotPC1, y = PlotPC2, color = Status), level = 0.95, linetype = 2)
    }

  ## 3D Plot Data
  UMAPPlotData <- base::data.frame(Sample = base::row.names(UMAP$layout),
                                   UMAP1 = UMAP$layout[, 1],
                                   UMAP2 = UMAP$layout[, 2],
                                   UMAP3 = UMAP$layout[, 3]) %>%
    base::merge(UmapDataClin, by = "Sample")

  # 3D Plot using plotly
  UMAPPlot3D <- plotly::plot_ly(UMAPPlotData, x = ~UMAP1, y = ~UMAP2, z = ~UMAP3, color = ~Status, type = "scatter3d", mode = "markers",colors = colors[unique(dataset$Status)]) %>%
    plotly::layout(title = base::paste(plotname),
                   scene = list(xaxis = list(title = "UMAP1"),
                                yaxis = list(title = "UMAP2"),
                                zaxis = list(title = "UMAP3")))

  # Creating Output Object
  Output <- base::list(UMAP = UMAP,
                       UMAPPlot2D = UMAPPlot2D,
                       UMAPPlot3D = UMAPPlot3D,
                       UMAPPlot_12 = UMAPPlot_12)

  return(Output)
}



## tSNE
## add roxygen comments
#' @title tSNE
#' @description t-distributed Stochastic Neighbor Embedding for the specified dataset.
#' @param dataset The dataset to be plotted
#' @param plotname The name to be displayed on created plots
#' @param show_ellipse Logical value indicating if the 95% confidence ellipse should be plotted
#' @return A list object containing the results of the tSNE calculations and the tSNE plot
#' @export
tSNE <- function(dataset, plotname = "", show_ellipse = F, ColPalette = "custom_vibrant") {

  ## Error if missing values in object
  if (any(is.na(dataset$Intensity))) {
    stop("Missing values in dataset, please impute")
  }

  if ("Protein" %in% colnames(dataset)) {

    ## Data management
    tSNEData <- dataset %>%
      tidyr::pivot_wider(names_from = Protein, values_from = Intensity) %>%
      dplyr::select(Sample, dplyr::contains("_")) %>%
      tibble::column_to_rownames(var = "Sample") %>%
      as.matrix()

    ClinicalData <- dataset %>%
      tidyr::pivot_wider(names_from = Protein, values_from = Intensity) %>%
      dplyr::select(Sample, !dplyr::contains("_"))
  }

  if ("Peptide" %in% colnames(dataset)) {

    ## Data management
    tSNEData <- dataset %>%
      tidyr::pivot_wider(names_from = Peptide, values_from = Intensity) %>%
      dplyr::select(Sample, dplyr::contains("_")) %>%
      tibble::column_to_rownames(var = "Sample") %>%
      as.matrix()

    ClinicalData <- dataset %>%
      tidyr::pivot_wider(names_from = Peptide, values_from = Intensity) %>%
      dplyr::select(Sample, !dplyr::contains("_"))
  }

  set.seed(2105)  # Use any fixed number for reproducability

  # Get the number of samples
  n_samples <- nrow(tSNEData)

  # Set perplexity to a suitable value (e.g., less than n_samples / 3)
  adjusted_perplexity <- min(30, floor(n_samples / 4))

  # Run t-SNE with adjusted perplexity
  tSNEResult <- Rtsne::Rtsne(tSNEData, dims = 3, perplexity = adjusted_perplexity)

  ## Make plot data
  tSNEPlotData <- data.frame(Sample = row.names(tSNEData),
                             Dim1 = tSNEResult$Y[, 1],
                             Dim2 = tSNEResult$Y[, 2],
                             Dim3 = tSNEResult$Y[, 3])

  ## Only Dim1 and Dim2
  tSNEPlotData12 <- tSNEPlotData %>%
    dplyr::select(-c("Dim3")) %>%
    dplyr::mutate(facet = "Dim1 (x) on Dim2 (y)") %>%
    dplyr::mutate(Sample = as.character(Sample))
  colnames(tSNEPlotData12)[2:3] <- c("PlotDim1", "PlotDim2")

  ## Only Dim1 and Dim3
  tSNEPlotData13 <- tSNEPlotData %>%
    dplyr::select(-c("Dim2")) %>%
    dplyr::mutate(facet = "Dim1 (x) on Dim3 (y)") %>%
    dplyr::mutate(Sample = as.character(Sample))
  colnames(tSNEPlotData13)[2:3] <- c("PlotDim1", "PlotDim2")

  ## Only Dim2 and Dim3
  tSNEPlotData23 <- tSNEPlotData %>%
    dplyr::select(-c("Dim1")) %>%
    dplyr::mutate(facet = "Dim2 (x) on Dim3 (y)") %>%
    dplyr::mutate(Sample = as.character(Sample))
  colnames(tSNEPlotData23)[2:3] <- c("PlotDim1", "PlotDim2")

  ## Recombine Dataframes in long format
  tSNEPlotData <- dplyr::bind_rows(tSNEPlotData12, tSNEPlotData13, tSNEPlotData23) %>%
    base::merge(ClinicalData, by = "Sample")

  ## assign colors to the status variable
  colors <- assign_colors(unique(dataset$Status),palette = ColPalette)

  ## Making 2D score plot
  tSNEPlot2D <- ggplot2::ggplot(tSNEPlotData, ggplot2::aes(x = PlotDim1, y = PlotDim2, colour = Status)) +
    ## use colors
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::geom_jitter() +
    ggplot2::facet_wrap(~facet, labeller = "label_value") +
    ggplot2::ylab("") +
    ggplot2::xlab("") +
    ggplot2::theme_light(base_size = 13) +
    ggplot2::theme(strip.background = ggplot2::element_rect(fill = "white"),
                   strip.text = ggplot2::element_text(colour = "black")) +
    ggplot2::ggtitle(plotname) +
    ## optional ellipse
    if(show_ellipse == T){
      ggplot2::geom_ellipse(ggplot2::aes(x = PlotDim1, y = PlotDim2, color = Status), level = 0.95, linetype = 2)
    }

  ## make plot with only the first 2 tSNE components
  tSNEPlotData_12 <- tSNEPlotData %>% dplyr::filter(facet == "Dim1 (x) on Dim2 (y)")

  tSNEPlot_12 <- ggplot(tSNEPlotData_12, aes(x =PlotDim1 , y = PlotDim2, colour = Status)) +
    ggplot2::scale_color_manual(values = colors) +
    geom_jitter() +
    ggtitle(plotname) +
    ylab("Dim2") +
    ## name x axis properly
    xlab("Dim1") +
    theme_light(base_size = 13) +
    ## make facet title background white and text black
    theme(strip.background = element_rect(fill = "white"),
          strip.text = element_text(colour = "black")) +
    if(show_ellipse == T){
      geom_ellipse(aes(x = PlotPC1, y = PlotPC2, color = Status), level = 0.95, linetype = 2)
    }



  ## Making 3D score plot
  PlotData3ED <- data.frame(Sample = row.names(tSNEData),
                            Dim1 = tSNEResult$Y[, 1],
                            Dim2 = tSNEResult$Y[, 2],
                            Dim3 = tSNEResult$Y[, 3]) %>%
    base::merge(ClinicalData, by = "Sample")

  tSNEPlot3D <- PlotData3ED %>%
    plotly::plot_ly(x = ~Dim1, y = ~Dim2, z = ~Dim3, color = ~Status, type = "scatter3d", mode = "markers",colors = colors[unique(dataset$Status)]) %>%
    plotly::layout(scene = list(xaxis = list(title = "Dim1"),
                                yaxis = list(title = "Dim2"),
                                zaxis = list(title = "Dim3")),
                   title = plotname)

  ## Make output list
  Output <- list(tSNE = tSNEResult,
                 ScorePlot2D = tSNEPlot2D,
                 ScorePlot3D = tSNEPlot3D,
                 tSNEPlot_12 = tSNEPlot_12)

  return(Output)
}



## EffectAnalysis
## Correlates clinical Variables with principal components to visualize the main effects in the Dataset
## add roxygen comments
#' @title EffectAnalysis
#' @description Calculates the correlation of clinical variables with principal components.
#' @param dataset The dataset to be tested
#' @return A list object containing the results of the effect analysis
#' @export
EffectAnalysis <- function(dataset, plotname = ""){

  PCA <- PCA(dataset)

  ## extract first 3 PCs
  PCs <- PCA$PCA@scores[,1:3] %>% data.frame

  ## Extract explained variances
  ExplainedVariances <- PCA$PCA@R2

  ##
 if ("Protein" %in% colnames(dataset)) {
   CorrelationData <- dataset  %>%
     tidyr::pivot_wider(names_from = Protein, values_from = Intensity) %>%
     dplyr::select(!dplyr::contains("_")) %>%
     dplyr::filter(Sample %in% dataset$Sample) %>%
     cbind(PCs)
  }

  if("Peptide" %in% colnames(dataset)){
    CorrelationData <- dataset  %>%
      tidyr::pivot_wider(names_from = Peptide, values_from = Intensity) %>%
      dplyr::select(!dplyr::contains("_")) %>%
      dplyr::filter(Sample %in% dataset$Sample) %>%
      cbind(PCs)
  }

  ## Replace " " with "." in colnames(dataset)
  colnames(CorrelationData) <- gsub(" ", ".", colnames(CorrelationData))

  ## make all columns numeric
  for(i in 1 : ncol(CorrelationData)){

    if(!is.numeric(CorrelationData[[i]])){
      ## get unique Entries of the column
      uniqueEntries <- unique(CorrelationData[[i]])
      ## make factor
      CorrelationData[[i]] <- as.numeric(factor(CorrelationData[[i]], levels = uniqueEntries))
    }
  }

  ## calculate correlations between PCs and clinical variables
  CorResults <- CorrelationData %>%
    ungroup() %>%
    dplyr::select(-Sample) %>%
    rstatix::cor_test(vars = dplyr::contains("PC"), vars2 = !dplyr::contains("PC")) %>%
    dplyr::mutate(var1 = ifelse(var1 == "PC1", paste("PC1 (", round(ExplainedVariances[1], digits = 2)*100, " %)", sep = ""), var1),
                  var1 = ifelse(var1 == "PC2", paste("PC2 (", round(ExplainedVariances[2], digits = 2)*100, " %)", sep = ""), var1),
                  var1 = ifelse(var1 == "PC3", paste("PC3 (", round(ExplainedVariances[3], digits = 2)*100, " %)", sep = ""), var1)) %>%
    rstatix::adjust_pvalue()

  ## plot results

  EffectPlot <- ggplot2::ggplot(CorResults) +
    ggplot2::geom_tile(ggplot2::aes(x = var1, y = var2, fill = cor)) +
    ## use color gradient
    ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                                  limits = c(-1, 1)) +
    ## highlight cells with significant p-values
    ggplot2::geom_text(ggplot2::aes(x = var1, y = var2, label = ifelse(p.adj < 0.05, "*", "")), size = 10) +
    ## add second star if p.adj < 0.01
    ggplot2::geom_text(ggplot2::aes(x = var1, y = var2, label = ifelse(p.adj < 0.01, "  *", "")), size = 10) +
    ## add grid of cells
    ggplot2::geom_vline(xintercept = seq(0.5, nrow(CorResults) + 0.5, 1), color = "white") +
    ggplot2::geom_hline(yintercept = seq(0.5, nrow(CorResults) + 0.5, 1), color = "white") +
    ## rename axis
    ggplot2::xlab("PCs") +
    ggplot2::ylab("Variables") +
    ggplot2::ggtitle(paste("Correlation between principal components and metadata", "\n *p.adj < 0.05 | **p.adj < 0.01", plotname)) +
    ggplot2::theme_minimal()

  ## Calculate correlations between clinical variables
  CorrResultsClin <- CorrelationData %>%
    ungroup() %>%
    dplyr::select(-Sample) %>%
    rstatix::cor_test(vars = !dplyr::contains("PC"), vars2 = !dplyr::contains("PC")) %>%
    rstatix::adjust_pvalue(method = "BH") %>%
    ## Filter out NA correlations
    dplyr::filter(!is.na(cor))

  ## make matrix of correlations
  CorrMatrix <- CorrResultsClin %>%
    dplyr::select(var1, var2, cor) %>%
    tidyr::pivot_wider(names_from = var2, values_from = cor) %>%
    tibble::column_to_rownames("var1") %>%
    as.matrix()

  ## plot CorrMatrix using complexheatmapl
  CorrMatrixPlot <- ComplexHeatmap::Heatmap(CorrMatrix, name = "correlation", column_title = paste("Correlation between clinical variables", plotname))


  output <- list(CorResults, EffectPlot = EffectPlot, PCA = PCA, EffectCorrelationsPlot = CorrMatrixPlot, CorrResultsClin = CorrResultsClin)

  return(output)
}


## Correlation Analysis
## This function takes  as input either
## Two vectors containing the names of Proteins of Interest
##(part of) Gene names (use = "gene" | "Gene)
## and calculates the correlation pearson coefficients between them
## It writes the CorResults to the global environment and plots a Heatmap as a Results
## add roxygen comments
#' @title CorrelateFeatures
#' @description Correlation analysis for the specified dataset.
#' @param dataset The dataset to be tested
#' @param vars1 The first list of features of interest
#' @param vars2 The second list of features of interest
#' @param use The type of features to be used (Protein, Peptide, Gene)
#' @param clustdist The distance metric to be used for clustering
#' @return A list object containing the correlation matrix, the correlation results and the heatmap
#' @export
CorrelateFeatures <- function(dataset, vars1 , vars2, use = "Protein", clustdist = "euclidean"){

  ## cleaning the InputData
  dataset <- dataset %>%
    ## disregarding Proteins with less than 10 Observations
    nObsPerFeature(n = 30) %>%
    ## Imputing the Intensities of Missing values to half the minimum observed (for each Protein)
    ImputeFeatureIntensity(method = "halfmin")

  ## creating vectors of Protein names containing the PoIs
  if (use == "gene" | use =="Gene"){

    PoIs1 <- dataset %>%
      mutate(gene = str_split_i(Protein, pattern = "_", 2)) %>%
      filter(grepl(gene, pattern = vars1)) %>% distinct(Protein)
    PoIs1 <- PoIs1$Protein

    PoIs2 <- dataset %>%
      mutate(gene = str_split_i(Protein, pattern = "_", 2)) %>%
      filter(grepl(gene, pattern = vars2)) %>% distinct(Protein)
    PoIs2 <- PoIs2$Protein

  }
  if (use == "Protein" | use =="protein"){

    PoIs1 <- dataset %>% select(Sample, Protein, Intensity) %>%
      filter(Protein %in% vars1)
    PoIs1 <- PoIs1$Protein

    PoIs2 <- dataset %>% select(Sample, Protein, Intensity) %>%
      filter(Protein %in% vars2)
    PoIs2 <- PoIs2$Protein

    ## Calculating correlation coefficients
    ## creating input Dataframe
    CorResults <<- dataset %>%
      pivot_wider(names_from = "Protein", values_from = "Intensity") %>%
      ## Imputing to half the minimum
      mutate(across(everything(), ~ifelse(is.na(.), (min(., na.rm = TRUE)/2), .))) %>%
      ## selecting proteins
      select(matches("_")) %>%
      ## calculating correaltions
      cor_test(vars = PoIs1, vars2 = PoIs2, method = "pearson", use = "pairwise.complete.obs")

  }
  if (use == "Peptide" | use =="peptide"){

    PoIs1 <- dataset %>% select(Sample, Peptide, Intensity) %>%
      filter(Peptide %in% vars1)
    PoIs1 <- PoIs1$Peptide

    PoIs2 <- dataset %>% select(Sample, Peptide, Intensity) %>%
      filter(Peptide %in% vars2)
    PoIs2 <- PoIs2$Peptide

    ## Calculating correlation coefficients
    ## creating input Dataframe
    CorResults <- dataset %>%
      pivot_wider(names_from = "Peptide", values_from = "Intensity") %>%
      ## Imputing to half the minimum
      mutate(across(everything(), ~ifelse(is.na(.), (min(., na.rm = TRUE)/2), .))) %>%
      ## selecting Peptides
      select(matches("_")) %>%
      ## calculating correaltions
      cor_test(vars = PoIs1, vars2 = PoIs2, method = "pearson", use = "pairwise.complete.obs")
  }


  #### Making correlation Matrix
  CorrelationMatirix <- CorResults %>%
    ## making Var1 and Var2 the gene names
    mutate(var1 = str_split_i(var1, pattern = "_", 2)) %>%
    mutate(var2 = str_split_i(var2, pattern = "_", 2)) %>%
    select("var1", "var2", "cor") %>%
    pivot_wider(names_from = "var2", values_from = "cor") %>%
    column_to_rownames("var1") %>%
    as.matrix()

  ## making custom colour pallate
  my_palette <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

  ## plotting Output heapmat  Output heatmap
  Heatmap <- ComplexHeatmap::Heatmap(CorrelationMatirix,

                                     ## clustering specifics
                                     ## Clustering columns
                                     cluster_columns = TRUE,
                                     clustering_distance_columns = clustdist,


                                     ## clustering Rows
                                     cluster_rows = TRUE,
                                     clustering_distance_rows = clustdist,


                                     ## colours and names
                                     name = "pearson correlation",
                                     col = my_palette)

  Output <- list()
  Output$Correlationmatrix <- CorrelationMatirix
  Output$Resutls <- CorResults
  Output$Heatmap <- Heatmap

  return(Output)
}


## BiomarkerPanel
## This function takes as input a dataset and a vector of Proteins of Interest
## It calculates the AUC of a logistic regression model for each possible combination of PoIs
## It returns the best combination of PoIs as a Biomarker Panel
## add roxygen comments
#' @title BiomarkerPanel
#' @description Biomarker panel for the specified dataset.
#' @param dataset The dataset to be tested
#' @param PoIs A vector containing the Proteins of interest
#' @param n The maximum number of PoIs to be included in the Biomarker Panel
#' @return A list object containing the results of the Biomarker Panel analysis
#' @export

BiomarkerPanel <- function(dataset, PoIs, n){

  ## make every permutation of PoIs up to length n
  PoIs_permutations <- lapply(1:n, function(i) utils::combn(PoIs, i, simplify = FALSE))

  ## how many entries are in PoIPermutarions
  nIterations <- sum(sapply(PoIs_permutations, length))

  print(paste("Calculating AUC for ", nIterations, "combinations of PoIs"))

  Combination_AUCs <- data.frame()

  ## add progress bar
  pb <- utils::txtProgressBar(min = 0, max = length(PoIs_permutations), style = 3)

  if("Protein" %in% colnames(dataset)){

    ## loop through every permutation
    for (i in seq_along(PoIs_permutations)) {
      ## get the current permutation
      current_permutation <- PoIs_permutations[[i]]

      ## loop through every combination in the current permutation
      for (j in seq_along(current_permutation)) {
        ## get the current combination
        current_combination <- current_permutation[[j]]

        ## model
        model <- BiomarkR::GLM(dataset, PoIs = current_combination, crossvalidation = T)
        AUC <- model$AUC

        ## Put the AUC and combination in a data frame
        Combination_AUCs <- rbind(Combination_AUCs, data.frame(Combination = paste(current_combination, collapse = ", "), AUC = AUC)) %>%
          arrange(-AUC)

      }

      ## update progress bar
      utils::setTxtProgressBar(pb, i)

    }


  }
  if("Peptide" %in% colnames(dataset)){

    ## loop through every permutation
    for (i in seq_along(PoIs_permutations)) {
      ## get the current permutation
      current_permutation <- PoIs_permutations[[i]]

      ## loop through every combination in the current permutation
      for (j in seq_along(current_permutation)) {
        ## get the current combination
        current_combination <- current_permutation[[j]]

        ## model
        model <- BiomarkR::GLM(dataset, PoIs = current_combination)

        AUC <- model$AUC

        ## Put the AUC and combination in a data frame
        Combination_AUCs <- rbind(Combination_AUCs, data.frame(Combination = paste(current_combination, collapse = ", "), AUC = AUC))

      }

      ## update progress bar
      utils::setTxtProgressBar(pb, i)

    }

  }

  BestCombination <- Combination_AUCs %>% dplyr::filter(AUC == max(AUC)) %>% dplyr::pull(Combination)

  ## separate item in best combination
  BestCombination <- strsplit(BestCombination, ", ")[[1]]

  ## sort Combination AUCs
  Combination_AUCs <- Combination_AUCs %>% dplyr::arrange(-AUC)

  return(list(BestCombination = BestCombination, Combination_AUCs = Combination_AUCs))
}

## Machine learning

## PLSDA
## add roxygen comments
#' @title PLSDA
#' @description Partial least squares discriminant analysis for the specified dataset.
#' @param dataset The dataset to be tested
#' @param PoIs A vector containing the Proteins of interest. Example: c(""Q8TF72_SHROOM3" , "Q9ULZ3_PYCARD") or (unique(dataset$Protein))
#' @param ncomp The number of components to be calculated
#' @param plotname The name to be displayed on created plots
#' @return A list object containing the results of the PLSDA model, 2D and 3D score and loading plots and the ROC plot if the dataset has 2 classes
#' @export
PLSDA <- function(dataset, plotname = "", PoIs){

  ## Partial least squares discriminant analysis using caret
  ## make machine learning data
  if("Protein" %in% colnames(dataset)){
    MLData <- dataset %>%
      dplyr::filter(Protein %in% PoIs) %>%
      tidyr::pivot_wider(names_from = Protein, values_from = Intensity) %>%
      dplyr::select(Status, dplyr::contains("_")) %>%
      ## remove spaces from Status
      dplyr::mutate(Status = gsub(" ", ".", Status))
  }

  if("Peptide" %in% colnames(dataset)){
    MLData <- dataset %>%
      dplyr::filter(Peptide %in% PoIs) %>%
      tidyr::pivot_wider(names_from = Peptide, values_from = Intensity) %>%
      dplyr::select(Status, dplyr::contains("_")) %>%
      ## remove spaces from Status
      dplyr::mutate(Status = gsub(" ", ".", Status))
  }

  ## Error message when missing values in MLData > 0
  if(sum(is.na(MLData)) > 0){
    base::message("Missing values in Data, please impute or remove them")
  }

  ## Create tune grid for plsda
  tune_grid <- base::expand.grid(ncomp = c(3, 4, 5, 6, 7, 8, 9, 10))

  ## make train control including recursive feature elimination
  if(length(base::unique(MLData$Status)) == 2){
    ## make train control for two class summary (10 x cross validation)
    train_control <- caret::trainControl(method = "cv", number = 10, classProbs = TRUE,
                                         summaryFunction = caret::twoClassSummary,
                                         savePredictions = TRUE,
                                         preProc = c("center", "scale"))
  } else {
    ## make train control for multi class summary (10 x cross validation)
    train_control <- caret::trainControl(method = "cv", number = 10, classProbs = TRUE,
                                         summaryFunction = caret::multiClassSummary,
                                         savePredictions = TRUE,
                                         preProc = c("center", "scale"))
  }

  ## train the model
  plsda_model <- caret::train(Status ~ ., data = MLData, method = "pls", trControl = train_control, tuneGrid = tune_grid)

  ## plot the tune results
  optimization_plot <- graphics::plot(plsda_model)

  ## extract confusion matrix
  Confusion_Matrix <- caret::confusionMatrix(plsda_model)

  ## extract ROC curve from the saved predictions if only two classes
  if(length(base::unique(MLData$Status)) == 2){
    # Extract predictions and actual outcomes
    predictions <- plsda_model$pred[4] %>% pull()   # Predicted probabilities
    actuals <- plsda_model$pred$obs        # Actual outcomes

    roc_curve <- pROC::roc(response = actuals, predictor = predictions, levels = base::levels(actuals))

    ## plot ROC curve using ggplot
    AUC <- base::round(pROC::auc(roc_curve), 2)

    ## get sensitivity and specificity from roc_curve
    ## plot ROC curve using ggplot

    roc_plot <- pROC::ggroc(roc_curve, legacy.axes = TRUE) +
      ggplot2::ggtitle("ROC Curve") +
      ggplot2::theme_minimal() +
      ## add 1:1 line (red, dashed)
      ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
      ggplot2::ggtitle("ROC Curve", subtitle = base::paste("AUC =", AUC))

  }

  ## extract scores
  scores1 <- plsda_model$finalModel$scores[,1]
  scores2 <- plsda_model$finalModel$scores[,2]
  scores3 <- plsda_model$finalModel$scores[,3]

  scores <- base::cbind(scores1, scores2, scores3) %>% data.frame()

  ##  combine scores with the original data
  scores <- base::cbind(MLData$Status, scores) %>% data.frame()

  ## rename columns
  colnames(scores) <- c("Status", "LV1", "LV2", "LV3")

  ## plot the scores of LV1 on LV2, LV1 on LV3 and LV2 on LV3 using ggplot
  scores_plot <- ggplot2::ggplot(scores, ggplot2::aes(x = LV1, y = LV2, color = Status)) +
    ggplot2::geom_point() +
    ggplot2::ggtitle("PLS-DA Scores Plot") +
    ggplot2::theme_minimal()


  ## make 3D score plot
  scores_plot_3D <- plotly::plot_ly(scores, x = ~LV1, y = ~LV2, z = ~LV3, color = ~Status, type = "scatter3d", mode = "markers") %>%
    plotly::layout(title = "PLS-DA Scores Plot")

  ## make a list of outputs
  output <- list(optimization_plot = optimization_plot,
                 Confusion_Matrix = Confusion_Matrix,
                 scores_plot = scores_plot,
                 scores_plot_3D = scores_plot_3D,
                 Model = plsda_model)

  ## Include ROC plot for 2 classes
  if(length(base::unique(MLData$Status)) == 2){
    output$roc_plot <- roc_plot
  }

  return(output)
}



## Gradient Boosting machine
## Gradient boosting machine with generic parameter grid settings
## add roxygen comments
#' @title GBM
#' @description Gradient boosting machine for the specified dataset.
#' @param dataset The dataset to be tested
#' @param PoIs A vector containing the Proteins of interest. Example: c(""Q8TF72_SHROOM3" , "Q9ULZ3_PYCARD") or (unique(dataset$Protein))
#' @param plotname The name to be displayed on created plots
#' @return A list object containing the results of the GBM model, the confusion matrix and the ROC plot and the AUC if the dataset has 2 classes
#' @export
GBM <- function(dataset, PoIs, plotname = "") {

  if ("Protein" %in% colnames(dataset)) {

    MLData <- dataset %>%
      normalizeIntensityOnSample(plot = FALSE) %>%
      dplyr::filter(Protein %in% PoIs) %>%

      # Convert to wide format
      tidyr::pivot_wider(names_from = Protein, values_from = Intensity) %>%
      dplyr::select(Status, Sample, dplyr::contains("_"))

  }

  if ("Peptide" %in% colnames(dataset)) {

    MLData <- dataset %>%
      normalizeIntensityOnSample(plot = FALSE) %>%
      dplyr::filter(Peptide %in% PoIs) %>%

      # Convert to wide format
      tidyr::pivot_wider(names_from = Peptide, values_from = Intensity) %>%
      dplyr::select(Status, Sample, dplyr::contains("_"))

  }


  # Create tune grid
  tuneGrid <- base::expand.grid(n.trees = c(100, 200, 300),
                                interaction.depth = c(1, 3, 5),
                                shrinkage = c(0.1, 0.3, 0.5),
                                n.minobsinnode = 10)

  # Gradient Boosting Machine for 2 classes
  if (length(base::unique(MLData$Status)) == 2) {

    train_control <- caret::trainControl(method = "cv",
                                         number = 10,
                                         classProbs = TRUE,
                                         summaryFunction = caret::twoClassSummary,
                                         savePredictions = TRUE)

    # Train the GBM model
    GBMModel <- caret::train(Status ~ .,
                             data = dplyr::select(MLData, -Sample),
                             method = "gbm",
                             preProc = c("center", "scale"),
                             trControl = train_control,
                             tuneGrid = tuneGrid,
                             metric = "ROC")

      # Generate ROC curve
    # Extract predictions and actual outcomes
    predictions <- GBMModel$pred[4] %>% pull()   # Predicted probabilities
    actuals <- GBMModel$pred$obs        # Actual outcomes

    roc_curve <- pROC::roc(response = actuals, predictor = predictions, levels = base::levels(actuals))

    ## plot ROC curve using ggplot
    AUC <- base::round(pROC::auc(roc_curve), 2)

    ## get sensitivity and specificity from roc_curve
    ## plot ROC curve using ggplot

    roc_plot <- pROC::ggroc(roc_curve, legacy.axes = TRUE) +
      ggplot2::ggtitle("ROC Curve") +
      ggplot2::theme_minimal() +
      ## add 1:1 line (red, dashed)
      ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
      ggplot2::ggtitle("ROC Curve", subtitle = base::paste("AUC =", AUC))

  }

  # Gradient Boosting Machine for 3 or more classes
  if (length(base::unique(MLData$Status)) > 2) {

    trainData$Status <- base::as.factor(trainData$Status)

    train_control <- caret::trainControl(method = "cv",
                                         number = 10,
                                         classProbs = TRUE,
                                         summaryFunction = caret::multiClassSummary)

    # Train the GBM model
    GBMModel <- caret::train(Status ~ .,
                             data = dplyr::select(trainData, -Sample),
                             method = "gbm",
                             preProc = c("center", "scale"),
                             trControl = train_control,
                             tuneGrid = tuneGrid,
                             metric = "Accuracy")

  }

  ## extract confusion matrix
  ConfusionMatrix <- caret::confusionMatrix(GBMModel)

  # Creating output list
  Output <- base::list(GBMModel = GBMModel,
                       ConfusionMatrix = ConfusionMatrix)

  # Add ROC and AUC to output list if there are 2 classes
  if (length(base::unique(MLData$Status)) == 2) {
    Output$ROC_Plot <- roc_plot
    Output$AUC <- AUC
  }

  return(Output)
}


## Support Vector Machine
## Support Vector Machine for 2 classes
## add roxygen comments
#' @title SVM
#' @description Support vector machine for the specified dataset.
#' @param dataset The dataset to be tested
#' @param PoIs A vector containing the Proteins of interest. Example: c(""Q8TF72_SHROOM3" , "Q9ULZ3_PYCARD") or (unique(dataset$Protein))
#' @param plotname The name to be displayed on created plots
#' @return A list object containing the results of the SVM model, the confusion matrix and the ROC plot and the AUC if the dataset has 2 classes
#' @export
SVM <- function(dataset, PoIs) {

  ## Support vector machine using caret
  ## make machine learning data
  if ("Protein" %in% colnames(dataset)) {
    MLData <- dataset %>%
      dplyr::filter(Protein %in% PoIs) %>%
      tidyr::pivot_wider(names_from = Protein, values_from = Intensity) %>%
      dplyr::select(Status, dplyr::contains("_")) %>%
      ## remove spaces from Status
      dplyr::mutate(Status = gsub(" ", ".", Status))
  }

  if ("Peptide" %in% colnames(dataset)) {
    MLData <- dataset %>%
      dplyr::filter(Peptide %in% PoIs) %>%
      tidyr::pivot_wider(names_from = Peptide, values_from = Intensity) %>%
      dplyr::select(Status, dplyr::contains("_")) %>%
      ## remove spaces from Status
      dplyr::mutate(Status = gsub(" ", ".", Status))
  }

  ## Error message when missing values in MLData > 0
  if (base::sum(base::is.na(MLData)) > 0) {
    base::message("Missing values in Data, please impute or remove them")
  }

  ## make tune grid for linear SVM
  tune_grid <- base::expand.grid(C = c(0.01, 0.1, 1, 10, 100, 1000))

  ## make train control including recursive feature elimination
  if (length(base::unique(MLData$Status)) == 2) {
    ## make train control for two class summary (10 x cross-validation)
    train_control <- caret::trainControl(method = "cv", number = 10, classProbs = TRUE,
                                         summaryFunction = caret::twoClassSummary,
                                         savePredictions = TRUE,
                                         selectionFunction = "best",
                                         preProc = c("center", "scale"))
  } else {
    ## make train control for multi-class summary (10 x cross-validation)
    train_control <- caret::trainControl(method = "cv", number = 10, classProbs = TRUE,
                                         summaryFunction = caret::multiClassSummary,
                                         savePredictions = TRUE,
                                         selectionFunction = "best",
                                         preProc = c("center", "scale"))
  }

  ## train the model
  svm_model <- caret::train(Status ~ ., data = MLData, method = "svmLinear", trControl = train_control, tuneGrid = tune_grid)

  ## plot the tune results
  optimization_plot <- graphics::plot(svm_model)

  ## extract confusion matrix
  Confusion_Matrix <- caret::confusionMatrix(svm_model)

  ## extract ROC curve from the saved predictions if only two classes
  if (length(base::unique(MLData$Status)) == 2) {
    # Extract predictions and actual outcomes
    predictions <- svm_model$pred[3] %>% pull()   # Predicted probabilities
    actuals <- svm_model$pred$obs        # Actual outcomes

    # Ensure predictions are numeric
    predictions <- as.numeric(predictions)

    # Ensure actuals are a factor
    actuals <- factor(actuals, levels = levels(actuals))

    # Generate ROC curve
    roc_curve <- pROC::roc(response = actuals, predictor = predictions, levels = base::levels(actuals))

    ## plot ROC curve using ggplot
    AUC <- base::round(pROC::auc(roc_curve), 2)

    ## get sensitivity and specificity from roc_curve
    ## plot ROC curve using ggplot

    roc_plot <- pROC::ggroc(roc_curve, legacy.axes = TRUE) +
      ggplot2::ggtitle("ROC Curve") +
      ggplot2::theme_minimal() +
      ## add 1:1 line (red, dashed)
      ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
      ggplot2::ggtitle("ROC Curve", subtitle = base::paste("AUC =", AUC))

  }



  ScorePlotData <- dataset %>%
    dplyr::filter(Protein %in% PoIs)

  PCA <- BiomarkR::PCA(ScorePlotData)

  ScorePlot2D <- PCA$ScorePlot2D
  ScorePlot3D <- PCA$ScorePlot3D

  ## Create a list of all the outputs
  output <- base::list(Optimization_Plot = optimization_plot,
                       Confusion_Matrix = Confusion_Matrix,
                       ScorePlot2D = ScorePlot2D,
                       ScorePlot3D = ScorePlot3D)

  ## add ROC plot if available
  if (length(base::unique(MLData$Status)) == 2) {
    output$ROC_Plot <- roc_plot
  }

  return(output)

}



## GLM
## add roxygen comments
#' @title GLM
#' @description Generalized linear model for the specified dataset.
#' @param dataset The dataset to be tested
#' @param PoIs A vector containing the Proteins of interest. Example: c(""Q8TF72_SHROOM3" , "Q9ULZ3_PYCARD") or (unique(dataset$Protein))
#' @param crossvalidation A boolean value indicating whether to use cross-validation> I recommend false for desriptive analysese and true for predictife tasks
#' @return A list object containing the results of the GLM model, the confusion matrix and the ROC plot and the AUC if the dataset has 2 classes
#' @param plotname The name to be displayed on created plots
#' @export
GLM <- function(dataset, PoIs, crossvalidation = F, plotname = "") {

  ## Error message if more than 2 classes, stop execution
  if (length(base::unique(dataset$Status)) > 2) {
    base::stop("Generalized Linear Model only works for two classes")
  }

  ## Generalized Linear Model using caret
  ## make machine learning data
  if ("Protein" %in% base::colnames(dataset)) {
    MLData <- dataset %>%
      dplyr::filter(Protein %in% PoIs) %>%
      tidyr::pivot_wider(names_from = Protein, values_from = Intensity) %>%
      dplyr::select(Status, dplyr::contains("_")) %>%
      ## remove spaces from Status
      dplyr::mutate(Status = base::gsub(" ", ".", Status))
  }

  if ("Peptide" %in% base::colnames(dataset)) {
    MLData <- dataset %>%
      dplyr::filter(Peptide %in% PoIs) %>%
      tidyr::pivot_wider(names_from = Peptide, values_from = Intensity) %>%
      dplyr::select(Status, dplyr::contains("_")) %>%
      ## remove spaces from Status
      dplyr::mutate(Status = base::gsub(" ", ".", Status))
  }

  ## Error message when missing values in MLData > 0
  if (base::sum(base::is.na(MLData)) > 0) {
    base::message("Missing values in Data, please impute or remove them")
  }

  ## train control
  train_control <- caret::trainControl(method = ifelse(crossvalidation == T, "cv", "none"),
                                       number = ifelse(crossvalidation == T, 10, NA) ,
                                       classProbs = TRUE, summaryFunction = caret::twoClassSummary,
                                       savePredictions = TRUE)

  ## train the model
  glm_model <- caret::train(Status ~ ., data = MLData, method = "glm", trControl = train_control)

  ## Compute confusion matrix
  if(crossvalidation == T){
    Confusion_Matrix <- caret::confusionMatrix(glm_model)
  } else {
    predictions <- stats::predict(glm_model, MLData) # Predicted probabilities for the positive class
    actuals <- factor(MLData$Status)
    Confusion_Matrix <-caret::confusionMatrix(predictions, actuals)
  }



  ## Use cross-validated predictions for ROC curve
  if (crossvalidation == TRUE) {
    # Extract the cross-validated predictions
    crossval_preds <- glm_model$pred %>% data.frame()

    # Make sure levels of predictions match the actual levels
    crossval_preds$obs <- factor(crossval_preds$obs, levels = unique(MLData$Status))

    predictor <- crossval_preds[[3]]

    # Generate ROC curve using cross-validated predictions
    suppressMessages(suppressWarnings(
      roc_curve <- pROC::roc(response = crossval_preds$obs,
                             predictor = predictor, # Predicted probabilities for positive class
                             ci = TRUE, boot.n = 100, conf.level = 0.95)
    ))

    # Calculate confidence intervals for specificities
    ci_band <- pROC::ci.se(roc_curve, specificities = seq(0, 1, by = 0.01), boot.n = 100) %>%
      data.frame() %>% rownames_to_column(var = "specificity")

    ci_band$specificity <- as.numeric(ci_band$specificity)

    colnames(ci_band)[2:4] <- c("lower", "mean", "upper")

    AUC_CI <- round(c(roc_curve$ci)[3],2)

  } else {
    # Use the non-CV predictions
    probabilities <- stats::predict(glm_model, MLData, type = "prob")[, 2]  # Predicted probabilities for the positive class
    actuals <- factor(MLData$Status)

    suppressMessages(suppressWarnings(
      roc_curve <- pROC::roc(response = actuals,
                             predictor = probabilities)
    ))
  }


  ## plot ROC curve using ggplot
  AUC <- base::round(pROC::auc(roc_curve), 2)


  roc_coords <- coords(roc_curve, x = "all", ret = c("1-specificity", "sensitivity"))
  roc_data <- data.frame(FPR = roc_coords$`1-specificity`, TPR = roc_coords$sensitivity)

  # Customizing aesthetics within ggroc
  # Create a manual ggplot
  roc_plot <- ggplot(roc_data, aes(x = FPR, y = TPR)) +
    geom_path(color = "black", size = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") + # Diagonal line
    labs(x = "False Positive Rate", y = "True Positive Rate", title = paste("AUC:", round(auc(roc_curve), 2))) +
    theme_minimal()

  ## add shaded confidence interval if cross validation == TRUE
  if(crossvalidation == T){
    roc_plot <- roc_plot +
      ggplot2::geom_ribbon(data = ci_band, aes(x = 1 - specificity, ymin = lower, ymax = upper), alpha = 0.2, linewidth = 1, inherit.aes = F)+
      ggplot2::ggtitle(plotname, paste("AUC:", AUC, "±" , round(AUC_CI - AUC, 2), "(95 % CI) ", "| 10 x cross validated"))
  }




  ## Create a list of all the output
  output <- base::list(Confusion_Matrix = Confusion_Matrix,
                       ROC_Plot = roc_plot,
                       AUC = c(AUC),
                       model = glm_model)

  return(output)
}

## Random forest
## add roxygen comments
#' @title Random Forest
#' @description Random forest for the specified dataset.
#' @param dataset The dataset to be tested
#' @param PoIs A vector containing the Proteins of interest. Example: c(""Q8TF72_SHROOM3" , "Q9ULZ3_PYCARD") or (unique(dataset$Protein))
#' @param mtry The number of variables randomly sampled as candidates at each split
#' @param plotname The name to be displayed on created plots
#' @return A list object containing the results of the RF model, the confusion matrix and the ROC plot and the AUC if the dataset has 2 classes
#' @export
RandomForest <- function(dataset, PoIs) {

  ## Random Forest using caret
  ## make machine learning data
  if ("Protein" %in% base::colnames(dataset)) {
    MLData <- dataset %>%
      dplyr::filter(Protein %in% PoIs) %>%
      tidyr::pivot_wider(names_from = Protein, values_from = Intensity) %>%
      dplyr::select(Status, dplyr::contains("_")) %>%
      ## remove spaces from Status
      dplyr::mutate(Status = base::gsub(" ", ".", Status))
  }

  if ("Peptide" %in% base::colnames(dataset)) {
    MLData <- dataset %>%
      dplyr::filter(Peptide %in% PoIs) %>%
      tidyr::pivot_wider(names_from = Peptide, values_from = Intensity) %>%
      dplyr::select(Status, dplyr::contains("_")) %>%
      ## remove spaces from Status
      dplyr::mutate(Status = base::gsub(" ", ".", Status))
  }

  ## Error message when missing values in MLData > 0
  if (base::sum(base::is.na(MLData)) > 0) {
    base::message("Missing values in Data, please impute or remove them")
  }

  ## make tune grid for random forest
  tune_grid <- base::expand.grid(mtry = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))

  ## make train control
  if (length(base::unique(MLData$Status)) == 2) {
    ## make train control for two class summary (10 x cross-validation)
    train_control <- caret::trainControl(method = "cv", number = 10, classProbs = TRUE, summaryFunction = caret::twoClassSummary, savePredictions = TRUE)
  } else {
    ## make train control for multi class summary (10 x cross-validation)
    train_control <- caret::trainControl(method = "cv", number = 10, classProbs = TRUE, summaryFunction = caret::multiClassSummary, savePredictions = TRUE)
  }

  ## train the model
  rf_model <- caret::train(Status ~ ., data = MLData, method = "rf", trControl = train_control, tuneGrid = tune_grid)

  ## plot the tune results
  optimization_plot <- graphics::plot(rf_model)

  ## extract confusion matrix
  Confusion_Matrix <- caret::confusionMatrix(rf_model)

  if (length(base::unique(MLData$Status)) == 2) {
    # Extract predictions and actual outcomes
    predictions <- rf_model$pred[3] %>% pull()   # Predicted probabilities
    actuals <- rf_model$pred$obs        # Actual outcomes
    # Ensure predictions are numeric
    predictions <- as.numeric(predictions)

    # Ensure actuals are a factor
    actuals <- factor(actuals, levels = levels(actuals))

    # Generate ROC curve
    roc_curve <- pROC::roc(response = actuals, predictor = predictions, levels = base::levels(actuals))

    ## plot ROC curve using ggplot
    AUC <- base::round(pROC::auc(roc_curve), 2)

    ## get sensitivity and specificity from roc_curve
    ## plot ROC curve using ggplot

    roc_plot <- pROC::ggroc(roc_curve, legacy.axes = TRUE) +
      ggplot2::ggtitle("ROC Curve") +
      ggplot2::theme_minimal() +
      ## add 1:1 line (red, dashed)
      ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
      ggplot2::ggtitle("ROC Curve", subtitle = base::paste("AUC =", AUC))

  }

  ## extract feature importance
  FeatureImportance <- caret::varImp(rf_model, scale = FALSE)
  FeatureImportance <- FeatureImportance$importance %>%
    base::data.frame() %>%
    tibble::rownames_to_column(var = "Protein") %>%
    dplyr::arrange(dplyr::desc(Overall)) %>%
    utils::head(20)

  ## plot feature importance using ggplot
  FeatureImportancePlot <- FeatureImportance %>%
    ggplot2::ggplot(ggplot2::aes(x = forcats::fct_reorder(Protein, Overall), y = Overall)) +
    ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal() +
    ggplot2::ggtitle("Feature Importance") +
    ggplot2::xlab("Protein") +
    ggplot2::ylab("Importance")

  ## visualize the top 20 features using PCA
  ## Top 20 features
  Features <- FeatureImportance$Protein %>%
    utils::head(20)

  PCA <- BiomarkR::PCA(dataset %>% dplyr::filter(Protein %in% Features), plotname = "Random Forest top 20 features")

  ScorePlot2D <- PCA$ScorePlot2D
  ScorePlot3D <- PCA$ScorePlot3D

  ## Create a list of all the output
  output <- base::list(Optimization_Plot = optimization_plot,
                       Confusion_Matrix = Confusion_Matrix,
                       Feature_Importance = FeatureImportancePlot,
                       ScorePlot2D = ScorePlot2D,
                       ScorePlot3D = ScorePlot3D,
                       model = rf_model)
  ## add ROC plot if available
  if (length(base::unique(MLData$Status)) == 2) {
    output$ROC_Plot <- roc_plot
    output$AUC <- AUC
  }

  return(output)
}



## K-Nearest Neighbors
## add roxygen comments
#' @title KNN
#' @description K-Nearest Neighbors calssification for the specified dataset.
#' @param dataset The dataset to be tested
#' @param PoIs A vector containing the Proteins of interest. Example: c(""Q8TF72_SHROOM3" , "Q9ULZ3_PYCARD") or (unique(dataset$Protein))
#' @return A list object containing the results of the KNN model, the confusion matrix and the ROC plot and the AUC if the dataset has 2 classes
#' @export
KNN <- function(dataset, PoIs) {

  ## K-Nearest Neighbors using caret
  ## make machine learning data
  if ("Protein" %in% base::colnames(dataset)) {
    MLData <- dataset %>%
      dplyr::filter(Protein %in% PoIs) %>%
      tidyr::pivot_wider(names_from = Protein, values_from = Intensity) %>%
      dplyr::select(Status, dplyr::contains("_")) %>%
      ## remove spaces from Status
      dplyr::mutate(Status = base::gsub(" ", ".", Status))
  }

  if ("Peptide" %in% base::colnames(dataset)) {
    MLData <- dataset %>%
      dplyr::filter(Peptide %in% PoIs) %>%
      tidyr::pivot_wider(names_from = Peptide, values_from = Intensity) %>%
      dplyr::select(Status, dplyr::contains("_")) %>%
      ## remove spaces from Status
      dplyr::mutate(Status = base::gsub(" ", ".", Status))
  }

  ## Error message when missing values in MLData > 0
  if (base::sum(base::is.na(MLData)) > 0) {
    base::message("Missing values in Data, please impute or remove them")
  }

  ## make tune grid for KNN
  tune_grid <- base::expand.grid(k = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))

  ## make train control
  if (length(base::unique(MLData$Status)) == 2) {
    ## make train control for two class summary (10 x cross-validation)
    train_control <- caret::trainControl(method = "cv", number = 10, classProbs = TRUE, summaryFunction = caret::twoClassSummary, savePredictions = TRUE)
  } else {
    ## make train control for multi class summary (10 x cross-validation)
    train_control <- caret::trainControl(method = "cv", number = 10, classProbs = TRUE, summaryFunction = caret::multiClassSummary, savePredictions = TRUE)
  }

  ## train the model
  knn_model <- caret::train(Status ~ ., data = MLData, method = "knn", trControl = train_control, tuneGrid = tune_grid)

  ## extract confusion matrix
  Confusion_Matrix <- caret::confusionMatrix(knn_model)

  if (length(base::unique(MLData$Status)) == 2) {
    # Extract predictions and actual outcomes
    predictions <- knn_model$pred[3] %>% pull()   # Predicted probabilities
    actuals <- knn_model$pred$obs        # Actual outcomes

    # Ensure predictions are numeric
    predictions <- as.numeric(predictions)

    # Ensure actuals are a factor
    actuals <- factor(actuals, levels = levels(actuals))

    # Generate ROC curve
    roc_curve <- pROC::roc(response = actuals, predictor = predictions, levels = base::levels(actuals))

    ## plot ROC curve using ggplot
    AUC <- base::round(pROC::auc(roc_curve), 2)

    ## get sensitivity and specificity from roc_curve
    ## plot ROC curve using ggplot

    roc_plot <- pROC::ggroc(roc_curve, legacy.axes = TRUE) +
      ggplot2::ggtitle("ROC Curve") +
      ggplot2::theme_minimal() +
      ## add 1:1 line (red, dashed)
      ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
      ggplot2::ggtitle("ROC Curve", subtitle = base::paste("AUC =", AUC))

  }

  ## Create a list of all the output
  output <- base::list(Confusion_Matrix = Confusion_Matrix,
                       ## ROC plot if available
                       if (length(base::unique(MLData$Status)) == 2) "ROC Plot" = roc_plot,
                       ## add AUC if available
                       if (length(base::unique(MLData$Status)) == 2) AUC = AUC,
                       model = knn_model)

  return(output)
}



## Naive Bayes
## add roxygen comments
#' @title Naive Bayes
#' @description Naive Bayes calssification for the specified dataset.
#' @param dataset The dataset to be tested
#' @param PoIs A vector containing the Proteins of interest. Example: c(""Q8TF72_SHROOM3" , "Q9ULZ3_PYCARD") or (unique(dataset$Protein))
#' @return A list object containing the results of the Naive Bayes model, the confusion matrix and the ROC plot and the AUC if the dataset has 2 classes
#' @export
NaiveBayes <- function(dataset, PoIs) {

  ## Naive Bayes using caret
  ## make machine learning data
  if ("Protein" %in% base::colnames(dataset)) {
    MLData <- dataset %>%
      dplyr::filter(Protein %in% PoIs) %>%
      tidyr::pivot_wider(names_from = Protein, values_from = Intensity) %>%
      dplyr::select(Status, dplyr::contains("_")) %>%
      ## remove spaces from Status
      dplyr::mutate(Status = base::gsub(" ", ".", Status))
  }

  if ("Peptide" %in% base::colnames(dataset)) {
    MLData <- dataset %>%
      dplyr::filter(Peptide %in% PoIs) %>%
      tidyr::pivot_wider(names_from = Peptide, values_from = Intensity) %>%
      dplyr::select(Status, dplyr::contains("_")) %>%
      ## remove spaces from Status
      dplyr::mutate(Status = base::gsub(" ", ".", Status))
  }

  ## Error message when missing values in MLData > 0
  if (base::sum(base::is.na(MLData)) > 0) {
    base::message("Missing values in Data, please impute or remove them")
  }

  ## make train control
  if (length(base::unique(MLData$Status)) == 2) {
    ## make train control for two class summary (10 x cross-validation)
    train_control <- caret::trainControl(method = "cv", number = 10, classProbs = TRUE, summaryFunction = caret::twoClassSummary, savePredictions = TRUE)
  } else {
    ## make train control for multi class summary (10 x cross-validation)
    train_control <- caret::trainControl(method = "cv", number = 10, classProbs = TRUE, summaryFunction = caret::multiClassSummary, savePredictions = TRUE)
  }

  ## train the model
  nb_model <- caret::train(Status ~ ., data = MLData, method = "nb", trControl = train_control)

  ## extract confusion matrix
  Confusion_Matrix <- caret::confusionMatrix(nb_model)

  if (length(base::unique(MLData$Status)) == 2) {
    # Extract predictions and actual outcomes
    predictions <- nb_model$pred[3] %>% pull()   # Predicted probabilities
    actuals <- nb_model$pred$obs        # Actual outcomes

    # Ensure predictions are numeric
    predictions <- as.numeric(predictions)

    # Ensure actuals are a factor
    actuals <- factor(actuals, levels = levels(actuals))

    # Generate ROC curve
    roc_curve <- pROC::roc(response = actuals, predictor = predictions, levels = base::levels(actuals))

    ## plot ROC curve using ggplot
    AUC <- base::round(pROC::auc(roc_curve), 2)

    ## get sensitivity and specificity from roc_curve
    ## plot ROC curve using ggplot

    roc_plot <- pROC::ggroc(roc_curve, legacy.axes = TRUE) +
      ggplot2::ggtitle("ROC Curve") +
      ggplot2::theme_minimal() +
      ## add 1:1 line (red, dashed)
      ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
      ggplot2::ggtitle("ROC Curve", subtitle = base::paste("AUC =", AUC))

  }

  ## Create a list of all the output
  output <- base::list(Confusion_Matrix = Confusion_Matrix,
                       ## ROC plot if available
                       if (length(base::unique(MLData$Status)) == 2) "ROC Plot" = roc_plot,
                       ## add AUC if available
                       if (length(base::unique(MLData$Status)) == 2) AUC = AUC,
                       model = nb_model)

  return(output)
}


## k fold feature selection
## add roxygen comments
#' @title kffs
#' @description k fold feature selection for the specified dataset.
#' @param dataset The dataset to be tested
#' @param method The method to be used for feature selection (SVM, RF, LDA, GLM)
#' @param kfolds The number of folds to be used for feature selection
#' @return A list object containing the results of the feature selection
#' @export
kffs <- function(dataset, method = "SVM", kfolds = 10){

  # Initialize progress bar
  pb <- utils::txtProgressBar(min = 0, max = kfolds, style = 3)

  if("Protein" %in% colnames(dataset)){

    dataset <- dataset %>%
      dplyr::select(Sample, Status, Protein, Intensity) %>%
      NaCutoff(50) %>%
      tidyr::pivot_wider(names_from = Protein, values_from = Intensity) %>%
      tidyr::pivot_longer(cols = contains("_"), names_to = "Protein", values_to = "Intensity") %>%
      ImputeFeatureIntensity() %>%
      tidyr::pivot_wider(names_from = Protein, values_from = Intensity)
  }

  if("Peptide" %in% colnames(dataset)){

    dataset <- dataset %>%
      dplyr::select(Sample, Status, Peptide, Intensity) %>%
      NaCutoff(50) %>%
      tidyr::pivot_wider(names_from = Peptide, values_from = Intensity) %>%
      tidyr::pivot_longer(cols = contains("_"), names_to = "Peptide", values_to = "Intensity") %>%
      ImputeFeatureIntensity() %>%
      tidyr::pivot_wider(names_from = Peptide, values_from = Intensity)
  }

  ## If k folds > Number of samples, set to number of samples
  if(kfolds > nrow(dataset)){
    kfolds <- nrow(dataset)
    message("Number of folds is greater than the number of samples, setting kfolds to the number of samples")
  }


  train <- dataset

  ## Train the models for feature selection

  ## split the data into 10 folds
  folds <- caret::createFolds(train$Status, k = kfolds, list = TRUE, returnTrain = TRUE)

  Features <- list()

  if (method == "SVM"){

    ## loop through the folds using a linear SVM model
    for (i in 1:kfolds){

      ## select the ith fold
      trainData <- train[folds[[i]],]

      trainx <- trainData %>% dplyr::select(-Status) %>% dplyr::select(-Sample)

      trainy <- trainData$Status

      ## create tune grid
      tuneGrid <- expand.grid(C = c(0.1, 1, 10, 100, 1000))

      ## create train control for 10 fold cross-validation
      trainControl <- caret::trainControl(method = "cv", number = 10)

      ## train the model
      model <- caret::train(x = trainx, y = trainy, method = "svmLinear", trControl = trainControl, tuneGrid = tuneGrid)

      ## get features used in the final model
      Features[[i]] <- model$finalModel@SVindex

      ## match the feature index to feature names
      Features[[i]] <- colnames(trainx)[Features[[i]]]

      # Update progress bar
      utils::setTxtProgressBar(pb, i)
    }
  }

  if (method == "RF"){

    ## loop through the folds using a random forest model
    for (i in 1:kfolds){

      ## select the ith fold
      trainData <- train[folds[[i]],]

      trainx <- trainData %>% dplyr::select(-Status) %>% dplyr::select(-Sample)

      trainy <- trainData$Status

      ## create tune grid
      tuneGrid <- expand.grid(mtry = c(1, 2, 3, 4, 5))

      ## create train control for 10 fold cross-validation
      trainControl <- caret::trainControl(method = "cv", number = 10)

      ## train the model
      model <- caret::train(x = trainx, y = trainy, method = "rf", trControl = trainControl, tuneGrid = tuneGrid)

      ## get features used in the final model
      Features[[i]] <- model$finalModel$importance %>% as.data.frame() %>% dplyr::arrange(desc(MeanDecreaseGini)) %>% head(20)

      ## match the feature index to feature names
      Features[[i]] <- rownames(Features[[i]])

      # Update progress bar
      utils::setTxtProgressBar(pb, i)
    }
  }

  if(method == "LDA"){

    ## loop through the folds using a linear discriminant analysis model
    for (i in 1:kfolds){

      ## select the ith fold
      trainData <- train[folds[[i]],]

      trainx <- trainData %>% dplyr::select(-Status) %>% dplyr::select(-Sample)

      trainy <- trainData$Status

      ## create train control for 10 fold cross-validation
      trainControl <- caret::trainControl(method = "cv", number = 10)

      ## train the model
      model <- caret::train(x = trainx, y = trainy, method = "lda", trControl = trainControl)

      ## get features used in the final model
      Features[[i]] <- model$finalModel$scaling %>% as.data.frame() %>% dplyr::arrange(desc(abs(LD1))) %>% head(20)

      ## match the feature index to feature names
      Features[[i]] <- rownames(Features[[i]])

      # Update progress bar
      utils::setTxtProgressBar(pb, i)
    }
  }

  if(method == "GLM"){

    ## loop through the folds using a generalized linear model
    for (i in 1:kfolds){

      ## select the ith fold
      trainData <- train[folds[[i]],]

      trainx <- trainData %>% dplyr::select(-Status) %>% dplyr::select(-Sample)

      trainy <- trainData$Status

      ## create train control for 10 fold cross-validation
      trainControl <- caret::trainControl(method = "cv", number = 10)

      ## train the model
      model <- caret::train(x = trainx, y = trainy, method = "glm", trControl = trainControl)

      ## get features used in the final model
      Features[[i]] <- model$finalModel$coefficients %>% as.data.frame() %>%
        tibble::rownames_to_column(var = "Feature") %>%
        dplyr::filter(Feature != "(Intercept)") %>%
        dplyr::filter(!is.na(.)) %>%
        dplyr::arrange(desc(abs(.)))

      ## match the feature index to feature names
      Features[[i]] <- Features[[i]]$Feature

      # Update progress bar
      utils::setTxtProgressBar(pb, i)
    }
  }

  ## Close progress bar
  close(pb)

  ## extract feature frequency
  Features <- unlist(Features)

  ## make sure Feature is a character
  Features <- as.character(Features)

  Features <- Features %>% table() %>% as.data.frame() %>% dplyr::arrange(desc(Freq)) %>% dplyr::filter(Freq > 1)

  ## rename the columns
  colnames(Features) <- c("Feature", "Frequency")

  output <- list(Features = Features)

  return(output)
}


## Network Analysis
## STRING performs pathway enrichment analysis of a list of Proteins
## returns a List object with plots and a restults tables
## add roxygen comments
#' @title STRING
#' @description STRING analysis for the specified dataset.
#' @param PoIs A vector containing the Proteins of interest. Example: c("Q8TF72_SHROOM3" , "Q9ULZ3_PYCARD")
#' @param BackgroundSet A vector containing the background Proteins. Example: c(""Q8TF72_SHROOM3" , "Q9ULZ3_PYCARD") or (unique(dataset$Protein))
#' @param plotname The name to be displayed on created plots
#' @return A list object containing the results of the STRING analysis, the STRING network plot and the STRING enrichment plots
#' @export
STRING <- function(PoIs, STRINGBackground ,plotname = "", colPellet = "Blues"){

  ## Creating Input Dataframe
  Input <- data.frame(PoIs)
  colnames(Input) <- "Protein"

  ## Splitting Protein Names and Genes
  Input <- Input %>% dplyr::mutate(Gene = stringr::str_split_i(Protein, pattern = "_", 2))

  ## Extracting Gene List
  gene_list <- Input$Gene
  gene_df <- data.frame(gene = gene_list)

  background_gene_list <- STRINGBackground  %>% data.frame()
  colnames(background_gene_list) <- "Protein"

  background_gene_df <- background_gene_list %>%
    dplyr::mutate(Gene = stringr::str_split_i(Protein, pattern = "_", 2)) %>%
    dplyr::select(Gene)

  ## initilizing string object
  string_db <- STRINGdb$new(version = "11.0", species = 9606, score_threshold = 400, input_directory = "")


  ## Map genes to theit Identifier
  mapped_genes <- string_db$map(gene_df, "gene", removeUnmappedRows = TRUE)

  ## Creating Output Object
  Output <- list(PoIs = PoIs)

  ## make sure mapped genes are not empty, skip if empty
  if (nrow(mapped_genes) > 0) {

    ## map background genes
    background_mapped_genes <- string_db$map(background_gene_df, "Gene", removeUnmappedRows = TRUE)

    ## set background
    string_db$set_background(background_mapped_genes$STRING_id)

    ## perform enrichemnt Analysis

    ## Top 10 Processes
    ProcessEnrichment <- string_db$get_enrichment(mapped_genes$STRING_id ,category = "Process", iea = T) %>%
      dplyr::arrange(p_value) %>% utils::head(10) %>% dplyr::select(category,number_of_genes, p_value, description) %>% dplyr::mutate(GeneRatio = number_of_genes / length(PoIs)) %>%
      dplyr::mutate(`Number of Genes` = number_of_genes)


    ## Top 10 Functional Pathways
    FunctionEnrichment <- string_db$get_enrichment(mapped_genes$STRING_id, category = "Function", iea = T) %>%
      dplyr::arrange(p_value) %>% utils::head(10) %>% dplyr::select(category ,number_of_genes, p_value, description) %>% dplyr::mutate(GeneRatio = number_of_genes / length(PoIs)) %>%
      dplyr::mutate(`Number of Genes` = number_of_genes)


    ## creating results dataframe
    Results <- rbind(ProcessEnrichment,FunctionEnrichment)

    ## creating results plots
    ProcessPlot <- ggplot2::ggplot(Results %>% dplyr::filter(category == "Process"), aes(x = GeneRatio, y = reorder(description, `Number of Genes`),  size = `Number of Genes`))+
      geom_point(aes(col = -log10(p_value)))+
      xlim(0,NA) +
      ylab("") +
      xlab("Gene ratio")+
      ggtitle(paste("Process enrichment", plotname))+
      theme_light(base_size = 13)+
      scale_color_distiller(palette = colPellet, direction = 1)  # Reverse the palette so higher values have darker colors

    FunctionPlot <- ggplot2::ggplot(Results %>% dplyr::filter(category == "Function"), aes(x = GeneRatio, y = reorder(description, `Number of Genes`),  size = `Number of Genes`))+
      geom_point(aes(col = -log10(p_value)))+
      ylab("") +
      xlab("Gene ratio")+
      ggtitle(paste("Function enrichment", plotname))+
      theme_light(base_size = 13) +
      scale_color_distiller(palette = colPellet, direction = 1)  # Reverse the palette so higher values have darker colors

    ## populating List
    Output$ProcessPlot <- ProcessPlot
    Output$FunctionPlot <- FunctionPlot
    Output$Table <- Results
  }
  else {
    print("No mapped genes for the input list")
  }

  ## returning the output object
  return(Output)
}

## KEGGEnrichment
## add roxygen comments
#' @title KEGGEnrichment
#' @description KEGG enrichment analysis for the specified dataset.
#' @param dataset The dataset to be tested
#' @param PoIs A vector containing the Proteins of interest. Example: c("Q8TF72_SHROOM3" , "Q9ULZ3_PYCARD")
#' @param plotname The name to be displayed on created plots
#' @param folder The folder to save the pathway plots
#' @return A list object containing the results of the KEGG enrichment analysis, the KEGG enrichment plot. The pictures of the pathways are saved in the specified folder
#' @export
KeggEnrichmetn <- function(dataset, PoIs, folder = NA, plotname = ""){

  GoIs <- stringr::str_split_i(PoIs, "_", 2)

  ## encode GoIs
  gene_ids <- clusterProfiler::bitr(GoIs, fromType = "SYMBOL",
                                    toType = "ENTREZID",
                                    OrgDb = org.Hs.eg.db::org.Hs.eg.db)

  # Extract Entrez IDs
  gene_ids <- gene_ids$ENTREZID

  KeggBackground <- unique(dataset$Protein)
  BackgroundGoIs <- stringr::str_split_i(KeggBackground, "_", 2)
  KeggBackground_ids <- clusterProfiler::bitr(BackgroundGoIs, fromType = "SYMBOL",
                                              toType = "ENTREZID",
                                              OrgDb = org.Hs.eg.db::org.Hs.eg.db)

  # Extract Entrez IDs
  KeggBackground_ids <- KeggBackground_ids$ENTREZID

  KEGG <- clusterProfiler::enrichKEGG(gene = gene_ids, organism = 'hsa', keyType = 'kegg',
                                      qvalueCutoff = 0.01,
                                      pAdjustMethod = "BH",
                                      universe = KeggBackground_ids)

  ## skip rest of the function if KEGG == NULL
  if(is.null(KEGG)){
    print("No enriched KEGG pathways")
    return(NULL)
  }

  ## make KEGG plot
  KEGGTable <- KEGG@result %>% dplyr::as_tibble() %>%
    dplyr::filter(!category == "Metabolism") %>%
    filter(qvalue <0.05)

  ## skip if nrow of KEGGTable == 0
  if(nrow(KEGGTable) == 0){
    print("No enriched KEGG pathways")
    return(NULL)
  }

  KeggPlotData <- KEGGTable %>% dplyr::select(Description, GeneRatio, qvalue) %>%
    dplyr::arrange(qvalue) %>% utils::head(10) %>%
    dplyr::mutate(GenesFound = stringr::str_split_i(GeneRatio, "/", 1)) %>%
    dplyr::mutate(GenesTotal = stringr::str_split_i(GeneRatio, "/", 2)) %>%
    ## make Genesfound and Genetotal numeric
    dplyr::mutate(GenesFound = as.numeric(GenesFound), GenesTotal = as.numeric(GenesTotal)) %>%
    ## calculate gene ratio
    dplyr::mutate(GeneRatio = GenesFound/GenesTotal)

  KEGGPlot <- ggplot2::ggplot(KeggPlotData, ggplot2::aes(x = GeneRatio, y = stats::reorder(Description, GeneRatio),  size = GeneRatio))+
    ggplot2::geom_point(ggplot2::aes(col = -log10(qvalue)))+
    ggplot2::xlim(0,NA) +
    ggplot2::ylab("") +
    ggplot2::xlab("Gene ratio")+
    ggplot2::ggtitle(paste("KEGG enrichment", plotname))+
    ggplot2::theme_light(base_size = 13)

  ## make pathway plots

  gene_ids <- clusterProfiler::bitr(GoIs, fromType = "SYMBOL",
                                    toType = "ENTREZID",
                                    OrgDb = org.Hs.eg.db::org.Hs.eg.db) %>%
    dplyr::rename(Gene = SYMBOL, ID = ENTREZID)

  gene_ids <- gene_ids %>% dplyr::distinct(Gene, .keep_all = TRUE)

  WtestResults <- WTest(dataset %>% dplyr::filter(Protein %in% PoIs))

  PathwayData <- WtestResults$raw %>%
    dplyr::mutate(Gene = stringr::str_split_i(Protein,"_", 2)) %>%
    dplyr::select(Gene,estimate) %>%
    base::merge(gene_ids, by = "Gene") %>%
    dplyr::select(ID, estimate) %>%
    tibble::column_to_rownames(var= "ID")

  Pathways <-  KEGGTable %>% dplyr::as_tibble() %>%
    dplyr::arrange(qvalue) %>% utils::head(10) %>%
    dplyr::select(ID)

  ## set wd to folder
  if(!is.na(folder)){
    base::setwd(folder)
  }

  ## make KEGG plots with colored fold changes and save them to the folder
  for(i in 1:base::nrow(Pathways)){

    pathview::pathview(PathwayData, pathway.id = Pathways$ID[i], kegg.native = TRUE,
                       limit = list(gene = base::max(base::abs(PathwayData)), cpd = 1),
                       low = list(gene = "blue", cpd = "blue"),
                       mid = list(gene = "gray", cpd = "gray"),
                       high = list(gene = "red", cpd = "red"))
  }

  ## delete everything in folder that does not have "pathview" in its name
  if(!is.na(folder)){
    files <- base::list.files()
    files <- files[!base::grepl("pathview", files)]
    base::file.remove(files)
  }

  ## write KEGGTable to the folder
  if(!is.na(folder)){
    utils::write.csv(KEGGTable, file = base::paste0("KEGGTable.csv"), row.names = FALSE)
  }

  return(base::list(KEGGPlot = KEGGPlot,
                    KeggData = KEGGTable))
}


## High order functions

## qdd roxygen comments
#' @title SplineRegression
#' @description Spline regression for the specified dataset.
#' @param dataset The dataset to be tested
#' @param PoIs A vector containing the Proteins of interest. Example: c(""Q8TF72_SHROOM3" , "Q9ULZ3_PYCARD") or (unique(dataset$Protein))
#' @param Timecol The column containing the time points. Example "GestationalAge"
#' @param alpha The confidence level for the confidence intervals
#' @param split_status A boolean indicating if the model should be split by status
#' @param plotname The name to be displayed on created plots
#' @param show_all_proteins A boolean indicating if all proteins should be shown on the plot
#' @param B The number of bootstrap samples to use for estimating confidence intervals
#' @return A list object containing the results of the spline regression, the confidence intervals and the plot
#' @export
SplineRegression <- function(dataset, PoIs, Timecol, alpha = 0.05, split_status = FALSE, plotname = "", B = 100, col_pellet = "Set1", show_all_proteins = T) {
  GlobalData <- dataset

  # Ensure Timecol is treated as a column name
  Timecol <- rlang::ensym(Timecol)

  # Step 1: Fit a global spline model across all proteins
  global_spline <- stats::lm(Intensity ~ splines::ns(as.numeric(as.factor(dataset[[Timecol]])),
                                                     knots = length(unique(dataset[[Timecol]])) - 1),
                             data = GlobalData)

  # Step 2: Predict confidence intervals for the global model
  global_preds <- stats::predict(global_spline, newdata = GlobalData, interval = "confidence", level = (1 - alpha))
  GlobalData <- dplyr::mutate(GlobalData, global_fit = global_preds[, 1], global_lower = global_preds[, 2], global_upper = global_preds[, 3])

  # Calculate splines for every PoI
  SplineResults <- data.frame()

  for (i in seq_along(PoIs)) {
    PoI <- PoIs[[i]]
    subset_data <- dplyr::filter(dataset, Protein == PoI)

    if (split_status) {
      protein_spline <- stats::lm(Intensity ~ splines::ns(as.numeric(as.factor(subset_data[[Timecol]])),
                                                          knots = length(unique(subset_data[[Timecol]])) - 1) * Status,
                                  data = subset_data)
    } else {
      protein_spline <- stats::lm(Intensity ~ splines::ns(as.numeric(as.factor(subset_data[[Timecol]])),
                                                          knots = length(unique(subset_data[[Timecol]])) - 1),
                                  data = subset_data)
    }

    # Predict confidence intervals for the protein model
    protein_preds <- stats::predict(protein_spline, newdata = subset_data, interval = "confidence", level = (1 - alpha))

    subset_data <- dplyr::mutate(subset_data, fit = protein_preds[, 1], lower = protein_preds[, 2], upper = protein_preds[, 3])

    SplineResults <- dplyr::bind_rows(SplineResults, subset_data)
  }

  # Make plotable Data
  # Prepare data for plotting
  if (split_status) {
    # Perform bootstrap resampling B times to estimate confidence intervals for the mean fit per Timecol and Status
    PlotData <- replicate(B, {

      # Step 1: Create a bootstrapped sample by randomly sampling rows from SplineResults with replacement
      boot_sample <- dplyr::slice_sample(SplineResults, n = nrow(SplineResults), replace = TRUE)

      # Step 2: For each bootstrap sample, compute the mean fit per Timecol and Status
      boot_sample %>%
        dplyr::group_by(!!Timecol, Status) %>%
        dplyr::summarise(meanFit = mean(fit), .groups = "drop")  # Calculate the mean fit for each group

    }, simplify = FALSE) %>%  # Return the results as a list of data frames

      # Step 3: Combine all bootstrap results into a single data frame
      dplyr::bind_rows() %>%

      # Step 4: Group by Timecol and Status to compute the confidence intervals
      dplyr::group_by(!!Timecol, Status) %>%
      dplyr::summarise(

        # Compute lower bound of confidence interval (alpha/2 quantile, e.g., 2.5% for 95% CI)
        lower = stats::quantile(meanFit, probs = alpha / 2),

        # Compute upper bound of confidence interval (1 - alpha/2 quantile, e.g., 97.5% for 95% CI)
        upper = stats::quantile(meanFit, probs = 1 - alpha / 2),

        # Compute the mean of the bootstrapped meanFit values (expected mean trend)
        meanFit = mean(meanFit),

        .groups = "drop"  # Prevent unwanted grouping in subsequent operations
      )


  } else {
    # Perform bootstrap resampling B times to estimate confidence intervals for the mean fit per Timecol and Status
    PlotData <- replicate(B, {

      # Step 1: Create a bootstrapped sample by randomly sampling rows from SplineResults with replacement
      boot_sample <- dplyr::slice_sample(SplineResults, n = nrow(SplineResults), replace = TRUE)

      # Step 2: For each bootstrap sample, compute the mean fit per Timecol
      boot_sample %>%
        dplyr::group_by(!!Timecol) %>%
        dplyr::summarise(meanFit = mean(fit), .groups = "drop")  # Calculate the mean fit for each group

    }, simplify = FALSE) %>%  # Return the results as a list of data frames

      # Step 3: Combine all bootstrap results into a single data frame
      dplyr::bind_rows() %>%

      # Step 4: Group by Timecol to compute the confidence intervals
      dplyr::group_by(!!Timecol) %>%
      dplyr::summarise(

        # Compute lower bound of confidence interval (alpha/2 quantile, e.g., 2.5% for 95% CI)
        lower = stats::quantile(meanFit, probs = alpha / 2),

        # Compute upper bound of confidence interval (1 - alpha/2 quantile, e.g., 97.5% for 95% CI)
        upper = stats::quantile(meanFit, probs = 1 - alpha / 2),

        # Compute the mean of the bootstrapped meanFit values (expected mean trend)
        meanFit = mean(meanFit),

        .groups = "drop"  # Prevent unwanted grouping in subsequent operations
      ) %>%
      mutate(data = "Protein")

  }

  ## make global plot data
  GlobalPlotData <- GlobalData %>%
    dplyr::group_by(!!Timecol) %>%
    dplyr::summarise(lower = stats::quantile(global_lower, probs = alpha / 2),
                     upper = stats::quantile(global_upper, probs = 1 - alpha / 2),
                     meanFit = mean(global_fit), .groups = "drop") %>%
    mutate(data = "Global")

  ## make one bit dataframe from plotdata and glabalplotdata
  PlotDataAll <- dplyr::bind_rows(PlotData, GlobalPlotData)

  # Generate output plot
  colors <- assign_colors(unique(dataset %>% arrange(Status) %>% pull(Status)),palette = col_pellet)

  ## generate output ploy

  Plot <- ggplot2::ggplot(data = PlotDataAll) +
    ## plot global spline
    ggplot2::geom_line(data = GlobalPlotData, ggplot2::aes(x = as.numeric(as.factor(!!Timecol)), y = meanFit), color = "black") +
    ## plot CI of global spline
    ggplot2::geom_ribbon(data = GlobalPlotData, ggplot2::aes(x = as.numeric(as.factor(!!Timecol)), ymin = lower, ymax = upper), fill = "grey", alpha = 0.5)


  if (split_status == T) {

    if (show_all_proteins == T){
      ## plot Ribbons for all the proteins
      for(i in 1:length(PoIs)){

        Plot <- Plot +

          ## plot CI of Protein splines
          ggplot2::geom_line(data = SplineResults %>% dplyr::filter(Protein == PoIs[i]), ggplot2::aes(x = as.numeric(as.factor(!!Timecol)), y = fit ,col = Status), alpha = 0.2)

      }

    }


    ## plot mean Protein splines for each status
    Plot <- Plot +
      ggplot2::geom_line(data = PlotData, ggplot2::aes(x = as.numeric(as.factor(!!Timecol)), y = meanFit, color = Status), linewidth = 2) +
      ## plot CI of Protein splines for each status
      ggplot2::geom_ribbon(data = PlotData, ggplot2::aes(x = as.numeric(as.factor(!!Timecol)), ymin = lower, ymax = upper, fill = Status), alpha = 0.5) +
      ## Apply pre-defined colors
      ggplot2::scale_color_manual(values = colors) +
      ggplot2::scale_fill_manual(values = colors)


  }

  if (split_status == F) {

    if (show_all_proteins == T){

      ## plot ribbons for every protein
      for(i in 1:length(PoIs)){

        Plot <- Plot +
          ## plot CI of Protein splines
          ggplot2::geom_line(data = SplineResults %>% dplyr::filter(Protein == PoIs[i]), ggplot2::aes(x = as.numeric(as.factor(!!Timecol)), y = fit),col = "grey", alpha = 0.5)

      }

    }

    ## plot mean Protein spline
    Plot <- Plot +
      ggplot2::geom_line(data = PlotData, ggplot2::aes(x = as.numeric(as.factor(!!Timecol)), y = meanFit), linewidth = 2, color = "black") +
      # plot CI for mean spline
      ggplot2::geom_ribbon(data = PlotData, ggplot2::aes(x = as.numeric(as.factor(!!Timecol)), ymin = lower, ymax = upper), fill = "grey", alpha = 0.5) +

      ## Apply colors
      ggplot2::scale_fill_manual(values = colors)


  }

  Plot <- Plot +
    ggplot2::ggtitle(paste(plotname)) +
    ggplot2::theme_minimal() +
    ## make the x axis the character entries in Timecol
    ggplot2::scale_x_continuous(breaks = as.numeric(unique(dataset[[Timecol]])), labels = unique(dataset[[Timecol]])) +
    ## rotate x-axis labels by 45 degrees
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ## make x axis label Timecol
    ggplot2::xlab(rlang::as_label(Timecol))


  # Create output list
  Output <- list(
    GlobalPreds = GlobalData,
    PoIPreds = SplineResults,
    Plot = Plot,
    PlotDataAll = PlotDataAll
  )

  return(Output)
}



## add roxygen comments
#' @title HeatMapProteinClusterAnalysis
#' @description Extension of the Heatmap function. This function takes a heatmap object and clusters the proteins into n clusters. It then performs gene ontology analysis on each cluster.
#' @param HeatMap The heatmap object to be analyzed
#' @param n The number of Protein lusters to be created
#' @return A list object containing the proteins in each cluster, the original heatmap, the modified heatmap and the gene ontology analysis for each cluster
#' @export
HeatMapProteinClusterAnalysis <- function(HeatMap, STRINGBackground ,n){

  ## safe original heat map
  OriginalHeatMap <- HeatMap

  ## set the row split of the heat map to n
  HeatMap@matrix_param$row_split <- n

  ## Save modified heat map
  ModifiedHeatMap <- HeatMap

  ## get List of all Proteins
  Proteins <- row.names(HeatMap@matrix)

  ## get Names of Proteins in each cluster
  RowOrder <- ComplexHeatmap::row_order(HeatMap)

  ## Creating a list to store proteins in each cluster
  ProteinsInCluster <- list()

  ## Looping over the clusters and extract protin names
  for(i in 1 :n){

    ## get cluster
    Cluster <- RowOrder[[i]]

    ## get proteins in cluster
    ProteinsInCluster[[i]] <- Proteins[Cluster]

  }

  ## Creating output object
  Output <- list(ProteinsInCluster = ProteinsInCluster,
                 OriginalHeatMap = OriginalHeatMap,
                 ModifiedHeatMap = ModifiedHeatMap)

  ## looping over the List and perform gene ontology analysis, store outputs separately
  for (i in 1:n){
    ## Getting Proteins
    Proteins <- ProteinsInCluster[[i]]

    ## Perform Gene Ontology Analysis
    GO <- STRING(Proteins, plotname = paste0("Cluster ", i), STRINGBackground = STRINGBackground)
    Sys.sleep(1)

    ## Store Output
    Output$STRINGResults[[paste0("Cluster_", i)]] <- GO

  }

  ## return output
  return(Output)

}

## TwoWayComparison
## Advanced function that combines the Wilcox test, PCA and STRING analysis
## add roxygen comments
#' @title TwoWayComaprison
#' @description Two way comparison of the specified dataset.
#' @param dataset The dataset to be tested
#' @param plotname The name to be displayed on the resulting plots
#' @return A list object containing the results of the Wilcox test, the PCA as well a UMAP analysis and the STRING analysis if there are significant results
#' @export
TwoWayComparison <- function(dataset, plotname = ""){

  ## Error when Data has more than 2 groups
  if(length(unique(dataset$Status)) > 2){
    stop("Data has more than 2 groups in Status")
  }

  ## Wilcox Test
  WilcoxResults <- WTest(dataset, plotname = plotname)

  ## PCA
  PCAResutls <- PCA(dataset %>% dplyr::filter(Protein %in% unique(WilcoxResults$raw$Protein)), plotname = plotname)

  ## UMAP
  UMAPResults <- UMAP(dataset %>% dplyr::filter(Protein %in% unique(WilcoxResults$raw$Protein)), plotname = plotname)


  ## Gene Ontology Analysis (if there are significant results)
  if(length(unique(WilcoxResults$Significant$Protein > 0))){



    ## Setting up Protein Lists
    PoIsUp <- WilcoxResults$Significant %>%
      dplyr::filter(estimate > 0) %>%
      dplyr::pull(Protein)

    PoIsDown <- WilcoxResults$Significant %>%
      dplyr::filter(estimate < 0) %>%
      dplyr::pull(Protein)

    ## String analysis
    StringResultsUP <- STRING(PoIsUp, plotname = paste("Upregulated Proteins", plotname), STRINGBackground = unique(dataset$Protein))
    StringResultsDown <- STRING(PoIsDown, plotname = paste("Downregulated Proteins", plotname), STRINGBackground = unique(dataset$Protein))
 }

  ## create Output List
  OutputList <- list(WilcoxResults = WilcoxResults)

  ## add PCA and umap if availabe
  if(exists("PCAResutls")){
    OutputList$PCAResutls <- PCAResutls
    OutputList$UMAPResults <- UMAPResults
  }

  ## add STRING results if available
  if(exists("StringResultsUP")){
    OutputList$StringResults$Up<- StringResultsUP
    OutputList$StringResults$Down <- StringResultsDown
  }
  ## return Output List
  return(OutputList)
}
## MultiWaycomparison
## A advanced function that combines the Kruskal test and PCA
## add roxygen comments
#' @title MultiWayComaprison
#' @description Multi way comparison of the specified dataset.
#' @param dataset The dataset to be tested
#' @param plotname The name to be displayed on the resulting plots
#' @return A list object containing the results of the Kruskal test and the PCA, as well a UMAP analysis
#' @export
MultiWayComaprison <- function(dataset, plotname = ""){

  ## Error when Data has  2 groups
  if(length(unique(dataset$Status)) == 2){
    stop("Data has 2 groups in Status, plase use TwoWayComparison() instead")
  }

  ## Error When Data has less than 3 groups
  if(length(unique(dataset$Status)) == 1){
    stop("Data needs at least 3 groups in Staus")
  }

  ## Wilcox Test
  KruskalResults <- KruskalTest(dataset, plotname = plotname)

  ## PCA
  PCAResutls <- PCA(dataset %>% dplyr::filter(Protein %in% unique(KruskalResults$raw$Protein)), plotname = plotname)

  ## UMAP
  UMAPResults <- UMAP(dataset %>% dplyr::filter(Protein %in% unique(KruskalResults$raw$Protein)), plotname = plotname)


  if(nrow(KruskalResults$Significant) > 0){

  }

  ## create Output List
  OutputList <- list(
    KruskalResults = KruskalResults
  )

  ## add PCA and umap if availabe
  if(exists("PCAResutls")){
    OutputList$PCAResutls <- PCAResutls
    OutputList$UMAPResults <- UMAPResults
  }

  ## return Output List
  return(OutputList)
}


## Multiscale Embedded Gene Co-expression Network Analysis (MEGENA)
## add roxygen comments
#' @title MEGENA
#' @description Multiscale Embedded Gene Co-expression Network Analysis (MEGENA) for the specified dataset.
#' @param dataset The dataset to be tested
#' @return A list object containing the results of the MEGENA analysis, the module table, the hierarchy plot, the STRING analysis and the correlation plot
#' @export
MEGENA <- function(dataset, plotname = ""){

  TestData <- dataset %>%
    dplyr::select(Sample, Protein, Intensity) %>%
    tidyr::pivot_wider(names_from = Protein, values_from = Intensity) %>%
    tidyr::pivot_longer(cols = -Sample, names_to = "Protein", values_to = "Intensity") %>%
    ImputeFeatureIntensity() %>%
    tidyr::pivot_wider(names_from = Protein, values_from = Intensity) %>%
    tibble::column_to_rownames(var = "Sample") %>%
    as.matrix() %>% t()

  # Calculate the correlation matrix
  cor_matrix <- MEGENA::calculate.correlation(TestData)

  # Construct the PFN
  pfn <- MEGENA::calculate.PFN(cor_matrix)

  ## Construct the graph
  g <- igraph::graph_from_data_frame(pfn, directed = FALSE)

  ## Calculate the modules
  MEGENA.output <- MEGENA::do.MEGENA(g)

  ## Summarize the modules
  summary.output <- MEGENA::MEGENA.ModuleSummary(MEGENA.output,
                                                 min.size = 50, max.size = igraph::vcount(g)/2,
                                                 output.sig = TRUE)

  ## Extract module table
  module.table <- summary.output$module.table
  colnames(module.table)[1] <- "id" # first column of module table must be labelled as "id".

  ## make hierarchy plot
  hierarchy.obj <- MEGENA::plot_module_hierarchy(module.table = module.table, label.scaleFactor = 0.15,
                                                 arrow.size = 0.03, node.label.color = "blue")
  # Hierarchy plot
  HirarchyPlot <- hierarchy.obj[[1]]

  ## modify Hierarchy plot
  HirarchyPlot <- HirarchyPlot +
    ggplot2::ggtitle(plotname) +
    ## delete legend
    ggplot2::theme(legend.position = "none")

  ## Correlate modules with Status

  ## Proteins in Modules
  ProteinsInModules <- data.frame()
  ## extract modules
  for (i in 1:length(summary.output$modules)){
    Proteins <- summary.output$modules[[i]]
    modules <- names(summary.output$modules)[i]
    ProteinsInModules <- rbind(ProteinsInModules, data.frame(Protein = Proteins, Module = modules))
  }

  CorrelationData <- base::merge(dataset, ProteinsInModules, by = "Protein")


  ## calculate correlation between the modules and meta variables
  vars <- colnames(dataset %>% dplyr::select(-c(Protein, Intensity, Sample)))

  CorrelationResultsList <- list()

  for(i in 1:length(vars)){
    ## define var
    var <- vars[i]

    ## make var a factor
    CorrelationData[[var]] <- base::as.numeric(base::as.factor(CorrelationData[[var]]))

    ## Calculate Correlations and p-values
    CorrelationResults <- CorrelationData %>%
      dplyr::group_by(Module) %>%
      dplyr::summarize(correlation = stats::cor.test(x = Intensity,y = .data[[var]], method = "spearman")$estimate,
                                               p.value = stats::cor.test(Intensity, .data[[var]])$p.value) %>%
      dplyr::mutate(p.value.adj = stats::p.adjust(p.value, method = "fdr")) %>%
      dplyr::arrange(correlation)

    ## Plot Correlation Results
    CorrelationPlot <- ggplot2::ggplot(CorrelationResults, ggplot2::aes(x = correlation, y = -log10(p.value.adj))) +
      ggplot2::geom_point() +
      ggplot2::geom_abline(intercept = -log10(0.05), slope = 0, color = "grey", linetype = 2) +
      ggplot2::geom_abline(intercept = -log10(0.01), slope = 0, color = "red", linetype = 2) +
      ggrepel::geom_text_repel(ggplot2::aes(label = Module), box.padding = 0.5) +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = "Correlation", y = "-log10(p.adj)") +
      ggplot2::ggtitle(paste("Correlation between", var, "and Modules"))

    ## put results into list
    CorrelationResultsList$Data[[var]] <- CorrelationResults
    CorrelationResultsList$Plot[[var]] <- CorrelationPlot

  }


  output <- list()
  ## populate output object
  output$MEGENA <- MEGENA.output
  output$Summary <- summary.output
  output$ModuleTable <- module.table
  output$HirarchyPlot <- HirarchyPlot
  output$CorrelationResults <- CorrelationResultsList
  output$CorrelationMatrix <- cor_matrix

  return(output)
}


## ToDo
## Data Manipulation
## Machine learning
## MISC
