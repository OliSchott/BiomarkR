.onLoad <- function(libname, pkgname) {
        install_if_missing <- function(pkg) {
                if (!requireNamespace(pkg, quietly = TRUE)) {
                        message(paste("Installing missing package:", pkg))
                        if (pkg %in% c("clusterProfiler", "org.Hs.eg.db", "STRINGdb", "ComplexHeatmap", "pathview", "sva")) {
                                BiocManager::install(pkg, ask = FALSE)
                        } else {
                                install.packages(pkg)
                        }
                }
        }
        
        # Ensure BiocManager is available
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
                install.packages("BiocManager")
        }
        
        # List of dependencies
        packages <- c(
                "stats", "tidyverse", "readr", "readxl", "rstatix", "ComplexHeatmap",
                "InteractiveComplexHeatmap", "pROC", "caret", "circlize", "umap",
                "stringi", "gtools", "seqinr", "sva", "dendextend", "STRINGdb",
                "MLmetrics", "plotly", "ggrepel", "magrittr", "igraph", "Rtsne",
                "grid", "MEGENA", "clusterProfiler", "org.Hs.eg.db", "pathview"
        )
        
        # Install missing dependencies
        lapply(packages, install_if_missing)
}
