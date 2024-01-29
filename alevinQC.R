#!/usr/bin/env Rscript

# Check if alevinQC is installed, and install it if it's not
if (!require("alevinQC", character.only = TRUE)) {
  install.packages("alevinQC", repos = "http://bioconductor.org/packages/release/bioc/")
  if (!require("alevinQC", character.only = TRUE)) stop("Package alevinQC could not be loaded or installed.")
}

library(alevinQC)


### Set the working directory
working_dir <- ''
setwd(working_dir)

# Now, rest of the script
dir <- getwd()

alevinQCReport(baseDir = dir,                      # Base directory
               sampleId = "test",         # Replace with a sample ID of your choice
               outputFile = "alevinReport.html",   # Name of the output file
               outputFormat = "html_document",     # Format of the output file
               outputDir = dir,                    # Directory where you want to save the report
               forceOverwrite = TRUE)
