#!/usr/bin/env Rscript
# Restore the recorded packages into an isolated library, without editing .Rprofile.
args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grepl("^--file=", args)])
if (length(script) != 1L) stop("Run with Rscript --vanilla environment/restore.R")
root <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)
library_path <- file.path(root, "environment/library")
dir.create(library_path, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(library_path, .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv", lib = library_path)
renv::restore(project = root, lockfile = file.path(root, "environment/renv.lock"),
              library = library_path, prompt = FALSE)
cat("Restored packages into ", library_path, "\n", sep = "")
cat("scripts/run_analysis.R uses this library when it is present.\n")
