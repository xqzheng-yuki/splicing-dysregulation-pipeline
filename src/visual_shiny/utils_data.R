# utils_data.R
# Utility functions for handling file paths, loading input data (e.g. BAM, BigWig, BED),
# managing gene ID/name parsing, and helper functions for building Gviz tracks.
# This script supports both plotting logic and server-side reactivity.

source("R/file_paths.R")
source("R/bed_utils.R")
source("R/track_helpers.R")
source("R/gene_utils.R")