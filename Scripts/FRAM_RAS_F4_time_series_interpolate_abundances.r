### Define working directory
setwd('XXXXX')

### Define directories
output_figures = ('results/figures/')
output_tables = ('results/tables/')
scripts = ("scripts/")

### Load necessary scripts
source(paste0(scripts,"FRAM_RAS_F4_calcTimeSeriesInfo.R"))
source(paste0(scripts,"FRAM_RAS_F4_daily_approx.R"))

### Load libraries
required_libraries <- c("segmenTier", "segmenTools", "igraph", "Hmisc", "reshape2", "ggplot2", "readr", "data.table", "dplyr", "tibble")
  for (lib in required_libraries) {
    if (!require(lib, character.only = TRUE)) {
      install.packages(lib, dependencies = TRUE)
      library(lib, character.only = TRUE)
    }
  }

### Function to generate interpolated relative abundance profiles

interpolate_abundances <- function(input_file, ras_id_file, output_file) {
  
  # Load the input files
  rel_profile <- read.csv(file=input_file, sep="\t", header=T)
  #rel_profile <- read_delim(input_file, delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  # Erste Spalte als rownames setzen (und aus den Daten entfernen, falls nötig)
  rn = rel_profile[[1]]
  rel_profile[,1] <- NULL
  rownames(rel_profile) <- rn # Zeilennamen zuweisen
  metadata <- read.delim2(ras_id_file)
  
  # Extract ASV names and adjust columns
  ASVName <- rownames(rel_profile)
  rel_profile <- rel_profile[, 2:ncol(rel_profile)]
  id <- colnames(rel_profile)
  id <- gsub("X([0-9])_", "0\1_", id)
  id <- gsub("X", "", id)
  
  # Match IDs and assign dates
  mi <- match(id, metadata$RAS_id)
  dt <- metadata$date[mi]
  colnames(rel_profile) <- dt
  
  # Order data by time
  oi <- order(dt)
  rel_profile <- rel_profile[, oi]
  rel_profile <- apply(rel_profile, 2, as.numeric)
  rownames(rel_profile) <- rn
  
  # Calculate NaNs per row and filter
  nans_per_row <- apply(rel_profile, 1, function(x) sum(is.na(x)))
  num_of_genes <- length(nans_per_row)
  wi <- which(nans_per_row < 44)  # Keep rows with fewer than 44 NaNs
  rel_profile <- rel_profile[wi, ]
  
  # Filter statistics
  num_genes_filtered <- length(wi)
  drop_out_rate <- 1 - (num_genes_filtered / num_of_genes)
  
  # Replace all NaNs with 0 and sort by date
  rel_profile[is.na(rel_profile)] <- 0
  colnames(rel_profile) <- as.character(as.Date(colnames(rel_profile), "%d.%m.%Y"))
  rel_profile <- rel_profile[, order(colnames(rel_profile))]
  #rel_profile <- rel_profile[, 1:46]
  
  # Time series analysis using SegmenTier
  tset <- processTimeseries(rel_profile, use.fft = TRUE, dft.range = 2:12, use.snr = TRUE, na2zero = TRUE)
  ee <- tset$dat
  
  # Remove rows with only NAs
  ee1 <- ee[rowSums(is.na(ee)) != ncol(ee), ]
  
  # Calculate time series signals
  allTS <- apply(ee1, 1, calcTimeSeriesInfo)
  
  # Extract sinus model data
  xx <- lapply(allTS, function(x) x$mod_dataTS)
  xx <- do.call("rbind", xx)
  xx <- apply(xx, 1, function(x) (x - min(x)) / max(x - min(x)))
  xx <- t(xx)
  
  # Define time variable (one sinus cycle = one year)
  ts <- seq(0, 2 * pi, length.out = ncol(xx))
  
  # Calculate phase and frequency
  yy <- lapply(allTS, function(x) x$phase_freq)
  yy <- do.call("rbind", yy)
  
  # Sort data by frequency
  yy <- as.data.frame(cbind(yy, ii = seq(nrow(yy))))
  yy_sorted <- yy[order(yy$freq), ]

  # Daily approximation
  dap_set2All <- daily_approx(as.Date(colnames(tset$ts)), tset$ts, fix = FALSE)
  colnames(dap_set2All$ee) <- as.character(dap_set2All$dtw)
  
  # Write the results to a CSV file
  write.table(dap_set2All$ee, file = output_file, quote=F, sep="\t")
  
  # User feedback
  cat("File processing completed. Results saved in:", output_file, "\n")
}

interpolate_abundances(paste0(output_tables,"FRAM_RAS_F4_MIC_ASV_filt_rel.txt"), "FRAM_RAS_F4_META.txt", paste0(output_tables,"FRAM_RAS_F4_MIC_ASV_interpolated_abundances.txt"))
interpolate_abundances(paste0(output_tables,"FRAM_RAS_F4_EUK_ASV_filt_rel.txt"), "FRAM_RAS_F4_META.txt", paste0(output_tables,"FRAM_RAS_F4_EUK_ASV_interpolated_abundances.txt"))
interpolate_abundances(paste0(output_tables,"FRAM_RAS_F4_GENE_CLUST_filt_rel.txt"), "FRAM_RAS_F4_META.txt", paste0(output_tables,"FRAM_RAS_F4_GENE_CLUST_interpolated_abundances.txt"))
