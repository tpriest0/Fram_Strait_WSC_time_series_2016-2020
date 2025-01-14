####################

#### Determine oscillation signals

####################

### Define working directory
setwd('XXXXX')

### Define directories
output_figures = ('results/figures/')
output_tables = ('results/tables/')
scripts = ("scripts/")

### Load libraries"
required_libraries <- c("readr", "dplyr", "data.table")
  for (lib in required_libraries) {
    if (!require(lib, character.only = TRUE)) {
      install.packages(lib, dependencies = TRUE)
      library(lib, character.only = TRUE)
    }
  }

determine_oscillations <- function(input_file, output_file) {
  # Function to split a row into four parts
  split_row <- function(row, N) {
    part_size <- length(row) %/% N  # Use integer division
    s1 <- row[1:part_size]                # First part
    s2 <- row[(part_size + 1):(2 * part_size)]  # Second part
    s3 <- row[(2 * part_size + 1):(3 * part_size)]  # Third part
    s4 <- row[(3 * part_size + 1):length(row)]   # Fourth part
    splits <- list(s1, s2, s3, s4)
    return(splits)
  }
  
  # Function to calculate the number of oscillations using Fourier transformation
  oscillations_from_fourier <- function(time_series) {
    time_series <- as.numeric(time_series)
    fourier_transform <- fft(time_series)
    num_points <- length(time_series)
    frequencies <- abs(fourier_transform) / num_points
    frequencies[1] <- frequencies[1] / 2
    num_components <- num_points / 2 + 1
    freq_index <- which.max(frequencies[2:num_components])
    number_of_oscillations_from_fourier <- freq_index
    return(number_of_oscillations_from_fourier)
  }
  
  # Read the data from the CSV file
  interpolated_abund <- read.csv(input_file, sep="\t", check.names=F)
  
  # Store the gene names separately
  ids <- rownames(interpolated_abund)

  # Convert column names to Date format
  colnames(interpolated_abund) <- as.Date(colnames(interpolated_abund), origin = "1900-01-01")
  
  # Define the number of parts (years) for splitting
  N <- 4
  
  # Initialize empty lists to store results
  id_list <- vector("list", length = nrow(interpolated_abund) * N)
  year_list <- rep(c("year1", "year2", "year3", "year4"), times = nrow(interpolated_abund))
  oscillations_list <- vector("list", length = nrow(interpolated_abund) * N)
  
  # Split each row, calculate oscillations, and store results in the lists
  for (i in 1:nrow(interpolated_abund)) {
    splits <- split_row(interpolated_abund[i, ], N)
    for (j in 1:N) {
      id_list[[N * (i - 1) + j]] <- ids[i]
      oscillations_list[[N * (i - 1) + j]] <- oscillations_from_fourier(splits[[j]])
    }
  }
  
  # Create a dataframe to store the results
  result_df <- data.frame(ID = unlist(id_list),
                          year = unlist(year_list),
                          oscillations = unlist(oscillations_list))
  
  # Save the dataframe to a CSV file
  write.csv(result_df, file = output_file, row.names = FALSE, sep="\t")
  
  # Print the first few rows of the result dataframe as feedback
  print(head(result_df))
  
  # Return the output file path
  cat("Processing completed. Results saved in:", output_file, "\n")
}

# Apply oscillation signal function to ASV and gene cluster interpolated abundance profiles
determine_oscillations(paste0(output_tables,"FRAM_RAS_F4_MIC_ASV_interpolated_abundances.txt"), paste0(output_tables,"FRAM_RAS_F4_MIC_ASV_oscillations_per_year.txt"))
determine_oscillations(paste0(output_tables,"FRAM_RAS_F4_EUK_ASV_interpolated_abund.txt"), paste0(output_tables,"FRAM_RAS_F4_EUK_ASV_oscillations_per_year.txt"))
determine_oscillations(paste0(output_tables,"FRAM_RAS_F4_GENE_CLUST_interpolated_abund.txt"), paste0(output_tables,"FRAM_RAS_F4_GENE_CLUST_oscillations_per_year.txt"))
