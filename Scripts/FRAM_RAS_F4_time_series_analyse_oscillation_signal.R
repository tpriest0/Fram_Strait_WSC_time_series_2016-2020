####################

#### Processing oscillation signal data

####################

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
required_libraries <- c("segmenTier", "segmenTools", "stringr", "dplyr", "reshape2", "ggplot2", "readr", "tibble", "patchwork", "data.table", "lubridate")
  for (lib in required_libraries) {
    if (!require(lib, character.only = TRUE)) {
      install.packages(lib, dependencies = TRUE)
      library(lib, character.only = TRUE)
    }
  }

### Import oscillation tables
gene_clust_oscillations = read.table(paste0(output_tables,"FRAM_RAS_F4_GENE_CLUST_oscillations_per_year.txt"), sep="\t", header=T)
euk_asv_oscillations = read.table(paste0(output_tables,"FRAM_RAS_F4_EUK_ASV_oscillations_per_year.txt"), sep="\t", header=T)
prok_asv_oscillations = read.table(paste0(output_tables,"FRAM_RAS_F4_MIC_ASV_oscillations_per_year.txt"), sep="\t", header=T)

### Import relative abudnance tables
mic_asv_rel=read.table(file=paste0(output_tables,"FRAM_RAS_F4_MIC_ASV_filt_rel.txt"), sep="\t",
                       check.names=F, header=T, row.names=1)
euk_asv_rel=read.table(file=paste0(output_tables,"FRAM_RAS_F4_EUK_ASV_filt_rel.txt"), sep="\t",
                       check.names=F, header=T, row.names=1)
gene_clust_rel=read.table(
  file=paste0(output_tables,"FRAM_RAS_F4_GENE_CLUST_filt_rel.txt"), 
  sep="\t", check.names=F, header=T, row.names=1)

### Combine oscillation from asvs and gene clusters together
asv_and_geneclust_oscillations = rbind(gene_clust_oscillations,euk_asv_oscillations,prok_asv_oscillations) %>%
    mutate(Type = case_when(
        grepl("prok",ID) ~ "Prok",
        grepl("euk",ID) ~ "Euk",
        grepl("gene",ID) ~ "Gene_clust"
    ))

### Assess the oscillation signal data retrieved from the time-series analysis
# Calculate the mean oscillations per year
asv_and_geneclust_oscillations_mean = 
    asv_and_geneclust_oscillations %>%
    aggregate(oscillations~ID+Type, data=., FUN=mean) %>%
    mutate(Type = factor(Type, levels=c("Prok","Euk","Gene_clust")))

# For how many ASVs and gene clusters could oscillations be determined? And
# How many exhibited a single oscillation per year?
dim(mic_asv_rel)
asv_and_geneclust_oscillations_mean %>%
    filter(., Type == "Prok") %>%
    dim()
asv_and_geneclust_oscillations_mean %>%
    filter(., Type == "Prok" & 
    oscillations == 1) %>%
    dim()
(658/3748)*100

dim(euk_asv_rel)
asv_and_geneclust_oscillations_mean %>%
    filter(., Type == "Euk") %>%
    dim()
asv_and_geneclust_oscillations_mean %>%
    filter(., Type == "Euk" &
    oscillations == 1) %>%
    dim()
(448/3019)*100

dim(gene_clust_rel)
asv_and_geneclust_oscillations_mean %>%
    filter(., Type == "Gene_clust") %>%
    dim()
asv_and_geneclust_oscillations_mean %>%
    filter(., Type == "Gene_clust" &
    oscillations == 1) %>%
    dim()
(178586/306088)*100

### Determine fraction of community represented by each oscillation type

### FUNCTION: calculate relative abund in each sample for each oscillation pattern
### for microbial and microeukaryotic ASVs and gene clusters
calc_rel_abund_per_osc <- function(infile){
  osc_abund_df = 
    infile %>% 
    tibble::rownames_to_column(., var="ID") %>%
    filter(ID %in% asv_and_geneclust_oscillations_mean$ID) %>%
    reshape2::melt(id.vars="ID", variable.name="RAS_id", value.name="Rel_abund") %>%
    filter(., Rel_abund > 0) %>%
    left_join(asv_and_geneclust_oscillations_mean, by="ID") %>%
    mutate(Rel_abund = as.numeric(Rel_abund)*100) %>%
    aggregate(Rel_abund~RAS_id+Type+oscillations, data=., FUN=sum)
  return(osc_abund_df)
}

mic_asv_osc_per_year_rel_prop <- calc_rel_abund_per_osc(mic_asv_rel)
euk_asv_osc_per_year_rel_prop <- calc_rel_abund_per_osc(euk_asv_rel)
gene_clust_osc_per_year_rel_prop <- calc_rel_abund_per_osc(gene_clust_rel)

### Combine dataframes 
asv_gene_rel_prop_per_oscillation <- 
  rbind(mic_asv_osc_per_year_rel_prop,
        euk_asv_osc_per_year_rel_prop,
        gene_clust_osc_per_year_rel_prop) %>%
  mutate(Type = factor(Type, levels=c("Prok","Euk","Gene_clust")))

# Calculate the combined, mean relative abundance of asvs and genes
# with different oscillation signals

### Visualise relative abundance of ASVs and GENEs with different 
### oscillation signals 
plot_rel_prop_per_oscillation = function(infile){
  ggplot(data=infile,
         aes(y = oscillations, x = Rel_abund)) +
    geom_boxplot(aes(group=oscillations)) + 
    facet_wrap(Type~., scales="free_x", nrow = 4,
               strip.position = "top") + 
    theme_bw() + 
    scale_y_reverse() +
    labs(y = "Mean number of oscillations per year", x = "Relative abundance (%)") + 
    theme(strip.background.x = element_rect(fill = "white", colour = "black"),
          strip.text.x = element_text(size = 12),
          axis.title.y = element_blank(),
          axis.text.y = element_text(size = 11),
          axis.title.x = element_text(size = 14),
          axis.text.x = element_text(size = 11),
          plot.background = element_blank())
} 

### Visualise number of oscillations per year for ASVs and GENEs
plot_num_oscillations_py = function(infile){
  ggplot(data=infile,
         aes(y = oscillations)) +
    geom_histogram() + 
    facet_wrap(Type~., scales = "free_x", nrow = 4) + 
    theme_bw() + 
    scale_y_reverse() + 
    labs(y = "Number of oscillations per year", x = "Count") + 
    theme(strip.background.x = element_rect(fill = "white", colour = "black"),
          strip.text.x = element_text(size = 12),
          axis.title.y = element_text(size = 14),
          axis.text.y = element_text(size = 11),
          axis.title.x = element_text(size = 14),
          axis.text.x = element_text(size = 11),
          plot.background = element_blank())
} 

### Create figure illustrating mean oscillations and rel prop per oscillation signal
plot_asv_gene_mean_osc_per_year <- plot_num_oscillations_py(asv_and_geneclust_oscillations_mean)
plot_asv_gene_rel_prop_per_osc <- plot_rel_prop_per_oscillation(asv_gene_rel_prop_per_oscillation)

# Export figure
pdf(file=paste0(output_figures,"FRAM_RAS_F4_MIC_EUK_ASV_GENE_OSC_per_year_and_rel_abund.pdf"),
    height=10, width=8)
plot_asv_gene_mean_osc_per_year|plot_asv_gene_rel_prop_per_osc
dev.off()

### Create a list of ASVs and genes that exhibit one oscillation per year year
mic_euk_asv_gene_clust_osc4_id_list = 
  asv_and_geneclust_oscillations_mean %>%
  filter(., oscillations == 1) %>%
  select(ID)

### Create relative abundance dataframes for those with one oscillation per year as these will be 
### used as the input for network analysis
mic_asv_osc4_rel_wide = 
  mic_asv_rel %>%
  tibble::rownames_to_column(., var="ID") %>%
  filter(ID %in% mic_euk_asv_gene_clust_osc4_id_list$ID)

euk_asv_osc4_rel_wide = 
  euk_asv_rel %>%
  tibble::rownames_to_column(., var="ID") %>%
  filter(ID %in% mic_euk_asv_gene_clust_osc4_id_list$ID)

gene_clust_osc4_rel_wide = 
  gene_clust_rel %>%
  tibble::rownames_to_column(., var="ID") %>%
  filter(ID %in% mic_euk_asv_gene_clust_osc4_id_list$ID)

# Export tables
write.table(mic_euk_asv_gene_clust_osc4_id_list,
            file=paste0(output_tables,"FRAM_RAS_F4_ASV_GENE_CLUST_OSC4_ID_list.txt"),
            sep="\t", quote = F, row.names = F)
write.table(mic_asv_osc4_rel_wide,
            file=paste0(output_tables,"FRAM_RAS_F4_MIC_ASV_OSC4_rel_abund_wide.txt"),
            sep="\t", quote = F, row.names = F)
write.table(euk_asv_osc4_rel_wide,
            file=paste0(output_tables,"FRAM_RAS_F4_EUK_ASV_OSC4_rel_abund_wide.txt"),
            sep="\t", quote = F, row.names = F)
write.table(gene_clust_osc4_rel_wide,
            file=paste0(output_tables,"FRAM_RAS_F4_GENE_CLUST_OSC4_rel_abund_wide.txt"),
            sep="\t", quote = F, row.names = F)

#####

### Group gene clusters with single annual oscillations into functional clusters

#####

### Import gene cluster to function information
gene_to_func_to_tax=read.table(
  file="FRAM_RAS_F4_GENE_CLUST_to_FUNC_to_taxonomy.txt",
  sep="\t",check.names=F, header=T)

# Filter the gene cluster to functional information to retain only gene clusters
# with a single oscillation per year, then combine with relative abundance matrix and
# sum the abundances of gene clusters for each function in each sample. This will
# generate an abundance profile for the functional clusters
func_osc4_rel_wide = gene_clust_osc4_rel_wide %>%
  filter(ID %in% gene_to_func_to_tax$GeneClustID) %>%
  left_join(gene_to_func_to_tax, by=c("ID" = "GeneClustID")) %>%
  select(-ID,-GTDB_taxonomy,-Domain_tiara,-FUNC,-Gene) %>%
  reshape2::melt(., id.vars=c("FuncClustID"), variable.name="RAS_id",
                 value.name="Rel_abund") %>%
  aggregate(Rel_abund~FuncClustID+RAS_id, data=., FUN=sum) %>%
  reshape2::dcast(FuncClustID~RAS_id, value.var="Rel_abund", data=.)

# The functional profile will now be subject to fourier transformation and then
# combined with the prokaryotic ASV oscillations in a network analysis
write.table(func_osc4_rel_wide,
            file=paste0(output_tables,"FRAM_RAS_F4_FUNC_osc4_rel_wide.txt"),
            sep="\t", quote = F, row.names = F)


### Now apply the same function used to interpolate the abundances of ASVs and gene clusters to the 
### functional relative abundance profile of functions, as this will be used in the network analysis

### Function to generate interpolated relative abundance profiles
interpolate_abundances <- function(input_file, ras_id_file, output_file) {
  
  # Load the input files
  rel_profile <- read_delim(input_file, delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  # Erste Spalte als rownames setzen (und aus den Daten entfernen, falls nötig)
  rn = rel_profile[[1]]
  rel_profile[,1] <- NULL
  rownames(rel_profile) <- rn # Zeilennamen zuweisen
  metadata <- read.delim2(ras_id_file)
  
  # Extract ASV names and adjust columns
  ID <- rownames(rel_profile)
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
  
  # User feedrel_profilek
  cat("File processing completed. Results saved in:", output_file, "\n")
}

interpolate_abundances(paste0(output_tables,"FRAM_RAS_F4_FUNC_osc4_rel_wide.txt"), "FRAM_RAS_F4_META.txt", paste0(output_tables,"FRAM_RAS_F4_FUNC_osc4_interpolated_abundances.txt"))
