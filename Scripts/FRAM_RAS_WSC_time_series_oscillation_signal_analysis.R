####################

#### Processing time-series analysis output

####################

### Define working directory
setwd('XXXXX')

### Define output directories
output_figures <- ('')
output_tables <- ('')

### Load libraries
library(stringr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(tibble)
### Import files
options(scipen=999)

# Import microbial ASV oscillation signal patterns
mic_asv_oscillation_py=read.table(file="RAS_F4_MIC_ASV_num_osc_per_year.csv", 
                                  sep=",", check.names=F, header=T) %>%
  mutate(type = "MIC_ASV") %>%
  dplyr::rename(., ID = asv)

# Import gene clusters (95% ANI clusters) oscillation signal patterns
mic_clust_oscillation_py=read.table(file="FRAM_RAS_F4_MIC_GENE_CLUST_num_osc_per_year.csv", 
                                   sep=",", check.names=F, row.names=NULL, header=T) %>%
  dplyr::select(asv,year,oscillations) %>%
  mutate(type = "MIC_GENE_CLUSTERS") %>%
  dplyr::rename(., ID = asv)

# Import eukaryotic ASV oscillation signal patterns
euk_asv_oscillation_py=read.table(file="RAS_F4_EUK_ASV_num_osc_per_year.csv", 
                                  sep=",", check.names=F, header=T) %>%
  mutate(type = "EUK_ASV") %>%
  dplyr::rename(., ID = asv)

# Import microbial ASV relative abundance data
mic_asv_rel=read.table(file="RAS_F4_MIC_ASV_filt_rare_rel.txt", sep="\t",
                       check.names=F, header=T, row.names=1)

# Import EUK ASV relative abundance data
euk_asv_rel=read.table(file="RAS_F4_EUK_ASV_filt_rare_rel.txt", sep="\t",
                       check.names=F, header=T, row.names=1)

# Import gene clusters relative abundance data
mic_clust_rel=read.table(
  file="FRAM_RAS_F4_GENE_CLUSTID_filt_rel_wide.txt", 
  sep="\t", check.names=F, header=T, row.names=NULL) %>%
  tibble::column_to_rownames(., var="clustID") %>%
  dplyr::select(-row.names)

#####

### Process and filter the time-series output

#####

### Assess the oscillation signal data retrieved from the time-series analysis
# How many microbial ASVs could an oscillation be determined?
mic_asv_oscillation_py %>%
  select(ID) %>%
  unique() %>%
  dim()
dim(mic_asv_rel)
(933/3748)*100

# How many exhibited one oscillation per annual cycle?
mic_asv_oscillation_py %>%
  filter(., oscillations == 1) %>%
  select(ID) %>%
  as.data.frame() %>%
  filter(., grepl("bac",ID)) %>%
  table() %>%
  as.data.frame() %>%
  filter(., Freq == 4) %>%
  dim()
(658/933)*100

# Repeat above steps for EUK
euk_asv_oscillation_py %>%
  select(ID) %>%
  unique() %>%
  dim()
dim(euk_asv_rel)
(1019/3019)*100

euk_asv_oscillation_py %>%
  filter(., oscillations == 1) %>%
  select(ID) %>%
  as.data.frame() %>%
  filter(., grepl("euk",ID)) %>%
  table() %>%
  as.data.frame() %>%
  filter(., Freq == 4) %>%
  dim()
(448/1019)*100

# Repeat above steps for gene clusters
mic_clust_oscillation_py %>%
  select(ID) %>%
  unique() %>%
  dim()

mic_clust_oscillation_py %>%
  filter(., oscillations == 1) %>%
  select(ID) %>%
  as.data.frame() %>%
  table() %>%
  as.data.frame() %>%
  filter(., Freq == 4) %>%
  dim()
(482923/704158)*100

# Combine ASV and GENE oscillation dataframes and calculate mean number of oscillations 
# per year
mic_euk_asv_gene_clust_mean_osc_per_year = rbind(mic_asv_oscillation_py,
                                       euk_asv_oscillation_py,
                                       mic_clust_oscillation_py) %>%
  aggregate(oscillations~ID+type, data=., FUN=mean) %>%
  mutate(oscillations = as.numeric(oscillations))

# How many ASVs and gene clusters show oscillation of 1 per year?
mic_euk_asv_gene_clust_mean_osc_per_year %>%
  filter(oscillations == 1) %>%
  mutate(count = 1) %>%
  aggregate(count~type+oscillations, data=., FUN=sum)

### Determine fraction of community represented by each oscillation type

### FUNCTION: calculate relative abund in each sample for each oscillation pattern
### for microbial and microeukaryotic ASVs and gene clusters
calc_rel_abund_per_osc <- function(infile){
  osc_abund_df = 
    infile %>% 
    tibble::rownames_to_column(., var="ID") %>%
    filter(ID %in% mic_euk_asv_gene_clust_mean_osc_per_year$ID) %>%
    reshape2::melt(id.vars="ID", variable.name="RAS_id", value.name="Rel_abund") %>%
    left_join(mic_euk_asv_gene_clust_mean_osc_per_year, by="ID") %>%
    mutate(Rel_abund = as.numeric(Rel_abund)) %>%
    aggregate(Rel_abund~RAS_id+type+oscillations, data=., FUN=sum) %>%
    mutate(oscillations = round(as.numeric(oscillations),3)) 
  return(osc_abund_df)
}

mic_asv_osc_per_year_rel_prop <- calc_rel_abund_per_osc(mic_asv_rel)
euk_asv_osc_per_year_rel_prop <- calc_rel_abund_per_osc(euk_asv_rel)
mic_clust_osc_per_year_rel_prop <- calc_rel_abund_per_osc(mic_clust_rel)

### Combine dataframes 
asv_gene_rel_prop_per_oscillation <- 
  rbind(mic_asv_osc_per_year_rel_prop,
        euk_asv_osc_per_year_rel_prop,
        mic_clust_osc_per_year_rel_prop) %>%
  mutate(Rel_abund = Rel_abund*100)

# Calculate the combined, mean relative abundance of asvs and genes
# with different oscillation signals
asv_gene_rel_prop_per_oscillation %>%
  group_by(type,oscillations) %>%
  summarise(mean_per_osc = mean(Rel_abund)) %>%
  as.data.frame()

### Visualise relative abundance of ASVs and GENEs with different 
### oscillation signals 
plot_rel_prop_per_oscillation = function(infile){
  ggplot(data=infile,
         aes(y = oscillations, x = Rel_abund)) +
    geom_boxplot(aes(group=oscillations)) + 
    facet_wrap(type~., scales="free_x", nrow = 4,
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
    facet_wrap(type~., scales = "free_x", nrow = 4) + 
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
plot_asv_gene_mean_osc_per_year <- plot_num_oscillations_py(mic_euk_asv_gene_clust_mean_osc_per_year)
plot_asv_gene_rel_prop_per_osc <- plot_rel_prop_per_oscillation(asv_gene_rel_prop_per_oscillation)

# Export figure
library(patchwork)
pdf(file="figures_output/FRAM_RAS_F4_MIC_EUK_ASV_GENE_OSC_per_year_and_rel_abund.pdf",
    height=10, width=8)
plot_asv_gene_mean_osc_per_year|plot_asv_gene_rel_prop_per_osc
dev.off()

#####

### Analyse the ASVs and GENEs with one oscillation per year

#####

### FUNCTION - Identify those with one oscillation each year
create_osc4_id_list <- function(infile){
  infile %>%
  reshape2::dcast(ID~year, value.var="oscillations", data=.,) %>%
    filter(ID == case_when(
      year1 == 1 &
        year2 == 1 &
        year3 == 1 &
        year4 == 1 ~ ID
    )) %>%
    dplyr::select(ID)
}

# Run function and combine ASV and GENE output dataframes
mic_asv_osc4_id_list <- create_osc4_id_list(mic_asv_oscillation_py)
euk_asv_osc4_id_list <- create_osc4_id_list(euk_asv_oscillation_py)
mic_clust_osc4_id_list <- create_osc4_id_list(mic_clust_oscillation_py)
mic_euk_asv_gene_clust_osc4_id_list <- rbind(mic_asv_osc4_id_list,
                                       euk_asv_osc4_id_list,
                                       mic_clust_osc4_id_list)

### Create relative abundance (ASV) and normalised count (GENEs)
### dataframes for those with one oscillation per year as these will be 
### used as the input for network analysis
### The EUK ASVs will not be used for network analysis, but we will still create
### the relevant dataframes for other comparative analysis
mic_asv_osc4_rel_wide = 
  mic_asv_rel %>%
  filter(row.names(.) %in% mic_asv_osc4_id_list$ID)

euk_asv_osc4_rel_wide = 
  euk_asv_rel %>%
  filter(row.names(.) %in% euk_asv_osc4_id_list$ID)

mic_gene_clust_osc4_rel_wide = 
  mic_clust_rel %>%
  filter(row.names(.) %in% mic_clust_osc4_id_list$ID)

# Export tables
write.table(mic_asv_osc4_id_list,
            file="time_series_analysis/FRAM_RAS_F4_MIC_ASVs_OSC4_ID_list.txt",
            sep="\t")
write.table(mic_euk_asv_gene_clust_osc4_id_list,
            file="time_series_analysis/FRAM_RAS_F4_MIC_EUK_ASVs_GENE_CLUST_OSC4_ID_list.txt",
            sep="\t")
write.table(mic_asv_osc4_rel_wide,
            file="time_series_analysis/FRAM_RAS_F4_MIC_ASV_OSC4_rel_abund_wide.txt",
            sep="\t")
write.table(euk_asv_osc4_rel_wide,
            file="time_series_analysis/FRAM_RAS_F4_EUK_ASV_OSC4_rel_abund_wide.txt",
            sep="\t")
write.table(mic_gene_clust_osc4_rel_wide,
            file="time_series_analysis/FRAM_RAS_F4_MIC_GENE_CLUST_OSC4_rel_abund_wide.txt",
            sep="\t")

#####

### Assessing the timing of peaks for annually oscillating ASVs and genes

#####

library(lubridate)

# Import table containing information about the sampling time point when
# each annually oscillating ASV and gene cluster reached peak abundance each year
mic_asv_gene_osc4_peaks_per_year_julian=read.table(file="time_series_analysis/RAS_F4_MIC_ASV_GENES_OSC4_peak_timings_julian.csv", 
                                  sep=";", check.names=F, header=T) 

# Using the information, convert the peak sample dates to julian day
# and compare the all years to year 1 to identify those that oscillated
# within a 30 day window consistently across years
mic_asv_gene_osc_within_30_days = 
  mic_asv_gene_osc4_peaks_per_year_julian %>%
  mutate(date1_julian = lubridate::yday(date1)) %>%
  mutate(date2_julian = lubridate::yday(date2)) %>%
  mutate(date3_julian = lubridate::yday(date3)) %>%
  select(ID,date1_julian,date2_julian,date3_julian) %>%
  mutate(
    peak_timings = case_when(
      ((date2_julian <= date1_julian+15 & date2_julian >= date1_julian-15) &
        (date3_julian <= date1_julian+15 & date3_julian >= date1_julian-15)) |
        ((date2_julian <= date1_julian+30 & date2_julian >= date1_julian) &
           (date3_julian <= date1_julian+30 & date3_julian >= date1_julian)) | 
        ((date2_julian <= date1_julian & date2_julian >= date1_julian-30) &
           (date3_julian <= date1_julian & date3_julian >= date1_julian-30)) ~ "recurred_in_window",
      TRUE ~ "No"
    )) %>%
  filter(peak_timings == "recurred_in_window")

### How many of the Microbial ASVs and gene clusters recurred within a 30day period
### each year?
mic_asv_gene_osc_within_30_days %>%
  mutate(type = case_when(
    ID = grepl("bac_asv",ID) ~ "BAC_ASV",
    ID = grepl("euk_asv",ID) ~ "EUK_ASV",
    ID = grepl("Cluster",ID) ~ "gene",
  )) %>%
  select(type) %>%
  table()

### What proportion of annually oscillating components is this?
mic_euk_asv_gene_clust_osc4_id_list %>%
  mutate(type = case_when(
    ID = grepl("bac_asv",ID) ~ "BAC_ASV",
    ID = grepl("euk_asv",ID) ~ "EUK_ASV",
    ID = grepl("Cluster",ID) ~ "gene",
  )) %>%
  select(type) %>%
  table()

### There were 658 microbial ASVs and 485075 microbial gene clusters with 
### annual oscillations
(336/658)*100
(110208/485075)*100

####################

#####
### Process microbial gene time-series data ready for network input
### A correlation network will be constructed based on fourier transformed
### information. We are specifically interested in ASVs/genes
### with a single oscillation per year. 
### However, because there are too many genes to include in a network,
### the orthologous genes will be grouped based on functional annotations
### We have used the EGGNOG seed ortholog annotations. Where KEGG was available
### this was the priority, followed by PFAM.
#####

####################

### Import files

# Import EGGNOG ID to geneID to function information
clust_id_func_id_to_func=read.csv(
  file="metagenomes/community_gene_profiles/FRAM_RAS_F4_GENES_CLUSTID_to_FUNCID_to_FUNC.txt",
  sep="\t",check.names=F, header=T)

# How many gene clusters were assigned a function based on eggnog database?
dim(clust_id_func_id_to_func)

# Move rownames to a column in relative abundance table
mic_gene_clust_osc4_rel_wide_mod = 
  mic_gene_clust_osc4_rel_wide %>%
  tibble::rownames_to_column(., var="clustID")

# Filter the gene cluster to functional information to retain only gene clusters
# with a single oscillation per year, then combine with relative abundance matrix and
# sum the abundances of gene clusters for each function in each sample. This will
# generate a functional profile

mic_clust_func_rel_abund_wide = 
  clust_id_func_id_to_func %>%
  filter(clustID %in% mic_gene_clust_osc4_rel_wide_mod$clustID) %>%
  left_join(mic_gene_clust_osc4_rel_wide_mod, by="clustID") %>%
  reshape2::melt(., id.vars=c("clustID","FUNC","funcID"), variable.name="RAS_id",
                 value.name="Rel_abund") %>%
  aggregate(Rel_abund~funcID+RAS_id, data=., FUN= sum) %>%
  reshape2::dcast(funcID~RAS_id, value.var="Rel_abund") %>%
  tibble::column_to_rownames(., var="funcID")

# The functional profile will now be subject to fourier transformation and then
# combined with the prokaryotic ASV oscillations in a network analysis
write.table(mic_clust_func_rel_abund_wide,
            file="time_series_analysis/FRAM_RAS_F4_MIC_GENE_OSC4_FUNC_rel_abund_wide.txt",
            sep="\t")

### What is the average relative abundance of functional gene groups across samples
mic_clust_func_rel_abund_wide %>%
  tibble::rownames_to_column(., var="funcID") %>%
  melt(id.vars="funcID", variable.name="RAS_id", value.name="Rel_abund") %>%
  aggregate(Rel_abund~RAS_id, data=., FUN=sum) %>%
  arrange(desc(Rel_abund))
