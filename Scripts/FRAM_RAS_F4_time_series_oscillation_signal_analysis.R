####################

#### Processing oscillation signal data

####################

### Define working directory
setwd('XXXXXX')

### Define output directories
output_figures = ('results/figures/')
output_tables = ('results/tables/')

### Load libraries
library(stringr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(tibble)
library(patchwork)
library(data.table)

### Import oscillation tables
gene_clust_oscillations = read.table(paste0(output_tables,"FRAM_RAS_F4_GENE_CLUST_oscillations_per_year.txt"), sep="\t", header=T)
euk_asv_oscillations = read.table(paste0(output_tables,"FRAM_RAS_F4_EUK_oscillations_per_year.txt"), sep="\t", header=T)
prok_asv_oscillations = read.table(paste0(output_tables,"FRAM_RAS_F4_MIC_oscillations_per_year.txt"), sep="\t", header=T)

### Import relative abudnance tables
mic_asv_rel=read.table(file=paste0(output_tables,"FRAM_RAS_F4_MIC_ASV_filt_rel.txt"), sep="\t",
                       check.names=F, header=T, row.names=1)
euk_asv_rel=read.table(file=paste0(output_tables,"FRAM_RAS_F4_EUK_ASV_filt_rel.txt"), sep="\t",
                       check.names=F, header=T, row.names=1)
gene_clust_rel=read.table(
  file=paste0(output_tables,"FRAM_RAS_F4_GENE_CLUST_filt_rel_wide.txt"), 
  sep="\t", check.names=F, header=T, row.names=1)

### Combine oscillation from asvs and gene clusters together
asv_and_geneclust_oscillations = rbind(gene_clust_oscillations,euk_asv_oscillations,prok_asv_oscillations) %>%
    mutate(Type = case_when(
        grepl("prok",nodeId) ~ "Prok",
        grepl("euk",nodeId) ~ "Euk",
        grepl("gene",nodeId) ~ "Gene_clust"
    ))

### Assess the oscillation signal data retrieved from the time-series analysis
# Calculate the mean oscillations per year
asv_and_geneclust_oscillations_mean = 
    asv_and_geneclust_oscillations %>%
    aggregate(oscillations~nodeId+Type, data=., FUN=mean) %>%
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
    tibble::rownames_to_column(., var="nodeId") %>%
    filter(nodeId %in% asv_and_geneclust_oscillations_mean$nodeId) %>%
    reshape2::melt(id.vars="nodeId", variable.name="RAS_id", value.name="Rel_abund") %>%
    filter(., Rel_abund > 0) %>%
    left_join(asv_and_geneclust_oscillations_mean, by="nodeId") %>%
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
  select(nodeId)

### Create relative abundance dataframes for those with one oscillation per year as these will be 
### used as the input for network analysis
mic_asv_osc4_rel_wide = 
  mic_asv_rel %>%
  tibble::rownames_to_column(., var="nodeId") %>%
  filter(nodeId %in% mic_euk_asv_gene_clust_osc4_id_list$nodeId)

euk_asv_osc4_rel_wide = 
  euk_asv_rel %>%
  tibble::rownames_to_column(., var="nodeId") %>%
  filter(nodeId %in% mic_euk_asv_gene_clust_osc4_id_list$nodeId)

gene_clust_osc4_rel_wide = 
  gene_clust_rel %>%
  tibble::rownames_to_column(., var="nodeId") %>%
  filter(nodeId %in% mic_euk_asv_gene_clust_osc4_id_list$nodeId)

# Export tables
write.table(mic_euk_asv_gene_clust_osc4_id_list,
            file=paste0(output_tables,"FRAM_RAS_F4_MIC_EUK_ASVs_GENE_CLUST_OSC4_ID_list.txt"),
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
  filter(nodeId %in% gene_to_func_to_tax$Gene_clustID) %>%
  left_join(gene_to_func_to_tax, by=c("nodeId" = "Gene_clustID")) %>%
  select(-nodeId,-GTDB_taxonomy,-Domain_tiara,-FUNC,-Gene) %>%
  reshape2::melt(., id.vars=c("funcID"), variable.name="RAS_id",
                 value.name="Rel_abund") %>%
  aggregate(Rel_abund~funcID+RAS_id, data=., FUN=sum) %>%
  reshape2::dcast(funcID~RAS_id, value.var="Rel_abund", data=.)

# The functional profile will now be subject to fourier transformation and then
# combined with the prokaryotic ASV oscillations in a network analysis
write.table(func_osc4_rel_wide,
            file=paste0(output_tables,"FRAM_RAS_F4_FUNC_osc4_rel_wide.txt"),
            sep="\t", quote = F, row.names = F)


#####

### Determine the annual peaks of oscillation signals for those with an oscillation of once per year

#####

            
#####

### Assessing the timing of peaks for annually oscillating ASVs and genes

#####

library(lubridate)

# Import table containing information about the sampling time point when
# each annually oscillating ASV and gene cluster reached peak abundance each year
mic_asv_gene_osc4_peaks_per_year_julian=read.table(file="RAS_F4_MIC_ASV_GENES_OSC4_peak_timings_julian.txt", 
                                  sep="\t", check.names=F, header=T) 

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
    ID = grepl("prok_asv",ID) ~ "PROK_ASV",
    ID = grepl("euk_asv",ID) ~ "EUK_ASV",
    ID = grepl("Cluster",ID) ~ "gene",
  )) %>%
  select(type) %>%
  table()
