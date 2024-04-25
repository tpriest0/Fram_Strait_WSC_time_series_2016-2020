####################

#### Investigating richness, evenness and compositional differences in
#### samples

####################

### Define working directory
setwd('') # Set the working directory as the directory where the data_files/ of this repository were downloaded to

### Load libraries
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(vegan)
library(Hmisc)
library(iNEXT)
library(forcats)
library(confintr)
library(splancs)

options(scipen=999)

### Import files
# Import eukaryotic filtered ASV raw and rel abund data
euk_asv_filt_raw=read.table(file="ASV/RAS_F4_EUK_ASV_filt_rare_raw.txt", sep="\t",
                            check.names=F, header=T, row.names=1)

euk_asv_filt_rel=read.table(file="ASV/RAS_F4_EUK_ASV_filt_rare_rel.txt", sep="\t",
                            check.names=F, header=T, row.names=1) %>%
  mutate(across(is.numeric, ~ .* 100))

# Import microbial filtered ASV raw and rel abund data
mic_asv_filt_raw=read.table(file="ASV/RAS_F4_MIC_ASV_filt_rare_raw.txt", sep="\t",
                            check.names=F, header=T, row.names=1)

mic_asv_filt_rel=read.table(file="ASV/RAS_F4_MIC_ASV_filt_rare_rel.txt", sep="\t",
                            check.names=F, header=T, row.names=1)  %>%
  mutate(across(is.numeric, ~ .* 100))

mic_asv_filt_taxa=read.table(file="ASV/RAS_F4_MIC_ASV_filt_taxa.txt", sep="\t",
                            check.names=F, header=T, row.names=1)  %>%
  mutate(across(is.numeric, ~ .* 100))

# Import sample metadata
sample_meta=read.table(file="RAS_F4_META_MIC_Amplicon.txt", sep="\t",
                       check.names=F, header=T) %>%
  mutate(date = as.Date(date, "%d.%m.%Y"))

#####

### How many bacterial and archaeal ASVs are there?

mic_asv_filt_taxa %>%
  filter(., Kingdom == "Bacteria") %>%
  dim()

mic_asv_filt_taxa %>%
  filter(., Kingdom == "Archaea") %>%
  dim()

#####

### Estimate species richness and shannon diversity after rarefying

#####

### Calculate alpha diversity metrics
calculate_alpha_div <- function(abund_table, sample_meta){
  R = specnumber(t(abund_table))
  H = diversity(t(abund_table), "shannon")
  E = H/log(R)
  
  alpha_div_df = 
    data.frame(cbind(Richness = R, Shannon_diversity = H, Evenness = E)) %>%
    tibble::rownames_to_column(., var = "RAS_id") %>%
    left_join(sample_meta, by = "RAS_id") %>%
    mutate(month = factor(month, levels = c("Jan","Feb","Mar","Apr",
                                        "May","Jun","Jul","Aug",
                                        "Sep","Oct","Nov","Dec")))
  return(alpha_div_df)
}

mic_alpha_with_meta = calculate_alpha_div(mic_asv_filt_raw, sample_meta)
euk_alpha_with_meta = calculate_alpha_div(euk_asv_filt_raw, sample_meta)

### FUNCTION: Visualise diversity over time-series
plot_diversity_metrics <- function(infile,metric,yaxislabel){
  infile %>%
    select(date,{{metric}}) %>%
    mutate(date = as.Date(date, "%d.%m.%Y")) %>%
    ggplot(., aes(x = as.Date(date), y={{metric}})) + 
    geom_line() + 
    scale_x_date(date_breaks = "6 months", date_labels =  "%b %Y") +
    labs(y = yaxislabel,
         x = "Date") + 
    theme_bw() + 
    theme(axis.text.y = element_text(size = 10),
          axis.title.y = element_text(size = 12),
          axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
          axis.title.x = element_text(size = 12),
          title = element_text(size = 10))
}

# specify code to remove x-axis labels for plots (prevents repeated axes when stacking plots)
remove_x_axis <- theme(
  axis.text.x = element_blank(),
  axis.title.x = element_blank()
)

euk_richness = plot_diversity_metrics(euk_alpha_with_meta, Richness, "Richness")
euk_shannon = plot_diversity_metrics(euk_alpha_with_meta, Shannon_diversity, "Shannon diversity")
euk_evenness = plot_diversity_metrics(euk_alpha_with_meta, Evenness, "Evenness")

mic_richness = plot_diversity_metrics(mic_alpha_with_meta, Richness, "Richness")
mic_shannon = plot_diversity_metrics(mic_alpha_with_meta, Shannon_diversity, "Shannon diversity")
mic_evenness = plot_diversity_metrics(mic_alpha_with_meta, Evenness, "Evenness")

#####

### Compare diversity metrics to environmental conditions

# First let's filter our metadata table to continuous variables and zero-mean variance
# standardize them
sample_meta_stand = 
  sample_meta %>%
  select(RAS_id,MLD,temp,Chl_a,sal,O2_conc,CO2,pH,PAR_satellite) %>%
  tibble::column_to_rownames(., var = "RAS_id") %>%
  decostand(., method="standardize", MARGIN=2, na.rm=T) %>%
  tibble::rownames_to_column(., var = "RAS_id")

### FUNCTION: Calculate pearson's correlation between diversity metrics and env variables
diversity_vs_env_correlation <- function(infile, sample_meta_stand){
  infile %>%
    select(RAS_id,Richness,Shannon_diversity,Evenness) %>%
    left_join(sample_meta_stand, by = "RAS_id") %>% 
    tibble::column_to_rownames(., var="RAS_id") %>%
    as.matrix() %>%
    rcorr(.)
}

### FUNCTION: Process correlation output into dataframe
process_correlation_matrix <- function(infile_cor, infile_p) {
  x <- upper.tri(infile_cor)
  y = data.frame(
    row = rownames(infile_cor)[row(infile_cor)[x]],
    column = rownames(infile_cor)[col(infile_cor)[x]],
    cor  =(infile_cor)[x],
    p = infile_p[x]
  )
  return(y)
}

### FUNCTION: multiple testing correction and significance filtering
correct_and_filter_cor_results <- function(cor_df){
  adjusted_p = p.adjust(cor_df$p, method = "bonferroni")
  cor_df %>%
    mutate(p_adj = adjusted_p) %>% 
    filter(p_adj < 0.05) %>%
    filter(row == "Shannon_diversity" |
             row == "Evenness" | 
             row == "Richness")
}

# Run above functions to create a dataframe containing correlation coefficients
# associated with adjusted p-values < 0.05
mic_alpha_vs_env_cor_results = diversity_vs_env_correlation(mic_alpha_with_meta, sample_meta_stand)
mic_alpha_vs_env_cor_df = process_correlation_matrix(mic_alpha_vs_env_cor_results$r,mic_alpha_vs_env_cor_results$P)
mic_alpha_vs_env_cor_sig_df = correct_and_filter_cor_results(mic_alpha_vs_env_cor_df)

euk_alpha_vs_env_cor_results = diversity_vs_env_correlation(euk_alpha_with_meta, sample_meta_stand)
euk_alpha_vs_env_cor_df = process_correlation_matrix(euk_alpha_vs_env_cor_results$r,euk_alpha_vs_env_cor_results$P)
euk_alpha_vs_env_cor_sig_df = correct_and_filter_cor_results(euk_alpha_vs_env_cor_df)

### Now we want to calculate and add in the confidence intervals for our
### coefficients with a p-value < 0.05

### FUNCTION: calculate confidence interval for correlation coefficients
calculate_ci <- function(correlation, n) {
  z <- 0.5 * log((1 + correlation) / (1 - correlation))
  se <- 1 / sqrt(n - 3)
  alpha <- 0.05  # 95% confidence interval
  z_critical <- qnorm(1 - alpha / 2)
  ci <- tanh(c(z - z_critical * se, z + z_critical * se))
  return(ci)
}

# Add in columns to store the lower and upper bound confidence intervals
ci_lower <- numeric(nrow(mic_alpha_vs_env_cor_sig_df))
ci_upper <- numeric(nrow(mic_alpha_vs_env_cor_sig_df))

# Loop through significant correlation dataframe and calculate confidence intervals
for (i in 1:nrow(mic_alpha_vs_env_cor_sig_df)) {
  correlation <- mic_alpha_vs_env_cor_sig_df$cor[i]
  n <- 97
  ci <- calculate_ci(correlation, n)
  ci_lower[i] <- ci[1]
  ci_upper[i] <- ci[2]
}

# Add the confidence intervals to the respective columns
mic_alpha_vs_env_cor_sig_df$ci_lower <- ci_lower
mic_alpha_vs_env_cor_sig_df$ci_upper <- ci_upper

### Repeat above steps for microeukaryotic correlations
ci_lower <- numeric(nrow(euk_alpha_vs_env_cor_sig_df))
ci_upper <- numeric(nrow(euk_alpha_vs_env_cor_sig_df))

# Loop through significant correlation dataframe and calculate confidence intervals
for (i in 1:nrow(euk_alpha_vs_env_cor_sig_df)) {
  correlation <- euk_alpha_vs_env_cor_sig_df$cor[i]
  n <- 96
  ci <- calculate_ci(correlation, n)
  ci_lower[i] <- ci[1]
  ci_upper[i] <- ci[2]
}

# Add the confidence intervals to the respective columns
euk_alpha_vs_env_cor_sig_df$ci_lower <- ci_lower
euk_alpha_vs_env_cor_sig_df$ci_upper <- ci_upper

### Export the dataframes containing information on significant correlation results
write.table(mic_alpha_vs_env_cor_sig_df,
            file = "figures_output/FRAM_RAS_F4_MIC_div_vs_env_correlation_df.txt", 
            sep = "\t")

write.table(euk_alpha_vs_env_cor_sig_df,
            file = "figures_output/FRAM_RAS_F4_EUK_div_vs_env_correlation_df.txt", 
            sep = "\t")

###
### FUNCTION: visualise significant correlations
plot_sig_cor_hits <- function(infile){
  ggplot(infile, aes(y=row, x=column)) +
    geom_tile(aes(fill=cor)) + 
    scale_fill_gradient2(limits=c(-1,1),
                         low = "#40004B", high = "#00441B", mid = "#ffffff") + 
    labs(fill = "Pearson's correlation") + 
    theme_bw() + 
    theme(axis.title.y = element_blank(),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size = 12, angle=45, hjust=1),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 11))
}

### Visualise correlation results
euk_div_vs_env_sig_cor_plot = plot_sig_cor_hits(euk_alpha_vs_env_cor_sig_df)
mic_div_vs_env_sig_cor_plot = plot_sig_cor_hits(mic_alpha_vs_env_cor_sig_df)

#####

### Combine all diversity and correlation plots into single figure

#####

### Combine plots and export
library(patchwork)
pdf(file="figures_output/RAS_F4_MIC_EUK_ASV_diversity_and_correlations.pdf",
    height=10, width = 14)
((mic_richness+remove_x_axis)/(mic_evenness+remove_x_axis)/mic_shannon/mic_div_vs_env_sig_cor_plot)|
  ((euk_richness+remove_x_axis)/(euk_evenness+remove_x_axis)/euk_shannon/euk_div_vs_env_sig_cor_plot)
dev.off()

#######################

##### Beta-diversity analysis

### Hellinger transform counts and calculate bray-curtis dissimilarity
mic_asv_hel_bray = mic_asv_filt_raw %>%
  decostand(., method="hellinger", MARGIN=2) %>%
  t() %>%
  vegdist(., method="bray")

euk_asv_hel_bray = euk_asv_filt_raw %>%
  decostand(., method="hellinger", MARGIN=2) %>%
  t() %>%
  vegdist(., method="bray")

### NMDS

## Prokaryotes
mic_asv_hel_bray_nmds = metaMDS(mic_asv_hel_bray, distance="bray", k=2, pc=T, autotransform=T, stepacross=T, center=T)

mic_asv_bray_nmds_with_meta <- 
  as.data.frame(mic_asv_hel_bray_nmds$points) %>%
  tibble::rownames_to_column(., var="RAS_id") %>%
  left_join(sample_meta, by="RAS_id") %>%
  mutate(month = factor(month, levels=c("Jan","Feb","Mar","Apr","May",
                                        "Jun","Jul","Aug","Sep","Oct",
                                        "Nov","Dec")))

# Calculate convex hulls for each month
mic_nmds_month_hulls <- mic_asv_bray_nmds_with_meta %>%
  group_by(month) %>%
  slice(chull(MDS1,MDS2))

month_colours <- colorRampPalette(brewer.pal(11, "Spectral"))(12)

# Plot NMDS with convex hulls and points coloured by month
mic_asv_bray_nmds_plot_with_hulls = 
  ggplot(data=mic_asv_bray_nmds_with_meta, aes(x = MDS1, y=MDS2)) +
  geom_point(aes(colour = month), position="identity") + 
  geom_polygon(data=mic_nmds_month_hulls, aes(fill = month, colour = month), alpha=0.2) + 
  scale_color_manual(values = month_colours) + 
  scale_fill_manual(values = month_colours) + 
  labs(y = "NMDS2", x = "NMDS1", title = round(mic_asv_hel_bray_nmds$stress,3)) + 
  theme_bw() + 
  theme(legend.text = element_text(size = 10),
        legend.title = element_text(size =12),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size =12),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size =12))

## Microeukaryotes
euk_asv_hel_bray_nmds = metaMDS(euk_asv_hel_bray, distance="bray", k=2, pc=T, autotransform=T, stepacross=T, center=T)

euk_asv_bray_nmds_with_meta <- 
  as.data.frame(euk_asv_hel_bray_nmds$points) %>%
  tibble::rownames_to_column(., var="RAS_id") %>%
  left_join(sample_meta, by="RAS_id") %>%
  mutate(month = factor(month, levels=c("Jan","Feb","Mar","Apr","May",
                                        "Jun","Jul","Aug","Sep","Oct",
                                        "Nov","Dec")))

# Calculate convex hulls for each month
euk_nmds_month_hulls <- euk_asv_bray_nmds_with_meta %>%
  group_by(month) %>%
  slice(chull(MDS1,MDS2))

# Plot NMDS with convex hulls and points coloured by month
euk_asv_bray_nmds_plot_with_hulls = 
  ggplot(data=euk_asv_bray_nmds_with_meta, aes(x = MDS1, y=MDS2)) +
  geom_point(aes(colour = month), position="identity") + 
  geom_polygon(data=euk_nmds_month_hulls, aes(fill = month, colour = month), alpha=0.2) + 
  scale_color_manual(values = month_colours) + 
  scale_fill_manual(values = month_colours) + 
  labs(y = "NMDS2", x = "NMDS1", title = round(euk_asv_hel_bray_nmds$stress,3)) + 
  theme_bw() + 
  theme(legend.text = element_text(size = 10),
        legend.title = element_text(size =12),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size =12),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size =12))

### combine plots and export
pdf(file="figures_output/FRAM_RAS_F4_MIC_and_EUK_NMDS_plots.pdf",
    height=12, width = 8)
mic_asv_bray_nmds_plot_with_hulls/euk_asv_bray_nmds_plot_with_hulls
dev.off()

### Calculate convex hull areas

euk_areas <- c()
prok_areas <- c()
months <- unique(sample_meta$month)

# Iterate over unique months
for (selectmonth in months) {
  
  euk_area <- euk_nmds_month_hulls %>%
    as.data.frame() %>%
    filter(., grepl(selectmonth, month)) %>%
    select(MDS1, MDS2) %>%
    as.matrix() %>%
    areapl(.)
  
  prok_area <- mic_nmds_month_hulls %>%
    as.data.frame() %>%
    filter(., grepl(selectmonth, month)) %>%
    select(MDS1, MDS2) %>%
    as.matrix() %>%
    areapl(.)
  <
  # Store the results in respective vectors
  euk_areas <- c(euk_areas, euk_area)
  prok_areas <- c(prok_areas, prok_area)
}

# Combine results into dataframe
mic_and_euk_convex_hull_areas <- data.frame(month = months, microeukaryotic_hull_area = euk_areas, prokaryotic_area = prok_areas)

# Export
write.table(mic_and_euk_convex_hull_areas, 
            file = "figures_output/FRAM_RAS_F4_MIC_EUK_NMDS_convex_hull_areas.txt",sep="\t")
