
####################

#### Processing network clusters and visualise cluster dynamics over time
#### series

####################

### Define working directory
setwd('D:/Own_cloud_updated/PhD/Arctic Bacterioplankton/Fram Strait/RAS_metagenomes/WSC')


### Load libraries
library(dplyr)
library(reshape2)
library(tibble)
library(ggplot2)
library(patchwork)



mic_net_edges_between=
  read.table(file="Network_analysis/RAS_F4_MIC_ASV_FUNC_Network_final_edges_between_0818.csv",
             sep=",",check.names=F, header=T)
head(mic_net_edges_between)



######

### Compare variations in diversity of functions over time based on gene clusters

library(iNEXT)
library(codyn)

# Import gene clust ID to FUNC ID to FUNCTION info
clustID_to_funcID_to_FUNC_to_NetMod=read.table(
  file="Network_analysis/RAS_F4_MIC_GENE_CLUSTID_to_FUNC_to_FUNCID_to_NetMod.txt", sep="\t",
  check.names=F, header=T)

# Import relative abundance of microbial genes with oscillation 4
mic_clust_osc4_rel = read.table(
  file="time_series_analysis/RAS_F4_MIC_GENE_CLUST_OSC4_rel_abund_wide.txt",
  sep="\t", check.names=F, header=T) %>%
  mutate(across(is.numeric, ~ .* 100)) %>%
  tibble::rownames_to_column(., var="clustID")

# Import metadata but only retain sample name to year information
sample_meta_metag=read.table(file="Network_analysis/RAS_F4_METAG_dates_and_mod_times.txt", sep="\t",
                             check.names=F, header=T) %>%
  mutate(date = format(as.Date(date, format="%m/%d/%Y")))
head(sample_meta_metag)

create_namelist <- function(infile, Mod){
  infile %>%
    filter(Module == Mod) %>%
    select(clustID)
}
filter_rel_abund_table <- function(infile, infile2){
  infile %>%
    filter(clustID %in% infile2$clustID) %>%
    tibble::column_to_rownames(., var="clustID")
}
estimate_richness_at_stand_cov <- function(infile){
  rich_df_temp = estimateD(infile, q=0, base="coverage", level=NULL, nboot=10)
  rich_df = rich_df_temp %>%
    as.data.frame() %>% 
    dplyr::rename(., RAS_id = Assemblage) %>%
    dplyr::rename(., Richness = qD) %>%
    select(RAS_id, Richness)
  return(rich_df)
}
estimate_diversity_at_stand_cov <- function(infile){
  div_df_temp = estimateD(infile, q=1, base="coverage", level=NULL, nboot=10)
  div_df = div_df_temp %>%
    as.data.frame() %>% 
    dplyr::rename(., RAS_id = Assemblage) %>%
    dplyr::rename(., Shannon_diversity = qD) %>%
    select(RAS_id, Shannon_diversity)
  return(div_df)
}
combine_div_and_plot <- function(infile,infile2,infile3){
  infile %>%
    left_join(infile2, by="RAS_id") %>%
    left_join(infile3, by="RAS_id") %>%
    ggplot(data=., aes(x = as.Date(date))) + 
    geom_line(aes(y = Richness, group = 1), colour = "blue") +
    geom_line(aes(y = Shannon_diversity, group = 1), colour = "black") +
    scale_x_date(date_breaks = "3 months", date_labels =  "%b %Y") +
    theme_bw() + 
    theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1))
}

reform_filt_rel_abund_table <- function(infile,infile2){
  infile %>%
    tibble::rownames_to_column(., var="clustID") %>%
    reshape2::melt(id.vars="clustID", variable.name="RAS_id",
                   value.name="Rel_abund") %>%
    left_join(infile2, by="RAS_id")
}
calc_turnover <- function(infile){
  turnover(df = infile,
           time.var = "M1_time", 
           species.var = "clustID", 
           abundance.var = "Rel_abund",
           replicate.var=NA)
}
calc_rank_shift <- function(infile){
  rank_shift(df = infile,
           time.var = "M1_time", 
           species.var = "clustID", 
           abundance.var = "Rel_abund",
           replicate.var=NA)
}

plot_turnover <- function(x){
  infile %>%
    left_join(infile2, by="M1_time") %>%
    ggplot(data=., aes(x = as.Date(date))) +
    geom_area(data=mic_func_module_rel_abund_long_split$M1,
              aes(x = as.Date(date), y=Rel_abund, fill = Module), alpha=0.5, show.legend=F) + 
    geom_line(aes(y=total*50, group = 1), colour = "blue") +
    scale_y_continuous(sec.axis = sec_axis(~./50, name="Turnover")) +
    theme_bw()
}


M1_names <- create_namelist(clustID_to_funcID_to_FUNC_to_NetMod,"M1")
M1_filt_rel_abund <- filter_rel_abund_table(mic_clust_osc4_rel, M1_names)
M1_rich <- estimate_richness_at_stand_cov(M1_filt_rel_abund)
M1_rich
M1_shann <- estimate_diversity_at_stand_cov(M1_filt_rel_abund)
M1_shann
M1_filt_rel_abund_with_date = M1_filt_rel_abund %>%
  tibble::rownames_to_column(., var="clustID") %>%
  reshape2::melt(id.vars="clustID", variable.name="RAS_id", value.name = "Rel_abund") %>%
  left_join(sample_meta_metag, by="RAS_id")

m1_turnover <- turnover(df = M1_filt_rel_abund_with_date,
                        time.var = "M1_time", 
                        species.var = "clustID", 
                        abundance.var = "Rel_abund",
                        replicate.var=NA)
m1_turnover_plot =
  m1_turnover %>%
  left_join(sample_meta_metag, by="M1_time") %>%
    ggplot(data=., aes(x = as.Date(date))) +
  geom_area(data=mic_func_module_rel_abund_long_split$M1,
            aes(x = as.Date(date), y=Rel_abund, fill = Module), alpha=0.5, show.legend=F) + 
  geom_line(aes(y=total*50, group = 1), colour = "blue") +
  scale_y_continuous(sec.axis = sec_axis(~./50, name="Turnover")) +
  theme_bw()
m1_turnover_plot

m1_rich_and_shan/m1_turnover_plot

m1_rank_shift <- rank_shift(df=filt_table_with_date,
                            time.var = "M1_time",
                            species.var = "clustID",
                            abundance.var = "Rel_abund")

m1_rank_shift %>%
  mutate(M1_time = gsub(".*-","",year_pair)) %>%
  mutate(M1_time = as.numeric(M1_time)) %>%
  left_join(sample_meta_metag, by="M1_time") %>%
  mutate(peak = case_when(
    RAS_id %in% M1_peak_sample_ref$RAS_id ~ 100,
    TRUE ~ 0 
  )) %>%
  ggplot(data=., aes(x = as.Date(date))) +
  geom_area(data=mic_func_module_rel_abund_long_split$M1,
            aes(x = as.Date(date), y=Rel_abund, fill = Module), alpha=0.5, show.legend=F) + 
  geom_line(aes(y=MRS/350, group = 1), colour = "blue") +
  #geom_vline(aes(xintercept = ifelse(peak == 100, date, NA))) + 
  scale_y_continuous(sec.axis = sec_axis(~.*350, name="Turnover")) +
  theme_bw()

m1_rank_shift_plot

m1_rate_change <- rank_shift(df=filt_table_with_date,
                            time.var = "M1_time",
                            species.var = "clustID",
                            abundance.var = "Rel_abund")
m1_rate_change_plot = 
  m1_rate_change %>%
  mutate(M1_time = gsub("-.*","",year_pair)) %>%
  mutate(M1_time = as.numeric(M1_time)) %>%
  left_join(sample_meta_metag, by="M1_time") %>%
  ggplot(data=., aes(x = as.Date(date), y = MRS)) +
  geom_line(aes(group = 1), colour = "black") +
  theme_bw()

m1_turnover_plot/m1_rank_shift_plot/m1_rate_change_plot/net_module_M1_dynamics



#####

### Analyse AAI variability in functions of modules during peak events across
### time-series

# Import sample_name to date
sample_meta_metag=read.table(file="Network_analysis/RAS_F4_METAG_dates_and_mod_times.txt", sep="\t",
                             check.names=F, header=T) %>%
  mutate(date = format(as.Date(date, format="%m/%d/%Y")))

# Import AAI table for module and remove self hits
input_dir <- ("Network_analysis/")

import_aai_and_reformat_single_peaks <- function(input_file){
  mod_aai_df <- read.table(file = paste0(input_dir, input_file), 
                           sep="\t", check.names=F, header=T) %>%
    group_by(grp = paste(pmax(Gene_1, Gene_2), pmin(Gene_1, Gene_2), sep = "_")) %>%
    filter(row_number() == 1) %>%
    ungroup() %>%
    select(-grp) %>%
    as.data.frame() %>%
    left_join(sample_meta_metag, by=c("Sample_1" = "RAS_id")) %>%
    mutate(year_1 = format(as.Date(date, format="%m/%d/%Y"),"%Y")) %>%
    left_join(sample_meta_metag, by=c("Sample_2" = "RAS_id")) %>%
    mutate(year_2 = format(as.Date(date.y, format="%m/%d/%Y"),"%Y")) %>%
    mutate(year_comparison = paste0(year_1,"_vs_",year_2))
  return(mod_aai_df)
}
import_aai_and_reformat_all <- function(input_file){
  mod_aai_df <- read.table(file = paste0(input_dir, input_file), 
                           sep="\t", check.names=F, header=T) %>%
    group_by(grp = paste(pmax(Gene_1, Gene_2), pmin(Gene_1, Gene_2), sep = "_")) %>%
    filter(row_number() == 1) %>%
    ungroup() %>%
    select(-grp) %>%
    as.data.frame() %>%
    left_join(sample_meta_metag, by=c("Sample_2" = "RAS_id"))
  return(mod_aai_df)
}

# Run import function
mic_func_net_M1_aai = import_aai_and_reformat_all("M1_func_all_vs_all_aai_hits_refined.txt")
mic_func_net_M2_aai = import_aai_and_reformat_all("M2_func_all_vs_all_aai_hits_refined.txt")
mic_func_net_M3_aai = import_aai_and_reformat_all("M3_func_all_vs_all_aai_hits_refined.txt")
mic_func_net_M4_aai = import_aai_and_reformat_all("M4_func_all_vs_all_aai_hits_refined.txt")

head(mic_func_net_M2_aai)

calc_mean_aai <- function(infile){
  infile %>%
    aggregate(AAI~Sample_2, data=., FUN=mean) %>%
    rename(mean_AAI = AAI)
}
calc_sd_aai <- function(infile){
  infile %>%
    aggregate(AAI~Sample_2, data=., FUN=sd) %>%
    rename(sd_AAI = AAI)
}
combine_data <- function(infile,infile2){
  infile %>%
    left_join(infile2, by="Sample_2") %>%
    left_join(sample_meta_metag, by=c("Sample_2" = "RAS_id")) %>%
    mutate(date = as.Date(date))
}
plot_aai_temporal <- function(infile){
  ggplot(data=infile,
         aes(y=mean_AAI, x = date)) + 
    geom_point() +
    geom_errorbar(aes(ymin=mean_AAI-sd_AAI, ymax=mean_AAI+sd_AAI), 
                  stat="identity",width=.2,
                  position=position_dodge(.9)) + 
    scale_x_date(date_breaks = "3 months", date_labels =  "%b %Y") +
    theme_bw() + 
    theme(axis.text.x = element_text(size = 8, angle=45, hjust=1))
}

m1_mean_aai <- calc_mean_aai(mic_func_net_M1_aai)
m1_sd_aai <- calc_sd_aai(mic_func_net_M1_aai)
m1_mean_and_sd_aai <- combine_data(m1_mean_aai,m1_sd_aai)
m1_aai_temporal_point<- plot_aai_temporal(m1_mean_and_sd_aai)

m2_mean_aai <- calc_mean_aai(mic_func_net_M2_aai)
m2_sd_aai <- calc_sd_aai(mic_func_net_M2_aai)
m2_mean_and_sd_aai <- combine_data(m2_mean_aai,m2_sd_aai)
m2_aai_temporal_point<- plot_aai_temporal(m2_mean_and_sd_aai)

m3_mean_aai <- calc_mean_aai(mic_func_net_M2_aai)
m3_sd_aai <- calc_sd_aai(mic_func_net_M2_aai)
m3_mean_and_sd_aai <- combine_data(m2_mean_aai,m2_sd_aai)
m3_aai_temporal_point<- plot_aai_temporal(m2_mean_and_sd_aai)

m1_aai_temporal_point/net_module_M1_dynamics
m2_aai_temporal_point/net_module_M2_dynamics
m3_aai_temporal_point/net_module_M3_dynamics

# FUNCTION: Visualise intra-year
# and inter-year AAI variability as boxplots
visualise_inter_year_aai_boxplot <- function(input_aai_df){
  input_aai_df %>%
    filter(., year_comparison == "2016_vs_2017" | 
             year_comparison == "2016_vs_2018" |
             year_comparison == "2016_vs_2019" |
             year_comparison == "2016_vs_2020" |
             year_comparison == "2017_vs_2018" |
             year_comparison == "2017_vs_2019" | 
             year_comparison == "2017_vs_2020" | 
             year_comparison == "2018_vs_2019" |
             year_comparison == "2019_vs_2020") %>%
    ggplot(data=., aes(x = year_comparison, y = AAI)) + 
    geom_boxplot() + 
    theme_bw() + 
    theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
}
visualise_intra_year_aai_boxplot <- function(input_aai_df){
  input_aai_df %>%
    filter(year_1 == year_1) %>%
    ggplot(data=., aes(x = year_1, y = AAI)) + 
    geom_boxplot() + 
    theme_bw()
}

mic_func_net_M1_aai_intra_year_boxplot = visualise_intra_year_aai_boxplot(mic_func_net_M1_aai)
mic_func_net_M2_aai_intra_year_boxplot = visualise_intra_year_aai_boxplot(mic_func_net_M2_aai)
mic_func_net_M3_aai_intra_year_boxplot = visualise_intra_year_aai_boxplot(mic_func_net_M3_aai)
mic_func_net_M4_aai_intra_year_boxplot = visualise_intra_year_aai_boxplot(mic_func_net_M4_aai)

mic_func_net_M1_aai_inter_year_boxplot = visualise_inter_year_aai_boxplot(mic_func_net_M1_aai)
mic_func_net_M2_aai_inter_year_boxplot = visualise_inter_year_aai_boxplot(mic_func_net_M2_aai)
mic_func_net_M3_aai_inter_year_boxplot = visualise_inter_year_aai_boxplot(mic_func_net_M3_aai)
mic_func_net_M4_aai_inter_year_boxplot = visualise_inter_year_aai_boxplot(mic_func_net_M4_aai)

# Combine intra- and inter-year boxplots for each module
library(patchwork)
mic_func_net_m1_boxplots <- mic_func_net_M1_aai_intra_year_boxplot/
  mic_func_net_M1_aai_inter_year_boxplot
mic_func_net_m2_boxplots <- mic_func_net_M2_aai_intra_year_boxplot/
  mic_func_net_M2_aai_inter_year_boxplot
mic_func_net_m3_boxplots <- mic_func_net_M3_aai_intra_year_boxplot/
  mic_func_net_M3_aai_inter_year_boxplot
mic_func_net_m4_boxplots <- mic_func_net_M4_aai_intra_year_boxplot/
  mic_func_net_M4_aai_inter_year_boxplot

# FUNCTION: Visualise intra-year
# and inter-year AAI variability as histograms
visualise_inter_year_aai_histogram <- function(input_aai_df){
  input_aai_df %>%
    filter(., year_comparison == "2016_vs_2017" | 
             year_comparison == "2016_vs_2018" |
             year_comparison == "2016_vs_2019" |
             year_comparison == "2016_vs_2020" |
             year_comparison == "2017_vs_2018" |
             year_comparison == "2017_vs_2019" | 
             year_comparison == "2017_vs_2020" | 
             year_comparison == "2018_vs_2019" |
             year_comparison == "2019_vs_2020") %>%
    ggplot(data=., aes(x = AAI, colour = year_comparison)) + 
    geom_freqpoly(alpha = 0.75, binwidth = 1, linewidth=1.5) +
    scale_x_continuous(limits=c(40,100)) + 
    scale_colour_brewer(palette = "Paired") + 
    labs(x = "Average amino acid identity (%)",
         y = "Count", colour = "Year") + 
    theme_bw() + 
    theme(axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          axis.title.x = element_text(size = 11),
          axis.title.y = element_text(size = 11),
          legend.title = element_text(size = 11),
          legend.text = element_text(size = 10)) + 
    guides(colour = guide_legend(override.aes = list(size = 3)))
}
visualise_intra_year_aai_histogram <- function(input_aai_df){
  input_aai_df %>%
    filter(year_1 == year_1) %>%
    ggplot(data=., aes(x = AAI, colour = year_1)) + 
    geom_freqpoly(alpha=0.75, binwidth = 1, linewidth=1) +
    scale_x_continuous(limits=c(40,100)) + 
    scale_colour_brewer(palette = "Set2") + 
    labs(x = "Average amino acid identity (%)",
         y = "Count", colour = "Year") + 
    theme_bw() + 
    theme(axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          axis.title.x = element_text(size = 11),
          axis.title.y = element_text(size = 11),
          legend.title = element_text(size = 11),
          legend.text = element_text(size = 10)) + 
    guides(colour = guide_legend(override.aes = list(size = 3)))
}

mic_func_net_M1_aai_intra_year_histogram = visualise_intra_year_aai_histogram(mic_func_net_M1_aai)
mic_func_net_M2_aai_intra_year_histogram = visualise_intra_year_aai_histogram(mic_func_net_M2_aai)
mic_func_net_M3_aai_intra_year_histogram = visualise_intra_year_aai_histogram(mic_func_net_M3_aai)
mic_func_net_M4_aai_intra_year_histogram = visualise_intra_year_aai_histogram(mic_func_net_M4_aai)

mic_func_net_M1_aai_inter_year_histogram = visualise_inter_year_aai_histogram(mic_func_net_M1_aai)
mic_func_net_M2_aai_inter_year_histogram = visualise_inter_year_aai_histogram(mic_func_net_M2_aai)
mic_func_net_M3_aai_inter_year_histogram = visualise_inter_year_aai_histogram(mic_func_net_M3_aai)
mic_func_net_M4_aai_inter_year_histogram = visualise_inter_year_aai_histogram(mic_func_net_M4_aai)

### Export plots
pdf(file="figures_output/FRAM_RAS_F4_NET_MIC_FUNC_MOD1_AAI_boxplot.pdf",
    height=8, width=8)
mic_func_net_M1_aai_intra_year_boxplot/
  mic_func_net_M1_aai_inter_year_boxplot
dev.off()

pdf(file="figures_output/FRAM_RAS_F4_NET_MIC_FUNC_MOD2_AAI_boxplot.pdf",
    height=8, width=8)
mic_func_net_M2_aai_intra_year_boxplot/
  mic_func_net_M2_aai_inter_year_boxplot
dev.off()

pdf(file="figures_output/FRAM_RAS_F4_NET_MIC_FUNC_MOD3_AAI_boxplot.pdf",
    height=8, width=8)
mic_func_net_M3_aai_intra_year_boxplot/
  mic_func_net_M3_aai_inter_year_boxplot
dev.off()

pdf(file="figures_output/FRAM_RAS_F4_NET_MIC_FUNC_MOD4_AAI_boxplot.pdf",
    height=8, width=8)
mic_func_net_M4_aai_intra_year_boxplot/
  mic_func_net_M4_aai_inter_year_boxplot
dev.off()

pdf(file="figures_output/FRAM_RAS_F4_NET_MIC_FUNC_MOD1_AAI_histogram.pdf",
    height=5, width=6)
mic_func_net_M1_aai_intra_year_histogram/
  mic_func_net_M1_aai_inter_year_histogram
dev.off()

pdf(file="figures_output/FRAM_RAS_F4_NET_MIC_FUNC_MOD2_AAI_histogram.pdf",
    height=5, width=6)
mic_func_net_M2_aai_intra_year_histogram/
  mic_func_net_M2_aai_inter_year_histogram
dev.off()

pdf(file="figures_output/FRAM_RAS_F4_NET_MIC_FUNC_MOD3_AAI_histogram.pdf",
    height=5, width=6)
mic_func_net_M3_aai_intra_year_histogram/
  mic_func_net_M3_aai_inter_year_histogram
dev.off()

pdf(file="figures_output/FRAM_RAS_F4_NET_MIC_FUNC_MOD4_AAI_histogram.pdf",
    height=5, width=6)
mic_func_net_M4_aai_intra_year_histogram/
  mic_func_net_M4_aai_inter_year_histogram
dev.off()

################################

mic_net_edges_between=
  read.table(file="Network_analysis/RAS_F4_MIC_ASV_FUNC_Network_final_type_to_type_0818.csv",
             sep=",",check.names=F, header=T)
head(mic_net_edges_between)


RAS_F4_MIC_ASV_FUNC_Network_final_edges_between_0818


