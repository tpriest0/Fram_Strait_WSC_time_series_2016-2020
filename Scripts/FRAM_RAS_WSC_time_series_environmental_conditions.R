####################

#### Plot environmental conditions and sampling time points

####################

### Define working directory
setwd('C:/Users/tpriest/OneDrive - ETH Zurich/MPI - FRAM/WSC')

### Load libraries
library(dplyr)
library(ggplot2)
library(patchwork)
library(tibble)
library(vegan)
library(stats)
library(Hmisc)
library(RcmdrMisc)
library(reshape2)
library(purrr)
library(tidyr)

### Import files
# Import metdata
ras_metadata=read.table(file="RAS_F4_META.txt", sep="\t",
                        check.names=F, header=T) %>%
  mutate(date = as.Date(date, "%d.%m.%Y"))

ras_16s_samples = ras_metadata %>%
  filter(SSU_16S == "16S") %>%
  select(RAS_id,date)

ras_mg_samples = ras_metadata %>%
  filter(MetaG == "MG") %>%
  select(RAS_id,date)

### Create sampling time point plot
ras_wsc_asv_and_mg_sampling_timepoints_plot = 
  ggplot(ras_16s_samples, aes(x=date)) + 
  theme_void() +
  ylim(0,0.25) + 
  geom_point(data=ras_16s_samples, aes(y=0.1),
             colour="#009292", size=2, shape=17) + 
  geom_point(data=ras_mg_samples, aes(y=0.2),
             colour="#490092", size=2, shape=17) 

### Define function to plot environmental variables
plot_env_variables <- function(infile, env_var, colour_var, labs_var){
  env_var <- enquo(env_var)
  ggplot(infile, aes(x=date, y=!!env_var)) + 
    geom_line(colour = colour_var, linewidth=2) + 
    labs(y = labs_var) + 
    scale_x_date(date_breaks = "6 months", date_labels =  "%b %Y") +
    theme_bw() + 
    theme(axis.text.y = element_text(size = 10),
          axis.text.x = element_blank(),
          axis.title.y = element_text(size = 11),
          axis.title.x = element_blank())
}

plot_env_variables_bottom <- function(infile, env_var, colour_var, labs_var){
  env_var <- enquo(env_var)
  ggplot(infile, aes(x=date, y=!!env_var)) + 
    geom_line(colour = colour_var, linewidth=2) + 
    labs(y = labs_var) + 
    scale_x_date(date_breaks = "6 months", date_labels =  "%b %Y") +
    theme_bw() + 
    theme(axis.text.y = element_text(size = 10),
          axis.text.x = element_text(size = 8),
          axis.title.y = element_text(size = 11),
          axis.title.x = element_blank())
}

### Execute function for each variable
time_series_temp_plot <- plot_env_variables(ras_metadata, temp, "black", "Temperature (?C)")
time_series_MLD_plot <- plot_env_variables(ras_metadata, MLD, "#924900", "Mixed Layer Depth (m)")
time_series_PAR_plot <- plot_env_variables_bottom(ras_metadata, PAR_satellite, "#d6b85a", 
                                           "Photosynthetically active radiation 
                   (?mol photons ^m-2 d-1)")
time_series_chla_plot <- plot_env_variables(ras_metadata, Chl_a, "#009292", 
                                           "Chlorophyll a (mg m^3)")
time_series_oxygen_plot <- plot_env_variables(ras_metadata, O2_sat, "#006ddb", 
                                            "Oxygen saturation")


### Combine all plots together into Figure and export
pdf(file="figures_output/FRAM_RAS_F4_env_conditions_and_sampling_timepoints.pdf",height=8, width=10)
time_series_temp_plot/time_series_oxygen_plot/time_series_MLD_plot/
  time_series_chla_plot/time_series_PAR_plot/
  ras_wsc_asv_and_mg_sampling_timepoints_plot+
  plot_layout(nrow=6,heights=c(3,3,3,3,3,0.5))
dev.off()

### Run correlation analysis between env variables
head(ras_metadata)
View(ras_metadata)
ras_meta_select_stand = ras_metadata %>%
  select(RAS_id,depth,temp,sal,O2_conc,Chl_a,PAR_satellite,MLD) %>%
  tibble::column_to_rownames(., var="RAS_id") %>%
  decostand(., method="standardize", MARGIN=2) %>%
  as.data.frame(.)
  
ras_meta_select_stand_corr = rcorr(as.matrix(ras_meta_select_stand))
ras_meta_select_stand_corr = rcorr.adjust(as.matrix(ras_meta_select_stand), type = "pearson")

ras_meta_select_stand_corr_sig_p = as.data.frame(ras_meta_select_stand_corr$P) %>%
  tibble::rownames_to_column(., var="var_one") %>%
  reshape2::melt(id.vars="var_one", variable.name="var_two", value.name="p_val") %>%
  filter(., p_val < 0.05 & p_val > 0)
ras_meta_select_stand_corr$R$r
ras_meta_select_stand_corr_coeff = as.data.frame(ras_meta_select_stand_corr$R$r) %>%
  tibble::rownames_to_column(., var="var_one") %>%
  reshape2::melt(id.vars="var_one", variable.name="var_two", value.name="coefficient")

ras_meta_select_stand_corr_sig_p %>%
  left_join(ras_meta_select_stand_corr_coeff, by=c("var_one", "var_two"))


### Can alternatively run all vs all pairwise linear regressions
out <- crossing(y = colnames(ras_meta_select_stand), x = colnames(ras_meta_select_stand)) %>%
  mutate(mod = map2(x, y,
                    ~ summary(lm(reformulate(.x, response = .y), data = ras_meta_select_stand))))
out$mod

