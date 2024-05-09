####################

#### Plot environmental conditions and sampling time points

####################

### Define working directory
setwd('XXXXX')

### Define output directories
output_figures <- ('')
output_tables <- ('')

### Load libraries
library(dplyr)
library(ggplot2)
library(patchwork)

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
pdf(file=paste0(output_figures,"FRAM_RAS_F4_env_conditions_and_sampling_timepoints.pdf"),height=8, width=10)
time_series_temp_plot/time_series_oxygen_plot/time_series_MLD_plot/
  time_series_chla_plot/time_series_PAR_plot/
  ras_wsc_asv_and_mg_sampling_timepoints_plot+
  plot_layout(nrow=6,heights=c(3,3,3,3,3,0.5))
dev.off()


