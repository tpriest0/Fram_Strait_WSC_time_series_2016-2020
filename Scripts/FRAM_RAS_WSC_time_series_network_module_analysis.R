
####################

#### Processing network clusters and visualise cluster dynamics over time
#### series

####################

### Define working directory
setwd('XXXXX')

output_tables <- ('output_tables')
output_figures <- ('output_figures')
dir.create(output_tables)
dir.create(output_figures)

### Load libraries
library(dplyr)
library(reshape2)
library(tibble)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(forcats)
library(tidyverse)
library(stringr)
library(vegan)
library(tibble)
library(Hmisc)

### Import data
# Import prokaryotic ASV relative abundance data
mic_asv_rel=read.table(file="RAS_F4_MIC_ASV_filt_rel.txt", sep="\t",
                       check.names=F, header=T, row.names=1) %>%
  tibble::rownames_to_column(., var="nodeId")

# Import ASV and gene network cluster assignments
net_mod_components=
  read.table(file="RAS_F4_NET_MOD_assignments.txt",
                    sep="\t", check.names=F, header=T)


# Import gene cluster ID to function to function ID
clust_id_func_id_to_func=read.csv(
  file="RAS_F4_GENES_CLUSTID_to_FUNCID_to_FUNC.txt",
  sep="\t",check.names=F, header=T)

# Import function abundance profile
mic_func_rel_wide = 
  read.table(file="RAS_F4_MIC_OSC4_FUNC_CLUST_rel_abund_wide.txt",sep="\t",
             check.names=F) %>%
  tibble::rownames_to_column(., var="funcID")

# Import metadata
sample_meta = 
  read.table(file="RAS_F4_META.txt",sep="\t",
             check.names=F, header=T)

#####
### Firstly, combine funcID to network module with information on gene clusters to funcID
clust_to_func_to_net_mod = 
  clust_id_func_id_to_func %>%
  left_join(net_mod_components, by=c("funcID" = "nodeId"))

write.table(clust_to_func_to_net_mod,
            file = "RAS_F4_GENES_CLUSTID_to_FUNCID_to_FUNC_to_MOD.txt",
            sep="\t", row.names=F)

#####

##### 

### Visualise the number of functional and taxonomic components in each 
### module
net_mod_components %>%
  select(Module,type) %>%
  table() %>%
  as.data.frame()

net_mod_component_counts_barplot = 
  net_mod_components %>%
  select(Module,type) %>%
  table() %>%
  as.data.frame() %>%
  ggplot(., aes(x = type, y = Freq)) + 
  geom_bar(aes(fill = type), stat = "identity", position = "stack") + 
  facet_grid(.~Module, scales = "free_x") + 
  scale_y_log10() + 
  labs(y = "Number of components") +
  scale_fill_manual(values = c("ASV" = "#616569",
                               "FUNC" = "#C5C6C7")) + 
  theme_bw() + 
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 12),
        strip.background.x = element_rect(colour = "black",
                                          fill = "white"),
        strip.text.x = element_text(size = 14))

# Export figure
pdf(file = paste0(output_figures,"/FRAM_RAS_F4_MOD_component_counts_barplot.pdf"),
    height = 5, width = 7)
net_mod_component_counts_barplot
dev.off()

######

### Determine the number of ASVs and Functions assigned to each module
mic_asv_func_network_module_counts = 
  net_mod_components %>%
  select(Module, type) %>%
  table() %>%
  as.data.frame()

### Determine the combined module relative abundance of functional components
mic_func_module_rel_abund_long = 
  mic_func_rel_wide %>%
  filter(funcID %in% net_mod_components$nodeId) %>%
  left_join(net_mod_components, by=c("funcID" = "nodeId")) %>%
  reshape2::melt(id.vars=c("funcID","Module","type"),
                 variable.name = "RAS_id", value.name="Rel_abund") %>%
  aggregate(Rel_abund~Module+RAS_id, data=., FUN=sum) %>%
  left_join(sample_meta, by="RAS_id") %>%
  mutate(date = as.Date(date, "%d.%m.%Y")) %>%
  mutate(Rel_abund = Rel_abund*100)

### Determine the combined module relative abundance of ASV components
mic_asv_module_rel_abund_long = 
  mic_asv_rel %>%
  filter(nodeId %in% net_mod_components$nodeId) %>%
  left_join(net_mod_components, by="nodeId") %>%
  reshape2::melt(id.vars=c("nodeId","Module","type"), 
                 variable.name="RAS_id", value.name="Rel_abund") %>%
  aggregate(Rel_abund~Module+type+RAS_id,
            data=., FUN=sum) %>%
  left_join(sample_meta, by="RAS_id") %>%
  mutate(date = as.Date(date, "%d.%m.%Y")) %>%
  mutate(Rel_abund = Rel_abund*100)

### Visualise relative abundance dynamics of modules
# set cluster colours for plots
module_colours = c("M1" = "#008066",
                    "M2" = "#55DDFF",
                    "M3" = "#FFB300",
                    "M4" = "#F55BF5",
                    "M5" = "#490092",
                    "M6" = "#920000")

# Define function for plotting modules asv and gene relative abundance
plot_module_abund = function(infile_asv, infile_gene, text_var){
  ggplot(data=infile_asv,
         aes(x=as.Date(date), y=Rel_abund)) + 
    geom_area(aes(fill = Module), alpha=0.5, show.legend=F) + 
    geom_line(data=infile_gene, 
              aes(x=as.Date(date), y=Rel_abund*2, colour = Module), 
              linewidth=1, show.legend=F) +
    geom_point(data=infile_gene, 
               aes(x=as.Date(date), y=Rel_abund*2, colour = Module), 
               show.legend=F) +
    scale_x_date(date_breaks = "6 months", date_labels =  "%b %Y") +
    scale_y_continuous(sec.axis = 
                         sec_axis(~./2, name="Cluster genes relative abundance (%)")) + 
    scale_fill_manual(values = module_colours) + 
    scale_colour_manual(values = module_colours) + 
    labs(y = "Relative abundance (%)",
         title = text_var) + 
    theme_bw() + 
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          plot.background = element_blank())
}

# specify code to remove x-axis labels for plots (prevents repeated axes when stacking plots)
remove_x_axis <- theme(
  axis.text.x = element_blank(),
  axis.title.x = element_blank()
)


### Split the gene and asv module relative abundance dataframe by cluster
mic_asv_module_rel_abund_long_split = 
  split(mic_asv_module_rel_abund_long,mic_asv_module_rel_abund_long$Module)
mic_func_module_rel_abund_long_split = 
  split(mic_func_module_rel_abund_long,mic_func_module_rel_abund_long$Module)

### Create the plots for each cluster
net_module_M1_dynamics = plot_module_abund(mic_asv_module_rel_abund_long_split$M1,
                                           mic_func_module_rel_abund_long_split$M1, 
                                           "M1 (early polar night): 153 ASVs, 2165 Functional genes")
net_module_M2_dynamics = plot_module_abund(mic_asv_module_rel_abund_long_split$M2,
                                           mic_func_module_rel_abund_long_split$M2,
                                           "M2 (late polar night): 244 ASVs, 1778 Functional genes")
net_module_M3_dynamics = plot_module_abund(mic_asv_module_rel_abund_long_split$M3,
                                           mic_func_module_rel_abund_long_split$M3,
                                           "M3 (spring): 66 ASVs, 3140 Functional genes")
net_module_M4_dynamics = plot_module_abund(mic_asv_module_rel_abund_long_split$M4,
                                           mic_func_module_rel_abund_long_split$M4,
                                           "M4 (summer): 135 ASVs, 3078 Functional genes")
net_module_M5_dynamics = plot_module_abund(mic_asv_module_rel_abund_long_split$M5,
                                           mic_func_module_rel_abund_long_split$M5,
                                           "M5 (autumn): 60 ASVs, 1159 Functional genes")

### Combine all plots and export
library(patchwork)
pdf(file=paste0(output_figures,"/FRAM_RAS_F4_MOD_temporal_dynamics.pdf"),
    height=10, width=8)
(net_module_M1_dynamics+remove_x_axis)/(net_module_M2_dynamics+remove_x_axis)/
  (net_module_M3_dynamics+remove_x_axis)/(net_module_M4_dynamics+remove_x_axis)/
  net_module_M5_dynamics
dev.off()

### Export module ASV and gene abundance tables
write.table(mic_func_module_rel_abund_long,
            file="RAS_F4_NET_MOD_FUNC_relative_abundance.txt", sep="\t",
            quote = F, row.names = F)

write.table(mic_asv_module_rel_abund_long,
            file="RAS_F4_NET_MOD_ASV_relative_abundance.txt", sep="\t",
            quote = F, row.names = F)


### ALTERNATIVE VISUALISATION:
### area graphs to better illustrate successional dynamics of module ASVs and module
### functions 


pdf(file=paste0(output_figures,"/FRAM_RAS_F4_MOD_ASV_only_temporal_dynamics_area_nonstacked.pdf"),
    height=6, width=8)
mic_asv_module_rel_abund_long %>%
  mutate(Module = factor(Module, levels = c("M3","M4","M2","M1","M6","M5"))) %>%
  ggplot(data=.,
         aes(x=as.Date(date), y=Rel_abund)) + 
  geom_area(aes(fill = Module), position="identity", alpha=1, show.legend=F) + 
  scale_x_date(date_breaks = "6 months", date_labels =  "%b %Y") + 
  scale_fill_manual(values=module_colours) + 
  theme_bw() + 
  labs(y = "Relative abundance of module functional genes (%)") + 
  theme(axis.text.x = element_text(size = 12, angle=45, hjust=1),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        plot.background = element_blank())
dev.off()  

pdf(file=paste0(output_figures,"/FRAM_RAS_F4_MOD_function_only_temporal_dynamics.pdf"),
    height=6, width=8)
mod_functional_temporal_area_plot = 
  mic_func_module_rel_abund_long %>%
  mutate(Module = factor(Module, levels = c("M3","M4","M2","M1","M6","M5"))) %>%
  ggplot(data=.,
              aes(x=as.Date(date), y=Rel_abund)) + 
         geom_area(aes(fill = Module), position="identity", alpha=1, show.legend=F) + 
         scale_x_date(date_breaks = "6 months", date_labels =  "%b %Y") + 
  scale_fill_manual(values=module_colours) + 
  theme_bw() + 
  labs(y = "Relative abundance of module functional genes (%)") + 
  theme(axis.text.x = element_text(size = 12, angle=45, hjust=1),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        plot.background = element_blank())
mod_functional_temporal_area_plot
dev.off()  

### Create dataframe with module ASV and Func abundance over time
mod_list <- c("M1","M2","M3","M4","M5")
mod_asv_list <-c("M1_asv","M2_asv","M3_asv","M4_asv","M5_asv")
names(mic_asv_module_rel_abund_long_split) <- mod_asv_list
list2env(mic_asv_module_rel_abund_long_split,envir=.GlobalEnv)

mod_gen_list <-c("M1_function","M2_function","M3_function","M4_function","M5_function")
names(mic_func_module_rel_abund_long_split) <- mod_gen_list
list2env(mic_func_module_rel_abund_long_split,envir=.GlobalEnv)

refine_mod_df <- function(x){
  mic_asv_module_rel_abund_long_split$x %>%
    filter(Module,Type,RAS_id,Rel_abund)
}

###############

##### What are the factors driving module dynamics?

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
    select(RAS_id,Rel_abund) %>%
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
correct_and_filter_cor_results <- function(cor_df, asv_or_func, modname){
  adjusted_p = p.adjust(cor_df$p, method = "bonferroni")
  cor_df %>%
    mutate(p_adj = adjusted_p) %>% 
    filter(p_adj < 0.05) %>%
    filter(row == "Rel_abund") %>%
    mutate(Type = asv_or_func) %>%
    mutate(Module = modname)
}

# Run above functions for the both module ASV abundance and function abundance
# and create a dataframe containing correlation coefficients associated with 
# adjusted p-values < 0.05

for (modname in c("M1","M2","M3","M4","M5")){
  input_file = mic_asv_module_rel_abund_long %>%
    filter(., Module == modname)
  
  cor_results = diversity_vs_env_correlation(input_file, sample_meta_stand)
  
  cor_df = process_correlation_matrix(cor_results$r, cor_results$P)
  
  df_output_name  = paste0(modname,"_asv_env_cor_sig_df")
  assign(df_output_name, correct_and_filter_cor_results(cor_df, "ASV", modname))
}

for (modname in c("M1","M2","M3","M4","M5")){
  input_file = mic_func_module_rel_abund_long %>%
    filter(., Module == modname)
  
  cor_results = diversity_vs_env_correlation(input_file, sample_meta_stand)
  
  cor_df = process_correlation_matrix(cor_results$r, cor_results$P)
  
  df_output_name  = paste0(modname,"_func_env_cor_sig_df")
  assign(df_output_name, correct_and_filter_cor_results(cor_df, "FUNC", modname))
}

# Combine the outputs into single dataframe
module_asv_func_vs_env_sig_cor_df = 
  rbind(M1_asv_env_cor_sig_df, M2_asv_env_cor_sig_df, M3_asv_env_cor_sig_df,
      M4_asv_env_cor_sig_df, M5_asv_env_cor_sig_df,
      M1_func_env_cor_sig_df, M2_func_env_cor_sig_df, M3_func_env_cor_sig_df,
      M4_func_env_cor_sig_df, M5_func_env_cor_sig_df)


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
ci_lower <- numeric(nrow(module_asv_func_vs_env_sig_cor_df))
ci_upper <- numeric(nrow(module_asv_func_vs_env_sig_cor_df))

# Loop through significant correlation dataframe and calculate confidence intervals
for (i in 1:nrow(module_asv_func_vs_env_sig_cor_df)) {
  correlation <- module_asv_func_vs_env_sig_cor_df$cor[i]
  n <- 97
  ci <- calculate_ci(correlation, n)
  ci_lower[i] <- ci[1]
  ci_upper[i] <- ci[2]
}

# Add the confidence intervals to the respective columns
module_asv_func_vs_env_sig_cor_df$ci_lower <- ci_lower
module_asv_func_vs_env_sig_cor_df$ci_upper <- ci_upper


### Visualise the significant correlations

net_mod_env_sig_corr_heatmap = 
  ggplot(module_asv_func_vs_env_sig_cor_df, aes(y=Type, x=column)) +
  geom_tile(aes(fill=cor)) + 
  scale_fill_gradient2(limits=c(-1,1),
                       low = "#40004B", high = "#00441B", mid = "#ffffff") + 
  labs(fill = "Pearson's correlation") + 
  facet_wrap(Module~., scales = "free_y", ncol = 1, nrow = 6) +
  scale_y_discrete(limits=rev) + 
  theme_bw() + 
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        strip.background = element_blank(),
        strip.text = element_blank())

pdf(file=paste0(output_figures,"/RAS_F4_MOD_env_corr_heatmap.pdf"),
    width = 8, height = 10)
net_mod_env_sig_corr_heatmap
dev.off()

### Export correlation table
write.table(module_asv_func_vs_env_sig_cor_df,
            file=paste0(output_tables,"/RAS_F4_MOD_ASV_FUNC_vs_env_cor_results.txt"), sep="\t")



################ 

#####

### This section is completely optional!
### It creates and outputs a table for each module that contains only the kegg ko
### information. These files can directly be uploaded to kegg mapper to explore
### metabolic maps

# Filter above table to only include functions with KEGG annotation in preparation for 
# gene enrichment analysis
kegg_func_modules = 
  clust_to_func_to_net_mod %>%
  select(Module,funcID,FUNC) %>%
  filter(., grepl("^ko",FUNC)) %>%
  mutate(FUNC = gsub("ko:","",FUNC)) %>%
  select(Module,FUNC) 

for (modulename in c("M1","M2","M3","M4","M5")){
  name_out <- paste(modulename, "_kegg_only", sep = "")
  
  assign(name_out, kegg_func_modules %>%
                      filter(Module == modulename) %>%
           unique)
}

write.table(M5_kegg_only,
            file="Network_analysis/FRAM_RAS_F4_NET_MOD_M5_kegg_only.txt",
            sep="\t")


#######################

###### Examine the content of modules
###### in terms of microbial and functional composition

######
library(tidytext)
library(ggdendro)
library(dendextend)
library(igraph)
library(tidyverse)
library(tidyr)
library(ggstatsplot)

### import necessary files
options(scipen=999)

### MODULE TAXONOMIC COMPOSITION

# Import module to asv info
net_mod_components=read.table(
  file="RAS_F4_NET_MOD_assignments.txt", sep="\t",
  check.names=F, header=T)

# Import microbial ASV relative abundance data
mic_asv_rel=read.table(file="RAS_F4_MIC_ASV_filt_rare_rel.txt", sep="\t",
                       check.names=F, header=T)

# Import microbial ASV taxa info
mic_asv_taxa=read.table(file="RAS_F4_MIC_ASV_filt_taxa.txt", sep="\t",
                       check.names=F, header=T)

# Import module node connections
net_mod_connections = read.table(file="RAS_F4_NET_MOD_node_weights_and_pvals.txt",
                                 sep="\t", check.names=F, header=T)

#####
  
### Identify most dominant ASVs in modules
# filter net module components table for ASVs
net_mods_to_asv <- net_mod_components %>%
  filter(., grepl("prok_", nodeId)) %>%
  select(Module,nodeId)

# Create long format relative abundance table
mic_asv_rel_long = mic_asv_rel %>%
  melt(., variable.name = "RAS_id", value.name = "Rel_abund")

# Identify top 10 most abundant ASVs in each module
net_mod_top10_abund_asvs =
  net_mods_to_asv %>%
  left_join(mic_asv_taxa, by=c("nodeId" = "ASV_name")) %>%
  left_join(mic_asv_rel_long, by=c("nodeId" = "ASV_name")) %>%
  aggregate(Rel_abund~nodeId+Genus+Module, data=., FUN=max) %>%
  group_by(Module) %>%
  arrange(desc(Rel_abund)) %>%
  slice_head(n=10) %>%
  as.data.frame() %>%
  mutate(Rel_abund = Rel_abund*100) %>%
  mutate(nodeId = factor(nodeId, levels = unique(nodeId))) %>%
  mutate(ASV_genus = paste0(Genus," - ",nodeId))

# Assess how many connections these ASVs have in the network to other ASVs and to
# functions within the same module. Importantly, we will divide the number of 
# edges by the number of nodes in each module

# First create table with number of nodes per module
net_mod_node_counts = 
  net_mod_components %>%
  select(Module,type) %>%
  table() %>%
  as.data.frame()

# Now create a table containing the top ten most abundant ASVs and the number
# of connections they have to other ASVs and functions. Also normalise the 
# number of connections based upon the number of components in the respective module
net_mod_top10_abund_asv_connections = 
  net_mod_connections %>%
  filter(., from %in% net_mod_top10_abund_asvs$nodeId) %>%
  left_join(net_mod_components, by=c("to" = "nodeId"))  %>%
  rename(., Module_to = Module) %>%
  select(from,Module_to,type) %>%
  mutate(count = 1) %>%
  aggregate(count~from+Module_to+type, data=., FUN=sum) %>%
  left_join(net_mod_top10_abund_asvs, by=c("from" = "nodeId")) %>%
  rename(., "Module_from" = "Module") %>%
  rename(., nodeId = from) %>%
  select(nodeId,Genus,ASV_genus,Rel_abund,Module_from,Module_to,type,count) %>%
  arrange(Module_from,desc(Rel_abund)) %>%
  left_join(net_mod_node_counts, by=c("Module_to" = "Module", "type")) %>%
  mutate(norm_count = count/Freq) %>%
  select(-Freq)

# Filter connections to include only those within the same module
net_mod_top10_abund_asv_connections_within_mod = 
  net_mod_top10_abund_asv_connections %>%
  filter(Module_from == Module_to) %>%
  mutate(ASV_genus = factor(ASV_genus, levels=unique(ASV_genus)))

### Visualise
module_colours = c("M1" = "#008066",
                   "M2" = "#55DDFF",
                   "M3" = "#FFB300",
                   "M4" = "#F55BF5",
                   "M5" = "#490092",
                   "M6" = "#920000")

net_mod_top10_asv_max_abund_barplot = 
  ggplot(data=net_mod_top10_abund_asv_connections_within_mod, 
       aes(x = Rel_abund, y = ASV_genus)) + 
  geom_bar(aes(fill = Module_from), stat="identity", position="identity", show.legend=F) +
  scale_fill_manual(values = module_colours) + 
  scale_y_discrete(limits=rev) + 
  facet_wrap(Module_from~., scales="free_y", ncol = 1, nrow = 6) +
  labs(x = "Maximum relative abundance (%)") + 
  theme_bw() + 
  theme(strip.background.x = element_blank(),
        strip.text.x = element_blank(),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12))

net_mod_top10_num_connections_with_asvs_barplot = 
  net_mod_top10_abund_asv_connections_within_mod %>%
  filter(., type == "ASV") %>%
  ggplot(data=., aes(x = norm_count, y = ASV_genus)) + 
  geom_bar(fill="#D9DDDC", stat="identity", position="identity", show.legend=F) + 
  facet_wrap(Module_from~., scales="free_y", ncol = 1, nrow = 6) + 
  theme_bw() + 
  labs(x = "Normalised count of connections with other ASVs in module") + 
  scale_y_discrete(limits=rev) + 
  theme(strip.background.x = element_blank(),
        strip.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())

net_mod_top10_num_connections_with_func_barplot = 
  net_mod_top10_abund_asv_connections_within_mod %>%
  filter(., type == "FUNC") %>%
  ggplot(data=., aes(x = norm_count, y = ASV_genus)) + 
  geom_bar(fill="#787276", stat="identity", position="identity", show.legend=F) + 
  facet_wrap(Module_from~., scales="free_y", ncol = 1, nrow = 6) + 
  theme_bw() + 
  labs(x = "Normalised count of connections with functions in module") + 
  scale_y_discrete(limits=rev) + 
  theme(strip.background.x = element_blank(),
        strip.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())

### Merge figures and export
pdf(file = paste0(output_figures,"/FRAM_RAS_F4_MOD_most_abund_asvs_barplot.pdf"),
    height = 10, width = 10)
net_mod_top10_asv_max_abund_barplot+net_mod_top10_num_connections_with_asvs_barplot+
  net_mod_top10_num_connections_with_func_barplot+plot_layout(widths=c(4,2,2))
dev.off()

###

### Identify the top 10 most abundant genera in each module

# Sum ASV abundances at genus level for each sample in each module
# Then identify the ten genera with the highest relative abundance values in 
# each module
net_mods_most_abund_genera = net_mods_to_asv %>%
  left_join(mic_asv_taxa, by=c("nodeId" = "ASV_name")) %>%
  left_join(mic_asv_rel, by=c("nodeId" = "ASV_name")) %>%
  melt(., variable.name="RAS_id", value.name="Rel_abund") %>%
  select(Module, Genus, RAS_id, Rel_abund) %>%
  aggregate(Rel_abund~Genus+Module+RAS_id, data=., FUN=sum) %>%
  aggregate(Rel_abund~Genus+Module, data=., FUN=max) %>%
  group_by(Module) %>%
  arrange(desc(Rel_abund)) %>%
  slice_head(n=10) %>%
  as.data.frame() %>%
  mutate(Rel_abund = Rel_abund*100) %>%
  arrange(Module,desc(Rel_abund)) %>%
  mutate(Genus = factor(Genus, levels = unique(Genus)))

# Plot maximum relative abundance of genera as barchart
# First need to use interaction between Genus and Module to enforce the correct
# ordering of Genus in the faceted plot
net_mods_most_abund_genera_new <- 
  net_mods_most_abund_genera %>%
  arrange(Module, desc(Rel_abund)) %>%
  mutate(Genus_Module = factor(interaction(Genus, Module), levels = unique(interaction(Genus, Module))))

net_mod_most_abund_genera_bar_facet = 
  ggplot(net_mods_most_abund_genera_new, 
       aes(x = Rel_abund, y = Genus_Module, fill = Module)) +
  geom_bar(stat="identity", position="identity") + 
  scale_fill_manual(values = module_colours) + 
  facet_grid(Module ~ ., scales = "free_y", switch = "y") +
  labs(x = "Maximum relative abundance (%)") + 
  theme_bw() + 
  theme(strip.background.y = element_blank(),
        strip.text.y = element_blank(),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12))

pdf(file = paste0(output_figures,"/FRAM_RAS_F4_MOD_most_abund_genera_bar.pdf"),
    height = 10, width = 6)
net_mod_most_abund_genera_bar_facet
dev.off()




################

### Figure 5b - module functional composition

###############

### KEGG metabolism composition

# The primary comparison between module functional composition is based on
# those functions annotated by the KEGG database. In particular, it is 
# focused on KEGG pathway maps.

# As we are leveraging the functional annotations from eggnog orthologs,
# there are instances where a functional cluster is associated with multiple KEGG KO annotations.
# For the coarse grained analysis on KEGG pathway maps, we will consider all of
# the possible KO's, however, ensuring that no single functional cluster is assigned
# to two different KEGG ko's within the same pathway.
# This will allow us to assess the differences in counts of KO's associated with
# a pathway across modules. But we can not use this approach to make quantitative assessments
# between pathways only within pathways across modules.

# The module kegg functional profiles have already been provided to you

### Import Module KEGG KO compositions

# Import FuncID to functional information
func_id_to_func_to_net_mod = 
  read.table(file="RAS_F4_GENES_CLUSTID_to_FUNCID_to_FUNC_to_MOD.txt",
             sep="\t", header=T, check.names=F)

# Import table containing kegg composition information
net_mod_kegg_comp=read.csv(
  file = "RAS_F4_NET_MOD_kegg_functional_composition.txt",
  sep="\t",check.names=F, header=T)

### To gain an overview on differences in metabolic features between modules,
### we will focus on energy metabolisms as well as a selection of those with high
### variance between modules

# Create table with list of kEGG pathways in Energy metabolism
kegg_energy_metabolism_names = 
  net_mod_kegg_comp %>%
  filter(., grepl("Energy metabolism",KEGG_pathway_map)) %>%
  filter(., !grepl("Oxidative",KEGG_module) & 
           !grepl("antenna",KEGG_module) & 
           !grepl("in photosynth",KEGG_module)) %>%
  select(KEGG_module) %>%
  unique()

# Now sum up the number of Functions attributed to all other KEGG pathways and identify those
# with the highest variance between modules. Then print list of those and combine
# with the list of energy metabolism pathways
kegg_mods_high_var = 
  net_mod_kegg_comp %>%
  filter(., !grepl("Energy metabolism",KEGG_pathway_map)) %>%
  filter(., !grepl("Glycolysis",KEGG_module)) %>%
  select(Module,KEGG_module) %>%
  table() %>%
  as.data.frame() %>%
  group_by(KEGG_module) %>%
  mutate(max = max(Freq)) %>%
  mutate(min = min(Freq)) %>%
  mutate(range = max-min) %>%
  as.data.frame() %>%
  select(KEGG_module,range) %>%
  unique() %>%
  arrange(desc(range)) %>%
  slice_head(n=9) %>%
  as.data.frame() %>%
  select(KEGG_module) %>%
  unique()

# Combine the two lists of metabolisms of interest
kegg_modules_of_interest = rbind(kegg_mods_high_var,kegg_energy_mods)

# Determine the number of functions in each metabolism
kegg_modules_of_interest_counts = 
  net_mod_kegg_comp %>%
  select(Module,KEGG_module) %>%
  table() %>%
  as.data.frame() %>%
  filter(KEGG_module %in% kegg_modules_of_interest$KEGG_module) %>%
  mutate(KEGG_module = factor(KEGG_module, levels=unique(kegg_modules_of_interest$KEGG_module)))

###

### Carbohydrate-active enzyme (CAZymes) composition

# We also supplement the functional comparison between modules based on 
# the composition of CAZymes

# Filter out CAZyme annotations from the functional profile
net_mod_cazymes = 
  func_id_to_func_to_net_mod %>%
  filter(., grepl("^GH",FUNC) |
           grepl("^PL",FUNC) |
           grepl("^CE",FUNC) |
           grepl("^CBM",FUNC)) %>%
  filter(., !grepl("kinases",FUNC) & 
           !grepl("hydro",FUNC) & 
           !grepl("PLD",FUNC) & 
           !grepl("deac",FUNC) & 
           !grepl("_",FUNC) & 
           !grepl("DUF",FUNC) & 
           !grepl("PSD",FUNC)) %>%
  unique()

# Determine total count of CAZyme genes in each module and then combine this
# with the counts of the selected KEGG metabolisms
metabolisms_of_interest_counts = net_mod_cazymes %>%
  select(Module) %>%
  table() %>%
  as.data.frame() %>%
  mutate(KEGG_module = "CAZyme_families") %>%
  rbind(., kegg_modules_of_interest_counts) %>%
  mutate(KEGG_module = factor(KEGG_module, levels=unique(KEGG_module)))

### Create a bar plot summarising functional composition of modules 

net_mod_kegg_top_metab_modules_plot = 
  metabolisms_of_interest_counts %>%
  mutate(Module = factor(Module, levels=c("M1","M2","M3","M4","M5"))) %>%
  ggplot(.) + 
  geom_bar(aes(x = KEGG_module, y = Freq, fill = Module),
           position="dodge", stat="identity") + 
  theme_bw() + 
  facet_wrap(KEGG_module~., scales="free", nrow=3) + 
  scale_fill_manual(values=module_colours) + 
  labs(y = "Number of unique functions") + 
  theme(strip.background.x = element_rect(fill = "white", colour = "black"),
        strip.text.x = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size = 14))
net_mod_kegg_top_metab_modules_plot

# Export plot
pdf(file = paste0(output_figures,"/FRAM_RAS_F4_MOD_KEGG_top_metab_modules_barchart.pdf"),
    height = 10, width = 14)
net_mod_kegg_top_metab_modules_plot
dev.off()


################

### Figure 6 - correlations between module prokaryotic ASVs and microeukaryotic ASVs

###############

library(circlize)

### import tables

mic_asv_to_euk_asv_cor = read.table(file = "RAS_F4_MIC_and_EUK_ASV_correlations.txt",
           sep="\t", header=T, check.names=F) 

mic_asv_taxa = read.table(file = "RAS_F4_MIC_ASV_taxa.txt",
                                    sep="\t", header=T, check.names=F)

euk_asv_taxa = read.table(file = "RAS_F4_EUK_ASV_taxa.txt",
                          sep="\t", header=T, check.names=F)

net_mod_to_asv=read.table(
  file="RAS_F4_NET_MOD_assignments.txt", sep="\t",
  check.names=F, header=T) %>%
  filter(., grepl("prok_",nodeId)) %>%
  select(nodeId,Module)

net_mod_mic_asv_to_euk_corr = 
  mic_asv_to_euk_asv_cor %>%
  left_join(net_mod_to_asv, by=c("From" = "nodeId")) %>%
  filter(., Module != "NA") %>%
  filter(., Adjusted_p_value < 0.05) %>%
  left_join(euk_asv_taxa, by=c("To" = "ASV_name")) %>%
  filter(., Pearson_R > 0.7) %>%
  arrange(Module,Phylum,Class,Family)

from = paste(net_mod_mic_asv_to_euk_corr[[5]])
to = paste(net_mod_mic_asv_to_euk_corr[[8]])
mat = matrix(0, nrow = length(unique(from)), ncol = length(unique(to)))
rownames(mat) = unique(from)
colnames(mat) = unique(to)
for(i in seq_along(from)) mat[from[i], to[i]] = 1

# set cluster colours for plots
module_colours = c("M1" = "#008066",
                   "M2" = "#55DDFF",
                   "M3" = "#FFB300",
                   "M4" = "#F55BF5",
                   "M5" = "#490092",
                   "M6" = "#920000")

genus_cols_tmp = as.data.frame(colnames(mat)) %>%
  mutate(col = "grey")
genus_cols <- genus_cols_tmp$col
names(genus_cols) <- genus_cols_tmp$`colnames(mat)`

chord_cols_all <- c(module_colours,genus_cols)

circos.clear()
circos.par(start.degree = 0)
pdf(file=paste0(output_figures,"/FRAM_RAS_F4_NET_MOD_bac_ASVs_euk_corr_chord_diagram.pdf"),
    height = 12, width = 12)
chordDiagram(mat, order = union(rownames(mat), colnames(mat)), 
             big.gap = 25, directional = TRUE, grid.col=chord_cols_all,
             annotationTrack = "grid",
             preAllocateTracks = list(track.height = 0.4))

circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important

dev.off()

################

### Figure 7 - genetic diversity within functions

###############

### Compare variations in diversity of functions over time based on gene clusters

library(data.table)
library(doParallel)
library(foreach)
library(dplyr)
library(broom)
library(furrr)
library(vegan)

# Import gene clust ID to FUNC ID to FUNCTION info
clustID_to_funcID_to_FUNC_to_NetMod = 
  read.table(file="RAS_F4_GENES_CLUSTID_to_FUNCID_to_FUNC_to_MOD.txt",
             sep="\t", header=T, check.names=F)

# Import gene cluster raw count data
ras_mg_geneclust_raw_abund=read.table(
  file="RAS_F4_GENE_CLUSTID_filt_raw_wide.txt",
  sep="\t",check.names=F, header=T)

# Import sample name to RAS_id mapping file
RAS_id_to_date = read.table(file="RAS_F4_META.txt", sep="\t",
                                    check.names=F, header=T) %>%
  select(RAS_id,date) %>%
  mutate(date = format(as.Date(date, format="%d.%m.%Y")))

# Import function relative abundance table
mic_clust_func_rel_abund_long = read.table(
  file="RAS_F4_MIC_OSC4_FUNC_CLUST_rel_abund_wide.txt",
  check.names=F, sep="\t") %>%
  tibble::rownames_to_column(., var="nodeId") %>%
  reshape2::melt(id.vars="nodeId", variable.name="RAS_id", value.name="Rel_abund")

### AIM: compare the diversity of genetic variants within functions over time
### This can provide insights into whether at times of functional gene peaks, there is
### a single population driving the pattern as well as the degree of functional redundancy

### Firstly, let's assess the number of functions comprised of single and multiple
### gene clusters
clustID_to_funcID_to_FUNC_to_NetMod %>%
  select(funcID) %>%
  table() %>%
  as.data.frame() %>%
  filter(Freq > 1) %>%
  dim()

clustID_to_funcID_to_FUNC_to_NetMod %>%
  select(funcID) %>%
  table() %>%
  as.data.frame() %>%
  filter(Freq == 1) %>%
  dim()

(2478/(2478+8842))*100

### Calculate per module proportions of singleton and multiple gene-cluster functions
mod_func_prop_singleton_and_multiple = 
  clustID_to_funcID_to_FUNC_to_NetMod %>%
  select(Module,funcID) %>%
  table() %>%
  as.data.frame() %>%
  filter(Freq > 0) %>%
  mutate(Category = case_when(
    Freq < 2 ~ "Single",
    TRUE ~ "Multiple"
  )) %>%
  select(Module,Category) %>%
  table() %>%
  as.data.frame() %>%
  reshape2::dcast(Module~Category, value.var = "Freq", data=.) %>%
  tibble::column_to_rownames(., var="Module") %>%
  decostand(., method = "total", MARGIN=1) %>%
  tibble::rownames_to_column (., var="Module") %>%
  reshape2::melt(., variable.name = "Category", value.name = "Proportion")


### Visualise proportions in pie charts
mod_func_prop_singleton_and_multiple_piecharts = 
  ggplot(mod_func_prop_singleton_and_multiple, aes(x="", y=Proportion, fill=Category)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) + 
  facet_wrap(Module~., ncol = 5) + 
  theme_void() + 
  scale_fill_manual(values=c("Multiple" = "#000000",
                             "Single" = "gray")) + 
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 11),
        strip.text.x = element_text(size = 11))

### export figure
pdf(file=paste0(output_figures,"/FRAM_RAS_F4_NET_MOD_func_singleton_vs_multi_piecharts.pdf"), height=4, width=8)
mod_func_prop_singleton_and_multiple_piecharts
dev.off()

# Create dataframes containing with singleton gene cluster functions
single_cluster_funcs = clustID_to_funcID_to_FUNC_to_NetMod %>%
  select(Module,funcID,FUNC) %>%
  table() %>%
  as.data.frame() %>%
  filter(Freq == 1)

# Create dataframe with multi gene cluster functions
multi_cluster_funcs = clustID_to_funcID_to_FUNC_to_NetMod %>%
  select(Module,funcID,FUNC) %>%
  table() %>%
  as.data.frame() %>%
  filter(Freq > 1)

### Now we will explore the diversity and relative abundance in functions
### that contain multiple gene clusters

# Filter dataframes to contain only information for functions with multiple gene clusters
clustid_to_funcid_to_net_mod_multi = 
  clustID_to_funcID_to_FUNC_to_NetMod %>%
  filter(., funcID %in% multi_cluster_funcs$funcID)

func_multi_to_mod = 
  clustid_to_funcid_to_net_mod_multi %>%
  select(funcID,Module) %>%
  unique()

mic_func_rel_multi = 
  mic_func_rel %>%
  filter(., funcID %in% multi_cluster_funcs$funcID)

# Subset gene cluster raw count table to those within multiple-gene cluster functions
net_mod_func_multi_gene_cluster_raw_abund = 
  ras_mg_geneclust_raw_abund %>%
  filter(., clustID %in% clustid_to_funcid_to_net_mod_multi$clustID) %>%
  tibble::column_to_rownames(., var="clustID") %>%
  mutate_all(., ~replace_na(.,0))

# Identify lowest count across samples
net_mod_func_multi_gene_cluster_raw_abund %>%
  colSums() %>%
  sort() %>%
  head(n=2)

# clearly one sample significantly lower than others, will remove this one sample
# and rarefy to the next samllest
net_mod_func_multi_gene_cluster_raw_abund_rarefied = 
  net_mod_func_multi_gene_cluster_raw_abund %>%
  select(-"03_2019_F4_1") %>%
  t() %>%
  rrarefy(., 23000) %>%
  t() %>%
  as.data.frame()

# FUNCTION: filter rarefied gene cluster count matrix for specific function,
# then calculate shannon diversity of gene clusters within that function across
# all samples. 
# Due to this process involving calculating Shannon diversity for a large number of 
# dataframes it will be parallelized.

# Create empty dataframe to store results
result_df <- data.frame()
funcID_list <- unique(clustid_to_funcid_to_net_mod_multi$funcID)

# Set up parallel processing
cl <- makeCluster(6)
registerDoParallel(cl)

# Load necessary packages within each parallel worker
clusterEvalQ(cl, {
  library(dplyr)
  library(vegan)
  library(tibble)
})

# Use foreach to loop through all functions
result_df <- foreach (genefunction = funcID_list, .combine = 'rbind') %dopar% {
  clusters_of_interest <- clustid_to_funcid_to_net_mod_multi[clustid_to_funcid_to_net_mod_multi$funcID == genefunction, ] %>%
    select(funcID,clustID)
  temp <- net_mod_func_multi_gene_cluster_raw_abund_rarefied[rownames(net_mod_func_multi_gene_cluster_raw_abund_rarefied) %in% 
                                                               clustid_to_funcid_to_net_mod_multi$clustID, , drop = FALSE] %>%
    vegan::diversity(., index="shannon", MARGIN=2) %>%
    as.data.frame(., nm="Diversity") %>%
    tibble::rownames_to_column(., var="RAS_id") %>%
    mutate(nodeId = genefunction)
  print(paste(genefunction, "finished"))
  return(temp)
}

# end parallelization
stopCluster(cl)

# Filter for those with Diversity values (after rarefaction, some functions
# won't have abundance information and thus no diversity)
func_with_div_values = result_df %>%
  aggregate(Diversity~nodeId, FUN=sum, data=.) %>%
  filter(., Diversity > 0) %>%
  select(nodeId)

result_df_filt = result_df[result_df$nodeId %in% func_with_div_values$nodeId, ]

### Export this dataframe
write.table(result_df_filt,
            file=paste0(output_tables,"/FRAM_RAS_F4_FUNC_shannon_diversity_per_func.txt"),
            sep="\t")

### Now investigate linear relationship between shannon diversity and module functions
mic_func_div_and_abund_long = mic_clust_func_rel_abund_long %>%
  filter(., nodeId %in% clustid_to_funcid_to_net_mod_multi$funcID) %>%
  left_join(result_df_filt, by=c("nodeId","RAS_id")) %>%
  filter(., Diversity != "NA")

# Create list of funcID's to loop over
cl <- makeCluster(4)
registerDoParallel(cl)

# Function to perform linear regression for each unique nodeId
perform_regression <- function(nodeId, data) {
  test <- data[data$nodeId == nodeId, ]
  
  result <- tryCatch({
    lm_result <- lm(Rel_abund ~ Diversity, data = test)
    
    if (!is.null(lm_result) && class(lm_result) == "lm") {
      p_value <- round(summary(lm_result)$coefficients[2, 4], 5)
      adj_rsquared <- round(summary(lm_result)$adj.r.squared, 3)
      
      return(data.frame(nodeId = nodeId, p_value = p_value, adj_rsquared = adj_rsquared))
    } else {
      cat("Skipping regression for nodeId:", nodeId, "due to an error in linear regression\n")
      return(NULL)
    }
  }, error = function(e) {
    cat("Skipping regression for nodeId:", nodeId, "due to an error in linear regression\n")
    return(NULL)
  })
  
  return(result)
}# Get unique nodeIds
unique_nodeIds <- unique(mic_func_div_and_abund_long$nodeId)

# Use foreach to apply the function in parallel
results_df <- future_map_dfr(unique_nodeIds, ~perform_regression(.x, mic_func_div_and_abund_long))
results_df

# Stop parallel processing
stopCluster(cl)

# Perform multiple testing correction with the Benjamini-Hochberg procedure
results_df$p_adjusted <- p.adjust(results_df$p_value, method = "BH")

# Filter to retain only significant hits and add in module information
net_mod_func_signif_lm_df = 
  results_df %>%
  filter(., p_adjusted < 0.05) %>%
  left_join(func_multi_to_mod, by=c("nodeId" = "funcID")) %>%
  unique()

# Create a dataframe that indicates the total number of multi-gene cluster functions
# per module 
net_mod_total_func_multi_count = 
  clustid_to_funcid_to_net_mod_multi %>%
  select(funcID,Module) %>%
  unique() %>%
  select(Module) %>%
  table() %>%
  as.data.frame() %>%
  rename(total_func_count = Freq)

net_mod_func_signif_lm_summary =
  net_mod_func_signif_lm_df %>%
  select(Module) %>%
  table() %>%
  as.data.frame() %>%
  rename(Significant = Freq) %>%
  left_join(net_mod_total_func_multi_count, by="Module") %>%
  mutate(Not_significant = total_func_count-Significant) %>%
  select(-total_func_count) %>%
  reshape2::melt(id.vars="Module", variable.name="Type", value.name="Count") 

net_mod_func_signif_lm_summary_prop = 
  net_mod_func_signif_lm_summary %>%
  reshape2::dcast(Module~Type, value.var = "Count", data=.) %>%
  tibble::column_to_rownames(., var="Module") %>%
  decostand(., method="total", MARGIN=1) %>%
  tibble::rownames_to_column(., var="Module") %>%
  reshape2::melt(id.vars="Module", variable.name="Type", value.name="Proportion")

# Create file with nodeId of those that showed significant linear relationship
# along with an identifier (will be used for colouring in figure)
signif_func_with_label = 
  net_mod_func_signif_lm_df %>%
  select(nodeId) %>%
  mutate(type = "Significant")

# Create dataframe with diveristy and relative abundance info for functions
# that showed significant linear relationship
mic_func_div_and_abund_with_labels = 
  mic_func_div_and_abund_long %>%
  left_join(func_multi_to_mod, by=c("nodeId" = "funcID")) %>%
  left_join(signif_func_with_label, by="nodeId") %>%
    mutate(type = case_when(
      type == "Significant" ~ "Significant",
      TRUE ~ "Not_significant"
    )) %>%
  arrange(type) %>%
  mutate(type = factor(type, levels=unique(type)))

### Create figure
module_colours = c("M1" = "#008066",
                   "M2" = "#55DDFF",
                   "M3" = "#FFB300",
                   "M4" = "#F55BF5",
                   "M5" = "#490092",
                   "M6" = "#920000")


library(patchwork)

create_abund_vs_div_scatter <- function(div_vs_abund_data, modulename){
    
    scatter_plot = div_vs_abund_data %>%
      filter(Module == modulename) %>%
      ggplot(.,
             aes(x = Rel_abund, y = Diversity, colour = Module)) + 
      geom_point(stat="identity", show.legend=F) + 
      scale_colour_manual(values = module_colours) + 
      labs(x = "Relative abundance (%)",
           y = "Shannon diversity") + 
      theme_bw() + 
      theme(axis.title.y = element_text(size = 12),
            axis.title.x = element_text(size = 12),
            axis.text.x = element_text(size = 11),
            axis.text.y = element_text(size = 11))
    return(scatter_plot)
      
}

create_abund_vs_div_sig_bar <- function(signif_data, modulename){
  
  signif_barplot = signif_data %>%
    filter(Module == modulename) %>%
    mutate(Type = factor(Type, levels=c("Not_significant","Significant"))) %>%
    ggplot(., aes(x = Module, y = Proportion, fill = Type)) + 
    scale_fill_manual(values=c("Not_significant" = "#D9DDDC",
                               "Significant" = "#48494B")) + 
    geom_bar(stat="identity", position="stack") + 
    labs(y = "Proportion of functions with multiple gene clusters") + 
    theme_bw() + 
    theme(axis.title.y = element_text(size = 12),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 11))
  
  return(signif_barplot)
  
}
      
for (modulename in c("M1","M2","M3","M4","M5")){
  scatter_fig_out <- paste(modulename, "_functions_rel_abund_and_diversity_scatter", sep = "")
  barplot_fig_out <- paste(modulename, "_functions_rel_abund_and_diversity_signif_bar", sep = "")
  
  assign(scatter_fig_out, create_abund_vs_div_scatter(mic_func_div_and_abund_with_labels,
                                                 modulename))
  assign(barplot_fig_out, create_abund_vs_div_sig_bar(net_mod_func_signif_lm_summary_prop,
                                                 modulename))
}

## Export plots (will export separately and join into figure later as the plots are 
## large and their composite nature renders combining them challenging)
pdf(file = paste0(output_figures,"/FRAM_RAS_F4_NET_MOD_ALL_func_abund_vs_diversity_scatter.pdf"),
    height = 8, width = 10)
M1_functions_rel_abund_and_diversity_scatter+M2_functions_rel_abund_and_diversity_scatter+M3_functions_rel_abund_and_diversity_scatter+
  M4_functions_rel_abund_and_diversity_scatter+M5_functions_rel_abund_and_diversity_scatter+plot_spacer()+plot_layout(ncol=3,nrow=2)
dev.off()

pdf(file = paste0(output_figures,"/FRAM_RAS_F4_NET_MOD_ALL_func_abund_vs_diversity_sig_bars.pdf"),
    height = 8, width = 9)
M1_functions_rel_abund_and_diversity_signif_bar+M2_functions_rel_abund_and_diversity_signif_bar+M3_functions_rel_abund_and_diversity_signif_bar+
  M4_functions_rel_abund_and_diversity_signif_bar+M5_functions_rel_abund_and_diversity_signif_bar+plot_spacer()+plot_layout(ncol=3,nrow=2)
dev.off()

## Export table with significant lm results after correction to prevent 
## recalculation later 
write.table(net_mod_func_signif_lm_df,
            file=paste0(output_tables,"/FRAM_RAS_F4_NET_MOD_FUNC_shannon_div_vs_abund_lm_sig_results.txt"),
            sep="\t")

### Beyond exploring whether the diversity shifts over time, we are also interested
### to find out whether the same gene variants dominate the abundance of each 
### function each year, i.e. is the same species/close relative dominating that function
### recurrently? 

### In order to not be too stringent and accommodate for the fact that modules
### typically peak over more than a single sample, we will first identify
### the top 3 samples for each module for each year. These are the samples in which
### the functional genes of the module reach highest abundances during the year 

# Import the combined functional gene abundance information for modules
net_mod_func_rel_abund = 
  read.table(file="RAS_F4_NET_MOD_FUNC_relative_abundance.txt", 
             sep="\t", header=T)

# Import the functional cluster abundance information for modules
mic_func_rel = 
  read.table(file="RAS_F4_MIC_OSC4_FUNC_CLUST_rel_abund_wide.txt", 
             sep="\t", header=T, check.names = F)

# Import network module assignments
net_mod_components = 
  read.table(file="RAS_F4_NET_MOD_assignments.txt", 
             sep="\t", header=T)

# Create funcID to FUNC mapping
funcID_to_func_mapping =
  clustID_to_funcID_to_FUNC_to_NetMod %>%
  select(funcID,FUNC)

# identify top three samples per year based on abundance - here we only focus
# on the complete julian years (2017,2018,2019) to ensure the peak each year is
# covered
mod_sample_peaks = net_mod_func_rel_abund %>%
  mutate(year = case_when(
    RAS_id = grepl("2016",RAS_id) ~ "2016",
    RAS_id = grepl("2017",RAS_id) ~ "2017",
    RAS_id = grepl("2018",RAS_id) ~ "2018",
    RAS_id = grepl("2019",RAS_id) ~ "2019",
    RAS_id = grepl("2020",RAS_id) ~ "2020",
  )) %>%
  filter(., grepl("2017|2018|2019", year)) %>%
  group_by(Module,year) %>%
  arrange(desc(Rel_abund)) %>%
  slice_head(n=3) %>%
  as.data.frame() %>%
  select(Module,RAS_id,year) %>%
  mutate(mod_and_peak = paste(Module,"-",RAS_id,sep=""))

# Import gene cluster abundance table
mic_gene_clust_rel = read.table(file="RAS_F4_MIC_GENE_CLUST_OSC4_rel_abund_wide.txt",
           sep="\t", check.names=F, header=T) %>%
  tibble::rownames_to_column(., var="clustID") %>%
  melt(., variable.name="RAS_id", value.name="gene_rel_abund")

# Import module KEGG annotation information
net_mod_kegg_composition_single_hits_with_annotations = 
  read.table(file = "Network_analysis/FRAM_RAS_F4_MOD_kegg_singlehit_composition.txt",sep="\t", check.names=F, header=T) 
net_mod_kegg_composition_all_hits_with_annotations = 
  read.table(file = "Network_analysis/FRAM_RAS_F4_MOD_kegg_allhits_composition.txt",
             sep="\t",check.names=F, header=T)

# FUNCTION: The function will combine the information from the functional cluster
# and gene cluster abundances, then filter for the three samples per year identified 
# as the peak oscillation period (above) and then retain the gene cluster that 
# represents the largest proportion of each functions abundance in those samples.
# The result is a single gene cluster for each function on each peak sample date in each
# annual cycle
identify_top_gene_clust_for_func_in_peak_samples <- 
  function(gene_to_func,mod_components,func_rel,gene_rel,mod_peaks,module_name) {
    func_names = gene_to_func %>%
      filter(., Module == module_name) %>%
      select(funcID) %>%
      table() %>%
      as.data.frame() %>%
      filter(., Freq > 1) %>%
      select(funcID)
    
    temp1 = gene_to_func %>%
      filter(., funcID %in% func_names$funcID) %>%
      select(funcID,clustID)
    
    temp2 = func_rel %>%
      filter(., funcID %in% func_names$funcID)
    
    top_gene_clust_hits = gene_rel %>%
      filter(., clustID %in% temp1$clustID) %>%
      left_join(temp1, by="clustID") %>%
      left_join(temp2, by=c("funcID","RAS_id")) %>%
      mutate(func_proportion = (gene_rel_abund/func_rel_abund)*100) %>%
      mutate(Module = module_name) %>%
      mutate(mod_and_peak = paste(Module,"-",RAS_id,sep="")) %>%
      filter(., mod_and_peak %in% mod_peaks$mod_and_peak) %>%
      left_join(mod_peaks, by=c("mod_and_peak","Module","RAS_id")) %>%
      select(-gene_rel_abund,-func_rel_abund,mod_and_peak) %>%
      filter(., func_proportion > 0) %>%
      group_by(funcID,RAS_id,Module) %>%
      arrange(desc(func_proportion)) %>%
      slice_head(n=1) %>%
      as.data.frame() %>%
      arrange(funcID,year)
    return(top_gene_clust_hits)
  }

# FUNCTION: How many times did the same gene cluster appear across the three years
# across peak samples?
calc_freq_of_top_gene_clust_recurring <- function(mod_peak_gene_clust,funcID_to_func_mapping,module_name){
  mod_peak_gene_clust %>%
    select(clustID,funcID,Module,year) %>%
    unique() %>%
    group_by(clustID,funcID,Module) %>%
    summarise(count = n()) %>%
    as.data.frame() %>%
    mutate(category = case_when(
      count > 2 ~ "Same gene cluster",
      TRUE ~ "Different gene cluster"
    )) %>%
    select(Module,funcID,category) %>%
    group_by(funcID,Module) %>%
    arrange(funcID,Module,desc(category)) %>%
    slice_head(n=1) %>%
    as.data.frame() %>%
    left_join(funcID_to_func_mapping, by="funcID") %>%
    unique()
}

tmp = clustID_to_funcID_to_FUNC_to_NetMod %>%
  filter(., Module == "M1") %>%
  select(funcID) %>%
  table() %>%
  as.data.frame() %>%
  filter(., Freq > 1) %>%
  select(funcID)
head(mic_func_rel)
clustID_to_funcID_to_FUNC_to_NetMod %>%
  filter(., funcID %in% tmp$funcID)
  
# Perform a loop to run over the above two functions for each module and produce
# output tables with information
for (modulename in c("M1","M2","M3","M4","M5")){
  mod_peaks = mod_sample_peaks %>%
    filter(., Module == modulename)
  out1 <- paste(modulename, "_top_gene_clust_per_peak", sep = "")
  assign(out1, identify_top_gene_clust_for_func_in_peak_samples(clustID_to_funcID_to_FUNC_to_NetMod,
                                                                net_mod_components,
                                                                mic_func_rel,
                                                                mic_gene_clust_rel,
                                                                mod_peaks,
                                                                modulename))
  out2 <- paste(modulename, "_freq_top_gene_clust_recurrence", sep = "")
  assign(out2, calc_freq_of_top_gene_clust_recurring(get(out1),funcID_to_func_mapping,modulename))
}

# Combine the output tables from all modules
module_top_gene_clust_per_peak = 
  rbind(M1_top_gene_clust_per_peak,M2_top_gene_clust_per_peak,
        M3_top_gene_clust_per_peak,M4_top_gene_clust_per_peak,
        M5_top_gene_clust_per_peak) %>%
  unique()

module_top_gene_cluster_per_function_recurrence = 
  rbind(M1_freq_top_gene_clust_recurrence,M2_freq_top_gene_clust_recurrence,
      M3_freq_top_gene_clust_recurrence,M4_freq_top_gene_clust_recurrence,
      M5_freq_top_gene_clust_recurrence) %>%
  unique()

write.table(module_top_gene_clust_per_peak, 
            file = "Network_analysis/FRAM_RAS_F4_NET_MOD_multi_func_top_gene_clusters_per_peak.txt",
            sep="\t")

write.table(module_top_gene_cluster_per_function_recurrence, 
            file = "Network_analysis/FRAM_RAS_F4_NET_MOD_multi_func_same_or_different_gene_cluster.txt",
            sep="\t")

### Combine information on different or same gene cluster with output from
### linear regression
temp = net_mod_func_signif_lm_df %>%
  select(nodeId,adj_rsquared) %>%
  left_join(module_top_gene_cluster_per_function_recurrence, by=c("nodeId" = "funcID"))
View(temp)
module_top_gene_cluster_per_function_recurrence %>%
  select(funcID,category)
module_top_gene_cluster_per_function_recurrence %>%
  left_join(net_mod_func_signif_lm_df, by=c("Module","funcID" = "nodeId")) %>%
  filter(., p_value > 0)

net_mod_func_signif_lm_df %>%
  left_join(module_top_gene_cluster_per_function_recurrence, by=c("Module","nodeId" = "funcID"))
module_top_gene_cluster_per_function_recurrence
module_top_gene_cluster_per_function_recurrence %>%
  select(category) %>%
  table() %>%
  as.data.frame()
1113+1325
(1113/(1113+1325))*100

# Calculate the relative proportion of multi-gene cluster functions for which
# a significant positive linear relationship between diversity and abundance was
# observed
module_top_gene_cluster_per_function_recurrence_proportion  = 
  module_top_gene_cluster_per_function_recurrence %>%
  select(Module,category) %>%
  table() %>%
  as.data.frame() %>%
  reshape2::dcast(Module~category, value.var = "Freq", data=.) %>%
  tibble::column_to_rownames(., var="Module") %>%
  decostand(., method = "total", MARGIN=1) %>%
  tibble::rownames_to_column (., var="Module") %>%
  reshape2::melt(., variable.name = "category", value.name = "Proportion")
module_top_gene_cluster_per_function_recurrence_proportion

# Visualise the results about dominant gene cluster recurrence across years
# set cluster colours for plots
module_colours = c("M1" = "#008066",
                   "M2" = "#55DDFF",
                   "M3" = "#FFB300",
                   "M4" = "#F55BF5",
                   "M5" = "#490092",
                   "M6" = "#920000")

module_top_gene_cluster_freq_per_function_per_year = 
  ggplot(module_top_gene_cluster_per_function_recurrence_proportion,
       aes(x =  Module, y = Proportion)) + 
  geom_bar(aes(fill = category), position="stack", stat="identity") + 
  scale_fill_manual(values = c("Different gene cluster" = "#BEBDB8",
                               "Same gene cluster" = "#490092")) + 
  facet_wrap(.~Module, scales="free_x", ncol=5) + 
  labs(y = "Proportion of gene functions") + 
  theme_bw() + 
  theme(strip.background.x = element_blank(),
        strip.text.x = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 12))

### export
pdf(file=paste0(output_figures,"/FRAM_RAS_F4_MOD_FUNC_gene_variant_recurrence_barplot.pdf"),
    height = 6, width = 8)
module_top_gene_cluster_freq_per_function_per_year
dev.off()









###
### Explore the data further with a focus on those functions that do have
### the same recurrent, dominant gene variant. What are these functions? Do they
### represent specialised traits?

module_top_gene_cluster_per_function_recurrence = read.table(
  file = "Network_analysis/FRAM_RAS_F4_NET_MOD_multi_func_same_or_different_gene_cluster.txt",
  sep="\t", check.names=F, header=T)

# Import module KEGG annotation information
net_mod_kegg_composition = 
  read.table(file = "Network_analysis/FRAM_RAS_F4_NET_MOD_kegg_functional_composition.txt",sep="\t", check.names=F, header=T)  %>%
  select(Module,FUNC,KEGG_category,KEGG_pathway_map,KEGG_module,Gene_description) %>%
  arrange(Module,KEGG_category,KEGG_pathway_map,KEGG_module,Gene_description)

# Extract functional clusters dominated by same gene cluster each year
# and that were annotated by KEGG database. Then combine with KEGG metabolic category
# information
same_clust_kegg_functions = 
  module_top_gene_cluster_per_function_recurrence %>%
  filter(., grepl("Same", category)) %>%
  mutate(FUNC = gsub("ko:","",FUNC)) %>%
  separate_rows(., FUNC, sep=",") %>%
  as.data.frame() %>%
  left_join(net_mod_kegg_composition, by=c("FUNC","Module")) %>%
  unique() %>%
  filter(., KEGG_category != "NA") %>%
  arrange(Module,KEGG_category,KEGG_pathway_map,KEGG_module,Gene_description,FUNC)
  
diff_clust_kegg_functions = 
  module_top_gene_cluster_per_function_recurrence %>%
  filter(., grepl("Different", category)) %>%
  mutate(FUNC = gsub("ko:","",FUNC)) %>%
  separate_rows(., FUNC, sep=",") %>%
  as.data.frame() %>%
  left_join(net_mod_kegg_composition, by=c("FUNC","Module")) %>%
  unique() %>%
  arrange(Module,KEGG_category,KEGG_pathway_map,KEGG_module,Gene_description,FUNC)
