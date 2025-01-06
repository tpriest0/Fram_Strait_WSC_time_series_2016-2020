
####################

#### Processing network clusters and visualise cluster dynamics over time
#### series

####################

### Define working directory
setwd('XXXXX')

### Define and create output directories
output_figures = ('results/figures/')
output_tables = ('results/tables/')

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
library(tidytext)
library(ggdendro)
library(dendextend)
library(igraph)
library(tidyr)
library(ggstatsplot)

options(scipen=999)

### Import data

# Import Prokaryotic ASV data
mic_asv_rel=read.table(file=paste0(output_tables,"FRAM_RAS_F4_PROK_ASV_filt_rel.txt"), sep="\t",
                       check.names=F, header=T, row.names=1) 

mic_asv_raw=read.table(file=paste0(output_tables,"FRAM_RAS_F4_PROK_ASV_filt_raw.txt"), sep="\t",
                       check.names=F, header=T, row.names=1) 

# Import Eukaryotic ASV data
euk_asv_rel=read.table(file=paste0(output_tables,"FRAM_RAS_F4_EUK_ASV_filt_rel.txt"), sep="\t",
                       check.names=F, header=T, row.names=1) 

euk_asv_raw=read.table(file=paste0(output_tables,"FRAM_RAS_F4_EUK_ASV_filt_raw.txt"), sep="\t",
                       check.names=F, header=T, row.names=1) 

# Import gene cluster ID to function to function ID
gene_to_func_to_tax=read.table(
  file="FRAM_RAS_F4_GENE_CLUST_to_FUNC_to_taxonomy.txt",
  sep="\t",check.names=F, header=T)

# Import gene cluster abundance profile
func_osc4_rel_wide = 
  read.table(file=paste0(output_tables,"FRAM_RAS_F4_FUNC_osc4_rel_wide.txt"),sep="\t",
             check.names=F, header=T) 

# Import function abundance profile
gene_osc4_raw_wide = 
  read.table(file=paste0(output_tables,"FRAM_RAS_F4_GENE_CLUST_filt_raw_wide.txt"),sep="\t",
             check.names=F, header=T)

# Import metadata
sample_meta = 
  read.table(file="FRAM_RAS_F4_META.txt",sep="\t",
             check.names=F, header=T)

# Import eukaryotic ASV relative abundance data
euk_asv_taxa=read.table(file=paste0(output_tables,"FRAM_RAS_F4_EUK_ASV_filt_taxa.txt"), sep="\t",
                       check.names=F, header=T)

# Import microbial ASV taxa info
mic_asv_taxa=read.table(file=paste0(output_tables,"FRAM_RAS_F4_PROK_ASV_filt_taxa.txt"), sep="\t",
                       check.names=F, header=T)

# Import ASV and gene network cluster assignments
net_mod_components=
  read.table(file="FRAM_RAS_F4_NET_mod_components.txt",
                    sep="\t", check.names=F, header=T)
                    
# Import module node connections
net_mod_connections = read.table(file="FRAM_RAS_F4_NET_mod_connections.txt",
                                 sep="\t", check.names=F, header=T) %>%
                                 mutate(from = gsub("bac_","prok_",from)) %>%
                                 mutate(to = gsub("bac_","prok_",to)) 

#####

### Define palettes

#####

# set cluster colours for plots
module_colours = c("M1" = "#008066",
                    "M2" = "#55DDFF",
                    "M3" = "#FFB300",
                    "M4" = "#F55BF5",
                    "M5" = "#490092",
                    "M6" = "#920000")


#####

### Reformat and generate necessary tables

#####

### Change module numbers to reflect temporal appearance in annual cycle
net_mod_components_mod = net_mod_components %>%
  mutate(Module = case_when(
    cluster == 3 ~ "M1",
    cluster == 2 ~ "M2",
    cluster == 1 ~ "M3",
    cluster == 5 ~ "M4",
    cluster == 4 ~ "M5"
  )) %>%
  select(nodeId,Module,type) %>%
  mutate(type = case_when(
        grepl("euk",type) ~ "Euk",
        grepl("bac",type) ~ "Prok",
        grepl("func",type) ~ "FUNC"
    ))

write.table(net_mod_components_mod,
            file = paste0(output_tables,"FRAM_RAS_F4_NET_mod_components-renamed.txt"),
            sep="\t", row.names=F, quote=F)

### Combine funcID to network module with information on gene clusters to funcID
net_mod_to_gene_to_func_to_tax = 
  gene_to_func_to_tax %>%
  left_join(net_mod_components_mod, by=c("funcID" = "nodeId"))

write.table(net_mod_to_gene_to_func_to_tax,
            file = paste0(output_tables,"FRAM_RAS_F4_NET_mod_to_gene_to_func_to_taxonomy.txt"),
            sep="\t", row.names=F, quote=F)

#####

### Assess if sequencing depth has an impact on the diversity of components in modules

#####

### FUNCTION: rarefy the ASV count table 100 times and each time assess 
### the number of ASVs detected in each module
calculate_alpha_div_asv <- function(abund_table, module_mapping){
  
  smallest_count <- min(colSums(abund_table))

  results_list <- list()
  
  for (i in 1:100) {
    rarefied_abund <- t(rrarefy(t(abund_table), smallest_count))
    
    module_richness_df = rarefied_abund %>%
      as.data.frame() %>%
      tibble::rownames_to_column(., var="nodeId") %>%
      left_join(module_mapping, by="nodeId") %>%
      select(Module) %>%
      filter(., Module != "NA") %>%
      mutate(Count = 1) %>%
      aggregate(Count~Module, data=., FUN=sum)
    
    results_list[[i]] <- module_richness_df
  }
  
  combined_results <- bind_rows(results_list)
  
  summary_stats <- combined_results %>%
    group_by(Module) %>%
    summarise(
      Richness_mean = mean(Count),
      Richness_SD = sd(Count)
    )
  
  return(summary_stats)
}


### Apply above function to prokaryotic ASV table
module_asv_richness_results_df = calculate_alpha_div_asv(mic_asv_raw, net_mod_components_mod) %>%
  as.data.frame() %>%
  mutate(type = "Prok_ASV")

### Apply above function to eukaryotic ASV table
module_euk_asv_richness_results_df = calculate_alpha_div_asv(euk_asv_raw, net_mod_components_mod) %>%
  as.data.frame() %>%
  mutate(type = "Euk_ASV")

### Repeat the same for gene clusters

### FUNCTION: rarefy the gene cluster count table 100 times and each time,
### reassign gene clusters to functions and assess the number of functions detected
### for each module
calculate_alpha_div_func <- function(abund_table, clust_to_func_mapping){
  
  smallest_count <- min(colSums(abund_table))

  results_list <- list()
  
  for (i in 1:100) {
    rarefied_abund <- t(rrarefy(t(abund_table), smallest_count))
    
    module_richness_df = rarefied_abund %>%
        as.data.frame() %>%
        tibble::rownames_to_column(., var="Gene_clustID") %>%
        reshape2::melt(id.vars="Gene_clustID", variable.name="RAS_id", value.name="Count") %>%
        left_join(clust_to_func_mapping, by="Gene_clustID") %>%
        select(funcID,Module,RAS_id,Count) %>%
        aggregate(Count~funcID+Module+RAS_id, FUN=sum, data=.) %>%
        select(funcID,Module) %>%
        unique() %>%
        mutate(Count = 1) %>%
        aggregate(Count~Module, data=., FUN=sum)
    
    results_list[[i]] <- module_richness_df
  }
  
  combined_results <- bind_rows(results_list)
  
  summary_stats <- combined_results %>%
    group_by(Module) %>%
    summarise(
      Richness_mean = mean(Count),
      Richness_SD = sd(Count)
    )
  
  return(summary_stats)
}

# Run above function
module_func_richness_results_df = calculate_alpha_div_func(gene_osc4_raw_wide, net_mod_to_gene_to_func_to_tax) %>%
  as.data.frame() %>%
  mutate(type = "FUNC")

### Combine richness results from ASVs and Functions
module_comp_rare_richness = rbind(module_asv_richness_results_df,module_euk_asv_richness_results_df,module_func_richness_results_df)
write.table(module_comp_rare_richness, file=paste0(output_tables,"FRAM_RAS_F4_NET_mod_rarefied_richness.txt"), sep="\t", quote=FALSE, row.names=FALSE)

#####

### Assess if the number of sampling time points has an impact on the diversity of components within modules

#####

# To do this, we will select one sample from each sampling month and then assess diversity of module components detected
subset_samples = sample_meta %>%
    arrange(month) %>%
    group_by(month) %>%
    slice_head(n=1) %>%
    as.data.frame()

subset_mic_asv_raw = mic_asv_raw %>%
    select(all_of(subset_samples$RAS_id)) %>%
    mutate(RowSum = ifelse(rowSums(.) != 0, 1, 0)) %>%
      filter(., RowSum != 0) %>%
      select(-RowSum)

subset_euk_asv_raw = euk_asv_raw %>%
    select(all_of(subset_samples$RAS_id)) %>%
    mutate(RowSum = ifelse(rowSums(.) != 0, 1, 0)) %>%
      filter(., RowSum != 0) %>%
      select(-RowSum)

subset_gene_osc4_raw_wide = gene_osc4_raw_wide %>%
    select(all_of(subset_samples$RAS_id)) %>%
    mutate(RowSum = ifelse(rowSums(.) != 0, 1, 0)) %>%
      filter(., RowSum != 0) %>%
      select(-RowSum)

### Reapply rarefying and richness function from before
module_subset_asv_richness_results_df = calculate_alpha_div_asv(subset_mic_asv_raw, net_mod_components_mod) %>%
  as.data.frame() %>%
  mutate(type = "Prok_ASV")

module_subset_euk_asv_richness_results_df = calculate_alpha_div_asv(subset_euk_asv_raw, net_mod_components_mod) %>%
  as.data.frame() %>%
  mutate(type = "Euk_ASV")

module_subset_func_richness_results_df = calculate_alpha_div_func(subset_gene_osc4_raw_wide, net_mod_to_gene_to_func_to_tax) %>%
  as.data.frame() %>%
  mutate(type = "FUNC")

### Combine richness results from ASVs and Functions
module_subset_comp_rare_richness = rbind(module_subset_asv_richness_results_df,module_subset_euk_asv_richness_results_df,module_subset_func_richness_results_df)
write.table(module_subset_comp_rare_richness, file=paste0(output_tables,"FRAM_RAS_F4_NET_MOD_subset_and_rarefied_richness.txt"), sep="\t", quote=FALSE, row.names=FALSE)

#####

### Compare diversity of modules components from original data, rarefied data and reduced sampling+rarefied data

#####

net_mod_components_counts = net_mod_components_mod %>%
    mutate(type = case_when(
        type == "Prok" ~ "Prok_ASV",
        type == "Euk" ~ "Euk_ASV",
        type == "FUNC" ~ "FUNC"
    )) %>%
    select(Module,type) %>%
    table() %>%
    as.data.frame() %>%
    rename(Richness_mean = Freq) %>%
    mutate(Category = "Original data") 

module_comp_rare_richness_mod = module_comp_rare_richness %>%
    select(-Richness_SD) %>%
    mutate(Category = "Rarefied data")

module_subset_comp_rare_richness_mod = module_subset_comp_rare_richness %>%
    select(-Richness_SD) %>%
    mutate(Category = "Reduced samples and rarefied data")
  
net_mod_components_counts_all_variations = rbind(net_mod_components_counts,module_comp_rare_richness_mod,module_subset_comp_rare_richness_mod) %>%
    mutate(type = case_when(
        type == "Prok_ASV" ~ "ASV",
        type == "Euk_ASV" ~ "ASV",
        TRUE ~ "FUNC"
    ))

net_mod_components_counts_all_variations_barplot = ggplot(net_mod_components_counts_all_variations, aes(x = Module, y = Richness_mean)) + 
    geom_bar(aes(fill = type), stat = "identity", position = "stack") + 
    facet_grid(Category~.) + 
    labs(y = "Number of components") +
    scale_fill_manual(values = c("ASV" = "#004949", "FUNC" = "#b66dff")) + 
    theme_bw() + 
    theme(axis.text.x = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 12),
        strip.background.y = element_rect(fill = "white", colour = "black"),
        strip.text.y = element_text(size = 10, angle=0),
        legend.position = "bottom",
        legend.title = element_blank())

### Export plot
pdf(file=paste0(output_figures,"FRAM_RAS_F4_NET_MOD_components_after_rarefying_and_reducing_samples_barplot.pdf"), height=8, width=8)
net_mod_components_counts_all_variations_barplot
dev.off()

#####

### Compare domain-level composition of module ASVs and functions

#####

### Domain-resolved module ASV counts

# Raw counts
net_mod_asv_counts = net_mod_components_mod %>%
    filter(., type == "Prok" | 
    type == "Euk") %>%
    left_join(mic_asv_taxa, by=c("nodeId" = "ASV_name")) %>%
    mutate(Domain = case_when(
        Kingdom == "Bacteria" ~ "Bacteria",
        Kingdom == "Archaea" ~ "Archaea",
        TRUE ~ "Eukarya"
    )) %>%
    select(Module,Domain) %>%
    mutate(count = 1) %>%
    aggregate(count~Domain+Module, data=., FUN=sum)

net_mod_asv_comp_raw_count_barplot = 
    ggplot(net_mod_asv_counts, aes(x = Module, y = count)) + 
    geom_bar(aes(fill = Domain), stat = "identity", position = "stack", show.legend=F) + 
    labs(y = "Number of components") +
    scale_fill_manual(values = c("Bacteria" = "#009292", "Archaea" = "#004949", "Eukarya" = "#ffb6db")) + 
    theme_bw() + 
    theme(axis.text.x = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 12))

# Relative proportion of total module ASV counts
net_mod_asv_comp_rel_long =
    net_mod_asv_counts %>%
    dcast(Module~Domain, value.var="count") %>%
    tibble::column_to_rownames(., var="Module") %>%
    mutate(across(where(is.numeric), ~replace_na(., 0))) %>%
    decostand(., method="total", MARGIN=1) %>%
    tibble::rownames_to_column(., var="Module") %>%
    melt(., id.vars=c("Module"), variable.name="Domain", value.name="Relative_proportion")

net_mod_asv_comp_rel_barplot = ggplot(net_mod_asv_comp_rel_long, aes(x = Module, y = Relative_proportion)) + 
    geom_bar(aes(fill = Domain), stat = "identity", position = "stack", show.legend=F) + 
    labs(y = "Relative proportion of Module ASVs") +
    scale_fill_manual(values = c("Bacteria" = "#009292", "Archaea" = "#004949", "Eukarya" = "#ffb6db")) + 
    theme_bw() + 
    theme(axis.text.x = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 12))

### Domain-level composition of module functions

# Raw counts
net_mod_func_domain_comp_counts_long = net_mod_to_gene_to_func_to_tax %>%
    select(Module,Domain_tiara) %>%
    mutate(Domain = case_when(
        Domain_tiara == "archaea" ~ "Archaea",
        Domain_tiara == "bacteria" ~ "Bacteria",
        Domain_tiara == "prokarya" ~ "Undefined",
        Domain_tiara == "organelle" ~ "Undefined",
        Domain_tiara == "unknown" ~ "Undefined",
        TRUE ~ "Eukarya"
    )) %>%
    mutate(count = 1) %>%
    aggregate(count~Domain+Module, data=., FUN=sum)

net_mod_func_domain_comp_raw_barplot = 
    ggplot(net_mod_func_domain_comp_counts_long, aes(x = Module, y = count)) + 
    geom_bar(aes(fill = Domain), stat = "identity", position = "stack", show.legend=F) + 
    labs(y = "Number of gene clusters within functions") +
    scale_fill_manual(values = c("Bacteria" = "#009292", "Archaea" = "#004949", "Eukarya" = "#ffb6db", "Undefined" = "#D9DDDC")) + 
    theme_bw() + 
    theme(axis.text.x = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 12))

# Relative proportion of total module Function counts
net_mod_func_domain_comp_rel_long = net_mod_func_domain_comp_counts_long %>%
    dcast(Module~Domain, value.var="count") %>%
    tibble::column_to_rownames(., var="Module") %>%
    decostand(., method="total", MARGIN=1) %>%
    tibble::rownames_to_column(., var="Module") %>%
    melt(., id.vars=c("Module"), variable.name="Domain", value.name="Relative_proportion")

net_mod_func_domain_comp_rel_barplot = 
    ggplot(net_mod_func_domain_comp_rel_long, aes(x = Module, y = Relative_proportion)) + 
    geom_bar(aes(fill = Domain), stat = "identity", position = "stack") + 
    labs(y = "Relative proportion of gene clusters in Module functions") +
    scale_fill_manual(values = c("Bacteria" = "#009292", "Archaea" = "#004949", "Eukarya" = "#ffb6db", "Undefined" = "#D9DDDC")) + 
    theme_bw() + 
    theme(axis.text.x = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 12))

### Combine plots and export
pdf(file=paste0(output_figures,"FRAM_RAS_F4_NET_MOD_ASV_and_FUNC_domain_composition_raw_count_barplot.pdf"), height=6, width=8)
net_mod_asv_comp_raw_count_barplot+net_mod_func_domain_comp_raw_barplot
dev.off()

pdf(file=paste0(output_figures,"FRAM_RAS_F4_NET_MOD_ASV_and_FUNC_domain_composition_relative_barplot.pdf"), height=6, width=8)
net_mod_asv_comp_rel_barplot+net_mod_func_domain_comp_rel_barplot
dev.off()

######

### Visualize the temporal dynamics of modules

#####

### Determine the combined relative abundance of module components
func_module_rel_abund_long = 
  func_osc4_rel_wide %>%
  filter(funcID %in% net_mod_components_mod$nodeId) %>%
  left_join(net_mod_components_mod, by=c("funcID" = "nodeId")) %>%
  reshape2::melt(id.vars=c("funcID","Module","type"),
                 variable.name = "Sample_name", value.name="Rel_abund") %>%
  aggregate(Rel_abund~Module+Sample_name, data=., FUN=sum) %>%
  left_join(sample_meta, by=c("Sample_name" = "MG_sample_names"))%>%
  mutate(date = as.Date(date, "%d.%m.%Y")) %>%
  mutate(Rel_abund = Rel_abund*100) %>%
  mutate(Type = "FUNC") %>%
  select(Module,RAS_id,Rel_abund,Type,date,year,month,depth,daylight,temp,sal,O2_conc,Chl_a,CO2,pH,PAR_satellite,MLD)

mic_asv_module_rel_abund_long = 
    mic_asv_rel %>%
    tibble::rownames_to_column(., var="nodeId") %>%
    filter(nodeId %in% net_mod_components_mod$nodeId) %>%
    left_join(net_mod_components_mod, by="nodeId") %>%
    reshape2::melt(id.vars=c("nodeId","Module","type"),
                 variable.name = "RAS_id", value.name="Rel_abund") %>%
    aggregate(Rel_abund~Module+RAS_id, data=., FUN=sum) %>%
    left_join(sample_meta, by="RAS_id")%>%
    mutate(date = as.Date(date, "%d.%m.%Y")) %>%
    mutate(Rel_abund = Rel_abund*100) %>%
    mutate(Type = "PROK_ASV") %>%
    select(Module,RAS_id,Rel_abund,Type,date,year,month,depth,daylight,temp,sal,O2_conc,Chl_a,CO2,pH,PAR_satellite,MLD)

euk_asv_module_rel_abund_long = 
    euk_asv_rel %>%
    tibble::rownames_to_column(., var="nodeId") %>%
    filter(nodeId %in% net_mod_components_mod$nodeId) %>%
    left_join(net_mod_components_mod, by="nodeId") %>%
    reshape2::melt(id.vars=c("nodeId","Module","type"),
                 variable.name = "RAS_id", value.name="Rel_abund") %>%
    aggregate(Rel_abund~Module+RAS_id, data=., FUN=sum) %>%
    left_join(sample_meta, by="RAS_id")%>%
    mutate(date = as.Date(date, "%d.%m.%Y")) %>%
    mutate(Rel_abund = Rel_abund*100) %>%
    mutate(Type = "EUK_ASV") %>%
    select(Module,RAS_id,Rel_abund,Type,date,year,month,depth,daylight,temp,sal,O2_conc,Chl_a,CO2,pH,PAR_satellite,MLD)

asv_module_rel_abund_long = rbind(mic_asv_module_rel_abund_long,euk_asv_module_rel_abund_long)

### Visualise relative abundance dynamics of modules

module_temporal_dynamics = ggplot(data=asv_module_rel_abund_long,
         aes(x=as.Date(date), y=Rel_abund, fill=Type)) + 
    geom_area(aes(fill = Type), position="stack", colour="black", alpha=0.5, show.legend=F) + 
    geom_line(data=func_module_rel_abund_long, aes(x=as.Date(date), y=Rel_abund*5, colour = Module), 
              linewidth=1, show.legend=F) +
    geom_point(data=func_module_rel_abund_long, aes(x=as.Date(date), y=Rel_abund*5), shape = 24, size = 2, 
               show.legend=F) +
    scale_x_date(date_breaks = "6 months", date_labels =  "%b %Y") +
    scale_y_continuous(sec.axis = 
                         sec_axis(~./5, name="Cluster genes relative abundance (%)")) + 
    #scale_fill_manual(values = module_colours) + 
    scale_colour_manual(values = module_colours) + 
    facet_grid(Module~., scales="free_y") + 
    labs(y = "Relative abundance (%)") + 
    theme_bw() + 
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          plot.background = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank())

### Combine all plots and export
library(patchwork)
pdf(file=paste0(output_figures,"FRAM_RAS_F4_NET_mod_temporal_dynamics.pdf"),
    height=10, width=8)
module_temporal_dynamics
dev.off()

### Export module ASV and gene abundance tables
write.table(func_module_rel_abund_long,
            file=paste0(output_tables,"FRAM_RAS_F4_NET_mod_FUNC_relative_abundance.txt"), sep="\t",
            quote = F, row.names = F)

write.table(mic_asv_module_rel_abund_long,
            file=paste0(output_tables,"FRAM_RAS_F4_NET_mod_ASV_relative_abundance.txt"), sep="\t",
            quote = F, row.names = F)


temp1= func_module_rel_abund_long %>%
  select(Module,RAS_id,Rel_abund) %>%
  arrange(Module,RAS_id)
temp1

### ALTERNATIVE VISUALISATION:
### area graphs to better illustrate successional dynamics of module ASVs and module
### functions 

pdf(file=paste0(output_figures,"FRAM_RAS_F4_NET_mod_ASV_only_temporal_dynamics_area_nonstacked.pdf"),
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

pdf(file=paste0(output_figures,"FRAM_RAS_F4_NET_mod_function_only_temporal_dynamics.pdf"),
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

###############

##### What are the factors driving module dynamics?

# First let's filter our metadata table to continuous variables and zero-mean variance
# standardize them

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

pdf(file=paste0(output_figures,"FRAM_RAS_F4_MOD_env_corr_heatmap.pdf"),
    width = 8, height = 10)
net_mod_env_sig_corr_heatmap
dev.off()

### Export correlation table
write.table(module_asv_func_vs_env_sig_cor_df,
            file=paste0(output_tables,"FRAM_RAS_F4_NET_mod_ASV_FUNC_vs_env_cor_results.txt"), sep="\t")


#######################

###### Examine the content of modules
###### in terms of microbial and functional composition

#######################

### MODULE TAXONOMIC COMPOSITION

### FUNCTION: Create dataframe with the 10 ASVs reaching the highest relative abundances in each module

identify_top10_asvs <- function(net_mod_components_mod, domain, asv_rel_abund, asv_taxa){
  net_mods_to_asv = net_mod_components_mod %>%
    filter(., grepl(domain, nodeId)) %>%
    select(Module,nodeId)
  
  asv_rel_long = asv_rel_abund %>%
    tibble::rownames_to_column(., var="ASV_name") %>%
    melt(., variable.name = "RAS_id", value.name = "Rel_abund")

  net_mod_top10_abund_asvs =
    net_mods_to_asv %>%
    left_join(asv_taxa, by=c("nodeId" = "ASV_name")) %>%
    left_join(asv_rel_long, by=c("nodeId" = "ASV_name")) %>%
    aggregate(Rel_abund~nodeId+Genus+Module, data=., FUN=max) %>%
    group_by(Module) %>%
    arrange(desc(Rel_abund)) %>%
    slice_head(n=10) %>%
    as.data.frame() %>%
    mutate(Rel_abund = Rel_abund*100) %>%
    mutate(ASV_genus = paste0(Genus," - ",nodeId)) %>%
    arrange(Module,desc(Rel_abund)) %>%
    mutate(ASV_genus = factor(ASV_genus, levels=unique(ASV_genus)))

  return(net_mod_top10_abund_asvs)
}

net_mod_top10_mic_asvs = identify_top10_asvs(net_mod_components_mod, "prok", mic_asv_rel, mic_asv_taxa)
net_mod_top10_euk_asvs = identify_top10_asvs(net_mod_components_mod, "euk", euk_asv_rel, euk_asv_taxa)

### FUNCTION: create a dataframe containing information on the number of connections that the top10 asvs have within the module, normalised by the 
### total number of connections in that module

identify_top10_asvs_connections <- function(net_mod_components_mod, net_mod_connections, top10_asvs){
  net_mod_node_counts = 
    net_mod_components_mod %>%
    select(Module,type) %>%
    table() %>%
    as.data.frame()
  
  p1 = net_mod_connections %>%
    filter(., from %in% top10_asvs$nodeId) %>%
    rename(node_one = from) %>%
    rename(node_two = to)

  p2 = net_mod_connections %>%
    filter(., to %in% top10_asvs$nodeId) %>%
    rename(node_one = to) %>%
    rename(node_two = from)
  
  net_mod_top10_abund_asv_connections = 
    rbind(p1,p2) %>%
    select(node_one,node_two) %>%
    unique() %>%
    left_join(net_mod_components_mod, by=c("node_two" = "nodeId")) %>%
    rename(., Module_to = Module) %>%
    select(node_one,Module_to,type) %>%
    mutate(count = 1) %>%
    aggregate(count~node_one+Module_to+type, data=., FUN=sum) %>%
    left_join(top10_asvs, by=c("node_one" = "nodeId")) %>%
    filter(., Module == Module_to) %>%
    rename(., nodeId = node_one) %>%
    select(nodeId,Module,ASV_genus,Rel_abund,type,count) %>%
    left_join(net_mod_node_counts, by=c("Module", "type")) %>%
    mutate(norm_count = count/Freq) %>%
    select(-Freq) %>%
    arrange(Module,desc(Rel_abund)) %>%
    mutate(ASV_genus = factor(ASV_genus, levels=unique(ASV_genus)))
  
  return(net_mod_top10_abund_asv_connections)
}

net_mod_top10_mic_asvs_connections <- identify_top10_asvs_connections(net_mod_components_mod, net_mod_connections, net_mod_top10_mic_asvs)
net_mod_top10_euk_asvs_connections <- identify_top10_asvs_connections(net_mod_components_mod, net_mod_connections, net_mod_top10_euk_asvs)

### FUNCTION: Create figure summarising most abundant asvs and their connections to other asvs and functions in the module
visualise_top10_asvs_and_connections <- function(top10_asvs, top10_asvs_connections){
     top10_asv_max_abund = 
          ggplot(data=top10_asvs, 
               aes(x = Rel_abund, y = ASV_genus)) + 
          geom_bar(aes(fill = Module), stat="identity", position="identity", show.legend=F) +
          scale_fill_manual(values = module_colours) + 
          scale_y_discrete(limits=rev) + 
          facet_wrap(Module~., scales="free_y", ncol = 1, nrow = 6) +
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
     
     top10_asv_prok_connections = 
          top10_asvs_connections %>%
          filter(., type == "Prok") %>%
          ggplot(data=., aes(x = norm_count, y = ASV_genus)) + 
          geom_bar(fill="#D9DDDC", stat="identity", position="identity", show.legend=F) + 
          facet_wrap(Module~., scales="free_y", ncol = 1, nrow = 6) + 
          theme_bw() + 
          labs(x = "Normalised num. of connections with prokaryotic ASVs in module") + 
          scale_y_discrete(limits=rev) + 
          theme(strip.background.x = element_blank(),
               strip.text.x = element_blank(),
               axis.text.y = element_blank(),
               axis.title.y = element_blank())

     top10_asv_euk_connections = 
          top10_asvs_connections %>%
          filter(., type == "Euk") %>%
          ggplot(data=., aes(x = norm_count, y = ASV_genus)) + 
          geom_bar(fill="#D9DDDC", stat="identity", position="identity", show.legend=F) + 
          facet_wrap(Module~., scales="free_y", ncol = 1, nrow = 6) + 
          theme_bw() + 
          labs(x = "Normalised num. of connections with eukaryotic ASVs in module") + 
          scale_y_discrete(limits=rev) + 
          theme(strip.background.x = element_blank(),
               strip.text.x = element_blank(),
               axis.text.y = element_blank(),
               axis.title.y = element_blank())
     
     top10_asv_func_connections = 
          top10_asvs_connections %>%
          filter(., type == "FUNC") %>%
          ggplot(data=., aes(x = norm_count, y = ASV_genus)) + 
          geom_bar(fill="#787276", stat="identity", position="identity", show.legend=F) + 
          facet_wrap(Module~., scales="free_y", ncol = 1, nrow = 6) + 
          theme_bw() + 
          labs(x = "Normalised num. of connections with functions in module") + 
          scale_y_discrete(limits=rev) + 
          theme(strip.background.x = element_blank(),
               strip.text.x = element_blank(),
               axis.text.y = element_blank(),
               axis.title.y = element_blank())
     
     return(list(
         max_abund_plot = top10_asv_max_abund,
         prok_connections_plot = top10_asv_prok_connections,
         euk_connections_plot = top10_asv_euk_connections,
         func_connections_plot = top10_asv_func_connections
     ))
}

### Run function to illustrate number of connectins for top10 most abundant prokaryotic ASVs in each module, then combine figures and export
net_mod_top10_asvs_prok_plots = visualise_top10_asvs_and_connections(net_mod_top10_mic_asvs, net_mod_top10_mic_asvs_connections)

pdf(file = paste0(output_figures,"FRAM_RAS_F4_NET_mod_prok_most_abund_asvs_barplot.pdf"),
    height = 8, width = 10)
net_mod_top10_asvs_prok_plots$max_abund_plot+net_mod_top10_asvs_prok_plots$prok_connections_plot+
  net_mod_top10_asvs_prok_plots$euk_connections_plot+net_mod_top10_asvs_prok_plots$func_connections_plot+
  plot_layout(widths=c(4,2,2,2))
dev.off()

### Run function to illustrate number of connectins for top10 most abundant eukaryotic ASVs in each module, then combine figures and export
net_mod_top10_asvs_euk_plots = visualise_top10_asvs_and_connections(net_mod_top10_euk_asvs, net_mod_top10_euk_asvs_connections)

pdf(file = paste0(output_figures,"FRAM_RAS_F4_NET_mod_euk_most_abund_asvs_barplot.pdf"),
    height = 8, width = 10)
net_mod_top10_asvs_euk_plots$max_abund_plot+net_mod_top10_asvs_euk_plots$prok_connections_plot+
  net_mod_top10_asvs_euk_plots$euk_connections_plot+net_mod_top10_asvs_euk_plots$func_connections_plot+
  plot_layout(widths=c(4,2,2,2))
dev.off()

### Create figure just illustrating relative abundance of top10 asvs of prokaryotes and eukaryotes
pdf(file = paste0(output_figures,"FRAM_RAS_F4_NET_mod_prok_and_euk_most_abund_asvs_barplot.pdf"),
    height = 8, width = 10)
net_mod_top10_asvs_prok_plots$max_abund_plot+net_mod_top10_asvs_euk_plots$max_abund_plot
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

### Import table containing kegg composition information
# The below file can not be provided as KEGG is a paid database. The file should contain information mapping KEGG KO's to metabolisms and pathways.
kegg_db=read.csv(
  file = "XXXXXX.txt",
  sep="\t",check.names=F, header=T)

### Add kegg metabolism information to module functions
net_mod_func_kegg_comp = 
    gene_to_func_to_tax %>%
    mutate(FUNC = gsub("ko:","",FUNC)) %>%
    left_join(kegg_db, by=c("FUNC" = "KEGG_ko")) %>%
    left_join(net_mod_components_mod, by=c("funcID" = "nodeId")) %>%
    select(-Gene,-funcID,-Gene_clustID,-GTDB_taxonomy,-type) %>%
    unique()

### To gain an overview on differences in metabolic features between modules,
### we will focus on energy metabolisms as well as a selection of those with high
### variance between modules

# Create table with list of kEGG pathways in Energy metabolism
kegg_energy_metabolism_names = 
  net_mod_func_kegg_comp %>%
  filter(., grepl("Energy metabolism",KEGG_pathway_map)) %>%
  filter(., !grepl("Oxidative",KEGG_module),
             !grepl("antenna",KEGG_module)) %>%
  select(KEGG_module) %>%
  unique()

# Now sum up the number of Functions attributed to all other KEGG pathways and identify those
# with the highest variance between modules. Then print list of those and combine
# with the list of energy metabolism pathways
kegg_mods_high_var = net_mod_func_kegg_comp %>%
  filter(., !grepl("Energy metabolism",KEGG_pathway_map)) %>%
  filter(., !grepl("Glycolysis",KEGG_module)) %>%
  select(Module,funcID,KEGG_category,KEGG_pathway_map,KEGG_module) %>%
  unique() %>%
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
  slice_head(n=10) %>%
  as.data.frame() %>%
  select(KEGG_module) %>%
  unique()

# Combine the two lists of metabolisms of interest
kegg_modules_of_interest = rbind(kegg_mods_high_var,kegg_energy_metabolism_names)

# Determine the number of functions in each metabolism
kegg_modules_of_interest_counts = 
  net_mod_func_kegg_comp %>%
  select(Module,funcID,KEGG_category,KEGG_pathway_map,KEGG_module) %>%
  unique() %>%
  select(Module,KEGG_module)
  table() %>%
  as.data.frame() %>%
  filter(KEGG_module %in% kegg_modules_of_interest$KEGG_module) %>%
  mutate(KEGG_module = factor(KEGG_module, levels=unique(kegg_modules_of_interest$KEGG_module)))

### Carbohydrate-active enzyme (CAZymes) composition

# We also supplement the functional comparison between modules based on 
# the composition of CAZymes

# Filter out CAZyme annotations from the functional profile
net_mod_cazymes = 
  net_mod_to_gene_to_func_to_tax %>%
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
  unique() %>%
  select(-Gene_clustID,-GTDB_taxonomy,-Gene,-type,-Domain_tiara)

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

# Export plot
pdf(file = paste0(output_figures,"FRAM_RAS_F4_NET_mod_KEGG_top_metab_modules_barchart.pdf"),
    height = 10, width = 14)
net_mod_kegg_top_metab_modules_plot
dev.off()

################

### Figure 6 and Supplementary Figure 4 - transcription of module functions across Arctic Ocean

################


###

### Import data

# Import metatranscriptomic abundances from Tara Oceans Arctic dataset
gene_clust_tara_abund=read.csv(
  file="FRAM_RAS_F4_GENE_CLUST_arctic_tara_abund.txt",
  sep="\t",check.names=F, header=T)

# Import TaraOcean metadata
tara_meta=read.csv(
  file="Tara_oceans_sampling_metadata_all.txt",
  sep="\t",check.names=F, header=T)

###

### Reformat and subset data

tara_meta_subset =
    tara_meta %>%
    mutate(Month = case_when(
        grepl("-05-",`Event_Date/Time_End`) ~ "May",
        grepl("-06-",`Event_Date/Time_End`) ~ "Jun",
        grepl("-07-",`Event_Date/Time_End`) ~ "Jul",
        grepl("-08-",`Event_Date/Time_End`) ~ "Aug",
        grepl("-09-",`Event_Date/Time_End`) ~ "Sep",
        grepl("-10-",`Event_Date/Time_End`) ~ "Oct",
    )) %>%
    filter(., grepl("0.22-3",Sample_name)) %>%
    select(Sample_name,Sample_station,Depth_layer,Latitude_Start,Longitude_Start,Month)

###

### Calculate mean depth of transcription at module level

gene_clust_module_tara_abund_mean_wide = gene_clust_tara_abund %>%
    filter(., TAD80 > 0) %>%
    left_join(net_mod_to_gene_to_func_to_tax, by=c("geneID" = "Gene_clustID")) %>%
    aggregate(TAD80~funcID+Tara_sample+Module, data=., FUN=sum) %>%
    aggregate(TAD80~Tara_sample+Module, data=., FUN=mean) %>%
    left_join(tara_meta_subset, by=c("Tara_sample" = "Sample_name")) %>%
    reshape2::dcast(Depth_layer+Tara_sample+Latitude_Start+Longitude_Start+Month+Sample_station~Module, value.var="TAD80") %>%
    filter(., !grepl("ZZZ",Tara_sample)) %>%
    filter(., !grepl("IZZ",Tara_sample))

###

### Assess transcription of module functions across Arctic

# Reformat coordinates of sample
coords_utm = gene_clust_module_tara_abund_mean_wide %>% 
  sf::st_as_sf( coords = c("Longitude_Start", "Latitude_Start") ) %>%
  sf::st_set_crs( 4326 ) %>%  
  sf::st_transform( 3995 ) %>%  
  sf::st_coordinates() %>%
  as.data.frame()

# Subset data to only include surface water samples for visualising on map
gene_clust_module_tara_abund_mean_wide_srf_utm = cbind(gene_clust_module_tara_abund_mean_wide,coords_utm) %>%
    filter(., Depth_layer == "SRF")

# use rnaturalearth package to get coastline data in the sf format
world_sf <- ne_countries(scale = "medium", returnclass = "sf")
map_crs <- st_crs(3995) # WGS 84 / Arctic Polar Stereographic

# Visualise mean module transcript depth on map
mod_func_transc_arctic_map = ggplot() + geom_sf(data = world_sf %>% st_transform(map_crs)) +
  coord_sf(datum = map_crs,  
           ylim = c(-3e6,3e6),  
           xlim = c(-3e6,3e6)) + 
  geom_scatterpie(data = gene_clust_module_tara_abund_mean_wide_srf_utm, aes(x = X, y = Y, group = Tara_sample), cols = c("M1","M2","M3","M4","M5"),
  pie_scale = 3) + 
  scale_fill_manual(values=module_colours) + 
  theme_bw()

# Visualise mean module transcript depth across surface, DCM and mesopelagic layers in a barplot
mod_func_transc_arctic_barplot = gene_clust_module_tara_abund_mean_wide %>%
    select(-Latitude_Start,-Longitude_Start) %>%
    reshape2::melt(., id.vars=c("Tara_sample","Month","Depth_layer","Sample_station"), variable.name = "Module", value.name = "TAD80") %>%
    aggregate(TAD80~Tara_sample+Module+Month+Depth_layer+Sample_station, data=., FUN=mean) %>%
    mutate(Depth_layer = factor(Depth_layer, levels=c("SRF","DCM","MES"))) %>%
    #mutate(Month = factor(Month, levels=c("May","Jun","Jul","Aug","Sep","Oct"))) %>%
    ggplot(., aes(x = Sample_station, y = TAD80)) +  
    geom_bar(aes(fill = Module), stat = "identity", position = "stack") + 
    facet_grid(Depth_layer~., scales="free_y") + 
    labs(y = "Truncated average depth of transcripts") +
    scale_fill_manual(values = module_colours) + 
    theme_bw() + 
    theme(axis.text.x = element_text(size = 12, angle=45, hjust=1),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 12),
        strip.background.y = element_rect(fill = "white", colour = "black"),
        strip.text.y = element_text(size = 12))

# Export figures
pdf(file=paste0(output_figures,"FRAM_RAS_F4_NET_mod_FUNC_transcript_tara_map.pdf"), height=10, width=12)
mod_func_transc_arctic_map
dev.off()
pdf(file=paste0(output_figures,"FRAM_RAS_F4_NET_mod_FUNC_transcript_tara_barplot.pdf"), height=8, width=10)
mod_func_transc_arctic_barplot
dev.off()

###

### Assess mean transcript depth of module functions within energy metabolism pathways


# Create long-format table of combined transcript depth at functional cluster level
gene_clust_tara_abund_long = gene_clust_tara_abund %>%
    filter(., TAD80 > 0) %>%
    left_join(net_mod_to_gene_to_func_to_tax, by=c("geneID" = "Gene_clustID")) %>%
    aggregate(TAD80~funcID+FUNC+Tara_sample+Module, data=., FUN=sum) %>%
    left_join(tara_meta_subset, by=c("Tara_sample" = "Sample_name")) %>%
    select(funcID,FUNC,Tara_sample,Module,Depth_layer,TAD80)

# Combine functional transcript depth information with KEGG energy metabolism pathways and calculate mean transcript depth per pathway
gene_clust_tara_abund_module_top_metab = gene_clust_tara_abund_long %>%
    filter(., grepl("ko:",FUNC)) %>%
    mutate(FUNC= gsub("ko:","",FUNC)) %>%
    left_join(kegg_db, by=c("FUNC" = "KEGG_ko")) %>%
    aggregate(TAD80~KEGG_module+KEGG_pathway_map+Module+Depth_layer, data=., FUN=mean)

# Visualise mean module transcript of energy metabolism pathways across depth layers
gene_clust_tara_abund_module_top_metab_barplot = gene_clust_tara_abund_module_top_metab %>%
    filter(., KEGG_pathway_map == "Energy metabolism") %>%
    filter(., Module == "M3" | Module == "M5") %>%
    mutate(Depth_layer = factor(Depth_layer, levels=c("SRF","DCM","MES"))) %>%
    #mutate(Month = factor(Month, levels=c("May","Jun","Jul","Aug","Sep","Oct"))) %>%
    ggplot(., aes(x = KEGG_module, y = TAD80)) +  
    geom_bar(aes(fill = Module), stat = "identity", position = "stack") + 
    facet_grid(Depth_layer~., scales="free_y") + 
    labs(y = "Truncated average depth of transcripts") +
    scale_fill_manual(values = module_colours) + 
    theme_bw() + 
    theme(axis.text.x = element_text(size = 12, angle=45, hjust=1),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 12),
        strip.background.y = element_rect(fill = "white", colour = "black"),
        strip.text.y = element_text(size = 12))

# Export figure
pdf(file=paste0(output_figures,"FRAM_RAS_F4_NET_mod_FUNC_transcript_tara_energy_metab_barplot.pdf"), height=8, width=10)
gene_clust_tara_abund_module_top_metab_barplot
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

### Set up tables
RAS_id_to_date = sample_meta %>%
  select(RAS_id,date) %>%
  mutate(date = format(as.Date(date, format="%d.%m.%Y")))

### AIM: compare the diversity of genetic variants within functions over time
### This can provide insights into whether at times of functional gene peaks, there is
### a single population driving the pattern as well as the degree of functional redundancy

### Firstly, let's assess the number of functions comprised of single and multiple
### gene clusters
net_mod_to_gene_to_func_to_tax %>%
  select(funcID) %>%
  table() %>%
  as.data.frame() %>%
  filter(Freq > 1) %>%
  dim()

net_mod_to_gene_to_func_to_tax %>%
  select(funcID) %>%
  table() %>%
  as.data.frame() %>%
  filter(Freq == 1) %>%
  dim()

(5554/(5554+2075))*100

### Calculate per module proportions of singleton and multiple gene-cluster functions
mod_func_prop_singleton_and_multiple = 
  net_mod_to_gene_to_func_to_tax %>%
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
pdf(file=paste0(output_figures,"FRAM_RAS_F4_NET_MOD_func_singleton_vs_multi_piecharts.pdf"), height=4, width=8)
mod_func_prop_singleton_and_multiple_piecharts
dev.off()

### Create dataframes containing single gene-cluster and multi gene-cluster functions
single_cluster_funcs = net_mod_to_gene_to_func_to_tax %>%
  select(Module,funcID,FUNC) %>%
  table() %>%
  as.data.frame() %>%
  filter(Freq == 1)

multi_cluster_funcs = net_mod_to_gene_to_func_to_tax %>%
  select(Module,funcID,FUNC) %>%
  table() %>%
  as.data.frame() %>%
  filter(Freq > 1)

### Now we will explore the diversity and relative abundance in functions
### that contain multiple gene clusters

# Filter dataframes to contain only information for functions with multiple gene clusters
clustid_to_funcid_to_net_mod_multi = 
  net_mod_to_gene_to_func_to_tax %>%
  filter(., funcID %in% multi_cluster_funcs$funcID)

func_multi_to_mod = 
  clustid_to_funcid_to_net_mod_multi %>%
  select(funcID,Module) %>%
  unique()

mic_func_rel_multi = 
  func_osc4_rel_wide %>%
  filter(., funcID %in% multi_cluster_funcs$funcID)


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
        select(funcID,Gene_clustID)
    
    temp <- gene_osc4_rare_rel[rownames(gene_osc4_rare_rel) %in% clusters_of_interest$Gene_clustID, ]
    
    # Check if temp has any rows
    if (nrow(temp) > 0) {
        temp <- temp %>%
            vegan::diversity(index = "shannon", MARGIN = 2) %>%  # Shannon diversity
            as.data.frame(nm = "Diversity") %>%  # Convert to dataframe
            tibble::rownames_to_column(var = "RAS_id") %>%  # Add rownames as a column
            mutate(funcID = genefunction)  # Add the genefunction (funcID)
        
        print(paste(genefunction, "finished"))
        
        return(temp)  # Return the result
    } else {
        # If empty, return NULL or an empty dataframe
        return(NULL)
    }
}

# end parallelization
stopCluster(cl)

# Filter for those with Diversity values (after rarefaction, some functions
# won't have abundance information and thus no diversity)
func_with_div_values = result_df %>%
  aggregate(Diversity~funcID, FUN=sum, data=.) %>%
  filter(., Diversity > 0) %>%
  select(funcID)

result_df_filt = result_df[result_df$funcID %in% func_with_div_values$funcID, ]

write.table(result_df_filt,
            file=paste0(output_figures,"FRAM_RAS_F4_FUNC_shannon_diversity_per_func.txt"),
            sep="\t", quote=F, row.names=F)

### Now investigate linear relationship between shannon diversity and module functions
mic_func_div_and_abund_long = func_osc4_rel_wide %>%
  reshape2::melt(id.vars=c("funcID"), variable.name="RAS_id", value.name="Rel_abund") %>%
  filter(., funcID %in% clustid_to_funcid_to_net_mod_multi$funcID) %>%
  left_join(result_df_filt, by=c("funcID","RAS_id")) %>%
  filter(., Diversity != "NA")

# Create list of funcID's to loop over
cl <- makeCluster(6)
registerDoParallel(cl)

# Function to perform linear regression for each unique nodeId
perform_regression <- function(funcID, data) {
  test <- data[data$funcID == funcID, ]
  
  result <- tryCatch({
    lm_result <- lm(Rel_abund ~ Diversity, data = test)
    
    if (!is.null(lm_result) && class(lm_result) == "lm") {
      p_value <- round(summary(lm_result)$coefficients[2, 4], 5)
      adj_rsquared <- round(summary(lm_result)$adj.r.squared, 3)
      
      return(data.frame(funcID = funcID, p_value = p_value, adj_rsquared = adj_rsquared))
    } else {
      cat("Skipping regression for funcId:", funcID, "due to an error in linear regression\n")
      return(NULL)
    }
  }, error = function(e) {
    cat("Skipping regression for funcId:", funcID, "due to an error in linear regression\n")
    return(NULL)
  })
  
  return(result)
}# Get unique nodeIds
unique_funcIDs <- unique(mic_func_div_and_abund_long$funcID)

# Use foreach to apply the function in parallel
lm_results_df <- future_map_dfr(unique_funcIDs, ~perform_regression(.x, mic_func_div_and_abund_long))

# Stop parallel processing
stopCluster(cl)

# Perform multiple testing correction with the Benjamini-Hochberg procedure
lm_results_df$p_adjusted <- p.adjust(lm_results_df$p_value, method = "BH")

# Filter to retain only significant hits and add in module information
net_mod_func_signif_lm_df = 
  lm_results_df %>%
  filter(., p_adjusted < 0.05) %>%
  left_join(func_multi_to_mod, by=c("funcID")) %>%
  unique()


write.table(net_mod_func_signif_lm_df,
            file=paste0(output_tables,"FRAM_RAS_F4_FUNC_diversity_vs_abund_lm_signif.txt"),
            sep="\t", quote=F, row.names=F)

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
  select(funcID) %>%
  mutate(type = "Significant")

# Create dataframe with diveristy and relative abundance info for functions
# that showed significant linear relationship
mic_func_div_and_abund_with_labels = 
  mic_func_div_and_abund_long %>%
  left_join(func_multi_to_mod, by=c("funcID")) %>%
  left_join(signif_func_with_label, by="funcID") %>%
    mutate(type = case_when(
      type == "Significant" ~ "Significant",
      TRUE ~ "Not_significant"
    )) %>%
  arrange(type) %>%
  mutate(type = factor(type, levels=unique(type)))

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
pdf(file = paste0(output_figures,"FRAM_RAS_F4_NET_mod_ALL_func_abund_vs_diversity_scatter.pdf"),
    height = 8, width = 10)
M1_functions_rel_abund_and_diversity_scatter+M2_functions_rel_abund_and_diversity_scatter+M3_functions_rel_abund_and_diversity_scatter+
  M4_functions_rel_abund_and_diversity_scatter+M5_functions_rel_abund_and_diversity_scatter+plot_spacer()+plot_layout(ncol=3,nrow=2)
dev.off()

pdf(file = paste0(output_figures,"FRAM_RAS_F4_NET_mod_ALL_func_abund_vs_diversity_sig_bars.pdf"),
    height = 8, width = 9)
M1_functions_rel_abund_and_diversity_signif_bar+M2_functions_rel_abund_and_diversity_signif_bar+M3_functions_rel_abund_and_diversity_signif_bar+
  M4_functions_rel_abund_and_diversity_signif_bar+M5_functions_rel_abund_and_diversity_signif_bar+plot_spacer()+plot_layout(ncol=3,nrow=2)
dev.off()

### Beyond exploring whether the diversity shifts over time, we are also interested
### to find out whether the same gene variants dominate the abundance of each 
### function each year, i.e. is the same species/close relative dominating that function
### recurrently? 

### In order to not be too stringent and accommodate for the fact that modules
### typically peak over more than a single sample, we will first identify
### the top 3 samples for each module for each year. These are the samples in which
### the functional genes of the module reach highest abundances during the year 

# Create funcID to FUNC mapping
funcID_to_func_mapping =
  net_mod_to_gene_to_func_to_tax %>%
  select(funcID,FUNC)

# identify top three samples per year based on abundance - here we only focus
# on the complete julian years (2017,2018,2019) to ensure the peak each year is
# covered
mod_sample_peaks = mic_func_module_rel_abund_long %>%
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


# FUNCTION: The function will combine the information from the functional cluster
# and gene cluster abundances, then filter for the three samples per year identified 
# as the peak oscillation period (above) and then retain the gene cluster that 
# represents the largest proportion of each functions abundance in those samples.
# The result is a single gene cluster for each function on each peak sample date in each
# annual cycle
identify_top_gene_clust_for_func_in_peak_samples <- 
  function(gene_to_func,func_rel,gene_rel,mod_peaks,module_name) {
    func_names = gene_to_func %>%
      filter(., Module == module_name) %>%
      select(funcID) %>%
      table() %>%
      as.data.frame() %>%
      filter(., Freq > 1) %>%
      select(funcID)
    
    temp1 = gene_to_func %>%
      filter(., funcID %in% func_names$funcID) %>%
      select(funcID,Gene_clustID,Module)
    
    temp2 = func_rel %>%
      reshape2::melt(., id.vars=c("funcID"), variable.name="RAS_id", value.name="func_rel_abund") %>%
      filter(., funcID %in% func_names$funcID)
    
    top_gene_clust_hits = gene_rel %>%
      tibble::rownames_to_column(., var="Gene_clustID") %>%
      filter(., Gene_clustID %in% temp1$Gene_clustID) %>%
      reshape2::melt(., id.vars=c("Gene_clustID"), variable.name="RAS_id", value.name="gene_rel_abund") %>%
      left_join(temp1, by="Gene_clustID") %>%
      left_join(temp2, by=c("funcID","RAS_id")) %>%
      mutate(func_proportion = (gene_rel_abund/func_rel_abund)*100) %>%
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
    select(Gene_clustID,funcID,Module,year) %>%
    unique() %>%
    group_by(Gene_clustID,funcID,Module) %>%
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
  
# Perform a loop to run over the above two functions for each module and produce
# output tables with information
for (modulename in c("M1","M2","M3","M4","M5")){
  mod_peaks = mod_sample_peaks %>%
    filter(., Module == modulename)
  out1 <- paste(modulename, "_top_gene_clust_per_peak", sep = "")
  assign(out1, identify_top_gene_clust_for_func_in_peak_samples(net_mod_to_gene_to_func_to_tax,
                                                                func_osc4_rel_wide,
                                                                gene_osc4_rel_wide,
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
            file = paste0(output_tables,"FRAM_RAS_F4_NET_mod_multi_func_top_gene_clusters_per_peak.txt"),
            sep="\t", quote=F, row.names=F)

write.table(module_top_gene_cluster_per_function_recurrence, 
            file = paste0(output_tables,"FRAM_RAS_F4_NET_mod_multi_func_same_or_different_gene_cluster.txt"),
            sep="\t", quote=F, row.names=F)

### Combine information on different or same gene cluster with output from
### linear regression

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

pdf(file=paste0(output_figures,"FRAM_RAS_F4_NET_mod_FUNC_gene_variant_recurrence_barplot.pdf"),
    height = 6, width = 8)
module_top_gene_cluster_freq_per_function_per_year
dev.off()
