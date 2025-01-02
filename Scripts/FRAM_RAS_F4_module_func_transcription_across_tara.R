### Define working directory
setwd('/nfs/cds-peta/exports/biol_micro_cds_gr_sunagawa/scratch/tpriest/projects/fram_wsc/')

### Load libraries
library(dplyr)
library(tidyr) 
library(ggplot2)
library(reshape2)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(rgdal)
library(sp)
library(scatterpie)

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

### Import data

#####

# Import module to gene cluster ID to function to function ID
net_mod_to_gene_to_func_to_tax=read.csv(
  file="data_files/FRAM_RAS_F4_NET_mod_to_gene_to_func_to_taxonomy.txt",
  sep="\t",check.names=F, header=T)

# Import metatranscriptomic abundances from Tara Oceans Arctic dataset
gene_clust_tara_abund=read.csv(
  file="output_files/TOPC_mapping/FRAM_RAS_F4_proteins_all_clust_rep_filt_TOPC_TAD80.txt",
  sep="\t",check.names=F, header=T)

# Import TaraOcean metadata
tara_meta=read.csv(
  file="data_files/Tara_oceans_sampling_metadata_all.txt",
  sep="\t",check.names=F, header=T)

# Import table containing kegg composition information
kegg_db=read.csv(
  file = "data_files/KEGG_database_complete.txt",
  sep="\t",check.names=F, header=T)

#####

### Reformat and subset data

#####

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

#####

### Calculate mean depth of transcription at module level

#####

gene_clust_module_tara_abund_mean_wide = gene_clust_tara_abund %>%
    filter(., TAD80 > 0) %>%
    left_join(net_mod_to_gene_to_func_to_tax, by=c("geneID" = "Gene_clustID")) %>%
    aggregate(TAD80~funcID+Tara_sample+Module, data=., FUN=sum) %>%
    aggregate(TAD80~Tara_sample+Module, data=., FUN=mean) %>%
    left_join(tara_meta_subset, by=c("Tara_sample" = "Sample_name")) %>%
    reshape2::dcast(Depth_layer+Tara_sample+Latitude_Start+Longitude_Start+Month+Sample_station~Module, value.var="TAD80") %>%
    filter(., !grepl("ZZZ",Tara_sample)) %>%
    filter(., !grepl("IZZ",Tara_sample))

#####

### Assess transcription of module functions across Arctic

#####

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
           ylim = c(-4e6,4e6),  
           xlim = c(-4e6,4e6)) + 
  geom_scatterpie(data = gene_clust_module_tara_abund_mean_wide_srf_utm, aes(x = X, y = Y, group = Tara_sample), cols = c("M1","M2","M3","M4","M5"),
  pie_scale = 3) + 
  geom_text(data = gene_clust_module_tara_abund_mean_wide_srf_utm, aes(label = Tara_sample, x = X, y = Y)) + 
  scale_fill_manual(values=module_colours) + 
  theme_bw()

# Visualise mean module transcript depth across surface, DCM and mesopelagic layers in a barplot
mod_func_transc_arctic_barplot = gene_clust_module_tara_abund_mean_wide %>%
    select(-Latitude_Start,-Longitude_Start) %>%
    reshape2::melt(., id.vars=c("Tara_sample","Month","Depth_layer","Sample_station"), variable.name = "Module", value.name = "TAD80") %>%
    aggregate(TAD80~Tara_sample+Module+Month+Depth_layer+Sample_station, data=., FUN=mean) %>%
    mutate(Depth_layer = factor(Depth_layer, levels=c("SRF","DCM","MES"))) %>%
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
pdf(file="output_files/FRAM_RAS_F4_NET_MOD_FUNC_transcript_tara_map.pdf", height=10, width=12)
mod_func_transc_arctic_map
dev.off()
pdf(file="output_files/FRAM_RAS_F4_NET_MOD_FUNC_transcript_tara_barplot.pdf", height=8, width=10)
mod_func_transc_arctic_barplot
dev.off()

#####

### Assess mean transcript depth of module functions within energy metabolism pathways

#####

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
pdf(file="output_files/FRAM_RAS_F4_NET_MOD_FUNC_transcript_tara_energy_metab_barplot.pdf", height=8, width=10)
gene_clust_tara_abund_module_top_metab_barplot
dev.off()