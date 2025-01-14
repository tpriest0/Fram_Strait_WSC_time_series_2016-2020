####################

#### Pre-processing ASV and gene data read for analysis

####################

### Define working directory
setwd('XXXXX')

### Define and create output directories
dir.create('results')
dir.create('results/figures')
dir.create('results/tables')
output_figures = ('results/figures/')
output_tables = ('results/tables/')

#####

### Load libaries

library(dplyr)
library(vegan)
library(reshape2)
library(ggplot2)
library(tibble)

#####

### ASV data

#####

# Import microbial ASV count data
mic_asv_raw=read.table(file="FRAM_RAS_F4_MIC_ASV_raw_counts.txt", sep="\t",
                       check.names=F, header=T)

# Import microbial ASV taxa information
mic_asv_taxa=read.table(file="FRAM_RAS_F4_MIC_ASV_taxa.txt", sep="\t",
                       check.names=F, header=T)

# Import eukaryotic ASV count data
euk_asv_raw=read.table(file="FRAM_RAS_F4_EUK_ASV_raw_counts.txt", sep="\t",
                            check.names=F, header=T)

# Import eukaryotic ASV taxa information
euk_asv_taxa=read.table(file="FRAM_RAS_F4_EUK_ASV_taxa.txt", sep="\t",
                        check.names=F, header=T)

##### 

### Remove low abundant ASVs

# Identify ASVs with less than 3 counts in less than 3 samples and remove them

mic_asv_filt_names <- mic_asv_raw %>%
  tibble::column_to_rownames(., var="ASV_name") %>%
  mutate(across(where(is.numeric), function(x) ifelse(x < 3, 0, x))) %>%
  mutate(across(where(is.numeric), function(x) ifelse(x >= 3, 1, x))) %>%
  mutate(sums = rowSums(.)) %>%
  filter(sums >= 3) %>%
  tibble::rownames_to_column(., var="ASV_name") %>%
  select(ASV_name)

mic_asv_filt_raw_wide = mic_asv_raw %>%
  filter(ASV_name %in% mic_asv_filt_names$ASV_name) 

euk_asv_filt_names <- euk_asv_raw %>%
  tibble::column_to_rownames(., var="ASV_name") %>%
  mutate(across(where(is.numeric), function(x) ifelse(x < 3, 0, x))) %>%
  mutate(across(where(is.numeric), function(x) ifelse(x >= 3, 1, x))) %>%
  mutate(sums = rowSums(.)) %>%
  filter(sums >= 3) %>%
  tibble::rownames_to_column(., var="ASV_name") %>%
  select(ASV_name)

euk_asv_filt_raw_wide = euk_asv_raw %>%
  filter(ASV_name %in% euk_asv_filt_names$ASV_name) 

### Generate relative abundance tables from filtered ASV data
mic_asv_filt_rel = 
  as.data.frame(mic_asv_filt_raw_wide) %>%
  tibble::column_to_rownames(., var="ASV_name") %>%
  decostand(., method="total", MARGIN=2) %>%
  tibble::rownames_to_column(., var="ASV_name")

euk_asv_filt_rel = 
  as.data.frame(euk_asv_filt_raw_wide) %>%
  tibble::column_to_rownames(., var="ASV_name") %>%
  decostand(., method="total", MARGIN=2) %>%
  tibble::rownames_to_column(., var="ASV_name")

### Rarefy data

### Prokaryotic ASVs

# Identify smallest read count across prokaryotic samples and rarefy the data
mic_asv_filt_raw_wide_mod = mic_asv_filt_raw_wide %>%
  tibble::column_to_rownames(., var = "ASV_name")
euk_asv_filt_raw_wide_mod = euk_asv_filt_raw_wide %>%
  tibble::column_to_rownames(., var = "ASV_name")

head(sort(colSums(mic_asv_filt_raw_wide_mod)))

mic_smallest_count <- min(colSums(mic_asv_filt_raw_wide_mod))

mic_asv_filt_rarefied = t(rrarefy(t(mic_asv_filt_raw_wide_mod), mic_smallest_count)) %>%
  as.data.frame() %>%
  tibble::rownames_to_column(., var="ASV_name")

### Microeukaryotic ASVs

# Identify smallest read count
head(sort(colSums(euk_asv_filt_raw_wide_mod)))

# We will remove 02_2018_F4_2 based on lack of data and rarefy samples to the
# next lowest read count
euk_asv_filt_raw_wide_mod_subset = subset(euk_asv_filt_raw_wide_mod, select=-c(`02_2018_F4_2`))
euk_smallest_count <- min(colSums(euk_asv_filt_raw_wide_mod_subset))

euk_asv_filt_rarefied = t(rrarefy(t(euk_asv_filt_raw_wide_mod_subset), euk_smallest_count)) %>%
  as.data.frame() %>%
  tibble::rownames_to_column(., var="ASV_name")

### Create relative abundance tables from rarefied data
mic_asv_filt_rarefied_rel = mic_asv_filt_rarefied %>%
  tibble::column_to_rownames(., var = "ASV_name") %>%
  decostand(., method="total", MARGIN=2) %>%
  tibble::rownames_to_column(., var="ASV_name")

euk_asv_filt_rarefied_rel = euk_asv_filt_rarefied %>%
  tibble::column_to_rownames(., var = "ASV_name") %>%
  decostand(., method="total", MARGIN=2) %>%
  tibble::rownames_to_column(., var="ASV_name")

euk_asv_filt_rel_mod = subset(euk_asv_filt_rel, select=-c(`02_2018_F4_2`))

### Based on the low counts observed in 02_2018_F4_2, also remove this 
### from the filtered relative abundance table as it will not be included
### in further analysis
euk_asv_filt_rel_mod = subset(euk_asv_filt_rel, select=-c(`02_2018_F4_2`))

### Export tables

write.table(mic_asv_filt_raw_wide,
            file=paste0(output_tables,"FRAM_RAS_F4_MIC_ASV_filt_raw.txt"), sep="\t", row.names = F, quote = F)
write.table(mic_asv_filt_rel,
            file=paste0(output_tables,"FRAM_RAS_F4_MIC_ASV_filt_rel.txt"), sep="\t", row.names = F, quote = F)
write.table(euk_asv_filt_raw_wide,
            file=paste0(output_tables,"FRAM_RAS_F4_EUK_ASV_filt_raw.txt"), sep="\t", row.names = F, quote = F)
write.table(euk_asv_filt_rel_mod,
            file=paste0(output_tables,"FRAM_RAS_F4_EUK_ASV_filt_rel.txt"), sep="\t", row.names = F, quote = F) 
write.table(mic_asv_filt_rarefied,
            file=paste0(output_tables,"FRAM_RAS_F4_MIC_ASV_filt_rare_raw.txt"), sep="\t", row.names = F, quote = F)
write.table(euk_asv_filt_rarefied,
            file=paste0(output_tables,"FRAM_RAS_F4_EUK_ASV_filt_rare_raw.txt"), sep="\t", row.names = F, quote = F)
write.table(mic_asv_filt_rarefied_rel,
            file=paste0(output_tables,"FRAM_RAS_F4_MIC_ASV_filt_rare_rel.txt"), sep="\t", row.names = F, quote = F)
write.table(euk_asv_filt_rarefied_rel,
            file=paste0(output_tables,"FRAM_RAS_F4_EUK_ASV_filt_rare_rel.txt"), sep="\t", row.names = F, quote = F)

### Filtering taxa tables based on subset of ASVs

mic_asv_taxa_filt = 
  mic_asv_taxa %>%
  filter(., ASV_name %in% rownames(mic_asv_filt_rarefied))

euk_asv_taxa_filt = 
  euk_asv_taxa %>%
  filter(., ASV_name %in% rownames(euk_asv_filt_rarefied))

# Export tables 
write.table(mic_asv_taxa_filt,
            file=paste0(output_tables,"FRAM_RAS_F4_MIC_ASV_filt_taxa.txt"), sep="\t", quote=F, row.names=F)
write.table(euk_asv_taxa_filt,
            file=paste0(output_tables,"FRAM_RAS_F4_EUK_ASV_filt_taxa.txt"), sep="\t",quote=F, row.names=F)

#####

### Metagenome-derived gene cluster data

#####


### Read input
gene_clust_abund = fread(file="FRAM_RAS_F4_proteins_all_clust_rep_abund.txt", 
    sep="\t", header=T)

mg_num_genomes = fread(file="FRAM_RAS_F4_num_genomes.txt", 
    sep="\t", header=T)

sample_meta = fread(file="FRAM_RAS_F4_META.txt", 
    sep="\t", header=T)

### Create list of nonsingletons (i.e. remove gene cluster representatives with only one count in one sample)
gene_clust_abund_filt = gene_clust_abund %>%
    aggregate(Count~Gene_rep, data=., FUN=sum) %>%
    filter(., Count >= 2)

### Filer out singletons from abundance table, then normalise the counts by the number of genomes estimated in each sample
abund_nonsing_raw_and_norm_long = gene_clust_abund %>%
    filter(., Gene_rep %in% gene_clust_abund_filt$Gene_rep) %>%
    mutate(Sample = gsub("-","_",Sample)) %>%
    left_join(mg_num_genomes, by=c("Sample" = "Sample_name")) %>%
    mutate(Norm_abund = Count/Num_genomes) %>%
    select(-Num_genomes) %>%
    left_join(sample_meta, by=c("Sample" = "RAS_sample_names")) %>%
    select(RAS_id,Gene_rep,Count,Norm_abund)

### Calculate relative abundance based on normalised counts and reformat to wide format relative abundance table
abund_filt_wide_rel = abund_nonsing_raw_and_norm_long %>%
    reshape2::dcast(Gene_rep~RAS_id, value.var="Norm_abund", data=.) %>%
    mutate(across(where(is.numeric), ~replace_na(., 0))) %>%
    tibble::column_to_rownames(., var="Gene_rep") %>%
    decostand(., method="total", MARGIN=2) %>%
    tibble::rownames_to_column(., var="Gene_rep")

### As we are focused on those gene clusters with clear patterns over multiple years, 
### we will apply a second filtering, to remove those gene clusters with <3 counts in <3 samples.
selected_genereps = gene_clust_abund %>%
    filter(., Count >= 3) %>%
    mutate(Count = 1) %>%
    aggregate(Count~Gene_rep, data=., FUN=sum) %>%
    filter(Count >= 3) %>%
    select(Gene_rep) %>%
    mutate(clustID = paste0("geneclust_",row_number()))

### Create wide format raw count table for selected gene clusters
gene_clust_filt_wide_raw = abund_nonsing_raw_and_norm_long %>%
    filter(., Gene_rep %in% selected_genereps$Gene_rep) %>%
    left_join(selected_genereps, by="Gene_rep") %>%
    reshape2::dcast(clustID~RAS_id, value.var="Count", data=.) %>%
    mutate(across(where(is.numeric), ~replace_na(., 0)))

### Create wide format normalised counts table for selected gene clusters
gene_clust_filt_wide_norm = abund_nonsing_raw_and_norm_long %>%
    filter(., Gene_rep %in% selected_genereps$Gene_rep) %>%
    left_join(selected_genereps, by="Gene_rep") %>%
    reshape2::dcast(clustID~RAS_id, value.var="Norm_abund", data=.) %>%
    mutate(across(where(is.numeric), ~replace_na(., 0)))

### Create wide format relative abundance table for selected gene clusters (this file will be the input for the oscillation signal determination)
gene_clust_filt_wide_rel = abund_filt_wide_rel %>%
    filter(., Gene_rep %in% selected_genereps$Gene_rep) %>%
    left_join(selected_genereps, by="Gene_rep") %>%
    select(-Gene_rep) %>%
    tibble::column_to_rownames(., var="clustID") %>%
    tibble::rownames_to_column(., var="clustID")

### Export tables
fwrite(gene_clust_filt_wide_raw, file=paste0(output_tables,"FRAM_RAS_F4_GENE_CLUST_filt_raw.txt"), sep="\t", quote=F, row.names=F)
fwrite(gene_clust_filt_wide_norm, file=paste0(output_tables,"FRAM_RAS_F4_GENE_CLUST_filt_norm.txt"), sep="\t",quote=F, row.names=F)
fwrite(gene_clust_filt_wide_rel, file=paste0(output_tables,"FRAM_RAS_F4_GENE_CLUST_filt_rel.txt"), sep="\t",quote=F, row.names=F)
fwrite(abund_nonsing_raw_and_norm_long, file=paste0(output_tables,"FRAM_RAS_F4_GENE_CLUST_nonsing_raw_and_norm.txt"), sep="\t",quote=F, row.names=F)
fwrite(selected_genereps, file=paste0(output_tables,"FRAM_RAS_F4_GENE_CLUST_filt_name_mapping.txt"), sep="\t",quote=F, row.names=F)
