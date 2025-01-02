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

### Import data

# Import microbial ASV count data
mic_asv_raw=read.table(file="FRAM_RAS_F4_PROK_ASV_raw_counts.txt", sep="\t",
                       check.names=F, header=T)

# Import microbial ASV taxa information
mic_asv_taxa=read.table(file="FRAM_RAS_F4_PROK_ASV_taxa.txt", sep="\t",
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
  decostand(., method="total", MARGIN=2) %>%
  tibble::rownames_to_column(., var="ASV_name")

euk_asv_filt_rel = 
  as.data.frame(euk_asv_filt_raw_wide) %>%
  decostand(., method="total", MARGIN=2) %>%
  tibble::rownames_to_column(., var="ASV_name")

##### 

### Rarefy data

### Prokaryotic ASVs

# Identify smallest read count across prokaryotic samples and rarefy the data
mic_asv_filt_raw_wide_mod = mic_asv_filt_raw_wide %>%
  tibble::column_to_rownames(., var = "ASV_name")
euk_asv_filt_raw_wide_mod = euk_asv_filt_raw_wide %>%
  tibble::column_to_rownames(., var = "ASV_name")

head(sort(colSums(mic_asv_filt_raw_wide_mod)))

mic_smallest_count <- min(colSums(mic_asv_filt_raw_wide_mod))

mic_asv_filt_rarefied = t(rrarefy(t(mic_asv_filt_raw_wide_mod), mic_smallest_count))

### Microeukaryotic ASVs

# Identify smallest read count
head(sort(colSums(euk_asv_filt_raw_wide_mod)))

# We will remove 02_2018_F4_2 based on lack of data and rarefy samples to the
# next lowest read count
euk_asv_filt_raw_wide_mod_subset = subset(euk_asv_filt_raw_wide_mod, select=-c(`02_2018_F4_2`))
euk_smallest_count <- min(colSums(euk_asv_filt_raw_wide_mod_subset))

euk_asv_filt_rarefied = t(rrarefy(t(euk_asv_filt_raw_wide_mod_subset), euk_smallest_count))

### Create relative abundance tables from rarefied data
mic_asv_filt_rarefied_rel = 
  as.data.frame(mic_asv_filt_rarefied) %>%
  decostand(., method="total", MARGIN=2) %>%
  tibble::rownames_to_column(., var="ASV_name")

euk_asv_filt_rarefied_rel = 
  as.data.frame(euk_asv_filt_rarefied) %>%
  decostand(., method="total", MARGIN=2) %>%
  tibble::rownames_to_column(., var="ASV_name")

### Based on the low counts observed in 02_2018_F4_2, also remove this 
### from the filtered relative abundance table as it will not be included
### in further analysis
euk_asv_filt_rel_mod = subset(euk_asv_filt_rel, select=-c(`02_2018_F4_2`))

### Export tables

write.table(mic_asv_filt_raw_wide,
            file=paste0(output_tables,"FRAM_RAS_F4_PROK_ASV_filt_raw.txt"), sep="\t", row.names = F, quote = F)
write.table(mic_asv_filt_rel,
            file=paste0(output_tables,"FRAM_RAS_F4_PROK_ASV_filt_rel.txt"), sep="\t", row.names = F, quote = F)
write.table(euk_asv_filt_raw_wide,
            file=paste0(output_tables,"FRAM_RAS_F4_EUK_ASV_filt_raw.txt"), sep="\t", row.names = F, quote = F)
write.table(euk_asv_filt_rel_mod,
            file=paste0(output_tables,"FRAM_RAS_F4_EUK_ASV_filt_rel.txt"), sep="\t", row.names = F, quote = F) 
write.table(mic_asv_filt_rarefied,
            file=paste0(output_tables,"FRAM_RAS_F4_PROK_ASV_filt_rare_raw.txt"), sep="\t", row.names = F, quote = F)
write.table(euk_asv_filt_rarefied,
            file=paste0(output_tables,"FRAM_RAS_F4_EUK_ASV_filt_rare_raw.txt"), sep="\t", row.names = F, quote = F)
write.table(mic_asv_filt_rarefied_rel,
            file=paste0(output_tables,"FRAM_RAS_F4_PROK_ASV_filt_rare_rel.txt"), sep="\t", row.names = F, quote = F)
write.table(euk_asv_filt_rarefied_rel,
            file=paste0(output_tables,"FRAM_RAS_F4_EUK_ASV_filt_rare_rel.txt"), sep="\t", row.names = F, quote = F)

#####

### Filtering taxa tables based on subset of ASVs

mic_asv_taxa_filt = 
  mic_asv_taxa %>%
  filter(., ASV_name %in% rownames(mic_asv_filt_rarefied))

euk_asv_taxa_filt = 
  euk_asv_taxa %>%
  filter(., ASV_name %in% rownames(euk_asv_filt_rarefied))

# Export tables 
write.table(mic_asv_taxa_filt,
            file=paste0(output_tables,"FRAM_RAS_F4_PROK_ASV_filt_taxa.txt"), sep="\t", quote=F, row.names=F)
write.table(euk_asv_taxa_filt,
            file=paste0(output_tables,"FRAM_RAS_F4_EUK_ASV_filt_taxa.txt"), sep="\t",quote=F, row.names=F)

#####

### Processing and reformatting metagenome gene cluster profiles

# Import gene cluster profile
gene_raw_abund=read.csv(
  file="FRAM_RAS_F4_GENE_CLUST_raw_counts.txt",
  sep="\t",check.names=F, header=T)

# Import information on the number of genomes sequenced in each metagenome
num_genomes=read.csv(
  file="FRAM_RAS_F4_num_genomes.txt",
  sep="\t",check.names=F, header=T)

# Remove genes detected in less than 3 samples
gene_names = gene_raw_abund %>%
  tibble::column_to_rownames(., var="clustID") %>%
  mutate(across(where(is.numeric), function(x) ifelse(x < 3, 0, x))) %>%
  mutate(across(where(is.numeric), function(x) ifelse(x >= 3, 1, x))) %>%
  mutate(sums = rowSums(.)) %>%
  filter(sums >= 3) %>%
  tibble::rownames_to_column(., var="Gene_clustID") %>%
  select(Gene_clustID)

# Filter gene count table 
gene_filt_raw_wide = gene_raw_abund %>%
  filter(clustID %in% gene_names$Gene_clustID) %>%
  rename(Gene_clustID = clustID)

# Normalise gene counts by the number of genomes sequenced in each sample
gene_filt_norm_wide = 
  gene_filt_raw_wide %>%
  reshape2::melt(., id.vars="Gene_clustID", variable.name="Sample_name", value.name="Count") %>%
  left_join(num_genomes, by="Sample_name") %>%
  mutate(Norm_count = Count/Num_genomes) %>%
  select(-Count,-Num_genomes) %>%
  reshape2::dcast(Gene_clustID~Sample_name, value.var="Norm_count", data=.) %>%
  replace(., is.na(.), 0)

# Create relative abundance table from normalised counts
gene_filt_rel_wide = 
    gene_filt_norm_wide %>%
    tibble::column_to_rownames(., var="Gene_clustID") %>%
    decostand(., method = "total", MARGIN=2) %>%
    tibble::rownames_to_column(., var="Gene_clustID")

### Export tables
write.table(gene_filt_raw_wide, file=paste0(output_tables,"FRAM_RAS_F4_GENE_CLUST_filt_raw_wide.txt"), sep="\t", row.names=F, quote=F)
write.table(gene_filt_norm_wide, file=paste0(output_tables,"FRAM_RAS_F4_GENE_CLUST_filt_norm_wide.txt"), sep="\t", row.names=F, quote=F)
write.table(gene_filt_rel_wide, file=paste0(output_tables,"FRAM_RAS_F4_GENE_CLUST_filt_rel_wide.txt"), sep="\t", row.names=F, quote=F)
