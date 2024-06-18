####################

#### Pre-processing ASV and gene data read for analysis

####################

### Define working directory
setwd('XXXXX')

### Define output directories
output_figures <- ('')
output_tables <- ('')

dir.create(output_tables)
dir.create(output_figures)

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
mic_asv_raw=read.table(file="RAS_F4_MIC_ASV_raw_counts.txt", sep="\t",
                       check.names=F, header=T)

# Import microbial ASV taxa information
mic_asv_taxa=read.table(file="RAS_F4_MIC_ASV_taxa.txt", sep="\t",
                       check.names=F, header=T)

# Import eukaryotic ASV count data
euk_asv_raw=read.table(file="RAS_F4_EUK_ASV_raw_counts.txt", sep="\t",
                            check.names=F, header=T)

# Import eukaryotic ASV taxa information
euk_asv_taxa=read.table(file="RAS_F4_EUK_ASV_taxa.txt", sep="\t",
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
  filter(ASV_name %in% mic_asv_filt_names$ASV_name) %>%
  tibble::column_to_rownames(., var = "ASV_name")

euk_asv_filt_names <- euk_asv_raw %>%
  tibble::column_to_rownames(., var="ASV_name") %>%
  mutate(across(where(is.numeric), function(x) ifelse(x < 3, 0, x))) %>%
  mutate(across(where(is.numeric), function(x) ifelse(x >= 3, 1, x))) %>%
  mutate(sums = rowSums(.)) %>%
  filter(sums >= 3) %>%
  tibble::rownames_to_column(., var="ASV_name") %>%
  select(ASV_name)

euk_asv_filt_raw_wide = euk_asv_raw %>%
  filter(ASV_name %in% euk_asv_filt_names$ASV_name) %>%
  tibble::column_to_rownames(., var = "ASV_name")

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
head(sort(colSums(mic_asv_filt_raw_wide)))

# We will remove 02_2018_F4_2 based on lack of data and rarefy samples to the
# next lowest read count
mic_smallest_count <- min(colSums(mic_asv_filt_raw_wide))

mic_asv_filt_rarefied = t(rrarefy(t(mic_asv_filt_raw_wide), mic_smallest_count))

### Microeukaryotic ASVs

# Identify smallest read count
head(sort(colSums(euk_asv_filt_raw_wide)))

# We will remove 02_2018_F4_2 based on lack of data and rarefy samples to the
# next lowest read count
euk_asv_filt_raw_mod = subset(euk_asv_filt_raw_wide, select=-c(`02_2018_F4_2`))
euk_smallest_count <- min(colSums(euk_asv_filt_raw_mod))

euk_asv_filt_rarefied = t(rrarefy(t(euk_asv_filt_raw_mod), euk_smallest_count))

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
            file="RAS_F4_MIC_ASV_filt_raw.txt", sep="\t", row.names = F, quote = F)
write.table(mic_asv_filt_rel,
            file="RAS_F4_MIC_ASV_filt_rel.txt", sep="\t", row.names = F, quote = F)
write.table(euk_asv_filt_raw_mod,
            file="RAS_F4_EUK_ASV_filt_raw.txt", sep="\t", row.names = F, quote = F)
write.table(euk_asv_filt_rel_mod,
            file="RAS_F4_EUK_ASV_filt_rel.txt", sep="\t", row.names = F, quote = F) 
write.table(mic_asv_filt_rarefied,
            file="RAS_F4_MIC_ASV_filt_rare_raw.txt", sep="\t", row.names = F, quote = F)
write.table(euk_asv_filt_rarefied,
            file="RAS_F4_EUK_ASV_filt_rare_raw.txt", sep="\t", row.names = F, quote = F)
write.table(mic_asv_filt_rarefied_rel,
            file="RAS_F4_MIC_ASV_filt_rare_rel.txt", sep="\t", row.names = F, quote = F)
write.table(euk_asv_filt_rarefied_rel,
            file="RAS_F4_EUK_ASV_filt_rare_rel.txt", sep="\t", row.names = F, quote = F)

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
            file="RAS_F4_MIC_ASV_filt_taxa.txt", sep="\t")
write.table(euk_asv_taxa_filt,
            file="RAS_F4_EUK_ASV_filt_taxa.txt", sep="\t")

#####

### Processing and reformatting metagenome gene cluster profiles

# Import gene cluster profile
ras_mg_geneclust_profile=read.table(
  file="RAS_F4_GENES_CLUSTID_abund_filt.txt",
  sep="\t",check.names=F, header=T)

# Create wide format raw count table with RAS_id headers
ras_mg_geneclust_profile_raw_wide = 
  ras_mg_geneclust_profile %>%
  select(RAS_id,clustID,Count) %>%
  dcast(clustID~RAS_id, value.var="Count", data=.)

# Create wide format normalised abundance table
ras_mg_geneclust_profile_norm_wide = 
  ras_mg_geneclust_profile %>%
  select(RAS_id,clustID,Norm_count) %>%
  reshape2::dcast(clustID~RAS_id, value.var="Norm_count", data=.) %>%
  replace(., is.na(.), 0)

# Create wide format relative abundance table
ras_mg_geneclust_profile_rel_prop_wide =
  ras_mg_geneclust_profile %>%
  select(RAS_id,clustID,Norm_count) %>%
  reshape2::dcast(clustID~RAS_id, value.var="Norm_count", data=.) %>%
  replace(., is.na(.), 0) %>%
  tibble::column_to_rownames(., var="clustID") %>%
  decostand(., method="total", MARGIN=2) %>%
  tibble::rownames_to_column(., var="clustID")


### Export tables

write.table(ras_mg_geneclust_profile_raw_wide,
            file="RAS_F4_GENE_CLUSTID_filt_raw_wide.txt",
            sep="\t", quote = F, row.names = F)

write.table(ras_mg_geneclust_profile_norm_wide,
            file="RAS_F4_GENE_CLUSTID_filt_norm_wide.txt",
            sep="\t", quote = F, row.names = F)

write.table(ras_mg_geneclust_profile_rel_prop_wide,
            file="RAS_F4_GENE_CLUSTID_filt_rel_wide.txt",
            sep="\t", quote = F, row.names = F)


