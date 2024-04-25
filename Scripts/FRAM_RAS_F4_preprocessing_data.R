####################

#### Pre-processing ASV and gene data read for analysis

####################


# Define working directory
setwd('C:/Users/tpriest/OneDrive - ETH Zurich/MPI - FRAM/WSC')

# Load libaries
library(dplyr)
library(vegan)
library(reshape2)
library(ggplot2)
library(tibble)

### Import data

# Import microbial ASV count data
mic_asv_raw=read.table(file="ASV/RAS_F4_MIC_ASV_raw_counts.txt", sep="\t",
                       check.names=F, header=T)

# Import microbial ASV taxa information
mic_asv_taxa=read.table(file="ASV/RAS_F4_MIC_ASV_taxa.txt", sep="\t",
                       check.names=F, header=T)

# Import eukaryotic ASV count data
euk_asv_raw=read.table(file="ASV/RAS_F4_EUK_ASV_raw_counts.txt", sep="\t",
                            check.names=F, header=T)

# Import eukaryotic ASV taxa information
euk_asv_taxa=read.table(file="ASV/RAS_F4_EUK_ASV_taxa.txt", sep="\t",
                        check.names=F, header=T)

### Identify ASVs with less than 3 counts in less than 3 samples and remove them
# Filter ASV dataset and create list of ASV names passing filtering
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

### Rarefy data

## Prokaryotic
# Identify smallest read count across prokaryotic samples and rarefy the data
head(sort(colSums(mic_asv_filt_raw_wide)))

# We will remove 02_2018_F4_2 based on lack of data and rarefy samples to the
# next lowest read count
mic_smallest_count <- min(colSums(mic_asv_filt_raw_wide))

mic_asv_filt_rarefied = t(rrarefy(t(mic_asv_filt_raw_wide), mic_smallest_count))

## Microeukaryotic
# Identify smallest read count
head(sort(colSums(euk_asv_filt_raw_wide)))

# We will remove 02_2018_F4_2 based on lack of data and rarefy samples to the
# next lowest read count
euk_asv_filt_raw_mod = subset(euk_asv_filt_raw_wide, select=-c(`02_2018_F4_2`))
euk_smallest_count <- min(colSums(euk_asv_filt_raw_mod))

euk_asv_filt_rarefied = t(rrarefy(t(euk_asv_filt_raw_mod), euk_smallest_count))

### Export rarefied count tables
write.table(mic_asv_filt_rarefied,
            file="ASV/RAS_F4_MIC_ASV_filt_rare_raw.txt", sep="\t")
write.table(euk_asv_filt_rarefied,
            file="ASV/RAS_F4_EUK_ASV_filt_rare_raw.txt", sep="\t")

### Create relative abundance tables
mic_asv_filt_rarefied_rel = 
  as.data.frame(mic_asv_filt_rarefied) %>%
  decostand(., method="total", MARGIN=2) %>%
  tibble::rownames_to_column(., var="ASV_name")

euk_asv_filt_rarefied_rel = 
  as.data.frame(euk_asv_filt_rarefied) %>%
  decostand(., method="total", MARGIN=2) %>%
  tibble::rownames_to_column(., var="ASV_name")

# Export tables 
write.table(mic_asv_filt_rarefied_rel,
            file="ASV/RAS_F4_MIC_ASV_filt_rare_rel.txt", sep="\t")
write.table(euk_asv_filt_rarefied_rel,
            file="ASV/RAS_F4_EUK_ASV_filt_rare_rel.txt", sep="\t")

### Filter taxa tables to include only those ASVs that passed filtering 
mic_asv_taxa_filt = 
  mic_asv_taxa %>%
  filter(., ASV_name %in% rownames(mic_asv_filt_rarefied))

euk_asv_taxa_filt = 
  euk_asv_taxa %>%
  filter(., ASV_name %in% rownames(euk_asv_filt_rarefied))

# Export tables 
write.table(mic_asv_taxa_filt,
            file="ASV/RAS_F4_MIC_ASV_filt_taxa.txt", sep="\t")
write.table(euk_asv_taxa_filt,
            file="ASV/RAS_F4_EUK_ASV_filt_taxa.txt", sep="\t")

###

### Processing and reformatting metagenome gene cluster profiles

# Import metdata
ras_metadata=read.table(file="RAS_F4_META_METAG.txt", sep="\t",
                        check.names=F, header=T)

# Import functional gene profile
ras_mg_geneclust_profile=read.table(
  file="metagenomes/community_gene_profiles/FRAM_RAS_F4_GENES_CLUSTID_abund_filt.txt",
  sep="\t",check.names=F, header=T)

# Create metagenome sample name to RAS_id mapping file
name_list = ras_metadata %>%
  select(metaG_sample_title,RAS_id)
View(name_list)

ras_mg_geneclust_profile %>%
  select(Sample) %>%
  unique() %>%
  dim()

# Create wide format raw count table with RAS_id headers
ras_mg_geneclust_profile_raw_wide = 
  ras_mg_geneclust_profile %>%
  left_join(name_list, by=c("Sample" = "metaG_sample_title")) %>%
  select(RAS_id,clustID,Count) %>%
  dcast(clustID~RAS_id, value.var="Count", data=.)

# Create wide format normalised abundance table
ras_mg_geneclust_profile_norm_wide = 
  ras_mg_geneclust_profile %>%
  left_join(name_list, by=c("Sample" = "metaG_sample_title")) %>%
  select(RAS_id,clustID,Norm_count) %>%
  reshape2::dcast(clustID~RAS_id, value.var="Norm_count", data=.) %>%
  replace(., is.na(.), 0)

# Create wide format relative abundance table
ras_mg_geneclust_profile_rel_prop_wide =
  ras_mg_func_profile %>%
  left_join(name_list, by=c("Sample" = "metaG_sample_title")) %>%
  select(RAS_id,clustID,Norm_count) %>%
  reshape2::dcast(clustID~RAS_id, value.var="Norm_count", data=.) %>%
  replace(., is.na(.), 0) %>%
  tibble::column_to_rownames(., var="clustID") %>%
  decostand(., method="total", MARGIN=2) %>%
  tibble::rownames_to_column(., var="clustID")


### Export tables
write.table(ras_mg_geneclust_profile_raw_wide,
            file="metagenomes/community_gene_profiles/FRAM_RAS_F4_GENE_CLUSTID_filt_raw_wide.txt",
            sep="\t")

write.table(ras_mg_geneclust_profile_norm_wide,
            file="metagenomes/community_gene_profiles/FRAM_RAS_F4_GENE_CLUSTID_filt_norm_wide.txt",
            sep="\t")

write.table(ras_mg_geneclust_profile_rel_prop_wide,
            file="metagenomes/community_gene_profiles/FRAM_RAS_F4_GENE_CLUSTID_filt_rel_wide.txt",
            sep="\t")

######

### What proportion of community gene content was captured by EGGNOG orthologous
### groups?
ras_mg_func_profile_rel_prop_wide %>%
  reshape2::melt(id.vars="geneID", variable.name="RAS_id", value.name="Rel_abund") %>%
  filter(geneID == "Unassigned") %>%
  arrange(desc(Rel_abund)) %>%
  mutate(prop = 1.00 - Rel_abund) %>%
  summarise(mean = mean(prop))


###

### What was the composition of the metagenomes at kingdom level? i.e. what 
### proportion of reads was bacteria, archaea or eukarya? Read assignments were
### performed using tiara

# Import metdata
mg_read_taxonomy=read.table(file="metagenomes/FRAM_RAS_F4_read_tiara_assignments.txt", 
                        sep="\t",check.names=F, header=T)

# Import metadata and create sample name to date mapping file
ras_metadata=read.table(file="RAS_F4_META.txt", sep="\t",
                        check.names=F, header=T) %>%
  select(RAS_id,date) %>%
  mutate(date = as.Date(date, "%d.%m.%Y"))

# Reformat names to RAS_id and calculate relative proportions
mg_read_taxonomy_rel_prop_wide = mg_read_taxonomy %>%
  select(RAS_id,Assignment) %>%
  table() %>%
  as.data.frame() %>%
  reshape2::dcast(RAS_id~Assignment, value.var="Freq", data=.) %>%
  tibble::column_to_rownames(., var="RAS_id") %>%
  decostand(., method="total", MARGIN=1)

# Reform relative proportion table to long format and join with metadata for plotting
mg_read_taxonomy_rel_prop_long_w_date = 
  mg_read_taxonomy_rel_prop_wide %>%
  tibble::rownames_to_column(., var="RAS_id") %>%
  reshape2::melt(variable.name = "Taxonomy", value.name = "Rel_prop", data=.) %>%
  left_join(ras_metadata, by="RAS_id") %>%
  mutate(Taxonomy = factor(Taxonomy, levels = c("archaea","bacteria",
                                                "prokarya","organelle",
                                                "eukarya","unknown")))

# Visualise in stacked bar plot
taxonomy_colours = c("archaea" = "#000000",
                     "bacteria" = "#009292",
                     "eukarya" = "#ffb6db",
                     "organelle" = "#db6d00",
                     "prokarya" = "#006ddb",
                     "unknown" = "gray90")

mg_read_comp_stacked_bar = ggplot(mg_read_taxonomy_rel_prop_long_w_date,
       aes(x = as.Date(date), y = Rel_prop)) +
  geom_bar(aes(fill = Taxonomy), position="stack", stat="identity") + 
  scale_x_date(date_breaks = "6 months", date_labels =  "%b %Y") +
  scale_fill_manual(values = taxonomy_colours) +
  labs(y = "Relative proportion (%)") + 
  theme_bw() + 
  theme(axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))

# Export figure
pdf("figures_output/FRAM_RAS_MG_read_composition_kingdom_stacked_bar.pdf",
    height=6, width=8)
mg_read_comp_stacked_bar
dev.off()
