
####################

#### Analyse the composition of metagenome-assembled genomes

####################

### Define working directory
setwd('XXX')

### Define and create output directories
output_figures = ('results/figures/')
output_tables = ('results/tables/')

### Load libraries
library(webr)
library(dplyr)
library(reshape2)
library(tibble)
library(ggplot2)
library(tidyr)

#####

### Import data

#####
# Import MAG abundance table
mag_rbp_abund = read.table(file="FRAM_RAS_F4_MAG_reps-rbp-abund.txt", sep="\t", header=T, check.names=F)

# Import MAG summary information table
mag_info = read.table(file="FRAM_RAS_F4_MAG_reps-info.txt", sep="\t", header=T, check.names=F)

# Import sample metadata table
meta = read.table(file="FRAM_RAS_F4_META.txt", sep="\t", header=T, check.names=F)

#####

### Determine fraction of sequence genomes covered by MAGs

#####

# Create sample to month mapping
sample_to_month = meta %>%
    select(MG_sample_names,month)

# Create a long-format table that includes information on the fraction of genomes covered by the MAGs
mag_rbp_abund_prop_long = mag_rbp_abund %>%
    aggregate(Rel_prop~Sample_name, data=., FUN=sum) %>%
    rename(MAG_covered = Rel_prop) %>%
    mutate(Not_covered = 100-MAG_covered) %>%
    melt(., id.vars=c("Sample_name"), variable.name="Category", value.name="Relative_proportion") %>%
    left_join(sample_to_month, by=c("Sample_name" = "MG_sample_names"))

# Visualise the fraction of genomes covered by MAGs across months
mag_proportion_boxplot = ggplot(mag_rbp_abund_prop_long, aes(x = month, y = Relative_proportion)) + 
    geom_boxplot(aes(colour = Category)) + 
    labs(y = "Relative proportion of genomes sequenced (%)") +
    theme_bw() +
    theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        legend.text = element_text(size = 12),
        legend.title = element_blank()) +
    scale_colour_manual(values = c("MAG_covered" = "#004949", "Not_covered" = "#6db6ff"))

pdf(file = paste0(output_figures,"FRAM_RAS_F4_MAG_prop_genomes_covered_boxplot.pdf"), height = 6, width = 6)
mag_proportion_boxplot
dev.off()

#####

### Visualise completeness and contamination statistics of MAGs

#####

mag_comp_vs_con_scatter = ggplot(mag_info, aes(x = Completeness, y = Contamination)) +
    geom_point() + 
    labs(y = "Contamination (%)", x = "Completeness (%)") +
    theme_bw() +
    theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        legend.text = element_text(size = 12))

pdf(file = paste0(output_figures,"FRAM_RAS_F4_MAG_prop_genomes_covered_boxplot.pdf"), height = 6, width = 6)
mag_proportion_boxplot
dev.off()

#####

### Create PieChart illustrating taxonomic composition of MAGs

#####

# Split GTDB taxonomy string into individual rank columns 
mag_info_mod = mag_info %>%
    separate(GTDB_taxonomy, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";", fill = "right") %>%
    mutate(across(c(Domain, Phylum, Class, Order, Family, Genus, Species), ~ gsub(".__", "", .))) %>%
    mutate(across(c(Domain, Phylum, Class, Order, Family, Genus, Species), ~ ifelse(. == "", "Unknown", .)))

# Summarise Class and Family counts
mags_class_family = mag_info_mod %>%
select(Class,Family) %>%
mutate(count = 1)

mags_class_family_counts = mags_class_family %>%
    group_by(Class,Family) %>%
    summarise(n = sum(count))   

# Visualise and export PieChart
pdf(file = paste0(output_figures,"FRAM_RAS_F4_MAG_class_family_piedonut.pdf"), height = 8, width = 8)
PieDonut(mags_class_family_counts, aes(Class,Family, count=n), explode = 2)
dev.off()
