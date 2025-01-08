####################

#### Analyse oscillations of ecosystem components

####################
library(devtools)
install_github("raim/segmenTools")

### Define working directory (where data files from repository and scripts are stored)
setwd('XXXXX')

### If you change output directories then it should be consistent across all scripts
output_figures <- ('output_figures')
output_tables <- ('output_tables')

### Load libraries

library(segmenTier)
library(segmenTools)
library(igraph)
library(Hmisc)
library(reshape2)
library(ggplot2)
library(readr)
source("FRAM_RAS_F4_calcTimeSeriesInfo.R")
source("FRAM_RAS_F4_daily_approx.R")

### Import files

# Import metadata
rasID=read.delim2("RAS_F4_META.txt")

# Import the dataset whose components you want to calculate oscillation signals for. E.g. the prokaryotic ASV or the gene data.
# IMPORTANT: you need to import the filtered, relative abundance tables (these are produced by the FRAM_RAS_F4_time_series_preprocessing_data.R script)

mg <- read_delim("RAS_F4_GENE_CLUSTID_filt_rel_wide.txt", delim = "\t", 
                 escape_double = FALSE, trim_ws = TRUE)

mg$Zahl <- NULL

geneName=mg$clusterID

### Reformat tables
mg=mg[,2:ncol(mg)]
id=colnames(mg)
id=gsub("X([0-9])_","0\\1_",id)
id=gsub("X",'',id)

mi=match(id,rasID$RAS_id)
#time info
dt=rasID$date[mi] 
colnames(mg)=dt
#order data by time
oi=order(dt)
mg=mg[,oi]
mg=apply(mg,2,as.numeric)
rownames(mg)=geneName
#check NaNs per row
nans_per_row=apply(mg,1,function(x)sum(is.na(x)))

### Filter and retain genes with < n-4 NaNs (n = number of samples)
num_of_genes = length(nans_per_row)

#filter NANs
wi=which(nans_per_row<44)
mg=mg[wi,]

num_genes_filtered = length(wi)
drop_out_rate = 1-(num_genes_filtered/num_of_genes)

#set all NaNs to zero
mg[is.na(mg)]=0
colnames(mg) = as.character(as.Date(colnames(mg), "%m/%d/%Y"))
mg = mg[ , order(colnames(mg))]

### Determine oscillation signals using SegmenTier
tset=processTimeseries(mg,use.fft = T,dft.range =2:12,use.snr = T,na2zero = T)

# calculate timeseries signals, phase, frequence etc.
# use the script calcTimeSeriesInfo
ee=tset$dat
allTS=apply(ee,1,calcTimeSeriesInfo)
## extract sinus model data (seasonality info)
xx=lapply(allTS,function(x)x$mod_dataTS)
xx=do.call("rbind",xx)
#normalise scale from 0 to 1
xx=apply(xx,1,function(x)(x-min(x))/max(x-min(x)))
xx=t(xx)

# define time variable one sinus cycle = one year  
ts=seq(0,2*pi,length.out = ncol(xx))

# get the phase and frequence for each organism
yy=lapply(allTS,function(x)x$phase_freq)
yy=do.call("rbind",yy)

### sort data by frequence and phase
yy=as.data.frame(cbind(yy,ii=seq(nrow(yy))))
yy_sorted = yy[order(yy$freq),]

dap_set2All=daily_approx(as.Date(colnames(tset$ts)),tset$ts, fix=F)
colnames(dap_set2All$ee)=as.character(dap_set2All$dtw)

write.csv(dap_set2All$ee, file="FRAM_RAS_F4_oscillation_signals.r")
