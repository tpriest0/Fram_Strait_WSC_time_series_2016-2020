####################

#### Determine oscillation signals

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
source("calcTimeSeriesInfo.R")
source("daily_approx.R")


### Import files

# Metadata
metadata=read.delim2("FRAM_RAS_F4_meta.txt")

# Relative abundance profiles of ASVs and gene clusters
mg <- read_delim(paste0(output_tables,"FRAM_RAS_F4_GENE_CLUST_filt_rel_wide.txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE)

Bac <- read_delim(paste0(output_tables,"FRAM_RAS_F4_MIC_ASV_filt_rel.txt"), 
                                           delim = "\t", escape_double = FALSE,
                                           trim_ws = TRUE)

Euk <- read_delim(paste0(output_tables,"FRAM_RAS_F4_EUK_ASV_filt_rel.txt"), 
                  delim = "\t", escape_double = FALSE,
                  trim_ws = TRUE)

#####

### Determine oscillation signals of gene clusters

#####

geneName=mg$clustID

#get similar colname
mg=mg[,2:ncol(mg)]
id=colnames(mg)
id=gsub("X([0-9])_","0\\1_",id)
id=gsub("X",'',id)

mi=match(id,metadata$RAS_id)
#time info
dt=metadata$date[mi] 
colnames(mg)=dt
#order data by time
oi=order(dt)
mg=mg[,oi]
mg=apply(mg,2,as.numeric)
rownames(mg)=geneName
#check NaNs per row
nans_per_row=apply(mg,1,function(x)sum(is.na(x)))

## keep at least the genes with < n-4 NaNs /47-4
## if you expect to observe ones per year  - > 4 values
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

### TimeSeries analysis using SegmenTier
tset=processTimeseries(mg,use.fft = T,dft.range =2:12,use.snr = T,na2zero = T)

### Calculate timeseries signals, phase, frequence etc.
# using the script calcTimeSeriesInfo
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

#colnames(tset$ts) = as.character(as.Date(colnames(tset$ts), "%m/%d/%Y"))

dap_set2All=daily_approx(as.Date(colnames(tset$ts)),tset$ts, fix=F)
colnames(dap_set2All$ee)=as.character(dap_set2All$dtw)

### Export table with information on oscillation signals per year
write.csv(dap_set2All$ee, file=paste0(output_tables,"FRAM_RAS_F4_GENE_CLUST_oscillations_per_year.txt"))


#####

### Determine oscillation signals of prokaryotic ASVs

#####

ASVName=Bac$ASV_name

#get similar colname
Bac=Bac[,2:ncol(Bac)]
id=colnames(Bac)
id=gsub("X([0-9])_","0\\1_",id)
id=gsub("X",'',id)

mi=match(id,metadata$RAS_id)
#time info
dt=metadata$date[mi] 
colnames(Bac)=dt
#order data by time
oi=order(dt)
Bac=Bac[,oi]
Bac=apply(Bac,2,as.numeric)
rownames(Bac)=ASVName
#check NaNs per row
nans_per_row=apply(Bac,1,function(x)sum(is.na(x)))

## keep at least the genes with < n-4 NaNs /47-4
## if you expect to observe ones per year  - > 4 values
num_of_genes = length(nans_per_row)

#filter NANs
wi=which(nans_per_row<44)
Bac=Bac[wi,]

num_genes_filtered = length(wi)
drop_out_rate = 1-(num_genes_filtered/num_of_genes)

#set all NaNs to zero
Bac[is.na(Bac)]=0
colnames(Bac) = as.character(as.Date(colnames(Bac), "%m/%d/%Y"))
Bac = Bac[ , order(colnames(Bac))]
Bac = Bac[,1:46]

### TimeSeries analysis using SegmenTier
tset=processTimeseries(Bac,use.fft = T,dft.range =2:12,use.snr = T,na2zero = T)

### Calculate timeseries signals, phase, frequence etc.
# use the script calcTimeSeriesInfo
ee=tset$dat

# Entferne alle Zeilen, in denen alle Werte NA sind
ee1 <- ee[rowSums(is.na(ee)) != ncol(ee), ]

allTS=apply(ee1,1,calcTimeSeriesInfo)
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

### Export table with information on oscillation signals per year
write.csv(dap_set2All$ee, file=paste0(output_tables,"FRAM_RAS_F4_MIC_ASV_oscillations_per_year.txt"))

#####

### Determine oscillation signals of Eukaryotic ASVs

#####

ASVName=Euk$ASV_name

#get similar colname
Euk=Euk[,2:ncol(Euk)]
id=colnames(Euk)
id=gsub("X([0-9])_","0\\1_",id)
id=gsub("X",'',id)

mi=match(id,metadata$RAS_id)
#time info
dt=metadata$date[mi] 
colnames(Euk)=dt
#order data by time
oi=order(dt)
Euk=Euk[,oi]
Euk=apply(Euk,2,as.numeric)
rownames(Euk)=ASVName
#check NaNs per row
nans_per_row=apply(Euk,1,function(x)sum(is.na(x)))

## keep at least the genes with < n-4 NaNs /47-4
## if you expect to observe ones per year  - > 4 values
num_of_genes = length(nans_per_row)

#filter NANs
wi=which(nans_per_row<44)
Euk=Euk[wi,]

num_genes_filtered = length(wi)
drop_out_rate = 1-(num_genes_filtered/num_of_genes)

#set all NaNs to zero
Euk[is.na(Euk)]=0
colnames(Euk) = as.character(as.Date(colnames(Euk), "%m/%d/%Y"))
Euk = Euk[ , order(colnames(Euk))]
Euk = Euk[,1:46]

### TimeSeries analysis using SegmenTier
tset=processTimeseries(Euk,use.fft = T,dft.range =2:12,use.snr = T,na2zero = T)

### Calculate timeseries signals, phase, frequence etc.
# use the script calcTimeSeriesInfo
ee=tset$dat

# Entferne alle Zeilen, in denen alle Werte NA sind
ee1 <- ee[rowSums(is.na(ee)) != ncol(ee), ]

allTS=apply(ee1,1,calcTimeSeriesInfo)
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

#### Export table with information on oscillation signals per year
write.csv(dap_set2All$ee, file=paste0(output_tables,"FRAM_RAS_F4_EUK_ASV_oscillations_per_year.txt"))

