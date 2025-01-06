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

rasID=read.delim2("RAS_F4_META.txt")

mg <- read_delim("RAS_F4_GENE_CLUSTID_filt_rel_wide.txt", delim = "\t", 
                 escape_double = FALSE, trim_ws = TRUE)
mg$Zahl <- NULL

geneName=mg$clusterID

#get similar colname
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

### Export table
write.csv(dap_set2All$ee, file="GeneTSInterCluster.csv")

###

tsetMG=tset
dap_set2AllGEN = dap_set2All
yy_sortedGEN = yy_sorted
#save(tsetMG, file = "eigenMG_EGGNOG_yearly_oscill_OSC4_FUNC26_1_24.Rdata")

mooring="F4"
org="bac"

#### load data
load(paste0(mooring,"-TimeSeriesData-",org,".Rdata"))

tset_dat = TimeSeriesData$pts_data$tset$dat

yy = TimeSeriesData$phaseFreqInfo
yy_sorted = yy[order(yy$freq), ]

colnames(TimeSeriesData$pts_data$tset$ts) = as.character(TimeSeriesData$date_vec)

# vc = BAC_GENE_names$ID[grep("bac",BAC_GENE_names$ID)]
# tset_dat = tset_dat[rownames(tset_dat) %in% vc, ]
# tset_ts = TimeSeriesData$pts_data$tset$ts[rownames(TimeSeriesData$pts_data$tset$ts) %in% vc, ]
# TimeS = TimeSeriesData$pts_data$tset$ts[rownames(TimeSeriesData$pts_data$tset$ts) %in% vc, ]
# yy_sorted = yy_sorted[rownames(yy_sorted) %in% vc, ]

colnames(TimeS) = as.character(TimeSeriesData$date_vec)

mi=match(as.character(colnames(tsetMG$ts)),as.character(TimeSeriesData$date_vec))
tset_ts=tset_ts[,mi]
colnames(tset_ts)=as.character(TimeSeriesData$date_vec[mi])

tsetBac=tset_ts
yy_sortedBAC = yy_sorted

write.csv(tsetBac, file="BacTSInterCluster.csv")

# ### create network combination of families and clusters
# FamBac=rbind(tsetBac,tsetMG$ts)
# # normalize data to make it comparable
# FamBac=apply(FamBac,1,function(x)scale(x))
# FamBac=t(FamBac)
# ## extract time Signals
# ### TimeSeries analysis using SegmenTier
# FBtset=processTimeseries(FamBac,use.fft = T,dft.range =2:12,use.snr = T,na2zero = T)
# 
# abundMat=FBtset$dat
# 
# rcc=rcorr(t(FBtset$dat),type = "pearson")
# rcc$P[lower.tri(rcc$P,diag = T)]=NaN
# rcc$r[lower.tri(rcc$r,diag = T)]=NaN
# mP=reshape2::melt(rcc$P,na.rm = T)
# mR=reshape2::melt(rcc$r,na.rm = T)
# mR$p.adj=p.adjust(mP[,3],method = "fdr")
# iii=which(mR$value>0 & mR$p.adj<0.05)
# net=mR[iii,]
# colnames(net)=c("from","to","weight","pval")
# filtered_net=net[net$weight>=0.7,]
# 
# inet=graph_from_data_frame(net,directed = F)
# icls=cluster_louvain(inet)
# 
# # get month for max abund value from each cluster
# cluMonth=list()
# for(i in 1:max(icls$membership)){
#   mm=FamBac[icls$names[which(icls$membership==i)],]
#   cn=colnames(tsetMG$ts)[apply(mm,1,which.max)]
#   mn=as.character(as.Date(cn),format="%m")
#   cluMonth[[i]]=mn
# }
# cluMonth=melt(lapply(cluMonth,table))
# colnames(cluMonth)=c("month","freq","cluster")
# ggplot(cluMonth,aes(x=as.factor(month),y=freq,fill=factor(cluster)))+geom_bar(stat="identity",color="black")+
#   facet_grid(cluster~.,scales = "free")+
#   theme_bw()+xlab("month")+
#   scale_x_discrete(labels=month.abb)+
#   scale_fill_manual(values = c("#FF9900","#FF0033","#cc33cc","#6666FF","#66CCFF","#33FF99", "lightgreen"))
# 
# 
# ## summarize data
# type=rep("bac",length(icls$names))
# type[grep("Func",icls$names)]="func"
# smNetData=data.frame(nodeId=icls$names,cluster=icls$membership,type=type)
# write.table(smNetData,"GeneFam_and_Bac_Osc4_EGGNOG_4_OSC4_FUNC26_1_24.txt",quote = F,row.names = F,sep = "\t")
# write.table(filtered_net,"GeneFam_and_Bac_Network_Osc4_EGGNOG_4_OSC4_FUNC26_1_24.txt",quote = F,row.names = F,sep = "\t")
# 
# filtered_net$from = smNetData$cluster[which(smNetData$nodeId == filtered_net$from)]
