
BAC_GENE_names <- read_delim("~/Desktop/F4Paper/dataEuk/input/RAS_F4_MIC_ASVs_OSC4_ID_list.txt", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)

EUK_GENE_names <- read_delim("~/Desktop/F4Paper/dataEuk/input/RAS_F4_EUK_ASV_osc4_rel.txt", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)


tsetMG=tset
dap_set2AllGEN = dap_set2All
yy_sortedGEN = yy_sorted
#save(tsetMG, file = "~/Desktop/F4Paper/dataEuk/output/FRAM_RAS_Genfunc_BacEuk03092024.Rdata")

#########################
###### Bacteria #########
#########################


mooring="F4"
org="bac"
#### load data ####
load(paste0("../../../MOSAiC_Data/4Y/",mooring,"-TimeSeriesData-",org,".Rdata"))

tset_dat = TimeSeriesData$pts_data$tset$dat

yy = TimeSeriesData$phaseFreqInfo
yy_sorted = yy[order(yy$freq), ]

colnames(TimeSeriesData$pts_data$tset$ts) = as.character(TimeSeriesData$date_vec)

vc = BAC_GENE_names$ID[grep("bac",BAC_GENE_names$ID)]
tset_dat = tset_dat[rownames(tset_dat) %in% vc, ]
tset_ts = TimeSeriesData$pts_data$tset$ts[rownames(TimeSeriesData$pts_data$tset$ts) %in% vc, ]
TimeS = TimeSeriesData$pts_data$tset$ts[rownames(TimeSeriesData$pts_data$tset$ts) %in% vc, ]
yy_sorted = yy_sorted[rownames(yy_sorted) %in% vc, ]

colnames(TimeS) = as.character(TimeSeriesData$date_vec)


mi=match(as.character(colnames(tsetMG$ts)),as.character(TimeSeriesData$date_vec))
tset_ts=tset_ts[,mi]
colnames(tset_ts)=as.character(TimeSeriesData$date_vec[mi])

tsetBac=tset_ts
yy_sortedBAC = yy_sorted

#save(tsetBac, file = "~/Desktop/F4Paper/dataEuk/output/FRAM_RAS_Bac_03092024.Rdata")

##########################
###### Eukaryota #########
##########################

mooring="F4"
org="euk"
#### load data ####
load(paste0("../../../MOSAiC_Data/4Y/",mooring,"-TimeSeriesData-",org,".Rdata"))

tset_dat = TimeSeriesData$pts_data$tset$dat

yy = TimeSeriesData$phaseFreqInfo
yy_sorted = yy[order(yy$freq), ]

colnames(TimeSeriesData$pts_data$tset$ts) = as.character(TimeSeriesData$date_vec)

vc = EUK_GENE_names$ASV_name[grep("euk",EUK_GENE_names$ASV_name)]
tset_dat = tset_dat[rownames(tset_dat) %in% vc, ]
tset_ts = TimeSeriesData$pts_data$tset$ts[rownames(TimeSeriesData$pts_data$tset$ts) %in% vc, ]
TimeS = TimeSeriesData$pts_data$tset$ts[rownames(TimeSeriesData$pts_data$tset$ts) %in% vc, ]
yy_sorted = yy_sorted[rownames(yy_sorted) %in% vc, ]

colnames(TimeS) = as.character(TimeSeriesData$date_vec)


mi=match(as.character(colnames(tsetMG$ts)),as.character(TimeSeriesData$date_vec))
tset_ts=tset_ts[,mi]
colnames(tset_ts)=as.character(TimeSeriesData$date_vec[mi])

tsetEuk=tset_ts
yy_sortedEUK = yy_sorted

#save(tsetEuk, file = "~/Desktop/F4Paper/dataEuk/output/FRAM_RAS_Euk_03092024.Rdata")

### create network combination of families and clusters
FamBacEuk=rbind(tsetBac,tsetEuk,tsetMG$ts)
# normalize data to make it comparable
FamBacEuk=apply(FamBacEuk,1,function(x)scale(x))
FamBacEuk=t(FamBacEuk)
## extract time Signals
### TimeSeries analysis using SegmenTier
FBtset=processTimeseries(FamBacEuk,use.fft = T,dft.range =2:12,use.snr = T,na2zero = T)

abundMat=FBtset$dat

rcc=rcorr(t(FBtset$dat),type = "pearson")
rcc$P[lower.tri(rcc$P,diag = T)]=NaN
rcc$r[lower.tri(rcc$r,diag = T)]=NaN
mP=reshape2::melt(rcc$P,na.rm = T)
mR=reshape2::melt(rcc$r,na.rm = T)
mR$p.adj=p.adjust(mP[,3],method = "fdr")
iii=which(mR$value>0 & mR$p.adj<0.05)
net=mR[iii,]
colnames(net)=c("from","to","weight","pval")
filtered_net=net[net$weight>=0.7,]

inet=graph_from_data_frame(net,directed = F)
icls=cluster_louvain(inet)

# get month for max abund value from each cluster
cluMonth=list()
for(i in 1:max(icls$membership)){
  mm=FamBacEuk[icls$names[which(icls$membership==i)],]
  cn=colnames(tsetMG$ts)[apply(mm,1,which.max)]
  mn=as.character(as.Date(cn),format="%m")
  cluMonth[[i]]=mn
}
cluMonth=melt(lapply(cluMonth,table))
colnames(cluMonth)=c("month","freq","cluster")
ggplot(cluMonth,aes(x=as.factor(month),y=freq,fill=factor(cluster)))+geom_bar(stat="identity",color="black")+
  facet_grid(cluster~.,scales = "free")+
  theme_bw()+xlab("month")+
  scale_x_discrete(labels=month.abb)+
  scale_fill_manual(values = c("#FF9900","#FF0033","#cc33cc","#6666FF","#66CCFF","#33FF99", "lightgreen"))


## summarize data
type=rep("bac",length(icls$names))
type[grep("func",icls$names)]="func"
type[grep("euk",icls$names)]="euk"
smNetData=data.frame(nodeId=icls$names,cluster=icls$membership,type=type)
write.table(smNetData,"~/Desktop/F4Paper/dataEuk/output/FRAM_RAS_F4_Networkmeta_OSC4_Bac_Euk_FUNC_03092024.txt",quote = F,row.names = F,sep = "\t")
write.table(filtered_net,"~/Desktop/F4Paper/dataEuk/output/FRAM_RAS_F4_Network_OSC4_Bac_Euk_FUNC_03092024.txt",quote = F,row.names = F,sep = "\t")



