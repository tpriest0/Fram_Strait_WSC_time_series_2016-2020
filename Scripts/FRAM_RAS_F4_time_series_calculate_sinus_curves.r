##################################
############## F4 ###############
#################################

#GeneData.csv is RAS_F4_GENE_CLUSTID_filt_rel_wide with dates a colmn names

gene <- read_delim("GeneData.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)#
gene <- as.data.frame(gene)
rownames(gene)=gene$Gene
gene$Gene = NULL
#bac <- read_delim("bac_oz1_abundance.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)

geneID <- read_delim("GeneTSInterCluster25_8_output_file_3years.csv", delim = ",", escape_double = FALSE, trim_ws = TRUE)
genID = geneID$asv[which(geneID$oscillations==1)]

gene = gene[which(rownames(gene) %in% genID),]

load("F4-TimeSeriesData-bac.Rdata")
Bac <- TimeSeriesData$abundMatRaw

BacID <- read_delim("New_BacALLTaylor_3years_output_file.csv", delim = ",", escape_double = FALSE, trim_ws = TRUE)
BacID = BacID$asv[which(BacID$oscillations==1)]

Bac = Bac[which(rownames(Bac) %in% BacID),]

sampelsgen = colnames(gene)

bac = Bac[,which(colnames(Bac) %in% sampelsgen)]

sampelsbac = colnames(bac)
gene = gene[,which(colnames(gene) %in% sampelsbac)]

both = as.data.frame(rbind(bac,gene))
# both$Zahl = NULL
# rownames(both) = both$Gene
# both$Gene = NULL

source("calcTimeSeriesInfo.R")
source("clusterEnrich.R")

library(MASS)
library(reshape2)
library(reshape)
library(ggplot2)


# year1 =[,4:13]
# year2 =[,14:23]
# year3 =[,24:36]
# year4 =[,33:45]

clusAVG=both[,4:13]
#phase shifts
pts=apply(clusAVG,1,calcTimeSeriesInfo)
allSins=sapply(pts,function(x) x$modData_sin)
ppf=sapply(pts,function(x)x$phase_freq,simplify = T)




ddf=data.frame(date=rep(colnames(clusAVG),ncol(allSins)),
               sin=melt(allSins)[,c(2,3)],
               #avg=melt(t(TimeSeriesData$cluster_average))[,3])
               avg=melt(t(clusAVG))[,3])
colnames(ddf)=c("date","module","sinus","avg")
#ddf$module=factor(ddf$module,levels=order(ppf[1,],decreasing = F),ordered=T)
#modSize=table(clusterVector)



ddf1 = stack(ddf[,c(2,4)])
ddf1$date = ddf$date
ddf1$module = ddf$module
ddf1$module = as.character(ddf1$module)
ddf1$module = as.factor(ddf1$module)


ddf2 = stack(ddf[,c(2,3)])
ddf2$date = ddf$date
ddf2$module = ddf$module
ddf2$module = as.character(ddf2$module)
ddf2$module = as.factor(ddf2$module)

#ddf2= ddf2[which(ddf2$module==c("bac_asv1","bac_asv2")),]
#ddf1= ddf1[which(ddf1$module==c("bac_asv1","bac_asv2")),]


# ggplot() +
#   geom_line(data=ddf1, aes(x=date, y=values, group=module, color=module)) +
#   geom_line(data=ddf2, aes(x=date, y=values, group=module, color=module)) +
#   scale_color_viridis(discrete = TRUE) +
#   theme_ipsum() +
#   ylab("abundance")

year1_sin = ddf2
year1_org = ddf1

year1_sin_w =reshape(year1_sin,v.names="values",timevar="date", idvar=c("module"),
                     direction="wide")


# peaks <- read_delim("~/Desktop/F4Paper/Tay_New_Data/RAS_F4_MIC_ASVs_OSC4_ID_list.txt", 
#                     delim = "\t", escape_double = FALSE, trim_ws = TRUE)
# 
# vec = peaks$ID
# 
# year1_sin_w = year1_sin_w[which(year1_sin_w$module %in% vec),]
# year2_sin_w = year2_sin_w[which(year2_sin_w$module %in% vec),]
# year3_sin_w = year3_sin_w[which(year3_sin_w$module %in% vec),]
# year4_sin_w = year4_sin_w[which(year4_sin_w$module %in% vec),]

colnames(year1_sin_w) <- gsub("^values\\.", "", colnames(year1_sin_w))
colnames(year2_sin_w) <- gsub("^values\\.", "", colnames(year2_sin_w))
colnames(year3_sin_w) <- gsub("^values\\.", "", colnames(year3_sin_w))
#colnames(year4_sin_w) <- gsub("^values\\.", "", colnames(year4_sin_w))

# Function to get the maximum value and corresponding column name(s) for each row
max_value_and_column <- function(row) {
  max_value <- max(row)
  max_columns <- names(row)[row == max_value]
  return(list(MaxValue = max_value, MaxColumns = paste(max_columns, collapse = ', ')))
}

rownames(year1_sin_w) = year1_sin_w$module
rownames(year2_sin_w) = year2_sin_w$module
rownames(year3_sin_w) = year3_sin_w$module
#rownames(year4_sin_w) = year4_sin_w$module


# Apply the function row-wise
result1 <- apply(year1_sin_w[3:ncol(year1_sin_w)], 1, max_value_and_column)
result1_df <- do.call(rbind.data.frame, result1)

result2 <- apply(year2_sin_w[3:ncol(year1_sin_w)], 1, max_value_and_column)
result2_df <- do.call(rbind.data.frame, result2)

result3 <- apply(year3_sin_w[3:ncol(year3_sin_w)], 1, max_value_and_column)
result3_df <- do.call(rbind.data.frame, result3)

#result4 <- apply(year4_sin_w[3:ncol(year4_sin_w)], 1, max_value_and_column)
#result4_df <- do.call(rbind.data.frame, result4)

all = cbind(result1_df,result2_df,result3_df)#,result4_df)
colnames(all)= c("maxvalue_1", "date1", "maxvalue_2", "date2", "maxvalue_3", "date3")#, "maxvalue_4", "date4")

#write.table(all, file = "../../../../MOSAiC/MicGen_Sinusodial_approximation_of_the_abundance_curve_for_peak_difference_calculation.csv", sep = ";")

all_ <- lapply(all, function(x) gsub(",.*", "", x))

#all1_ <- lapply(all, function(x) gsub("^(.*?),.*$", "\\1", x))


test = as.data.frame(all_)
rownames(test)=rownames(all)
#all$max <- apply(all[,2:ncol(all)], 1, max, na.rm=TRUE)

#test <- read_delim("~/Desktop/test1.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)

# Convert character dates to Date objects
for (col in c("date1", "date2", "date3")) {
  test[[col]] <- as.Date(test[[col]])
}

# Calculate the differences between dates
test$Diff1_2 <- test$date2 - test$date1
test$Diff2_3 <- test$date3 - test$date2
test$Diff3_4 <- test$date4 - test$date3

# Print the resulting dataframe with date differences
print(test)


write.table(test, file = "RAS_F4_MIC_ASV_GENES_OSC4_peak_timings_julian.txt", sep = ";")


