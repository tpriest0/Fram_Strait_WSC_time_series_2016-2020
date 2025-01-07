#### CLUSTER ENRICHMENT ANALYSIS ####
# load data
clusterEnrich=function(TimeSeriesData,nclu=4,cluster_vec=NULL,p_tol=0.05,tax_level="Class"){
    ## raw abundance matrix
    xx=TimeSeriesData$abundMatRaw
    tax=TimeSeriesData$taxa_info
    xx=apply(xx,2,function(x)x/sum(x)*mean(colMeans(xx)))
    #cluster numbers
    if(is.null(cluster_vec)){
      cluster_vec=TimeSeriesData$pts_data$cset$clusters[,paste0("K:",nclu)]
    }
      uclu=unique(cluster_vec)
    allEn=list()
    for(i in 1:length(uclu)){
        jj=which(cluster_vec==uclu[i])
        clu_xx=rowSums(xx[jj,])
        clu_tax=tax[jj,tax_level]  
        CluTaxSum=rowsum(clu_xx,clu_tax)
        AllTaxSum=rowsum(rowSums(xx),tax[,tax_level])
        
        #calculate for each cluster taxa enrichment by glm
        res=matrix(0,ncol=2,nrow=length(CluTaxSum))
        rownames(res)=rownames(CluTaxSum)
        for (j in 1:length(CluTaxSum)){
            #cluster taxa name
            tx=rownames(CluTaxSum)[j]
            #contingency 2x2 matrix
            # cluster tax o.i. abundance, other cluster tax abundance
            # all tax o.i. abundance, all other tax abundance
            mm=matrix(c(CluTaxSum[tx,],AllTaxSum[tx,], 
                        sum(CluTaxSum)-CluTaxSum[tx,],sum(AllTaxSum)-AllTaxSum[tx,]),
                      nrow=2)
            
            d <- data.frame(g=factor(1:2),s=mm[2,],f=mm[1,])
            g <- glm(s/(s+f)~g,weights=colSums(mm),data=d,family="binomial")
            cf=coef(summary(g))["g2",c("Estimate","Pr(>|z|)")]
            res[j,]=cf
            
        }
        res=cbind(res,p.adjust(res[,2],method = "BH"),rep(uclu[i],nrow(res)),rownames(res))
        colnames(res)=c("logFold","pval","p.adjust","module","taxa")
        rownames(res)=NULL
        res=res[which(as.numeric(res[,3]) < p_tol ),]
        allEn[[as.character(uclu[i])]]=res
  
}

return(allEn)
}

