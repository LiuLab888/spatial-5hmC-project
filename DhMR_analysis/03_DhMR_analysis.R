####************************************************************************####
####************************************************************************####
#### R packages ####
library(DSS)
library(bsseq)
library(genomation)


#### Sample name ####

#### 数据导入 ####


#### raw data ####

rm(list=ls())

setwd(".../results")


Cluster1 <- read.table("A_preWork_outputdir_Sample_name/bismark/matrix/all_Clusters/Cluster1/merge_sorted_filter_end_new_chr_dss.bed",
                       header=T)

Cluster2 <- read.table("A_preWork_outputdir_Sample_name/bismark/matrix/all_Clusters/Cluster2/merge_sorted_filter_end_new_chr_dss.bed",
                       header=T)


#### Cluster2 vs cluster1 ####

BSobj <- makeBSseqData( list(Cluster2, Cluster1),
                        c("C","N") )


index <- c(paste0("chr",1:19),"chrX")

for(i in 1:20){
  
  BSobj_subset <- BSobj[which(BSobj@rowRanges@seqnames == index[i]),]
  
  
  
  #without smothing
  
  #dmlTest = DMLtest(BSobj, group1=c("C"), group2=c("N"))
  
  
  #with smothing
  
  dmlTest.sm <- DMLtest(BSobj_subset, group1=c("C"), group2=c("N"), 
                        smoothing=TRUE)
  
  #DML
  
  # dmlTest.sm <- DMLtest(bs, group1 = c("Abnormal"), group2 = c("Normal"), 
  #                       smoothing = TRUE, BPPARAM = MulticoreParam())
  
  #dmls <- callDML(dmlTest.sm, delta=0.05, p.threshold = 0.001)
  
  #detal参数确定原假设的差异程度
  
  
  #DMR
  
  #detal参数确定原假设的差异程度
  
  dmrs <- callDMR(dmlTest.sm, delta = 0 , minlen=50, minCG=3, 
                  dis.merge=100, p.threshold = 1e-5)
  
  dmrs_up <- dmrs[which(dmrs$areaStat > 0),]
  
  dmrs_down <- dmrs[which(dmrs$areaStat < 0),]
  
  
  
  setwd(".../results")
  
  write.table(dmrs_up, paste0("DMR_analysis/sample_name/CpG/Cluster1_vs_2/DMRs_dss_up_",index[i],".txt"),
              row.names = F, quote = F, sep = "\t")
  
  write.table(dmrs_down, paste0("DMR_analysis/sample_name/CpG/Cluster1_vs_2/DMRs_dss_down_",index[i],".txt"),
              row.names = F, quote = F, sep = "\t")
  
  
}


showOneDMR(dmrs[1,], BSobj)

showOneDMR(dmrs[2,], BSobj)

