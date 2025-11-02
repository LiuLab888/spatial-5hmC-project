#### R packages ####
library(DOSE)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)
library(ggnewscale)
library(ggupset)
library(europepmc)
library(enrichplot)
library(biomaRt)
library(RColorBrewer)
library(openxlsx)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ReactomePA)
library(openxlsx)


library(GenomicFeatures)
library(ChIPseeker)
library(ggplot2)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene


####************************************************************************####
#### 5hmC -- Apical_progenitor_vs_MigratingNeuron ####
####************************************************************************####

library(GenomicFeatures)
library(ChIPseeker)
library(ggplot2)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

rm(list=ls())

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene


#### raw data ####

region1 <- read.table( 
  ".../merge_three_stages/Apical_progenitor_vs_MigratingNeuron/DMRs_dss_down_merge_filter_end.txt",
  sep = "\t", header = T)


region2 <- read.table( 
  ".../merge_three_stages/Apical_progenitor_vs_MigratingNeuron/DMRs_dss_up_merge_filter_end.txt",
  sep = "\t", header = T)


#### Peak提取 #### 

peak1 <- data.frame(chr=region1$chr,
                    start=region1$start,
                    end=region1$end)

write.table(peak1, 
            ".../merge_three_stages/Apical_progenitor_vs_MigratingNeuron/DMRs_dss_down_merge_filter_bed_end.txt",
            row.names = F, col.names = T, quote = F, sep = "\t")



peak2 <- data.frame(chr=region2$chr,
                    start=region2$start,
                    end=region2$end)

write.table(peak2,
            ".../merge_three_stages/Apical_progenitor_vs_MigratingNeuron/DMRs_dss_up_merge_filter_bed_end.txt",
            row.names = F, col.names = T, quote = F, sep = "\t")



#### 注释 ####

file_path1 <- ".../merge_three_stages/Apical_progenitor_vs_MigratingNeuron/DMRs_dss_down_merge_filter_bed_end.txt"

file_path2 <- ".../merge_three_stages/Apical_progenitor_vs_MigratingNeuron/DMRs_dss_up_merge_filter_bed_end.txt"


peakAnno1 <- annotatePeak(file_path1, tssRegion=c(-1000, 0), TxDb=txdb, 
                          annoDb = 'org.Mm.eg.db', addFlankGeneInfo=TRUE, flankDistance=5000)
plotAnnoBar(peakAnno1)

genelist1 <- peakAnno1@anno$geneId[which(peakAnno1@anno$annotation == "Promoter (<=1kb)") ]



peakAnno2 <- annotatePeak(file_path2, tssRegion=c(-1000, 0), TxDb=txdb, 
                          annoDb = 'org.Mm.eg.db', addFlankGeneInfo=TRUE, flankDistance=5000)
plotAnnoBar(peakAnno2)

genelist2 <- peakAnno2@anno$geneId[which(peakAnno2@anno$annotation == "Promoter (<=1kb)") ]


# DMR-Down N=1765
# DMR-Up N=40

peakAnnoList <- lapply(list("DMR-Down"=file_path1, "DMR-Up"=file_path2), 
                       annotatePeak, 
                       TxDb =txdb,
                       tssRegion=c(-1000, 0))

plotAnnoBar(peakAnnoList)


setwd(".../merge_three_stages/Apical_progenitor_vs_MigratingNeuron")


png(file="peak_distribution.png",width=500,height=300,res=100)

plotAnnoBar(peakAnnoList)

dev.off()


pdf(file="peak_distribution.pdf", width=5, height=3)

plotAnnoBar(peakAnnoList)

dev.off()



#### peak 注释 -- Down ####

library(GenomicFeatures)
library(ChIPseeker)
library(ggplot2)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

rm(list=ls())

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene


file_path <- ".../merge_three_stages/Apical_progenitor_vs_MigratingNeuron/DMRs_dss_down_merge_filter_bed_end.txt"


peakAnno <- annotatePeak(file_path, tssRegion=c(-1000, 0), TxDb=txdb, 
                         annoDb = 'org.Mm.eg.db', addFlankGeneInfo=TRUE, flankDistance=5000)

gene_list <- peakAnno@anno$geneId[which(peakAnno@anno$distanceToTSS <= 5000)]


#### Go ####

setwd(".../merge_three_stages/Apical_progenitor_vs_MigratingNeuron")


geneGO_enrich_All_results <- enrichGO(gene = gene_list, minGSSize = 5, 
                                      OrgDb = org.Mm.eg.db, ont = "ALL", pAdjustMethod = "BH",
                                      keyType = 'ENTREZID', pvalueCutoff = 0.01, qvalueCutoff = 0.2, 
                                      readable = TRUE)

View(geneGO_enrich_All_results@result)


geneGO_enrich_All_results_sorted <- geneGO_enrich_All_results@result[order(geneGO_enrich_All_results@result$pvalue),]

write.xlsx(geneGO_enrich_All_results_sorted, "GO_down.xlsx")


geneGO_enrich_All_results_sorted$Description[1:20]



#### 使用ggplot2进行更加自由的可视化呈现 ####

#先提取富集结果表
geneGO_enrich_All_results_subset <- geneGO_enrich_All_results_sorted[1:20,]


#指定绘图顺序（转换为因子）：
geneGO_enrich_All_results_subset$terms <- factor(geneGO_enrich_All_results_subset$Description,
                                                 levels = rev(geneGO_enrich_All_results_subset$Description))

#Top20富集条形图：
mytheme <-  theme(legend.title = element_text(size=14, family="sans", face="plain"))+
  theme(legend.text = element_text(size=14, family="sans", face="plain"))+
  theme(plot.title = element_text(size = 14,family="sans", face="plain", hjust = 0.5))+
  theme(axis.title.x = element_text(size = 14, family="sans", face="plain")) +
  theme(axis.title.y = element_text(size = 14, family="sans", face="plain"))+
  theme(axis.text.x = element_text(size = 14, color="black", face="plain",
                                   family="sans", colour="black")) +
  theme(axis.text.y = element_text(size = 14, color="black", face="plain",
                                   family="sans", colour="black")) +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, 
                                    linetype = "solid")) 

geneGO_enrich_All_results_subset_plot <- ggplot(data = geneGO_enrich_All_results_subset,
                                                aes(x = Count, y = terms, fill = -log10(pvalue)))+
  scale_fill_distiller(palette = "RdBu",direction = -1) +
  geom_bar(stat = "identity", width = 0.75) +
  theme_bw() +
  labs(x = "Number of Gene",
       y = "Gene terms",
       title = NULL) + mytheme

geneGO_enrich_All_results_subset_plot



png(file="GO_down.png", width=800, height=600)

geneGO_enrich_All_results_subset_plot

dev.off()


pdf(file="GO_down.pdf", width=8, height=6)
geneGO_enrich_All_results_subset_plot

dev.off()




#### peak 注释 -- Up ####


file_path <- ".../merge_three_stages/Apical_progenitor_vs_MigratingNeuron/DMRs_dss_up_merge_filter_bed_end.txt"


peakAnno <- annotatePeak(file_path, tssRegion=c(-1000, 0), TxDb=txdb, 
                         annoDb = 'org.Mm.eg.db', addFlankGeneInfo=TRUE, flankDistance=5000)

gene_list <- peakAnno@anno$geneId[which(peakAnno@anno$distanceToTSS <= 5000)]


#### Go ####

setwd(".../merge_three_stages/Apical_progenitor_vs_MigratingNeuron")


geneGO_enrich_All_results <- enrichGO(gene = gene_list, minGSSize = 5, 
                                      OrgDb = org.Mm.eg.db, ont = "ALL", pAdjustMethod = "BH",
                                      keyType = 'ENTREZID', pvalueCutoff = 0.01, qvalueCutoff = 0.2, 
                                      readable = TRUE)

View(geneGO_enrich_All_results@result)


geneGO_enrich_All_results_sorted <- geneGO_enrich_All_results@result[order(geneGO_enrich_All_results@result$pvalue),]

write.xlsx(geneGO_enrich_All_results_sorted, "GO_up.xlsx")


geneGO_enrich_All_results_sorted$Description[1:20]



#### 使用ggplot2进行更加自由的可视化呈现 ####

#先提取富集结果表
geneGO_enrich_All_results_subset <- geneGO_enrich_All_results_sorted[1:5,]


#指定绘图顺序（转换为因子）：
geneGO_enrich_All_results_subset$terms <- factor(geneGO_enrich_All_results_subset$Description,
                                                 levels = rev(geneGO_enrich_All_results_subset$Description))

#Top20富集条形图：
mytheme <-  theme(legend.title = element_text(size=14, family="sans", face="plain"))+
  theme(legend.text = element_text(size=14, family="sans", face="plain"))+
  theme(plot.title = element_text(size = 14,family="sans", face="plain", hjust = 0.5))+
  theme(axis.title.x = element_text(size = 14, family="sans", face="plain")) +
  theme(axis.title.y = element_text(size = 14, family="sans", face="plain"))+
  theme(axis.text.x = element_text(size = 14, color="black", face="plain",
                                   family="sans", colour="black")) +
  theme(axis.text.y = element_text(size = 14, color="black", face="plain",
                                   family="sans", colour="black")) +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, 
                                    linetype = "solid")) 

geneGO_enrich_All_results_subset_plot <- ggplot(data = geneGO_enrich_All_results_subset,
                                                aes(x = Count, y = terms, fill = -log10(pvalue)))+
  scale_fill_distiller(palette = "PuOr",direction = -1) +
  geom_bar(stat = "identity", width = 0.75) +
  theme_bw() +
  labs(x = "Number of Gene",
       y = "Gene terms",
       title = NULL) + mytheme

geneGO_enrich_All_results_subset_plot



png(file="GO_up.png", width=800, height=600)

geneGO_enrich_All_results_subset_plot

dev.off()


pdf(file="GO_up.pdf", width=8, height=6)
geneGO_enrich_All_results_subset_plot

dev.off()



####************************************************************************####
#### Heatmap -- Down ####
#### data #### 

rm(list=ls())

DhMR_level <- read.table( 
  ".../merge_three_stages/Apical_progenitor_vs_MigratingNeuron/DMRs_dss_down_merge_filter_end.txt",
  sep = "\t", header = T)

matrix_score <- DhMR_level[,c("meanMethy1","meanMethy2")]


#### plot ####

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("ComplexHeatmap")

library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggsci)
library(scales)
library(openxlsx)

matrix_score_plot <- as.matrix(matrix_score)

rownames(matrix_score_plot) <- rownames(matrix_score)

colnames(matrix_score_plot) <- c("Apical progenitor", "Migrating neuron")

#col <- colorRamp2(seq(0, 0.1, length.out=5), rev(brewer.pal(5, "PuOr")))


c("#8C510A", "#BF812D", "#DFC27D", "#F6E8C3", "#F5F5F5", "#C7EAE5", "#80CDC1", "#35978F", "#01665E")


col <- colorRamp2(seq(0, 0.2, length.out=11), rev(brewer.pal(11, "BrBG")))

setwd(".../merge_three_stages/Apical_progenitor_vs_MigratingNeuron/DhMR")

#"ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC


pdf("DhMR_5hmC.pdf", width = 4, height = 8)

set.seed(1)

Heatmap(matrix_score_plot, cluster_rows = T, show_row_dend = F, cluster_columns = F,na_col = "grey90", 
        rect_gp = gpar(col = NA),  col = col, border_gp = gpar(col = "black",lwd=0.25),row_title ="Hypo DhMR (N=1765)",
        show_row_names = F, show_column_names = T, #clustering_method_rows="ward.D",
        column_names_rot = 60, column_names_gp = gpar(fontsize = 16),
        heatmap_legend_param = list(legend_height = unit(4, "cm"),
                                    grid_width = unit(0.5, "cm"),
                                    title = "5hmC",
                                    title_gp = gpar(fontsize = 16, 
                                                    fontface = "plain"),
                                    labels_gp = gpar(fontsize = 16))) 

dev.off()



####************************************************************************####
#### 5hmC -- MigratingNeuron_vs_ProjectionNeuron ####
####************************************************************************####

library(GenomicFeatures)
library(ChIPseeker)
library(ggplot2)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

rm(list=ls())

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene


#### raw data ####

region1 <- read.table( 
  ".../merge_three_stages/MigratingNeuron_vs_ProjectionNeuron/DMRs_dss_down_merge_filter_end.txt",
  sep = "\t", header = T)


region2 <- read.table( 
  ".../merge_three_stages/MigratingNeuron_vs_ProjectionNeuron/DMRs_dss_up_merge_filter_end.txt",
  sep = "\t", header = T)


#### Peak提取 #### 

peak1 <- data.frame(chr=region1$chr,
                    start=region1$start,
                    end=region1$end)

write.table(peak1, 
            ".../merge_three_stages/MigratingNeuron_vs_ProjectionNeuron/DMRs_dss_down_merge_filter_bed_end.txt",
            row.names = F, col.names = T, quote = F, sep = "\t")



peak2 <- data.frame(chr=region2$chr,
                    start=region2$start,
                    end=region2$end)

write.table(peak2,
            ".../merge_three_stages/MigratingNeuron_vs_ProjectionNeuron/DMRs_dss_up_merge_filter_bed_end.txt",
            row.names = F, col.names = T, quote = F, sep = "\t")



#### 注释 ####

file_path1 <- ".../merge_three_stages/MigratingNeuron_vs_ProjectionNeuron/DMRs_dss_down_merge_filter_bed_end.txt"

file_path2 <- ".../merge_three_stages/MigratingNeuron_vs_ProjectionNeuron/DMRs_dss_up_merge_filter_bed_end.txt"


peakAnno1 <- annotatePeak(file_path1, tssRegion=c(-1000, 0), TxDb=txdb, 
                          annoDb = 'org.Mm.eg.db', addFlankGeneInfo=TRUE, flankDistance=5000)
plotAnnoBar(peakAnno1)

genelist1 <- peakAnno1@anno$geneId[which(peakAnno1@anno$annotation == "Promoter (<=1kb)") ]



peakAnno2 <- annotatePeak(file_path2, tssRegion=c(-1000, 0), TxDb=txdb, 
                          annoDb = 'org.Mm.eg.db', addFlankGeneInfo=TRUE, flankDistance=5000)
plotAnnoBar(peakAnno2)

genelist2 <- peakAnno2@anno$geneId[which(peakAnno2@anno$annotation == "Promoter (<=1kb)") ]


# DMR-Down N=1432
# DMR-Up N=49

peakAnnoList <- lapply(list("DMR-Down"=file_path1, "DMR-Up"=file_path2), 
                       annotatePeak, 
                       TxDb =txdb,
                       tssRegion=c(-1000, 0))

plotAnnoBar(peakAnnoList)


setwd(".../merge_three_stages/MigratingNeuron_vs_ProjectionNeuron")


png(file="peak_distribution.png",width=500,height=300,res=100)

plotAnnoBar(peakAnnoList)

dev.off()


pdf(file="peak_distribution.pdf", width=5, height=3)

plotAnnoBar(peakAnnoList)

dev.off()



#### peak 注释 -- Down ####

library(GenomicFeatures)
library(ChIPseeker)
library(ggplot2)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

rm(list=ls())

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene


file_path <- ".../merge_three_stages/MigratingNeuron_vs_ProjectionNeuron/DMRs_dss_down_merge_filter_bed_end.txt"


peakAnno <- annotatePeak(file_path, tssRegion=c(-1000, 0), TxDb=txdb, 
                         annoDb = 'org.Mm.eg.db', addFlankGeneInfo=TRUE, flankDistance=5000)

gene_list <- peakAnno@anno$geneId[which(peakAnno@anno$distanceToTSS <= 5000)]


#### Go ####

setwd(".../merge_three_stages/MigratingNeuron_vs_ProjectionNeuron")


geneGO_enrich_All_results <- enrichGO(gene = gene_list, minGSSize = 5, 
                                      OrgDb = org.Mm.eg.db, ont = "ALL", pAdjustMethod = "BH",
                                      keyType = 'ENTREZID', pvalueCutoff = 0.01, qvalueCutoff = 0.2, 
                                      readable = TRUE)

View(geneGO_enrich_All_results@result)


geneGO_enrich_All_results_sorted <- geneGO_enrich_All_results@result[order(geneGO_enrich_All_results@result$pvalue),]

write.xlsx(geneGO_enrich_All_results_sorted, "GO_down.xlsx")


geneGO_enrich_All_results_sorted$Description[1:20]



#### 使用ggplot2进行更加自由的可视化呈现 ####

#先提取富集结果表
geneGO_enrich_All_results_subset <- geneGO_enrich_All_results_sorted[1:20,]


#指定绘图顺序（转换为因子）：
geneGO_enrich_All_results_subset$terms <- factor(geneGO_enrich_All_results_subset$Description,
                                                 levels = rev(geneGO_enrich_All_results_subset$Description))

#Top20富集条形图：
mytheme <-  theme(legend.title = element_text(size=14, family="sans", face="plain"))+
  theme(legend.text = element_text(size=14, family="sans", face="plain"))+
  theme(plot.title = element_text(size = 14,family="sans", face="plain", hjust = 0.5))+
  theme(axis.title.x = element_text(size = 14, family="sans", face="plain")) +
  theme(axis.title.y = element_text(size = 14, family="sans", face="plain"))+
  theme(axis.text.x = element_text(size = 14, color="black", face="plain",
                                   family="sans", colour="black")) +
  theme(axis.text.y = element_text(size = 14, color="black", face="plain",
                                   family="sans", colour="black")) +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, 
                                    linetype = "solid")) 

geneGO_enrich_All_results_subset_plot <- ggplot(data = geneGO_enrich_All_results_subset,
                                                aes(x = Count, y = terms, fill = -log10(pvalue)))+
  scale_fill_distiller(palette = "RdBu",direction = -1) +
  geom_bar(stat = "identity", width = 0.75) +
  theme_bw() +
  labs(x = "Number of Gene",
       y = "Gene terms",
       title = NULL) + mytheme

geneGO_enrich_All_results_subset_plot



png(file="GO_down.png", width=1000, height=600)

geneGO_enrich_All_results_subset_plot

dev.off()


pdf(file="GO_down.pdf", width=10, height=6)
geneGO_enrich_All_results_subset_plot

dev.off()




#### peak 注释 -- Up ####


file_path <- ".../merge_three_stages/MigratingNeuron_vs_ProjectionNeuron/DMRs_dss_up_merge_filter_bed_end.txt"


peakAnno <- annotatePeak(file_path, tssRegion=c(-1000, 0), TxDb=txdb, 
                         annoDb = 'org.Mm.eg.db', addFlankGeneInfo=TRUE, flankDistance=5000)

gene_list <- peakAnno@anno$geneId[which(peakAnno@anno$distanceToTSS <= 5000)]


#### Go ####

setwd(".../merge_three_stages/MigratingNeuron_vs_ProjectionNeuron")


geneGO_enrich_All_results <- enrichGO(gene = gene_list, minGSSize = 5, 
                                      OrgDb = org.Mm.eg.db, ont = "ALL", pAdjustMethod = "BH",
                                      keyType = 'ENTREZID', pvalueCutoff = 0.01, qvalueCutoff = 0.2, 
                                      readable = TRUE)

View(geneGO_enrich_All_results@result)


geneGO_enrich_All_results_sorted <- geneGO_enrich_All_results@result[order(geneGO_enrich_All_results@result$pvalue),]

write.xlsx(geneGO_enrich_All_results_sorted, "GO_up.xlsx")


geneGO_enrich_All_results_sorted$Description[1:20]



#### 使用ggplot2进行更加自由的可视化呈现 ####

#先提取富集结果表
geneGO_enrich_All_results_subset <- geneGO_enrich_All_results_sorted[1:10,]


#指定绘图顺序（转换为因子）：
geneGO_enrich_All_results_subset$terms <- factor(geneGO_enrich_All_results_subset$Description,
                                                 levels = rev(geneGO_enrich_All_results_subset$Description))

#Top20富集条形图：
mytheme <-  theme(legend.title = element_text(size=14, family="sans", face="plain"))+
  theme(legend.text = element_text(size=14, family="sans", face="plain"))+
  theme(plot.title = element_text(size = 14,family="sans", face="plain", hjust = 0.5))+
  theme(axis.title.x = element_text(size = 14, family="sans", face="plain")) +
  theme(axis.title.y = element_text(size = 14, family="sans", face="plain"))+
  theme(axis.text.x = element_text(size = 14, color="black", face="plain",
                                   family="sans", colour="black")) +
  theme(axis.text.y = element_text(size = 14, color="black", face="plain",
                                   family="sans", colour="black")) +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, 
                                    linetype = "solid")) 

geneGO_enrich_All_results_subset_plot <- ggplot(data = geneGO_enrich_All_results_subset,
                                                aes(x = Count, y = terms, fill = -log10(pvalue)))+
  scale_fill_distiller(palette = "PuOr",direction = -1) +
  geom_bar(stat = "identity", width = 0.75) +
  theme_bw() +
  labs(x = "Number of Gene",
       y = "Gene terms",
       title = NULL) + mytheme

geneGO_enrich_All_results_subset_plot



png(file="GO_up.png", width=800, height=600)

geneGO_enrich_All_results_subset_plot

dev.off()


pdf(file="GO_up.pdf", width=8, height=6)
geneGO_enrich_All_results_subset_plot

dev.off()

####************************************************************************####
#### Heatmap -- Down ####
#### data #### 

rm(list=ls())

DhMR_level <- read.table( 
  ".../merge_three_stages/MigratingNeuron_vs_ProjectionNeuron/DMRs_dss_down_merge_filter_end.txt",
  sep = "\t", header = T)

matrix_score <- DhMR_level[,c("meanMethy1","meanMethy2")]


#### plot ####

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("ComplexHeatmap")

library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggsci)
library(scales)
library(openxlsx)

matrix_score_plot <- as.matrix(matrix_score)

rownames(matrix_score_plot) <- rownames(matrix_score)

colnames(matrix_score_plot) <- c("Migrating neuron", "Projection neuron")

#col <- colorRamp2(seq(0, 0.1, length.out=5), rev(brewer.pal(5, "PuOr")))


c("#8C510A", "#BF812D", "#DFC27D", "#F6E8C3", "#F5F5F5", "#C7EAE5", "#80CDC1", "#35978F", "#01665E")


col <- colorRamp2(seq(0, 0.2, length.out=11), rev(brewer.pal(11, "BrBG")))

setwd(".../merge_three_stages/MigratingNeuron_vs_ProjectionNeuron/DhMR")

#"ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC


pdf("DhMR_5hmC.pdf", width = 4, height = 8)

set.seed(1)

Heatmap(matrix_score_plot, cluster_rows = T, show_row_dend = F, cluster_columns = F,na_col = "grey90", 
        rect_gp = gpar(col = NA),  col = col, border_gp = gpar(col = "black",lwd=0.25),row_title ="Hypo DhMR (N=1432)",
        show_row_names = F, show_column_names = T, #clustering_method_rows="ward.D",
        column_names_rot = 60, column_names_gp = gpar(fontsize = 16),
        heatmap_legend_param = list(legend_height = unit(4, "cm"),
                                    grid_width = unit(0.5, "cm"),
                                    title = "5hmC",
                                    title_gp = gpar(fontsize = 16, 
                                                    fontface = "plain"),
                                    labels_gp = gpar(fontsize = 16))) 

dev.off()

