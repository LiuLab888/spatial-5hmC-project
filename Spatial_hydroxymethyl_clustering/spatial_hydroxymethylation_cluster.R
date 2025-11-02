####*************************************************************************####
####*************************************************************************####
#### mE155S1-16H_split  ####

rm(list=ls())

setwd("/media/tangym/ETD/Spatial_Methylation/data_analysis/mE155S1-16H_split")

#### 空间位置和甲基化数据 ####

input_name1 <- ".../mE155S1-16H_split/coverage_matrix_total.txt"

input_name2 <- ".../mE155S1-16H_split/methylation_matrix_total.txt"

refer_path <- ".../mE155S1-16H_split/mE155S1-16H_ref_barcode.txt"


# 1
ref_barcode <- read.table(refer_path,stringsAsFactors = F)

our_barcode1 <- read.table(input_name1,stringsAsFactors = F)

our_barcode2 <- read.table(input_name2,stringsAsFactors = F)


# 2
merge_barcode0 <- merge(our_barcode1[which(our_barcode1$V1 %in% ref_barcode$V1),], 
                        our_barcode2[which(our_barcode2$V1 %in% ref_barcode$V1),], 
                        by="V1", all=TRUE)

merge_barcode1 <- merge(ref_barcode, merge_barcode0, by="V1", all=TRUE)

colnames(merge_barcode1) <- c("V1","V2","V3","coverage","methylation")

merge_barcode1$coord <- paste0(merge_barcode1$V2,"x", merge_barcode1$V3)


#### 数值过滤 ####

merge_barcode_filter <- merge_barcode1[which(merge_barcode1$coverage >= 50000 & 
                                               merge_barcode1$methylation >= 0),]


#### 边界过滤 ####

if(FALSE){
  
  source(".../data_analysis/0_rm_outliner.shaorui.new.R")
  
  coord <- merge_barcode_filter$coord
  
  dat <- coord_to_dat(coord)
  
  value <- merge_barcode_filter$methylation
  
  dat$value <- value
  
  tissue_selector(dat, 800, 800)
  
  
}

index_filtered_EM <- read.table(".../mE155S1-16H_split/pixels_delete.txt", header=T)

merge_barcode_filter1 <- merge_barcode_filter[-which(merge_barcode_filter$coord %in% index_filtered_EM$coord ),]


dim1 <-  merge_barcode_filter1$V2

dim2 <-  merge_barcode_filter1$V3



P1 <- c(72,17); P2 <- c(43,8)

A1 <- P2[2]-P1[2]; B1 <- P1[1]-P2[1]; C1 <- P2[1]*P1[2] - P1[1]*P2[2]

value1 <-scale((A1*dim1 + B1*dim2 + C1)/sqrt(A1^2 + B1^2))


P3 <- c(72,33); P4 <- c(43,22)

D1 <- P4[2]-P3[2]; E1 <- P3[1]-P4[1]; F1 <- P4[1]*P3[2] - P3[1]*P4[2]

value2 <-scale((D1*dim1 + E1*dim2 + F1)/sqrt(D1^2 + E1^2))


P5 <- c(72,24); P6 <- c(43,12)

G1 <- P6[2]-P5[2]; H1 <- P5[1]-P6[1]; I1 <- P6[1]*P5[2] - P5[1]*P6[2]

value3 <-scale((G1*dim1 + H1*dim2 + I1)/sqrt(G1^2 + H1^2))


P7 <- c(72,12); P8 <- c(43,3)

J1 <- P8[2]-P7[2]; K1 <- P7[1]-P8[1]; L1 <- P8[1]*P7[2] - P7[1]*P8[2]

value4 <-scale((J1*dim1 + K1*dim2 + L1)/sqrt(J1^2 + K1^2))



value5 <- scale(sqrt((dim1-57)^2 + (dim2-36)^2 ))

value6 <- scale(sqrt((dim1-57)^2 + (dim2-1)^2 ))



methylation <- (merge_barcode_filter1$methylation-mean(merge_barcode_filter1$methylation))/
  (sd(merge_barcode_filter1$methylation)/2)


data_merge <- data.frame(PCA1=value1,
                         PCA2=value2,
                         PCA3=value3,
                         PCA4=value4,
                         PCA5=value5,
                         PCA6=value6,
                         PCA7=1*methylation)


PCA_selection <- data_merge


#### 聚类分析 Kmeans ####

if(FALSE){
  
  ## 数据群体性估计，计算hopkins统计量 ##
  
  index <- get_clust_tendency(PCA_selection, 50, graph=TRUE)
  
  
  ## 认为Hopkins统计量的值<0.5，表明数据是高度可聚合的 ## 
  
  index$hopkins_stat
  
  index$plot
  
  
  ## 估计聚类簇
  
  set.seed(100)
  
  gap_state <- clusGap(PCA1_10, FUN = kmeans, nstart = 25, K.max = 10, B = 50)
  
  fviz_gap_stat(gap_state)   
  
}


#### 利用Kmeans进行聚类 ####

set.seed(100)

km_results <- kmeans(PCA_selection, 5, nstart = 25)

#fviz_cluster(km_results, PCA_selection) #显示聚类分布情况

PCA_selection_cluster <- as.data.frame(km_results$cluster)

colnames(PCA_selection_cluster)[1] <- "cluster"

data_cluster <- data.frame(V1 = merge_barcode_filter1$coord,
                           V2 = PCA_selection_cluster$cluster)

summary(factor(data_cluster$V2))

# 1   2   3   4   5 
# 135 189 168 258 268 

write.table(data_cluster,"data_cluster_EM_part.txt",
            row.names = FALSE,
            col.names = TRUE)



#### umap ####
library(umap)

## 利用PCA处理后的数据umap处理 ##
## 参数选择1：100, 0.8, 8
## 参数选择1：100, 0.7, 8年

config.params <- umap.defaults

config.params$random_state=100

config.params$min_dist=0.8 #控制局部距离，值越大距离越离散

config.params$n_neighbors=10#控制整体，值越大总体越分散

config.params$n_components=2

#config.params$set_op_mix_ratio=0.5

#config.params$local_connectivity=2

umap_result <- umap(PCA_selection, config = config.params)

umap1 <- umap_result$layout[,1]

umap2 <- umap_result$layout[,2]

#umap3 <- umap_result$layout[,3]


#reticulate::py_install("umap-learn")


data_merge_all <- cbind(merge_barcode_filter1, cluster=data_cluster$V2, umap1, umap2)

write.table(data_merge_all,"data_merge_all_EM_part.txt",
            row.names = FALSE,
            col.names = TRUE)


#### 原始聚类 ####

#### plot UMAP #### 

data_group <- factor(data_merge_all$cluster)

summary(data_group)

umap1 <- data_merge_all$umap1

umap2 <- data_merge_all$umap2

data_plot <- data.frame(umap1, umap2, data_group)

cols <- c("#699ECA", "#FF8C00", "#F898CB", "#4DAF4A", "#D65190", "#731A73", "#FFCB5B", "#e66f51", "#0076B9", 
          "#3D505A", "#0098B2","#FBEA2E", "#F8B072", "#8582BD", "#4F99C9", "#A8D3A0", "#A6D0E6", "#EC3E31")     

plot_scatter <- ggplot(data=data_plot, aes(x=umap1, y=umap2, 
                                           color = data_group)) + 
  geom_point(size=3, alpha=1,shape=19 ) +
  
  scale_color_manual(values = cols) 

# stat_ellipse(aes(fill=data_group2), geom='polygon', type="norm",
#              level=0.68, alpha=0.5, show.legend = F) 

plot_scatter


#消除背景和网格
plot_scatter <- plot_scatter + theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())+
  theme(axis.ticks = element_line(size = 0.75, color ="black"),
        axis.ticks.length=unit(1.25, 'mm'))


#调节文本大小，face取值包含，plain普通，bold加粗，italic斜体，bold.italic斜体加粗            
plot_scatter <-plot_scatter + theme(legend.title = element_text(size=18, family = "sans", face = "plain"))+
  theme(legend.text = element_text(size=16, family = "sans", face = "plain")) +
  theme(plot.title = element_text(size = 18, family = "sans", face = "plain", hjust = 0.5)) +
  theme(axis.title.x = element_text(size = 18, family = "sans", face = "plain")) +
  theme(axis.title.y = element_text(size = 18, family = "sans", face = "plain"))+
  theme(axis.text.x = element_text(size = 18, color = "black", face = "plain",
                                   family = "sans", colour = "black")) +
  theme(axis.text.y = element_text(size = 18, color = "black", face = "plain",
                                   family = "sans", colour = "black")) +  
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1.25, 
                                    linetype = "solid")) + #调节图形的四周边框
  labs(title = "",x = "UMAP1", y = "UMAP2") +
  theme(plot.margin = unit(c(0.5,0.5,0.2,0.2),"cm"))+
  
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05)), limits = c(-12, 12),
                     breaks = c(-10, -5, 0, 5, 10),
                     labels = c( "-10", "-5", "0", "5", "10")) +  
  
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)), limits = c(-13, 14),
                     breaks = c(-10, -5, 0, 5),
                     labels = c("-10", "-5", "0", "5")) +
  
  # theme(legend.position=c(0.97, 0.97),
  #       legend.justification = c(0.97, 0.97)) +
  # 
  # theme(legend.background = element_rect(fill="white",colour="black")) +
  
  guides(color=guide_legend(title = "Cluster"))

plot_scatter


png(file="Cluster_UMAP_part.png", width=600, height=600)

plot_scatter

dev.off()


pdf(file="Cluster_UMAP_part.pdf", width=6, height=6)

plot_scatter

dev.off()


#### 去掉非cell部分 ####

#### plot UMAP #### 

data_merge_all_filter <- data_merge_all[-which(data_merge_all$cluster==9),]

data_merge_all_filter$cluster[data_merge_all_filter$cluster == 10] <- 9


data_group <- factor(data_merge_all_filter$cluster)

summary(data_group)

umap1 <- data_merge_all_filter$umap1

umap2 <- data_merge_all_filter$umap2

data_plot <- data.frame(umap1, umap2, data_group)

cols0 <- c("#699ECA", "#FF8C00", "#F898CB", "#4DAF4A", "#D65190", "#731A73", "#FFCB5B", "#e66f51", "#0076B9", 
           "#3D505A", "#0098B2","#FBEA2E", "#F8B072", "#8582BD", "#4F99C9", "#A8D3A0", "#A6D0E6", "#EC3E31")

cols <- c("#699ECA", "#FF8C00", "#F898CB", "#4DAF4A", "#D65190", "#731A73", "#FFCB5B", "#e66f51", 
          "#3D505A", "#0098B2","#FBEA2E", "#F8B072", "#8582BD", "#4F99C9", "#A8D3A0", "#A6D0E6", "#EC3E31")      

plot_scatter <- ggplot(data=data_plot, aes(x=umap1, y=umap2, 
                                           color = data_group)) + 
  geom_point(size=2, alpha=1,shape=19 ) +
  
  scale_color_manual(values = cols) 

# stat_ellipse(aes(fill=data_group2), geom='polygon', type="norm",
#              level=0.68, alpha=0.5, show.legend = F) 

plot_scatter


#消除背景和网格
plot_scatter <- plot_scatter + theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())+
  theme(axis.ticks = element_line(size = 0.75, color ="black"),
        axis.ticks.length=unit(1.25, 'mm'))


#调节文本大小，face取值包含，plain普通，bold加粗，italic斜体，bold.italic斜体加粗            
plot_scatter <-plot_scatter + theme(legend.title = element_text(size=18, family = "sans", face = "plain"))+
  theme(legend.text = element_text(size=16, family = "sans", face = "plain")) +
  theme(plot.title = element_text(size = 18, family = "sans", face = "plain", hjust = 0.5)) +
  theme(axis.title.x = element_text(size = 18, family = "sans", face = "plain")) +
  theme(axis.title.y = element_text(size = 18, family = "sans", face = "plain"))+
  theme(axis.text.x = element_text(size = 18, color = "black", face = "plain",
                                   family = "sans", colour = "black")) +
  theme(axis.text.y = element_text(size = 18, color = "black", face = "plain",
                                   family = "sans", colour = "black")) +  
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1.25, 
                                    linetype = "solid")) + #调节图形的四周边框
  labs(title = "",x = "UMAP1", y = "UMAP2") +
  theme(plot.margin = unit(c(0.5,0.5,0.2,0.2),"cm"))+
  
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05)), limits = c(-15, 17),
                     breaks = c( -10, 0, 10),
                     labels = c("-10", "0", "10")) +  
  
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)), limits = c(-17, 17),
                     breaks = c(-10, 0, 10),
                     labels = c("-10", "0", "10")) +
  
  # theme(legend.position=c(0.97, 0.97),
  #       legend.justification = c(0.97, 0.97)) +
  # 
  # theme(legend.background = element_rect(fill="white",colour="black")) +
  
  guides(color=guide_legend(title = "Cluster"))

plot_scatter


png(file="Cluster_UMAP.png", width=600, height=600)

plot_scatter

dev.off()


pdf(file="Cluster_UMAP.pdf", width=8, height=8)

plot_scatter

dev.off()


#### 原始聚类 ####

#### plot location #### 

data_merge_all <- read.table(".../data_analysis/mE155S1-16H_split/data_merge_all_EM_part.txt",
                             header=T)

data_group <- factor(data_merge_all$cluster)

summary(data_group)

nA <- data_merge_all$V2

nB <- data_merge_all$V3

data_plot <- data.frame(nA, nB, data_group)

cols <- c("#699ECA", "#FF8C00", "#F898CB", "#4DAF4A", "#D65190", "#731A73", "#FFCB5B", "#e66f51", "#0076B9", 
          "#3D505A", "#0098B2","#FBEA2E", "#F8B072", "#8582BD", "#4F99C9", "#A8D3A0", "#A6D0E6", "#EC3E31")   


spatialPlot <- ggplot(data=data_plot, aes(x=nA,y=nB,color=data_group))+
  geom_point(shape=16,size=3,show.legend=FALSE)+
  scale_color_manual(values =cols) +
  scale_y_continuous(trans = "reverse",breaks = seq(0,45,5)) +
  scale_x_continuous(trans = "reverse",breaks = seq(0,12,5)) +
  expand_limits(y=c(0,45),x=c(0,12))+ labs(color="Cluster") + theme_void() 

spatialPlot


setwd(".../data_analysis/mE155S1-16H_split")

png(file="Cluster_location_part.png",width=120/1.2,height=450/1.2)

spatialPlot 

dev.off()


pdf(file="Cluster_location_part.pdf",width=1.2/0.9,height=4.5/0.9)

spatialPlot 

dev.off()


#### plot location -- 与RNA对齐画法 ####

data_merge_all <- read.table(".../data_analysis/mE155S1-16H_split/data_merge_all_EM_part.txt",
                             header=T)

data_group <- factor(data_merge_all$cluster)

summary(data_group)

nA <- as.numeric(data_merge_all$V2)

nB <- as.numeric(data_merge_all$V3)

data_plot <- data.frame(nA, nB, data_group)

cols <- c("#699ECA", "#FF8C00", "#F898CB", "#4DAF4A", "#D65190", "#731A73", "#FFCB5B", "#e66f51", "#0076B9", 
          "#3D505A", "#0098B2","#FBEA2E", "#F8B072", "#8582BD", "#4F99C9", "#A8D3A0", "#A6D0E6", "#EC3E31")   


spatialPlot <- ggplot(data=data_plot, aes(x=nA,y=nB,color=data_group))+
  geom_point(shape=15,size=1,show.legend=FALSE)+
  scale_color_manual(values =cols) +
  scale_y_continuous(trans = "reverse",breaks = seq(1,36,5)) +
  scale_x_continuous(trans = "reverse",breaks = seq(43,72,5)) +
  expand_limits(x=c(43,72),y=c(1,36))+ labs(color="Cluster") + theme_void()+
  theme(plot.margin = margin(0,0,0,0,unit="pt"))

spatialPlot


setwd(".../data_analysis/mE155S1-16H_split")

png(file="Cluster_location_part_new.png",width=300/2,height=360/2)

spatialPlot 

dev.off()


pdf(file="Cluster_location_part_new.pdf",width=3/2.35,height=3.6/2.2)

spatialPlot 

dev.off()



#### 原始聚类 ####

#### plot location -- methylation level #### 

library(circlize)
library(RColorBrewer)
library(viridis)
library(ggsci)

data_merge_all <- read.table(".../data_analysis/mE155S1-16H_split/data_merge_all_EM_part.txt",
                             header=T)

data_level <- data_merge_all$methylation

nA <- data_merge_all$V2

nB <- data_merge_all$V3

data_plot <- data.frame(nA, nB, data_level)

spatialPlot <- ggplot(data=data_plot, aes(x=nB,y=nA,color=data_level))+
  geom_point(shape=16,size=3,show.legend=T)+
  #scale_color_continuous(type = "viridis") +
  scale_color_gradientn(colors=viridis(11), 
                        breaks=c(0.25, 0.5, 0.75), 
                        label=c("0.25","0.50","0.75") ,
                        limits=c(0.2, 0.85)) +
  scale_y_continuous(breaks = seq(0,20,5)) +
  scale_x_continuous(trans = "reverse", breaks = seq(0,15,5)) +
  expand_limits(y=c(1,20),x=c(2,12))+ labs(color="ML") + theme_void() 

spatialPlot


setwd(".../data_analysis/mE155S1-16H_split")

png(file="Cluster_location_raw_EM_ML_part_new_scaled.png",width=200/1.6,height=250/1.65)

spatialPlot 

dev.off()


pdf(file="Cluster_location_raw_EM_ML_part_new_scaled.pdf",width=2/1.2,height=2.5/1.25)

spatialPlot 

dev.off()


#### 去掉非cell部分 ####

#### plot location #### 

data_group <- factor(data_merge_all_filter$cluster)

summary(data_group)

nA <- data_merge_all_filter$V2

nB <- data_merge_all_filter$V3

data_plot <- data.frame(nA, nB, data_group)


cols <- c("#699ECA", "#FF8C00", "#F898CB", "#4DAF4A", "#D65190", "#731A73", "#FFCB5B", "#e66f51", 
          "#3D505A", "#0098B2","#FBEA2E", "#F8B072", "#8582BD", "#4F99C9", "#A8D3A0", "#A6D0E6", "#EC3E31")                   

spatialPlot <- ggplot(data=data_plot, aes(x=nB,y=nA,color=data_group))+
  geom_point(shape=16,size=2,show.legend=F)+
  scale_color_manual(values =cols) +
  scale_y_continuous(breaks = seq(0,95,5)) +
  scale_x_continuous(trans = "reverse",breaks = seq(0,95,5)) +
  expand_limits(y=c(0,96),x=c(0,96))+ labs(color="Cluster") + theme_void() 

spatialPlot


png(file="Cluster_location.png",width=500,height=500)

spatialPlot 

dev.off()


pdf(file="Cluster_location.pdf",width=8,height=8)

spatialPlot 

dev.off()



#### 提取特定Cluster ####

rm(list=ls())

setwd(".../data_analysis/mE155S1-16H_split")

data_all <- read.table("data_merge_all_EM_part.txt", header = T)

cell_index  <-  c(data_all$V1[which(data_all$cluster == 1)])

write.table(cell_index, "Cluster_1_cell_index.txt", quote = FALSE, 
            sep = "\t", row.names = FALSE, col.names = FALSE)


data_all <- read.table("data_merge_all_EM_part.txt", header = T)

cell_index  <-  c(data_all$V1[which(data_all$cluster == 2)])

write.table(cell_index, "Cluster_2_cell_index.txt", quote = FALSE, 
            sep = "\t", row.names = FALSE, col.names = FALSE)


data_all <- read.table("data_merge_all_EM_part.txt", header = T)

cell_index  <-  c(data_all$V1[which(data_all$cluster == 3)])

write.table(cell_index, "Cluster_3_cell_index.txt", quote = FALSE, 
            sep = "\t", row.names = FALSE, col.names = FALSE)


data_all <- read.table("data_merge_all_EM_part.txt", header = T)

cell_index  <-  c(data_all$V1[which(data_all$cluster == 4)])

write.table(cell_index, "Cluster_4_cell_index.txt", quote = FALSE, 
            sep = "\t", row.names = FALSE, col.names = FALSE)


data_all <- read.table("data_merge_all_EM_part.txt", header = T)

cell_index  <-  c(data_all$V1[which(data_all$cluster == 5)])

write.table(cell_index, "Cluster_5_cell_index.txt", quote = FALSE, 
            sep = "\t", row.names = FALSE, col.names = FALSE)


####*************************************************************************####
####*************************************************************************####