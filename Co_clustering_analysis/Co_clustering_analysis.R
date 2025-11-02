library(openxlsx)
library(ggplot2)

rm(list=ls())


#### genebody ####

if(FALSE){

index1 <- read.table(".../data/100kb_bin/mE75spM_6_6_ref_barcode_filter.txt", 
                     header=F)


for(i in 1:2400){
  
  if(file.exists(paste0(".../data/feature_level/genebody/mE75spM_6_6_split/",index1[i,1],
                        "_genebody_hydroxymethylation_level_new.txt"))){
    
    spatial_methyl <- read.table(paste0(".../data/feature_level/genebody/mE75spM_6_6_split/",index1[i,1],
                                        "_genebody_hydroxymethylation_level_new.txt"), header=T)
    
    
    spatial_methyl_filter <- spatial_methyl[]
    
    assign(paste0("spatial_methyl_sample", i), spatial_methyl)
    
  }
  
}



index2 <- read.table(".../data/multi_omics_gastrulation/multi_omics_gastrulation_met_dir.txt",
                     header=F)

for(j in 1:390){
  
  if(file.exists(paste0(".../data/multi_omics_gastrulation/feature_level/genebody/mE75/",
                        index2[j,1], "_genebody_hydroxymethylation_level_new.txt"))){
    
    multi_omics <- read.table(paste0(".../data/multi_omics_gastrulation/feature_level/genebody/mE75/",
                                     index2[j,1], "_genebody_hydroxymethylation_level_new.txt"), header=T)
    
    assign(paste0("multi_omics_sample", j), multi_omics)
    
    
  }
  
}



#### Co-clustering  genebody ####

## spatial hydroxymethylation 


spatial_methyl_sample_level_merge <- data.frame(gene_name = spatial_methyl_sample1$gene_name )

for(i in 1:2400){
  
  if(exists(paste0("spatial_methyl_sample",i))){
    
    spatial_methyl_sample_filter <- get(paste0("spatial_methyl_sample",i))[,c(1,5)]
    
    spatial_methyl_sample_level_merge <- merge(spatial_methyl_sample_level_merge, spatial_methyl_sample_filter, by="gene_name" )
    
  }
  
}



spatial_methyl_sample_CpG_merge <- data.frame(gene_name = spatial_methyl_sample1$gene_name )

for(i in 1:2400){
  
  if(exists(paste0("spatial_methyl_sample",i))){
    
    spatial_methyl_sample_filter <- get(paste0("spatial_methyl_sample",i))[,c(2,5)]
    
    spatial_methyl_sample_CpG_merge <- merge(spatial_methyl_sample_CpG_merge, spatial_methyl_sample_filter, by="gene_name" )
    
  }
  
}



order <- c()

for(i in 1:2400){
  
  order[i] <- exists(paste0("spatial_methyl_sample",i))
  
}

colnames(spatial_methyl_sample_level_merge) <- c("gene_name",index1$V1[order])


colnames(spatial_methyl_sample_CpG_merge) <- c("gene_name",index1$V1[order])


setwd(".../data/genebody/results")


write.table(spatial_methyl_sample_level_merge, "spatial_methyl_sample_level_merge.txt")

write.table(spatial_methyl_sample_CpG_merge, "spatial_methyl_sample_CpG_merge.txt")

}


#### used  01 -- spatial methyl genebody extra-embryonic and embryonic ####

rm(list=ls())

setwd(".../data/genebody/results")

spatial_methyl_sample_level_merge <- read.table("spatial_methyl_sample_level_merge.txt",header = TRUE)

spatial_methyl_sample_CpG_merge <- read.table("spatial_methyl_sample_CpG_merge.txt",header = TRUE)


spatial_methyl_sample_level_merge[spatial_methyl_sample_CpG_merge < 1] <- -1


#### filter cells ####

data_T_and_F <- spatial_methyl_sample_level_merge[,2:length(colnames(spatial_methyl_sample_level_merge))] == -1 

data_T_and_F[which(data_T_and_F=="TRUE")] <- 1

spatial_methyl_sample_level_merge_filter <-  spatial_methyl_sample_level_merge[,c(1,which(colSums(data_T_and_F[,1:length(colnames(data_T_and_F))]) <= round(56847*0.9)) +1)]



#### filter feature ####

data_T_and_F1 <- spatial_methyl_sample_level_merge_filter[,2:length(colnames(spatial_methyl_sample_level_merge_filter))] == -1 

data_T_and_F1[which(data_T_and_F1=="TRUE")] <- 1

spatial_methyl_sample_level_merge_filter1 <-  spatial_methyl_sample_level_merge_filter[which(rowSums(data_T_and_F1[,1:length(colnames(data_T_and_F1))]) <= round(1337*0.9)),]


#### filter by embryo cells ####

ref_index <- read.table(".../data/100kb_bin/mE75spM_6_6_split_ref_barcode.txt", header=F)

EM_index <- read.table(".../data/100kb_bin/results/cell_filter_EM_new.txt", header=T)

other_index <- read.table(".../data/100kb_bin/results/cell_filter_other.txt", header=T)

cavity_index <- read.table(".../data/100kb_bin/results/mE75_6_6_cell_filter_4th_rev.txt", header=T)


ref_index$coord <- paste0(ref_index$V3, "x", ref_index$V2)

cell_index <- ref_index$V1[-which(ref_index$coord %in% c(other_index$coord, cavity_index$coord))]


spatial_methyl_sample_level_merge_filter1_EM <- spatial_methyl_sample_level_merge_filter1[,
                                                                                          c(1,which(colnames(spatial_methyl_sample_level_merge_filter1) %in% cell_index ))]
#cell number : 130

gene_index <- read.table("/mnt/server1/data2/tangym/p300s/tangym/mouse_embryo/RNA_results/reference/mm39/GRCm39_gencode_vM33_annotation_filter.txt",
                         header=F)

spatial_methyl_sample_level_merge_filter1_EM$gene_name <- gene_index$V4[spatial_methyl_sample_level_merge_filter1_EM$gene_name]


#### save data ####

write.table(spatial_methyl_sample_level_merge_filter1, "spatial_methyl_sample_level_merge_filter1.txt")

write.table(spatial_methyl_sample_level_merge_filter1_EM, "spatial_methyl_sample_level_merge_filter1_EM.txt")




#### used genebody -- multi omics gastrulation ####

if(FALSE){

multi_omics_sample_level_merge <- data.frame(gene_name = multi_omics_sample1$gene_name)

for(i in 1:390){
  
  if(exists(paste0("multi_omics_sample",i))){
    
    multi_omics_sample_filter <- get(paste0("multi_omics_sample",i))[,c(1,5)]
    
    multi_omics_sample_level_merge <- merge(multi_omics_sample_level_merge, multi_omics_sample_filter, by="gene_name" )
    
  }
  
}



multi_omics_sample_CpG_merge <- data.frame(gene_name = multi_omics_sample1$gene_name)

for(i in 1:390){
  
  if(exists(paste0("multi_omics_sample",i))){
    
    multi_omics_sample_filter <- get(paste0("multi_omics_sample",i))[,c(2,5)]
    
    multi_omics_sample_CpG_merge <- merge(multi_omics_sample_CpG_merge, multi_omics_sample_filter, by="gene_name" )
    
  }
  
}


order1 <- c()

for(i in 1:390){
  
  order1[i] <- exists(paste0("multi_omics_sample",i))
  
}

colnames(multi_omics_sample_level_merge) <- c("gene_name",substr(index2$V1, 1, nchar(index2$V1)-7)[order1])

colnames(multi_omics_sample_CpG_merge) <- c("gene_name",substr(index2$V1, 1, nchar(index2$V1)-7)[order1])


setwd(".../data/genebody/results")

write.table(multi_omics_sample_level_merge, "multi_omics_sample_level_merge.txt")

write.table(multi_omics_sample_CpG_merge, "multi_omics_sample_CpG_merge.txt")

}


#### used 02 -- multi omics genebody ####

setwd(".../data/genebody/results")

multi_omics_sample_level_merge <- read.table("multi_omics_sample_level_merge.txt",header = TRUE)

multi_omics_sample_CpG_merge <- read.table("multi_omics_sample_CpG_merge.txt",header = TRUE)

multi_omics_sample_level_merge[multi_omics_sample_CpG_merge < 1] <- -1



#### filter cells ####

data_T_and_F <- multi_omics_sample_level_merge[,2:length(colnames(multi_omics_sample_level_merge))] == -1 

data_T_and_F[which(data_T_and_F=="TRUE")] <- 1

multi_omics_sample_level_merge_filter <-  multi_omics_sample_level_merge[,c(1,which(colSums(data_T_and_F[,1:length(colnames(data_T_and_F))]) <= round(56847*0.9)) +1)]



#### filter feature ####

data_T_and_F1 <- multi_omics_sample_level_merge_filter[,2:length(colnames(multi_omics_sample_level_merge_filter))] == -1 

data_T_and_F1[which(data_T_and_F1=="TRUE")] <- 1

multi_omics_sample_level_merge_filter1 <-  multi_omics_sample_level_merge_filter[which(rowSums(data_T_and_F1[,1:length(colnames(data_T_and_F1))]) <= round(320*0.9)),]
#cell number : 267

#### save data ####

write.table(multi_omics_sample_level_merge_filter1, "multi_omics_sample_level_merge_filter1.txt")



#### used merge two sets genebody ####

#merge_two_sets <- merge(spatial_methyl_sample_level_merge_filter1, multi_omics_sample_level_merge_filter1, by="bin_name")

merge_two_sets <- merge(spatial_methyl_sample_level_merge_filter1_EM, multi_omics_sample_level_merge_filter1, by="gene_name")

merge_two_sets[merge_two_sets == -1] <- 0.5



index <- c()

for (i in 1:15709){
  
  value1 <-  mean(as.numeric(merge_two_sets[i, which(colnames(merge_two_sets) %in% ref_index$V1[ref_index$coord %in% EM_index$coord])]))
  
  value2 <- mean(as.numeric(merge_two_sets[i, grep("ES1", colnames(merge_two_sets))]))
  
  gap <- abs(value1-value2)
  
  if(gap <= 0.05){
    
    index[i] <- "TRUE"
    
  } else {
    
    index[i] <- "FALSE"
  }
  
}



col_index1 <- which(colnames(merge_two_sets) %in% ref_index$V1[ref_index$coord %in% EM_index$coord])

col_index2 <- which( !(colnames(merge_two_sets) %in% c(colnames(merge_two_sets)[grep("ES1", colnames(merge_two_sets))],
                                                      ref_index$V1[ref_index$coord %in% EM_index$coord] ) ) )

index1 <- c()

for (i in 1:15709){
  
  value1 <-  mean(as.numeric(merge_two_sets[i, col_index1]))
  
  value2 <- mean(as.numeric(merge_two_sets[i, col_index2[2:length(col_index2)] ]))
  
  gap <- abs(value1-value2)
  
  if(gap > 0.2){
    
    index1[i] <- "TRUE"
    
  } else {
    
    index1[i] <- "FALSE"
  }
  
}



#### used merge two sets genebody new ####

#merge_two_sets <- merge(spatial_methyl_sample_level_merge_filter1, multi_omics_sample_level_merge_filter1, by="bin_name")

merge_two_sets <- merge(spatial_methyl_sample_level_merge_filter1_EM, multi_omics_sample_level_merge_filter1, by="gene_name")

merge_two_sets[merge_two_sets == -1] <- 0.5


index_spatial <- read.table(".../data/100kb_bin/mE75spM_6_6_ref_barcode_filter.txt", 
                     header=F)


cell_type <- read.table(".../data/100kb_bin/multi_omics_sample_metadata.txt",
                        header=T)

cell_type_filter <- cell_type[match(cell_type$sample[grep("ES1",cell_type$sample)], cell_type$sample), ] 


c("Ectoderm","Endoderm","Epiblast",
  "Mesoderm","Primitive_Streak","None")


index <- c()

for (i in 1:8690){
  
  value1 <-  mean(as.numeric(merge_two_sets[i, which(colnames(merge_two_sets) %in%
                                                       cell_type_filter$sample[cell_type_filter$lineage10x_2=="Ectoderm"]) ]))
  
  value2 <- mean(as.numeric(merge_two_sets[i, which(colnames(merge_two_sets) %in%
                                                      cell_type_filter$sample[cell_type_filter$lineage10x_2=="Endoderm"]) ]))
  
  gap <- abs(value1-value2)
  
  if(gap > 0.1){
    
    index[i] <- "TRUE"
    
  } else {
    
    index[i] <- "FALSE"
  }
  
}



index1 <- c()

for (i in 1:8690){
  
  value1 <-  mean(as.numeric(merge_two_sets[i, which(colnames(merge_two_sets) %in%
                                                       cell_type_filter$sample[cell_type_filter$lineage10x_2=="Ectoderm"]) ]))
  
  value2 <- mean(as.numeric(merge_two_sets[i, which(colnames(merge_two_sets) %in%
                                                      cell_type_filter$sample[cell_type_filter$lineage10x_2=="Mesoderm"])]))
  
  gap <- abs(value1-value2)
  
  if(gap > 0.1){
    
    index1[i] <- "TRUE"
    
  } else {
    
    index1[i] <- "FALSE"
  }
  
}



index2 <- c()

for (i in 1:8690){
  
  value1 <-  mean(as.numeric(merge_two_sets[i, which(colnames(merge_two_sets) %in%
                                                       cell_type_filter$sample[cell_type_filter$lineage10x_2=="Mesoderm"]) ]))
  
  value2 <- mean(as.numeric(merge_two_sets[i, which(colnames(merge_two_sets) %in%
                                                      cell_type_filter$sample[cell_type_filter$lineage10x_2=="Endoderm"]) ]))
  
  gap <- abs(value1-value2)
  
  if(gap > 0.1){
    
    index2[i] <- "TRUE"
    
  } else {
    
    index2[i] <- "FALSE"
  }
  
}



#### Only spatial hydroxymethylation genebody ####

merge_two_sets <- spatial_methyl_sample_level_merge_filter_EM1

write.table(spatial_methyl_sample_level_merge_filter_EM1, "spatial_methyl_sample_level_merge_filter_EM1.txt")

merge_two_sets[merge_two_sets == -1] <- 0.5



#### R packages ####
library(ggsci)
library(ggplot2)
library(scales)
library(RColorBrewer)
library(umap)
library(devtools)
library(cluster)


#### PCA ####

data_merge_end <- merge_two_sets

data_merge_end_filter <-  data_merge_end[which(index=="TRUE" & index1=="TRUE"),2:length(colnames(data_merge_end))]


data_matrix <- data_merge_end_filter

data_raw_normlized <- scale(t(data_matrix))

#data_raw_normlized_filter <-  data_raw_normlized[-which(is.nan(data_raw_normlized[,1])),]

PCA_result <- prcomp(data_raw_normlized) # Occupy much resource of memory

#5733*337


#### Kmeans ####


PCA_selection <- PCA_result$x[,1:20]

set.seed(100)

km_results <- kmeans(PCA_selection, 3, nstart = 25)

#fviz_cluster(km_results, PCA_selection) #显示聚类分布情况

PCA_selection_cluster <- as.data.frame(km_results$cluster)

colnames(PCA_selection_cluster)[1] <- "cluster"

data_cluster <- data.frame(V1 = rownames(PCA_selection_cluster),
                           V2 = PCA_selection_cluster$cluster)

summary(factor(data_cluster$V2))


write.table(data_cluster,"data_cluster_merge_all_two_sets.txt",
            row.names = FALSE,
            col.names = TRUE)

#### umap ####


config.params <- umap.defaults

config.params$random_state=100

config.params$min_dist=0.8 

config.params$n_neighbors=8 

umap_result <- umap(PCA_selection, config = config.params)

PCA1 <- umap_result$layout[,1]

PCA2 <- umap_result$layout[,2]


PCA1_output <- as.data.frame(PCA1)

PCA2_output <- as.data.frame(PCA2)


PCA1_output$PCA2 <- PCA2_output$PCA2

PCA1_output$PCA1 <- PCA1_output$PCA1


PCA_output <- PCA1_output


PCA_output$barcode <- rownames(PCA1_output)


## merge cluster

data_cluster_used <-  data_cluster

colnames(data_cluster_used) <- c("barcode", "cluster")

data_all <- merge(PCA_output, data_cluster_used, by="barcode")



## merge part

data_all$part <- "None"

data_all$part[which(data_all$barcode %in% ref_index$V1[ref_index$coord %in% EM_index$coord])] <- "Embryonic"

data_all$part[-which(data_all$barcode %in% ref_index$V1[ref_index$coord %in% EM_index$coord])] <- "Extra-embryonic"

data_all$part[grep("ES1",data_all$barcode)] <- "Multi-omics embryonic cells"


write.table(data_all,"data_cluster_PCA_all_two_sets.txt",
            row.names = FALSE,
            col.names = TRUE)




#### plot with cluster ####

col <- c("#00468BB2", "#ED0000B2", "#42B540B2", "#0099B4B2", "#925E9FB2", "#FDAF91B2", "#AD002AB2", "#ADB6B6B2")

data_plot <- data.frame(PCA1=data_all$PCA1, 
                        PCA2=data_all$PCA2, 
                        data_cluster=factor(data_all$cluster,level=c(1:8)))

plot_scatter <- ggplot(data=data_plot, aes(x=PCA1, y=PCA2, 
                                            fill = data_cluster)) + 
                geom_point(size=4, alpha=1,shape=21, stroke=0.1) +
                scale_fill_manual(values = col) 

# stat_ellipse(aes(fill=data_group2), geom='polygon', type="norm",
#              level=0.68, alpha=0.5, show.legend = F) 

plot_scatter


plot_scatter <- plot_scatter + theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())+
  theme(axis.ticks = element_line(size = 0.75, color ="black"),
        axis.ticks.length=unit(1.25, 'mm'))


plot_scatter <- plot_scatter + theme(legend.title = element_text(size=18, family = "sans", face = "plain"))+
  theme(legend.text = element_text(size=16, family = "sans", face = "plain"))+
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
  
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05)), limits = c(-5.5, 9),
                     breaks = c( -4, 0, 4, 8),
                     labels = c("-4", "0", "4","8")) +  
  
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)), limits = c(-6, 8),
                     breaks = c( -5, 0, 5),
                     labels = c("-5", "0", "5")) +
  
  # theme(legend.position=c(0.97, 0.97),
  #       legend.justification = c(0.97, 0.97)) +
  # 
  # theme(legend.background = element_rect(fill="white",colour="black")) +
  
  guides(fill=guide_legend(title = "Cluster"))

plot_scatter



png(file="UMAP_cluster_genebody.png",width=600,height=600,res=100)

plot_scatter 

dev.off()


pdf(file="UMAP_cluster_genebody.pdf",width=6,height=6)

plot_scatter 

dev.off()



#### used plot with group ####

col <- c("#BC3C29B2", "#0072B5B2")

data_plot <- data.frame(PCA1=data_all$PCA1, 
                        PCA2=data_all$PCA2, 
                        data_group=factor(data_all$group))

plot_scatter <- ggplot(data=data_plot, aes(x=PCA1, y=PCA2, 
                                           fill = data_group)) + 
  geom_point(size=4, alpha=1,shape=21, stroke=0.05) +
  scale_fill_manual(values = col) 

# stat_ellipse(aes(fill=data_group2), geom='polygon', type="norm",
#              level=0.68, alpha=0.5, show.legend = F) 

plot_scatter


plot_scatter <- plot_scatter + theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())+
  theme(axis.ticks = element_line(size = 0.75, color ="black"),
        axis.ticks.length=unit(1.25, 'mm'))


plot_scatter <- plot_scatter + theme(legend.title = element_text(size=18, family = "sans", face = "plain"))+
  theme(legend.text = element_text(size=16, family = "sans", face = "plain"))+
  theme(plot.title = element_text(size = 18, family = "sans", face = "plain", hjust = 0.5)) +
  theme(axis.title.x = element_text(size = 18, family = "sans", face = "plain")) +
  theme(axis.title.y = element_text(size = 18, family = "sans", face = "plain"))+
  theme(axis.text.x = element_text(size = 18, color = "black", face = "plain",
                                   family = "sans", colour = "black")) +
  theme(axis.text.y = element_text(size = 18, color = "black", face = "plain",
                                   family = "sans", colour = "black")) +  
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1.25, 
                                    linetype = "solid")) + 
  labs(title = "",x = "UMAP1", y = "UMAP2") +
  theme(plot.margin = unit(c(0.5,0.5,0.2,0.2),"cm"))+
  
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05)), limits = c(-5.1, 5.3),
                     breaks = c( -4, 0, 4),
                     labels = c("-4", "0", "4")) +  
  
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)), limits = c(-8.5, 8),
                     breaks = c( -6, 0, 6),
                     labels = c("-6", "0", "6")) +
  
  # theme(legend.position=c(0.97, 0.97),
  #       legend.justification = c(0.97, 0.97)) +
  # 
  # theme(legend.background = element_rect(fill="white",colour="black")) +
  
  guides(fill=guide_legend(title = "Group"))

plot_scatter



png(file="UMAP_group_genebody_new.png",width=900,height=600,res=100)

plot_scatter 

dev.off()


pdf(file="UMAP_group_genebody_new.pdf",width=9,height=6)

plot_scatter 

dev.off()



#### used plot with part ####

col <- c("#E18727B2", "#20854EB2", "#7876B1B2")

data_plot <- data.frame(PCA1=data_all$PCA1, 
                        PCA2=data_all$PCA2, 
                        data_group=factor(data_all$part))

plot_scatter <- ggplot(data=data_plot, aes(x=PCA1, y=PCA2, 
                                           fill = data_group)) + 
  geom_point(size=4, alpha=1,shape=21, stroke=0.05) +
  scale_fill_manual(values = col) 

# stat_ellipse(aes(fill=data_group2), geom='polygon', type="norm",
#              level=0.68, alpha=0.5, show.legend = F) 

plot_scatter



plot_scatter <- plot_scatter + theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())+
  theme(axis.ticks = element_line(size = 0.75, color ="black"),
        axis.ticks.length=unit(1.25, 'mm'))


plot_scatter <- plot_scatter + theme(legend.title = element_text(size=18, family = "sans", face = "plain"))+
  theme(legend.text = element_text(size=16, family = "sans", face = "plain"))+
  theme(plot.title = element_text(size = 18, family = "sans", face = "plain", hjust = 0.5)) +
  theme(axis.title.x = element_text(size = 18, family = "sans", face = "plain")) +
  theme(axis.title.y = element_text(size = 18, family = "sans", face = "plain"))+
  theme(axis.text.x = element_text(size = 18, color = "black", face = "plain",
                                   family = "sans", colour = "black")) +
  theme(axis.text.y = element_text(size = 18, color = "black", face = "plain",
                                   family = "sans", colour = "black")) +  
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1.25, 
                                    linetype = "solid")) + 
  labs(title = "",x = "UMAP1", y = "UMAP2") +
  theme(plot.margin = unit(c(0.5,0.5,0.2,0.2),"cm"))+
  
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05)), limits = c(-5.1, 5.3),
                     breaks = c( -4, 0, 4),
                     labels = c("-4", "0", "4")) +  
  
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)), limits = c(-8.5, 8),
                     breaks = c( -6, 0, 6),
                     labels = c("-6", "0", "6")) +
  
  # theme(legend.position=c(0.97, 0.97),
  #       legend.justification = c(0.97, 0.97)) +
  # 
  # theme(legend.background = element_rect(fill="white",colour="black")) +
  
  guides(fill=guide_legend(title = "Group"))

plot_scatter



png(file="UMAP_part_genebody_new.png",width=1000,height=600,res=100)

plot_scatter 

dev.off()


pdf(file="UMAP_part_genebody_new.pdf",width=10,height=6)

plot_scatter 

dev.off()


#### used  plot with location two parts ####

setwd(".../data/genebody/results")

data_all <- read.table("data_cluster_PCA_all_two_sets.txt",header=T)


ref_barcode <- read.table(".../data/genebody/results/mE75spM_6_6_split_ref_barcode.txt",
                          header=F)


ref_barcode$coord <- paste0(ref_barcode$V3, "x", ref_barcode$V2)


# EM_and_TPs_index <- read.table(".../data/100kb_bin/results/cell_filtered_EM_and_TPs.txt",
#                                header=T)
# 
# 
# ref_barcode <- ref_barcode[which(ref_barcode$coord %in% EM_and_TPs_index$coord ),]



colnames(ref_barcode) <- c("barcode", "nA", "nB","coord")

data_merge <- merge(ref_barcode, data_all, by="barcode", all=F )

data_merge_filter <- data_merge

data_merge_filter$part <- factor(data_merge_filter$part )


data_merge_filter$nA <- as.numeric(factor( data_merge_filter$nA, label=1:length(levels(factor(data_merge_filter$nA))) ))

data_merge_filter$nB <- as.numeric(factor( data_merge_filter$nB, label=1:length(levels(factor(data_merge_filter$nB))) ))



col <- c("#E18727B2", "#20854EB2")

info_plot <- ggplot(data=data_merge_filter, 
                    aes(x=nB,y=nA,color=part))+
  geom_point(shape=16,size=4, show.legend=T,stroke=0)+
  scale_color_manual(values = col) +
  scale_y_continuous(breaks = seq(1,47,5)) +
  scale_x_continuous(transform = "reverse", breaks = seq(1,18,5)) +
  expand_limits(y=c(min(data_merge_filter$nA)-2,
                    max(data_merge_filter$nA)+2),
                x=c(min(data_merge_filter$nB-2),
                    max(data_merge_filter$nB))+2)+ 
  labs(fill="Group") + theme_void()+
  theme(legend.title = element_text(size=18, family="sans", face="plain",hjust=0))+
  theme(legend.text=element_text(size=16, family="sans", face="plain"))

info_plot


png(file="UMAP_location_EM_and_TPs_genebody_All.png",width=500,height=1000,res=100)

info_plot 

dev.off()


pdf(file="UMAP_location_EM_and_TPs_genebody_All.pdf",width=5,height=10)

info_plot 

dev.off()



