library(getopt)
command=matrix(c("combine_barcode","b",1,"character",
                 "work_path","p",1,"character",
                 "file_name","n",1,"character",
                 "help","h",0,"logical"),byrow=T,ncol=4)

args=getopt(command)
if (!is.null(args$help) || is.null(args$combine_barcode) || is.null(args$work_path) || is.null(args$file_name)) {
  cat(paste(getopt(command, usage = T), "\n"))
  q()
}

work_path = args$work_path 
input_name = args$file_name 
refe_path = args$combine_barcode 

setwd(work_path)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
# 1
ref_barcode = read.table(refe_path,stringsAsFactors = F)
our_barcode = read.table(input_name,stringsAsFactors = F)


# 2
merge_barcode <- merge(ref_barcode, our_barcode, by="V1", all=TRUE)
ref_barcode$V2 <- merge_barcode[,2]
ref_barcode$V3 <- merge_barcode[,3]
ref_barcode$count <- merge_barcode[,4]


nchannel_A = max(ref_barcode$V2)
nchannel_B=max(ref_barcode$V3)
all_combi.m = matrix(0,nchannel_A,nchannel_B) # barcodeA*barcodeB
rownames(all_combi.m) = paste0("A",1:nchannel_A)
colnames(all_combi.m) = paste0("B",1:nchannel_B)
tmp <- lapply(1:dim(ref_barcode)[1], function(i) all_combi.m[ref_barcode$V2[i],ref_barcode$V3[i]] <<- ref_barcode$count[i])

# from B96 --> B1
all_combi.m = all_combi.m[,nchannel_B:1]
# form A96 --> A1
all_combi.m = all_combi.m[nchannel_A:1,]


all_combi.m[all_combi.m < 0.32] <- NA

pdf("hydroxymethylation_level.pdf",width = round(nchannel_B*12/96)+1,height = round(nchannel_A*12/96))
Heatmap(all_combi.m,cluster_rows = F,cluster_columns = F,
        rect_gp = gpar(col = "grey", lwd = 0.5), 
        show_row_names = F, show_column_names = F,
        col = colorRamp2(c(0.32, 0.49, 0.66, 0.83, 1.00),
                         c("#0D36A6", "#5992F2", "#FFC125","#CD0000","#FF0000")),
        heatmap_legend_param = list(legend_height = unit(3, "cm"),
                                    grid_width = unit(0.5, "cm"),
                                    title = "Methyl_level",
				    at = c(0.4, 0.6, 0.8, 1),
				    labels = c(0.4, 0.6, 0.8, 1),	
                                    title_gp = gpar(fontsize = 10, 
                                                    fontface = "bold"),
                                    labels_gp = gpar(fontsize = 10))) 
dev.off()

library(circlize)
svg("hydroxymethylation_level.svg",width = round(nchannel_B*12/96)+1,height = round(nchannel_A*12/96))
Heatmap(all_combi.m,cluster_rows = F,cluster_columns = F,
        rect_gp = gpar(col = "grey", lwd = 0.5),
        show_row_names = F, show_column_names = F,
        heatmap_legend_param = list(legend_height = unit(3, "cm"),
                                    grid_width = unit(0.5, "cm"),
                                    title = "Methyl_level",
                                    title_gp = gpar(fontsize = 10,
                                                    fontface = "bold"),
                                    labels_gp = gpar(fontsize = 10)))

dev.off()
# https://zhuanlan.zhihu.com/p/96444730

