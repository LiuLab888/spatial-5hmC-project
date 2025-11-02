# List of sample names, stages, and resolution values
samples <- c("Sample1", "Sample2", "Sample3")  # Add all sample names
stage <- c("E6.5", "E7.5", "E8.5")  # Define corresponding stages
reso <- c(0.5, 0.7, 1)  # Define corresponding resolution values for clusters

# Initialize an empty list to hold Seurat objects
seurat_objects <- list()

# Loop through each sample
for (n in 1:length(samples)) {
  # Construct file path dynamically
  file_path <- paste0("/path/to/data/", samples[n], "/1_stpipeline/data.tsv")
  
  # Read the data and create Seurat object
  sample_data <- read.table(file_path)
  sample_data <- as.matrix(sample_data)
  
  # Create Seurat object and process metadata
  seurat_obj <- CreateSeuratObject(counts = t(sample_data), project = samples[n], assay = "RNA")
  seurat_obj$coord <- paste0(substr(samples[n], 1, 3), "_", row.names(seurat_obj@meta.data))
  seurat_obj$nA <- as.integer(str_split_fixed(seurat_obj$coord, "_|x", 4)[, 2])
  seurat_obj$nB <- as.integer(str_split_fixed(seurat_obj$coord, "_|x", 4)[, 3])
  seurat_obj$stage <- stage[n]
  
  # Filter low-count cells
  seurat_obj <- seurat_obj[, seurat_obj$nCount_RNA > 300]
  
  # Apply SCTransform
  seurat_obj <- SCTransform(seurat_obj, assay = "RNA", verbose = FALSE, variable.features.n = 3000)
  
  # Add to the list of Seurat objects
  seurat_objects[[samples[n]]] <- seurat_obj
}

# Integration steps
Object.list <- seurat_objects

# Select integration features
The.features <- SelectIntegrationFeatures(object.list = Object.list, nfeatures = 3000)

# Prepare for integration
The.list <- PrepSCTIntegration(object.list = Object.list, anchor.features = The.features, verbose = FALSE)

# Find anchors
The.anchors <- FindIntegrationAnchors(object.list = The.list, normalization.method = "SCT", anchor.features = The.features)

# Integrate data
Combined_data.seurat <- IntegrateData(anchorset = The.anchors, normalization.method = "SCT")
DefaultAssay(Combined_data.seurat) <- "integrated"

# Run PCA
Combined_data.pca <- RunPCA(Combined_data.seurat, npcs = 30, verbose = FALSE)

# Run UMAP
Combined_data.umap <- RunUMAP(Combined_data.pca, reduction = "pca", dims = 1:30)
Combined_data.umap <- FindNeighbors(Combined_data.umap, reduction = "pca", dims = 1:30)

# Find clusters with the dynamic resolution
Combined_data.umap <- FindClusters(Combined_data.umap, verbose = FALSE, resolution = reso[n])

# Create UMAP plots
pdf(file="Combined_data.umap.bySource.pdf", width=7, height=5)
DimPlot(Combined_data.umap, reduction = "umap", label = TRUE, shuffle = TRUE, group.by = "orig.ident")
dev.off()

pdf(file="Combined_data.umap.bySeuratcluster.pdf", width=7, height=5)
DimPlot(Combined_data.umap, reduction = "umap", label = TRUE, shuffle = TRUE, group.by = "seurat_clusters")
dev.off()

pdf(file="Combined_data.umap.stage.pdf", width=7, height=5)
DimPlot(Combined_data.umap, reduction = "umap", label = TRUE, shuffle = TRUE, group.by = "stage")
dev.off()

# Set the identity of the Seurat object to Subcluster_L3
Idents(Combined_data.umap) <- "seurat_clusters"

# Find all markers
AllMarkers <- FindAllMarkers(Combined_data.umap, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.25)
write.csv(AllMarkers, "AllMarkers.csv", quote = FALSE)

# Extract markers for each specific cluster
celltype1.marker <- AllMarkers$gene[AllMarkers$cluster == "celltype1"]
celltype2.marker <- AllMarkers$gene[AllMarkers$cluster == "celltype2"]
celltype3.marker <- AllMarkers$gene[AllMarkers$cluster == "celltype3"]

# Collect all celltype markers
CelltypeMarkers <- list(celltype1.marker, celltype2.marker, celltype3.marker)

# Perform GO enrichment analysis for each set of markers
for (n in seq_along(CelltypeMarkers)) {
  print(paste0(Sys.time(), " ::: Start GO enrichment for Cluster ", n))
  
  MARKER <- CelltypeMarkers[[n]]
  GO <- enrichGO(gene = MARKER, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "All", pAdjustMethod = "BH", qvalueCutoff = 0.05)
  
  assign(x = paste0("Celltype", n, ".GO"), value = GO)
  
  # Plot GO results
  PLOT <- dotplot(GO, showCategory = 10, split = "ONTOLOGY") +
    labs(title = paste0("Celltype ", n, " GO Results")) +
    facet_grid(ONTOLOGY ~ ., scales = "free") +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) +
    theme(plot.margin = margin(10, 20, 10, 20))
  
  # Print plot
  print(PLOT)
}
