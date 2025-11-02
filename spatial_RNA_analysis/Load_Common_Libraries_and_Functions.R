# This script serves as a default setup and pre-processing step before conducting any analysis. 
# It includes:
# 1. Loading all necessary libraries for the analysis.
# 2. Defining commonly used functions such as:
#    - Retrieving Ensembl IDs from gene symbols.
#    - Generating UMAP feature plots for multiple genes.
#    - Creating spatial feature plots for gene expression.
#    - Comparing cell type compositions.
#    - Plotting spatial distributions of clusters in PDF or PNG formats.
# 3. The functions and setup are intended to streamline the analysis process and can be 
#    run before any other specific analysis steps in your project.

#--------------------------------------------------------
# 1. Load Necessary Libraries
#--------------------------------------------------------
library(Seurat)            # For handling single-cell RNA-seq data
library(SeuratWrappers)    # Additional Seurat functionality
library(dplyr)             # For data manipulation
library(tidyr)             # For data tidying
library(org.Mm.eg.db)      # For mouse gene annotations
library(AnnotationDbi)     # For annotation database interface
library(stringr)           # For string manipulations
library(ggplot2)           # For plotting
library(patchwork)         # For combining ggplot2 plots
library(viridis)           # For color scales
library(RColorBrewer)      # For color palettes
library(ggrepel)           # For adding labels to ggplot
library(clusterProfiler)  # For GO enrichment analysis
library(ReactomePA)       # For Reactome pathway analysis
library(ComplexHeatmap)   # For creating complex heatmaps
library(pheatmap)         # For heatmaps
library(reshape2)          # For reshaping data
library(forcats)           # For factor manipulation
library(ggpubr)            # For publication-ready plots

#--------------------------------------------------------
# 2. Functions Definition
#--------------------------------------------------------

# 2.1 Retrieve Ensembl IDs from Gene Symbols
get_ensembl_id <- function(gene_symbols) {
  ensembl_ids <- mapIds(
    x = org.Mm.eg.db,
    keys = gene_symbols,
    column = "ENTREZID",
    keytype = "SYMBOL",
    multiVals = "first"
  )
  result <- data.frame(
    gene_symbol = names(ensembl_ids),
    ensembl_id = ensembl_ids,
    stringsAsFactors = FALSE
  )
  return(result)
}

# 2.2 Create UMAP Feature Plots for Multiple Gene Symbols
plot_gene_expression_umap_multiple <- function(seurat_object, gene_symbols_string, assay = "SCT", n_col = 3, use_ensembl = TRUE) {
  plots <- list()
  gene_symbols <- unlist(strsplit(gene_symbols_string, ","))
  
  # Get Ensembl IDs if requested
  ensembl_ids <- if (use_ensembl) get_ensembl_id(gene_symbols) else NULL
  
  for (gene_symbol in gene_symbols) {
    ensembl_id <- if (use_ensembl) {
      ensembl_ids$ensembl_id[ensembl_ids$gene_symbol == gene_symbol]
    } else {
      gene_symbol
    }
    
    if (is.na(ensembl_id) || !(ensembl_id %in% rownames(seurat_object))) {
      warning(paste("Gene not found in the Seurat object:", gene_symbol))
      next
    }
    
    # Fetch and combine data
    expression_data <- FetchData(seurat_object, vars = ensembl_id, assay = assay)
    metadata <- seurat_object@meta.data
    metadata$umap_1 <- seurat_object@reductions$umap@cell.embeddings[, 1]
    metadata$umap_2 <- seurat_object@reductions$umap@cell.embeddings[, 2]
    plot_data <- cbind(metadata, expression_data)
    
    # Create the plot
    p <- ggplot(plot_data, aes_string(x = "umap_1", y = "umap_2", color = ensembl_id)) +
      geom_point(size = 0.5) +
      scale_color_viridis(option = "C", direction = 1, na.value = "lightgrey") +
      theme_minimal() +
      labs(title = paste("Expression of", gene_symbol), x = "UMAP 1", y = "UMAP 2", color = "Expression") +
      theme(legend.position = "right", legend.key.size = unit(0.5, "cm"), plot.title = element_text(hjust = 0.5))
    
    plots[[gene_symbol]] <- p
  }
  
  # Combine and return the plot
  combined_plot <- wrap_plots(plots, ncol = n_col) + plot_layout(guide = "keep")
  return(combined_plot)
}

# 2.3 Create Spatial Feature Plots for Multiple Gene Symbols
plot_gene_expression_spatial_multiple <- function(seurat_object, gene_symbols_string, n_col = 3, samplename = "", use_ensembl = TRUE) {
  plots <- list()
  gene_symbols <- unlist(strsplit(gene_symbols_string, ","))
  
  # Get Ensembl IDs if requested
  ensembl_ids <- if (use_ensembl) get_ensembl_id(gene_symbols) else NULL
  
  for (gene_symbol in gene_symbols) {
    ensembl_id <- if (use_ensembl) {
      ensembl_ids$ensembl_id[ensembl_ids$gene_symbol == gene_symbol]
    } else {
      gene_symbol
    }
    
    if (is.na(ensembl_id) || !(ensembl_id %in% rownames(seurat_object))) {
      warning(paste("Gene not found in the Seurat object:", gene_symbol))
      next
    }
    
    # Extract spatial data and gene expression
    spatial_data <- seurat_object@meta.data[, c("nA", "nB")]
    expression_data <- FetchData(seurat_object, vars = ensembl_id)
    
    # Combine data for plotting
    plot_data <- cbind(spatial_data, expression_data)
    colnames(plot_data)[3] <- "Expression"
    
    # Create the plot
    p <- ggplot(plot_data, aes(x = nA, y = nB, color = Expression)) +
      geom_point(size = 1.5, shape = 15) +
      scale_color_viridis(option = "C", direction = 1, na.value = "lightgrey") +
      theme_bw() +
      labs(title = paste(gene_symbol, "Expression"), x = "nA", y = "nB", color = "Expression") +
      xlim(0, 96) + ylim(0, 96)
    
    plots[[gene_symbol]] <- p
  }
  
  # Combine and return the plot
  combined_plot <- wrap_plots(plots, ncol = n_col) + 
    plot_layout(guides = 'keep') +
    plot_annotation(title = samplename, theme = theme(plot.title = element_text(hjust = 0.5)))
  
  return(combined_plot)
}

# 2.4 Compare Cell Type Composition
compare_cell_composition <- function(seurat_object, mycluster = "seurat_clusters", refcluster = "") {
  DF <- data.frame(MyCluster = seurat_object@meta.data[[mycluster]], RefCluster = seurat_object@meta.data[[refcluster]])
  DF <- DF[!is.na(DF$RefCluster), ]
  
  composition <- DF %>%
    group_by(MyCluster, RefCluster) %>%
    summarise(total_cells = n(), .groups = 'drop') %>%
    group_by(MyCluster) %>%
    mutate(composition = total_cells / sum(total_cells) * 100) %>%
    ungroup() %>%
    dplyr::select(MyCluster, RefCluster, composition) %>%
    arrange(MyCluster)
  
  composition_output <- composition %>%
    filter(composition > 10) %>%
    mutate(composition = paste0(round(composition, 0), "%"))
  
  return(composition_output)
}

# 2.5 Plot Spatial Distribution of Clusters (PDF Format)
plot_spatial_cluster_pdf <- function(seurat_object, filename, cluster_colors, orig_seuobject,
                                     cluster_col_by = "seurat_clusters", save_pdf = TRUE) {
  # Extract cluster levels and map colors
  all_clusters <- levels(orig_seuobject@meta.data[[cluster_col_by]])
  cluster_color_mapping <- setNames(cluster_colors[1:length(all_clusters)], all_clusters)
  
  # Prepare data for plotting
  plot_data <- data.frame(
    nB = seurat_object$nB,
    nA = seurat_object$nA,
    clusters = seurat_object@meta.data[[cluster_col_by]]
  )
  
  # Map colors to clusters
  subset_clusters <- levels(seurat_object@meta.data[[cluster_col_by]])
  subset_colors <- cluster_color_mapping[subset_clusters]
  
  # Create the ggplot
  SeuPlot <- ggplot(plot_data, aes(x = nA, y = nB, color = clusters)) +
    geom_point(size = 2, shape = 15) +
    scale_color_manual(values = subset_colors) +
    xlim(0, 96) +
    ylim(0, 96) +
    theme_bw()
  
  # Save the plot to a PDF if save_pdf is TRUE
  if (save_pdf) {
    pdf(file = filename, width = 11, height = 8.5)
    print(SeuPlot)
    dev.off()
  } else {
    print(SeuPlot)  # Just display the plot if not saving
  }
}

# 2.6 Plot Spatial Distribution of Clusters (PNG Format)
plot_spatial_cluster_png <- function(seurat_object, filename, cluster_colors = fancy_palette, orig_seuobject,
                                     cluster_col_by = "seurat_clusters", save_png = TRUE) {
  # Extract cluster levels and map colors
  all_clusters <- levels(orig_seuobject@meta.data[[cluster_col_by]])
  cluster_color_mapping <- setNames(cluster_colors[1:length(all_clusters)], all_clusters)
  
  # Prepare data for plotting
  plot_data <- data.frame(
    nB = seurat_object$nB,
    nA = seurat_object$nA,
    clusters = seurat_object@meta.data[[cluster_col_by]]
  )
  
  # Map colors to clusters
  subset_clusters <- levels(seurat_object@meta.data[[cluster_col_by]])
  subset_colors <- cluster_color_mapping[subset_clusters]
  
  # Create the ggplot
  SeuPlot <- ggplot(plot_data, aes(x = nA, y = nB, color = clusters)) +
    geom_point(size = 2, shape = 15) +
    scale_color_manual(values = subset_colors) +
    xlim(0, 96) +
    ylim(0, 96) +
    theme_bw()
  
  # Save the plot to a PNG if save_png is TRUE
  if (save_png) {
    png(filename = filename, width = 900, height = 800)
    print(SeuPlot)
    dev.off()
  } else {
    print(SeuPlot)  # Just display the plot if not saving
  }
}

# 2.7 Load Spatial Data and Create Seurat Object
LoadSpatialData <- function(i) {
  count.m = read.table(my_files[i])  
  count.m = as.matrix(count.m)
  tmp = CreateSeuratObject(counts = t(count.m), project = Sample.ls[i], assay = "RNA")
  
  # Assign spatial coordinates
  tmp$coord = paste0(Sample.ls[i], "_", row.names(tmp@meta.data))
  tmp$nA = as.integer(str_split_fixed(tmp$coord, "_|x", n = Inf)[, 4])
  tmp$nB = as.integer(str_split_fixed(tmp$coord, "_|x", n = Inf)[, 5])
  
  # Subset and process data
  tmp.subset = tmp[, tmp$nCount_RNA > Filter.m[i]]
  tmp.subset <- SCTransform(tmp.subset, assay = "RNA", verbose = FALSE, 
                            variable.features.n = 3000, return.only.var.genes = FALSE)
  
  return(tmp.subset)
}

# 2.8 Create Volcano Plot
create_volcano_plot <- function(de_data, dataset_name, adjusted_p_value_cutoff, log2FC_cutoff) {
  ggplot(de_data, aes(x = avg_log2FC, y = -log10(p_val_adj), colour = Gene_Type)) + 
    geom_point(alpha = 0.65, size = 1) + 
    scale_color_manual(values = c("#546de5", "#d2dae2", "#ff4757")) + 
    geom_vline(xintercept = c(-log2FC_cutoff, log2FC_cutoff), lty = 4, col = "black", lwd = 0.8, linetype = "dashed") + 
    geom_hline(yintercept = -log10(adjusted_p_value_cutoff), lty = 4, col = "black", lwd = 0.8, linetype = "dashed") + 
    labs(x = "Log2 Fold Change", y = "-log10 Adjusted P-value", title = paste("Volcano Plot -", dataset_name)) + 
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), 
          legend.position = "right", 
          legend.title = element_blank()) + 
    geom_text_repel(
      data = subset(de_data, p_val_adj < adjusted_p_value_cutoff & abs(avg_log2FC) >= log2FC_cutoff),
      aes(label = Gene), size = 3, 
      box.padding = unit(0.3, "lines"), point.padding = unit(0.1, "lines"),
      segment.color = "black", show.legend = FALSE, max.overlaps = 5, segment.size = 0.2
    )
}

# 2.9 Function to add 'Gene_Type' column based on p-value and log2 fold change
assign_gene_type <- function(de_data, p_val_adj_cutoff, log2FC_cutoff) {
  Data <- de_data %>%
    mutate(Gene_Type = case_when(
      p_val_adj < p_val_adj_cutoff & avg_log2FC >= log2FC_cutoff ~ "Down",
      p_val_adj < p_val_adj_cutoff & avg_log2FC <= -log2FC_cutoff ~ "Up",
      TRUE ~ "no"
    )) %>%
    mutate(Gene = rownames(de_data))
  Data$p_val[Data$p_val == 0.000000e+00] <- 1e-300
  Data$p_val_adj[Data$p_val_adj == 0.000000e+00] <- 1e-300
  
  return(Data)
}


# 2.9 Define function to create a feature plot for multiple gene symbols
plot_gene_expression_umap_multiple <- function(seurat_object, gene_symbols_string, assay = "SCT", n_col = 4, use_ensembl = TRUE, point_size=0.5) {
  plots <- list()
  
  # Split the input string into a vector of gene symbols
  gene_symbols <- unlist(strsplit(gene_symbols_string, ","))
  
  # Map gene symbols to Ensembl IDs if use_ensembl is TRUE
  if (use_ensembl) {
    ensembl_ids <- get_ensembl_id(gene_symbols)
  } else {
    ensembl_ids <- data.frame(gene_symbol = gene_symbols, ensembl_id = gene_symbols)
  }
  
  # Loop through each gene symbol
  for (gene_symbol in gene_symbols) {
    ensembl_id <- ensembl_ids$ensembl_id[ensembl_ids$gene_symbol == gene_symbol]
    
    # If use_ensembl is FALSE, use the gene symbol directly
    if (!use_ensembl) {
      ensembl_id <- gene_symbol
    }
    
    # Check if the gene exists in the dataset
    if (!(ensembl_id %in% rownames(seurat_object))) {
      warning(paste("Gene not found in the Seurat object:", gene_symbol))
      next
    }
    
    # Get the expression data for the gene using its Ensembl ID or symbol
    expression_data <- FetchData(seurat_object, vars = ensembl_id, assay = assay)
    
    # Combine expression data with metadata
    metadata <- seurat_object@meta.data
    metadata$umap_1 <- seurat_object@reductions$umap@cell.embeddings[, 1]
    metadata$umap_2 <- seurat_object@reductions$umap@cell.embeddings[, 2]
    plot_data <- cbind(metadata, expression_data)
    
    # Create the ggplot for the current gene
    p <- ggplot(plot_data, aes_string(x = "umap_1", y = "umap_2", color = ensembl_id)) +
      geom_point(size = point_size) +
      scale_color_viridis(option = "C", direction = 1, na.value = "lightgrey") +  # Use viridis color scale
      theme_minimal() +
      labs(title = paste("Expression of", gene_symbol),
           x = "UMAP 1",
           y = "UMAP 2",
           color = "Expression") +
      theme(legend.position = "right",
            legend.key.size = unit(0.5, "cm"),  # Adjust size of the legend key
            plot.title = element_text(hjust = 0.5))  # Center the title
    
    # Add the plot to the list
    plots[[gene_symbol]] <- p
  }
  
  # Combine all plots into one with specified number of columns
  combined_plot <- wrap_plots(plots, ncol = n_col) + plot_layout(guide="keep")  # Adjust layout
  return(combined_plot)
}

# 2.10 Define function to create a feature spatial plot for multiple gene symbols
plot_gene_expression_spatial_multiple <- function(seurat_object, gene_symbols_string, n_col = 3, samplename = "", use_ensembl = TRUE, point_size=1.5) {
  plots <- list()
  
  # Split the input string into a vector of gene symbols
  gene_symbols <- unlist(strsplit(gene_symbols_string, ","))
  
  # Map gene symbols to Ensembl IDs if use_ensembl is TRUE
  if (use_ensembl) {
    ensembl_ids <- get_ensembl_id(gene_symbols)  # Ensure this function is defined elsewhere
  } else {
    ensembl_ids <- data.frame(gene_symbol = gene_symbols, ensembl_id = gene_symbols)
  }
  
  # Loop through each gene symbol
  for (gene_symbol in gene_symbols) {
    ensembl_id <- ensembl_ids$ensembl_id[ensembl_ids$gene_symbol == gene_symbol]
    
    # If use_ensembl is FALSE, use the gene symbol directly
    if (!use_ensembl) {
      ensembl_id <- gene_symbol
    }
    
    if (is.na(ensembl_id) || !(ensembl_id %in% rownames(seurat_object))) {
      warning(paste("Gene not found in the Seurat object:", gene_symbol))
      next
    }
    
    # Extract spatial coordinates (nA, nB) and gene expression
    spatial_data <- seurat_object@meta.data[, c("nA", "nB", "orig.ident")]
    expression_data <- FetchData(seurat_object, vars = ensembl_id)
    
    # Combine data for plotting
    plot_data <- cbind(spatial_data, expression_data)
    colnames(plot_data)[4] <- "Expression"
    
    # Create spatial plot
    p <- plot_data %>%
      filter(orig.ident == samplename) %>%
      ggplot(aes(x = nA, y = nB, color = Expression)) +
      geom_point(size = point_size, shape = 15) +
      scale_color_gradientn(colors = c("#f5f5f4", "#fdf7f9", "#fae6ec","#f0b6c6","#e27594", "#9a2043")) +
      theme_bw() +
      labs(title = paste(gene_symbol, "Expression"),
           x = "nA",
           y = "nB",
           color = "Expression") +
      xlim(0, 96) +
      ylim(0, 96)
    
    # Add the plot to the list
    plots[[gene_symbol]] <- p
  }
  
  # Combine all plots into one with specified number of columns
  combined_plot <- wrap_plots(plots, ncol = n_col) + 
    plot_layout(guides = 'keep') +  # Adjust layout
    plot_annotation(title = samplename, 
                    theme = theme(plot.title = element_text(hjust = 0.5)))
  
  return(combined_plot)
}
