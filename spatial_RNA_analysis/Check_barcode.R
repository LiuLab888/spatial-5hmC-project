# Set R library path
.libPaths("/path/to/Rlib/4.3")

# Load required libraries
library(getopt)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

# Define command-line arguments
command <- matrix(c(
  "combine_barcode", "b", 1, "character",
  "work_path", "p", 1, "character",
  "file_name", "n", 1, "character",
  "help", "h", 0, "logical"
), byrow = TRUE, ncol = 4)

# Parse arguments
args <- getopt(command)
if (!is.null(args$help) || is.null(args$combine_barcode) || is.null(args$work_path) || is.null(args$file_name)) {
  cat(paste(getopt(command, usage = TRUE), "\n"))
  q()
}

# Assign arguments to variables
work_path <- args$work_path        # e.g., "/path/to/output_directory"
input_name <- args$file_name       # e.g., "sample_barcode.txt"
ref_path <- args$combine_barcode   # e.g., "/path/to/combine_barcode.txt"

# Set working directory
setwd(work_path)

# Load reference barcode data
ref_barcode <- read.table(ref_path, stringsAsFactors = FALSE)
our_barcode <- read.table(input_name, stringsAsFactors = FALSE)$V1

# Print total barcode count
cat("Total barcode count:", length(our_barcode), "\n")

# Filter barcodes that match reference list
our_barcode <- our_barcode[our_barcode %in% ref_barcode$V1]
cat("Barcode count in reference list:", length(our_barcode), "\n")

# Convert barcode to factor and count occurrences
our_barcode <- factor(our_barcode, levels = ref_barcode$V1)
ref_barcode$count <- table(our_barcode)

# Determine number of barcode channels
nchannel <- max(ref_barcode$V2)
all_combi.m <- matrix(0, nchannel, nchannel)
row.names(all_combi.m) <- paste0("A", 1:nchannel)
colnames(all_combi.m) <- paste0("B", 1:nchannel)

# Populate matrix with barcode pair counts
lapply(1:nrow(ref_barcode), function(i) {
  all_combi.m[ref_barcode$V2[i], ref_barcode$V3[i]] <<- ref_barcode$count[i]
})

# Reorder columns (from B96 to B1)
all_combi.m <- all_combi.m[, nchannel:1]

# Determine heatmap size
nwid <- ifelse(nchannel > 90, 15, 10)

# Generate PDF heatmaps
pdf("0_barcode_combination.pdf", width = nwid, height = nwid)
Heatmap(all_combi.m, cluster_rows = FALSE, cluster_columns = FALSE)
Heatmap(log2(all_combi.m + 1), cluster_rows = FALSE, cluster_columns = FALSE)
dev.off()

# Generate SVG heatmap
nlwd <- ifelse(nchannel > 90, 2, 5)
svg("1_barcode_combination.svg", width = 6, height = 6)
Heatmap(all_combi.m, cluster_rows = FALSE, cluster_columns = FALSE, 
        rect_gp = gpar(col = "white", lwd = nlwd),
        col = colorRamp2(c(0, quantile(all_combi.m, 0.9) / 4, max(all_combi.m)), c("blue", "green", "red")),
        show_heatmap_legend = FALSE, show_row_names = FALSE, show_column_names = FALSE)
dev.off()

cat("Processing completed. Heatmaps saved.\n")
