library(Seurat)
library(tximport)
library(SeuratData)
library(SeuratWrappers)
library(Azimuth)
library(ggplot2)
library(patchwork)
library(reticulate)

# Change the following line to include your own directory
setwd('/Users/juanjovel/jj/data_analysis/jeffBiernaskie/RangiferTarandus/transcriptomeAssembly/compare_jeff_against_combined/evaluation_hybrid_transcriptome/hybrid_transcriptome/alevin_results/')
main_dir <- getwd()

# List all items in the specified path that match the pattern "output"
all_items <- list.files(main_dir, pattern = "output", full.names = TRUE)

# Filter to include only directories
directories <- all_items[file.info(all_items)$isdir]

# Initialize counter
i <- 1
seurat_list <- list()

# Loop through each directory and load the RDS file
for (sub_dir in directories) {
  subdir_name <- basename(sub_dir)
  rds_name <- paste0(subdir_name, ".RDS")
  rds_path <- file.path(sub_dir, rds_name)
  
  if (file.exists(rds_path)) {
    obj_name <- paste0("obj", i)
    seurat_list[[obj_name]] <- readRDS(rds_path)
    i <- i + 1
  } else {
    message("File not found: ", rds_path)
  }
}

# Merge Seurat objects
obj <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = c("BM1","BM2","BM3","BM4","skin1", "skin2"), project = "obj_merged")

obj2 <- merge(seurat_list[[1]], y = c(seurat_list[2],seurat_list[3],seurat_list[4]), add.cell.ids = c("BM1","BM2","BM3","BM4"), project = "obj_merged")

# Subset to select only cells with at least 1000 features
obj <- subset(obj2, nFeature_RNA > 1000)

# Join layers
obj <- JoinLayers(obj)

# Normalize data
obj <- NormalizeData(obj)

# Find variable features
obj <- FindVariableFeatures(obj)

# Scale data
obj <- ScaleData(obj)

# Run azimuth
obj <- RunAzimuth(obj, reference = "bonemarrowref")

# Run PCA
obj <- RunPCA(obj, features = VariableFeatures(object = obj))


# Run UMAP
obj <- RunUMAP(obj, dims = 1:30)

# Find clusters
obj <- FindNeighbors(obj, dims = 1:30)
obj <- FindClusters(obj, resolution = 0.5)

p1 <- DimPlot(obj, reduction = "umap", group.by = "orig.ident") + 
  ggtitle('Unintegrated datasets') 

p2 <- DimPlot(obj, reduction = "umap", group.by = "predicted.celltype.l2", label = TRUE, repel = TRUE) + 
  ggtitle('Unintegrated datasets') + NoLegend()

merged_plots <- p1 + p2
print(merged_plots)
# Save the plot
umap_plot <- "unintegrated_data_UMAP.png"
ggsave(filename = umap_plot, plot = merged_plots, width = 12, height = 6)

feature <- c("CD79A")
features <- c("Bend3", "IGHG2","Aatk")
RidgePlot(obj, features = features, ncol = 2)
VlnPlot(obj, features = features)

FeaturePlot(obj, features = features)

find_markers <- function(obj, output_file){
  # Join layers of obj
  obj <- JoinLayers(obj)
  
  # Open a connection to the file
  sink(output_file)
  
  # Find markers for each cluster compared to all remaining cells
  clusters <- sort(unique(Idents(obj)))
  n_clusters <- length(clusters)
  for (i in seq_along(clusters)) {
    cluster <- clusters[i]
    cat("Cluster:", as.numeric(cluster) - 1, "\n__________________________\n")
    cluster.de.markers <- FindMarkers(obj, ident.1 = cluster, ident.2 = NULL, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
    cluster.de.markers.sorted <- cluster.de.markers[order(cluster.de.markers$avg_log2FC, decreasing = TRUE),]
    print(head(cluster.de.markers.sorted, n = 20))
    cat("__________________________\n")
  }
  
  # Close the connection to the file
  sink()
}

# CCA integration
cca_obj <- IntegrateLayers(
  object = obj, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  dims = 1:10, # Reduce the maximum dimensions to 10
  verbose = FALSE
)

# RPCA integration
rpca_obj <- IntegrateLayers(
  object = obj, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  dims = 1:10, # Reduce the maximum dimensions to 10
  k.anchor = 30, # Increasing k.anchor
  verbose = FALSE
)

# Harmony integration
harmony_obj <- IntegrateLayers(
  object = obj, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

# FastMNN integration 
mnn_obj <- IntegrateLayers(
  object = obj, method = FastMNNIntegration,
  new.reduction = "integrated.mnn",
  verbose = FALSE
)

# scVI integration 
# https://scvi-tools.org/
scvi_params <- list(
  n_latent = 30, 
  n_epochs = 400,
  batch_size = 128,
  learning_rate = 1e-3
)

scvi_obj <- IntegrateLayers(
  object = obj, method = scVIIntegration,
  new.reduction = "integrated.scvi",
  conda_env = "~/miniforge3/envs/scvi", 
  verbose = FALSE,
  integration.params = scvi_params
)

# Find markers per cluster
find_markers(cca_obj, "cca_cluster_markers.txt")
find_markers(rpca_obj, "rpca_cluster_markers.txt")
find_markers(harmony_obj, "harmony_cluster_markers.txt")
find_markers(mnn_obj, "mnn_cluster_markers.txt")
find_markers(scvi_obj, "scvi_cluster_markers.txt")

### Create plots
plotIntegratedData <- function(obj, keyword_string){
  if (keyword_string == "harmony"){
    red_method <- keyword_string
  } else {
    red_method <- paste0("integrated.", keyword_string)
  }
  
  clust_name <- paste0(keyword_string, "_clusters")
  red_name   <- paste0("umap.", keyword_string)
  
  obj <- FindNeighbors(obj, reduction = red_method, dims = 1:30)
  obj <- FindClusters(obj, resolution = 0.5, cluster.name = clust_name)
  
  obj <- RunUMAP(obj, reduction = red_method, dims = 1:30, reduction.name = red_name)
  
  p1 <- DimPlot(obj, reduction = red_name, group.by = c(clust_name),combine = FALSE, label.size = 2)
  p2 <- DimPlot(obj, reduction = red_name, group.by = c("orig.ident"),combine = FALSE, label.size = 2)
  
  combined_plot <- p1[[1]] + p2[[1]]
  
  umapplot_filename <- paste0(keyword_string, "_UMAPplot.png")
  png(umapplot_filename, width = 1600, height = 800)
  print(combined_plot)
  dev.off()
  cat("UMAP plot was saved to file: ", umapplot_filename)
}

plotIntegratedData(cca_obj, "cca")
plotIntegratedData(rpca_obj, "rpca")
plotIntegratedData(harmony_obj, "harmony")
plotIntegratedData(mnn_obj, "mnn")
plotIntegratedData(scvi_obj, "scvi")

features <- c("Cdan1", "IL7R", "SCML4", "TCF7", "KDR")
RidgePlot(cca_obj, features = features, ncol = 2)



### Not working yet
# Function to classify cells and visualize results using Azimuth
classify_cells <- function(obj, keyword_string) {
  obj <- JoinLayers(obj)
  
  # Normalize the data if not already done
  if (!"SCT" %in% Assays(obj)) {
    obj <- SCTransform(obj, assay = "RNA")
  }
  
  # Determine the reduction method
  if (keyword_string == "harmony") {
    red_method <- keyword_string
  } else {
    red_method <- paste0("integrated.", keyword_string)
  }
  
  clust_name <- paste0(keyword_string, "_clusters")
  red_name   <- paste0("umap.", keyword_string)
  
  # Run Azimuth with the correct reference and integrated reduction
  obj <- RunAzimuth(obj, reference = "bonemarrowref", reduction = red_method)
  obj <- FindNeighbors(obj, reduction = red_method, dims = 1:30)
  obj <- FindClusters(obj, resolution = 0.5, graph.name = paste0(red_method, "_snn"))
  #obj <- FindClusters(obj, resolution = 0.5)  # Removed the group.by argument
  obj <- RunUMAP(obj, reduction = red_method, dims = 1:30, reduction.name = red_name)
  
  # Create and combine DimPlots
  p1 <- DimPlot(obj, reduction = red_name, group.by = clust_name, combine = FALSE, label.size = 2)
  p2 <- DimPlot(obj, reduction = red_name, group.by = "predicted.celltype.l2", combine = FALSE, label.size = 2)
  p3 <- DimPlot(obj, reduction = red_name, group.by = "orig.ident", combine = FALSE, label.size = 2)
  
  combined_plot <- p1[[1]] + p2[[1]] + p3[[1]]
  umap_plot <- paste0(keyword_string, "integrated_UMAP_compare.png")
  ggsave(filename = umap_plot, plot = combined_plot, width = 18, height = 6)
  cat("UMAP plot saved to: ", umap_plot, "\n")
}


# Example usage of the function with a CCA-integrated object
classify_cells(cca_obj, "cca")
classify_cells(rpca_obj, "rpca")
classify_cells(harmony_obj, "harmony")
classify_cells(mnn_obj, "mnn")
classify_cells(scvi_obj, "scvi")

