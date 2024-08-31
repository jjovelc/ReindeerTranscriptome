library(Seurat)
library(tximport)
library(SeuratData)
library(SeuratWrappers)
library(Azimuth)
library(ggplot2)
library(patchwork)
library(reticulate)

setwd('/Users/juanjovel/jj/data_analysis/jeffBiernaskie/RangiferTarandus/transcriptomeAssembly/compare_jeff_against_combined/evaluation_hybrid_transcriptome/alevin_results/')
main_dir <- getwd()

# List all items in the specified path that match the pattern "output"
all_items <- list.files(main_dir, pattern = "output", full.names = TRUE)

# Filter to include only directories
directories <- all_items[file.info(all_items)$isdir]

# Initialize counter
i <- 1

# Loop through each directory and load the RDS file
for (sub_dir in directories) {
  subdir_name <- basename(sub_dir)
  rds_name <- paste0(subdir_name, ".RDS")
  rds_path <- file.path(sub_dir, rds_name)
  
  if (file.exists(rds_path)) {
    obj_name <- paste0("obj", i)
    assign(obj_name, readRDS(rds_path))
    # Load the Seurat object
    obj <- get(obj_name)
    
    # Normalize the query dataset using SCTransform
    #obj <- SCTransform(obj, assay = "RNA")
    
    # Run Azimuth with the correct reference
    obj <- tryCatch({
      RunAzimuth(obj, reference = "bonemarrowref")
    }, error = function(e) {
      message("Error running RunAzimuth: ", e$message)
      NULL
    })
    
    if (!is.null(obj)) {
      # Save the modified object back
      assign(obj_name, obj)
    }
    i <- i + 1
  } else {
    message("File not found: ", rds_path)
  }
}

annotate_ind_datasets <- function(obj, prefix){
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj)
  obj <- FindNeighbors(obj, dims = 1:30, reduction = "pca")
  obj <- FindClusters(obj, resolution = 2, cluster.name = "unintegrated_clusters")
  obj <- RunUMAP(obj, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
  plotName <- paste0(prefix, "_plusAzimuth.png")
  png(plotName)
  p <- DimPlot(obj, reduction = "umap.unintegrated",
          group.by = "predicted.celltype.l2")
  print(p)
  dev.off()
}

annotate_ind_datasets(obj1, "R_BM_1")
annotate_ind_datasets(obj2, "R_BM_2")
annotate_ind_datasets(obj3, "R_BM_A")
annotate_ind_datasets(obj4, "R_BM_S")
annotate_ind_datasets(obj5, "R_skin-N")




obj <- merge(obj1, y = c(obj2, obj3, obj4, obj5), add.cell.ids = c("BM1","BM2","BM3","BM4","skin1"), 
             project = "obj_merged")
# Subset to select only cells with at least 1000 features
obj <- subset(obj, nFeature_RNA > 1000)
obj <- SCTransform(obj, assay = "RNA", verbose = FALSE)
obj <- RunAzimuth(obj, reference = "bonemarrowref")
obj <- NormalizeData(obj)

obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)





# Visualize a standard analysis without integration
# While a UMAP analysis is just a visualization of this, clustering 
# this dataset would return predominantly batch-specific clusters. 
# Especially if previous cell-type annotations were not available, 
# this would make downstream analysis extremely challenging.
obj <- FindNeighbors(obj, dims = 1:30, reduction = "pca")
obj <- FindClusters(obj, resolution = 2, cluster.name = "unintegrated_clusters")
obj <- RunUMAP(obj, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

# Visualize by batch and cell type annotation
# cell type annotation were previously added by Azimuth
png("UMAPplot_by_experiment.png")
DimPlot(obj, reduction = "umap.unintegrated",
        group.by = "orig.ident")
dev.off()


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
  #png(umapplot_filename, width = 1600, height = 800)
  print(combined_plot)
  #dev.off()
  cat("UMAP plot was saved to file: ", umapplot_filename)
}

plotIntegratedData(cca_obj, "cca")
plotIntegratedData(rpca_obj, "rpca")
plotIntegratedData(harmony_obj, "harmony")
plotIntegratedData(mnn_obj, "mnn")
plotIntegratedData(scvi_obj, "scvi")





# Function to classify cells and visualize results using Azimuth
classify_cells <- function(obj, keyword_string) {
  # Normalize the data if not already done
  if (!"SCT" %in% Assays(obj)) {
    obj <- SCTransform(obj, assay = "RNA")
  }
  
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
  obj <- FindClusters(obj, resolution = 0.5, group.by = clust_name)
  obj <- RunUMAP(obj, reduction = red_method, dims = 1:30, reduction.name = red_name)
  
  # Create and combine DimPlots
  p1 <- DimPlot(obj, reduction = red_name, group.by = clust_name, combine = FALSE, label.size = 2)
  p2 <- DimPlot(obj, reduction = red_name, group.by = "orig.ident", combine = FALSE, label.size = 2)
  
  combined_plot <- p1[[1]] + p2[[1]]
  print(combined_plot)
}

# Example usage of the function with a CCA-integrated object
classify_cells(cca_obj, "cca")

#################
# Run Azimuth with the correct reference and integrated reduction
obj <- NormalizeData(cca_obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
obj <- SCTransform(obj, assay = "RNA")
obj <- RunAzimuth(obj, reference = "bonemarrowref")
obj <- FindNeighbors(cca_obj, reduction = "integrated.cca", dims = 1:30)
obj <- FindClusters(obj, resolution = 0.5, group.by = "cca_clusters")
obj <- RunUMAP(obj, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")
p1 <- DimPlot(obj, reduction = "integrated.cca", group.by = "cca_clusters", combine = FALSE, label.size = 2)
p2 <- DimPlot(obj, reduction = "umap.cca", group.by = "orig.ident", combine = FALSE, label.size = 2)

combined_plot <- p1[[1]] + p2[[1]]
print(combined_plot)



saveRDS(obj, "BM_and_skin_wIntegration.rds")

