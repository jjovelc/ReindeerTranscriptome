# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(Azimuth)
library(pals)

setwd('/Users/juanjovel/jj/data_analysis/jeffBiernaskie/RangiferTarandus/transcriptomeAssembly/compare_jeff_against_combined/evaluation_hybrid_transcriptome/hybrid_transcriptome/alevin_results')

# 1. Load data
R_BM_A <- readRDS("R_BM_A_output_bone_marrow.RDS")
R_BM_1 <- readRDS("R_BM_1_output_bone_marrow.RDS")
R_BM_S <- readRDS("R_BM_S_output_bone_marrow.RDS")
R_BM_2 <- readRDS("R_BM_2_output_bone_marrow.RDS")

# 2. Merge datasets
R_BM_A$origin <- "Antler"
R_BM_1$origin <- "Antler"
R_BM_S$origin <- "Sternum"
R_BM_2$origin <- "Sternum"

merged_seurat <- merge(R_BM_A, y = c(R_BM_1, R_BM_S, R_BM_2), 
                       add.cell.ids = c("A", "1", "S", "2"), 
                       project = "BoneMarrow")

# 3. Quality control and filtering
VlnPlot(merged_seurat, 
        features = c("nFeature_RNA", "nCount_RNA"), 
        group.by = "origin", 
        ncol = 3)

merged_seurat <- subset(merged_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)

VlnPlot(merged_seurat, 
        features = c("nFeature_RNA", "nCount_RNA"), 
        group.by = "origin", 
        ncol = 3)

# 4. Normalize and scale data
merged_seurat <- NormalizeData(merged_seurat)
merged_seurat <- FindVariableFeatures(merged_seurat)
merged_seurat <- ScaleData(merged_seurat)

# 5. Perform integration
seurat_list <- SplitObject(merged_seurat, split.by = "orig.ident")
features <- SelectIntegrationFeatures(object.list = seurat_list)
anchors <- FindIntegrationAnchors(object.list = seurat_list, anchor.features = features)
integrated_seurat <- IntegrateData(anchorset = anchors)

# 6. Dimension reduction and clustering
# Find variable features on the integrated data
DefaultAssay(integrated_seurat) <- "integrated"
all_genes <- rownames(integrated_seurat)
integrated_seurat <- ScaleData(integrated_seurat, features = all_genes)
integrated_seurat <- RunPCA(integrated_seurat, features = VariableFeatures(object = integrated_seurat))
integrated_seurat <- RunUMAP(integrated_seurat, dims = 1:30)
integrated_seurat <- FindNeighbors(integrated_seurat, dims = 1:30)
integrated_seurat <- FindClusters(integrated_seurat, resolution = 0.5)

# Visualize UMAP
DimPlot(integrated_seurat, reduction = "umap", group.by = "orig.ident")
DimPlot(integrated_seurat, reduction = "umap", label = TRUE)
DimPlot(integrated_seurat, reduction = "umap", group.by = "predicted.celltype.l2", label = T, repel = T)

# 7. Identify cell types
# This step requires manual inspection of marker genes
markers <- FindAllMarkers(integrated_seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top_markers <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

# Visualize top markers
DoHeatmap(integrated_seurat, features = top_markers$gene)

# Assign cell types based on marker genes (example)
new.cluster.ids <- c("CD8 memory", "CD4 memory", "Plasma", "Late Eryth", "CD14 Mono",
                     "HSC", "Memory B", "Early Eryth","CD4 Naive", "CD16 Mono",
                     "LMPP", "T proliferating", "pDC", "unk1", "BaEoMa", "unk2",
                     "unk3")
names(new.cluster.ids) <- levels(integrated_seurat)
integrated_seurat <- RenameIdents(integrated_seurat, new.cluster.ids)


# 8. Compare cell type proportions
cell_proportions <- prop.table(table(Idents(integrated_seurat), integrated_seurat$origin), margin = 2)
cell_proportions_df <- as.data.frame(cell_proportions)
colnames(cell_proportions_df) <- c("CellType", "Origin", "Proportion")

# Visualize cell type proportions

# Create the plot
ggplot(cell_proportions_df, aes(x = Origin, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=as.vector(polychrome(26))) +
  theme_minimal() +
  labs(title = "Cell Type Proportions in Antler vs Sternum Bone Marrow")

# Statistical test for differences in proportions
cell_counts <- table(Idents(integrated_seurat), integrated_seurat$origin)
chi_sq_test <- chisq.test(cell_counts)
print(chi_sq_test)

# Differential expression analysis

# Join layers
merged_seurat <- JoinLayers(merged_seurat)

# Set the tissue type as identity class
Idents(merged_seurat) <- merged_seurat$origin

# Extract the count data
counts <- as.matrix(merged_seurat[["RNA"]]$counts)

# Ensure the data is integer for DESeq2
counts <- round(counts)

# Create a DESeqDataSet object
colData <- data.frame(row.names = colnames(counts), tissue = merged_seurat$origin)
dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ tissue)
# Custom function to calculate geometric means, ignoring zeros
gm_mean <- function(x, na.rm=TRUE) {
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

# Estimate size factors using the custom geometric mean function
geoMeans <- apply(counts(dds), 1, gm_mean)
dds <- estimateSizeFactors(dds, geoMeans=geoMeans)

# Define the tissue types for comparison
group1 <- "Antler"
group2 <- "Sternum"

# List of DE methods available in Seurat
de_methods <- c("wilcox", "bimod", "roc", "t", "LR", "MAST", "DESeq2")

# Initialize a list to store DE results
de_results <- list()

# Perform DE analysis using each method
for (method in de_methods) {
  if (method == "DESeq2") {
    # Run DESeq2 analysis
    dds <- DESeq(dds)
    res <- results(dds, contrast = c("tissue", group1, group2))
    de_results[[method]] <- as.data.frame(res)
    
    # Save DESeq2 results to file in TSV format
    write.table(de_results[[method]], file = paste0("DESeq2_results_", group1, "_vs_", group2, ".tsv"), sep = "\t")
  } else {
    de_results[[method]] <- FindMarkers(merged_seurat, ident.1 = group1, ident.2 = group2, test.use = method)
    
    # Save Seurat FindMarkers results to file in TSV format
    write.table(de_results[[method]], file = paste0("FindMarkers_", method, "_results_", group1, "_vs_", group2, ".tsv"), sep = "\t")
  }
}
