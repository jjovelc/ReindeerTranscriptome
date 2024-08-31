library(Seurat)
library(tximport)

#################################################################
# NOTE: This script requires that ONLY the alevin outputs
#       are present in the pwd, each of which should have
#       an 'alevin' directory with the following files:
#       alevin.log, predictions.txt, quants_mat_cols.txt	
#       raw_cb_frequency.txt, cell_umi_graphs.gz, quants_mat.gz
#       quants_mat_rows.txt, whitelist.txt, featureDump.txt
#       quants_mat.mtx.gz, quants_tier_mat.gz
#################################################################

setwd('/Users/juanjovel/jj/data_analysis/jeffBiernaskie/RangiferTarandus/transcriptomeAssembly/compare_jeff_against_combined/evaluation_hybrid_transcriptome/hybrid_transcriptome/alevin_results')
db <- 'Hybrid'
dirs <- list.dirs(getwd(), full.names = T, recursive = F)
n_comp = 20
for (dir in dirs){
  prefix <- basename(dir)
  cat("Processing sample ", prefix, "\n")
  setwd(dir)

  # Path to the Alevin output directory
  alevin_dir <- "alevin"

  # File paths
  # path to the output directory of Alevin run
  files <- file.path("alevin/")
  file.exists(files)

  # Reading in the alevin quants quants
  txi <- tximport(files, type="alevin")

  obj <- CreateSeuratObject(counts = txi$counts , min.cells = 3, 
                          min.features = 200, project = prefix)

  obj <- NormalizeData(obj, normalization.method = "LogNormalize", 
                     scale.factor = 10000)
  obj <- FindVariableFeatures(obj, selection.method = "vst", 
                            nfeatures = 2000)

  all.genes <- rownames(obj)

  obj <- ScaleData(obj, features = all.genes)

  obj <- RunAzimuth(obj, reference = "bonemarrowref")
  obj <- RunPCA(obj, features = VariableFeatures(object = obj))

  obj <- FindNeighbors(obj, dims = 1:n_comp)
  obj <- FindClusters(obj, resolution = 0.5)

  obj <- RunUMAP(obj, dims = 1:n_comp)

  umap_plot <- paste0(prefix, "_umap_wAzimuth", "_", n_comp, "comp.png")
  azimuth_labels <- obj[["predicted.celltype.l2"]]
  p <- DimPlot(obj, reduction = "umap", group.by = "predicted.celltype.l2", label = TRUE, repel = TRUE) + 
    ggtitle(paste0(db, ' ', prefix)) + NoLegend()
  # Save the plot
  ggsave(filename = umap_plot, plot = p, width = 8, height = 6)
    
  rds_object <- paste0(prefix, ".RDS")
  saveRDS(obj, rds_object)
}
