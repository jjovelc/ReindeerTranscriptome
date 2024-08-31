library(cowplot)
library(Seurat)
library(tximport)
library(SeuratData)
library(SeuratWrappers)
library(Azimuth)
library(ggplot2)
library(patchwork)
library(reticulate)

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



# Function to generate descriptive statistics plots for a single Seurat object
generate_stats_plots <- function(seurat_obj, obj_name) {
  # Number of cells
  cell_count <- ncol(seurat_obj)
  
  # Get the count matrix
  count_matrix <- GetAssayData(seurat_obj, slot = "counts")
  
  # Number of features (genes) detected in each cell
  feature_counts <- colSums(count_matrix > 0)
  
  # Total UMI counts per cell
  umi_counts <- colSums(count_matrix)
  
  # Percentage of mitochondrial genes
  mt_genes <- grep("^MT-", rownames(count_matrix), value = TRUE)
  percent_mt <- colSums(count_matrix[mt_genes, , drop = FALSE]) / colSums(count_matrix) * 100
  
  # Create plots
  p1 <- ggplot(data.frame(features = feature_counts), aes(x = features)) +
    geom_histogram(bins = 50, fill = "blue", alpha = 0.7) +
    theme_minimal() +
    labs(title = paste("Number of Genes per Cell -", obj_name),
         x = "Number of Genes", y = "Cell Count")
  
  p2 <- ggplot(data.frame(umi = umi_counts), aes(x = umi)) +
    geom_histogram(bins = 50, fill = "red", alpha = 0.7) +
    theme_minimal() +
    labs(title = paste("UMI Counts per Cell -", obj_name),
         x = "UMI Count", y = "Cell Count")
  
  p3 <- ggplot(data.frame(percent_mt = percent_mt), aes(x = percent_mt)) +
    geom_histogram(bins = 50, fill = "green", alpha = 0.7) +
    theme_minimal() +
    labs(title = paste("Percentage of Mitochondrial Genes -", obj_name),
         x = "Percentage of Mitochondrial Genes", y = "Cell Count")
  
  p4 <- ggplot(data.frame(features = feature_counts, umi = umi_counts), aes(x = features, y = umi)) +
    geom_point(alpha = 0.5) +
    theme_minimal() +
    labs(title = paste("Features vs UMI Counts -", obj_name),
         x = "Number of Genes", y = "UMI Count")
  
  # Combine plots
  combined_plot <- plot_grid(p1, p2, p3, p4, ncol = 2)
  
  # Add a title to the combined plot
  title <- ggdraw() + 
    draw_label(paste("Descriptive Statistics for", obj_name), fontface = "bold", size = 16)
  
  final_plot <- plot_grid(title, combined_plot, ncol = 1, rel_heights = c(0.1, 1))
  
  # Print summary statistics
  cat("\nSummary Statistics for", obj_name, ":\n")
  cat("Number of Cells:", cell_count, "\n")
  cat("Median Genes per Cell:", median(feature_counts), "\n")
  cat("Median UMI per Cell:", median(umi_counts), "\n")
  cat("Median Percentage of Mitochondrial Genes:", median(percent_mt), "\n\n")
  
  return(final_plot)
}

# Generate plots for each Seurat object in the list
plot_list <- list()
for (i in seq_along(seurat_list)) {
  obj_name <- names(seurat_list)[i]
  plot_list[[obj_name]] <- generate_stats_plots(seurat_list[[i]], obj_name)
}

# Save all plots to a PDF file
pdf("descriptive_statistics_plots.pdf", width = 12, height = 10)
for (plot in plot_list) {
  print(plot)
}
dev.off()

# Generate a summary table
summary_table <- data.frame(
  Library = names(seurat_list),
  Cells = sapply(seurat_list, ncol),
  Median_Genes = sapply(seurat_list, function(x) {
    count_matrix <- GetAssayData(x, slot = "counts")
    median(colSums(count_matrix > 0))
  }),
  Median_UMI = sapply(seurat_list, function(x) {
    count_matrix <- GetAssayData(x, slot = "counts")
    median(colSums(count_matrix))
  }),
  Median_Percent_MT = sapply(seurat_list, function(x) {
    count_matrix <- GetAssayData(x, slot = "counts")
    mt_genes <- grep("^MT-", rownames(count_matrix), value = TRUE)
    percent_mt <- colSums(count_matrix[mt_genes, , drop = FALSE]) / colSums(count_matrix) * 100
    median(percent_mt)
  })
)

# Print the summary table
print(summary_table)

# Save the summary table to a CSV file
write.csv(summary_table, "library_summary_statistics.csv", row.names = FALSE)

########################
library(ggplot2)
library(tidyr)
library(dplyr)

# Assuming your data is stored in a dataframe called 'summary_table'
# If not, create it from the data you provided:
summary_table <- data.frame(
  Library = c("BM_Antler1", "BM_Antler2", "BM_Stern1", "BM_Stern2", "Skin_N", "Skin_P"),
  Cells = c(1324, 2982, 1266, 2950, 392, 392),
  Median_Genes = c(1172.0, 875.5, 1441.0, 936.5, 3673.0, 3673.0)
)

# Reshape the data for plotting
plot_data <- summary_table %>%
  select(Library, Cells, Median_Genes) %>%
  pivot_longer(cols = c(Cells, Median_Genes), names_to = "Metric", values_to = "Value")

# Create the paired bar plot
ggplot(plot_data, aes(x = Library, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8) +
  scale_fill_manual(values = c("Cells" = "#1f77b4", "Median_Genes" = "#ff7f0e")) +
  labs(title = "Number of Cells and Median Genes per Library",
       x = "Library",
       y = "Count",
       fill = "Metric") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top") +
  scale_y_continuous(labels = scales::comma) +
  geom_text(aes(label = scales::comma(Value)), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5, 
            size = 3)

# Save the plot
ggsave("cells_and_median_genes_plot.png", width = 12, height = 8, dpi = 300)


