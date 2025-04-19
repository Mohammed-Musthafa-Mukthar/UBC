# === Set custom library path ===
.libPaths("/gpfs/fs7/aafc/phenocart/PhenomicsProjects/UFPSGPSCProject/4_Assets/R/custom_lib")

# trajectory_analysis_slingshot.R
# Perform single-cell RNA-seq trajectory analysis using slingshot

# 1. Install and load necessary packages
#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#BiocManager::install(c("SingleCellExperiment", "slingshot", "scater", "scran", "RColorBrewer"), ask = FALSE)
#install.packages("ggplot2", repos = "https://cloud.r-project.org")

library(SingleCellExperiment)
library(slingshot)
library(scater)
library(scran)
library(RColorBrewer)
library(ggplot2)
library(ragg)

# 2. Load expression count matrix
counts_path <- "/gpfs/fs7/aafc/phenocart/PhenomicsProjects/RNASeq/musthafa_project_trajectory/cleaned_counts.csv"
counts_df   <- read.csv(counts_path, row.names = 1, check.names = FALSE)
expr_mat    <- as.matrix(counts_df)

# 3. Create SingleCellExperiment and normalize
sce <- SingleCellExperiment(assays = list(counts = expr_mat))
sce <- logNormCounts(sce)

# 4. Dimensionality reduction: PCA then UMAP
sce <- runPCA(sce, ncomponents = 50)
sce <- runUMAP(sce, dimred = "PCA")

# 5. Clustering
cluster_labels <- clusterCells(sce, use.dimred = "PCA")
colData(sce)$cluster <- factor(cluster_labels)

# 6. Trajectory inference with Slingshot
sce <- slingshot(sce, clusterLabels = 'cluster', reducedDim = 'PCA')

# 7. Extract and save pseudotime
pseudotime_values <- slingPseudotime(sce)[,1]
colData(sce)$pseudotime <- pseudotime_values
write.csv(
  data.frame(cell = colnames(expr_mat), pseudotime = pseudotime_values),
  file = "pseudotime_values_slingshot.csv", row.names = FALSE
)

# 8. PCA trajectory plot (PNG via ragg)
ragg::agg_png("trajectory_pca.png", width = 800, height = 600, units = "px")
colors <- brewer.pal(n = length(unique(cluster_labels)), name = "Set3")
plot(
  reducedDim(sce, "PCA")[,1:2],
  col   = colors[as.integer(cluster_labels)],
  pch   = 16,
  asp   = 1,
  xlab  = "PC1",
  ylab  = "PC2",
  main  = "Slingshot Trajectory on PCA"
)
lines(SlingshotDataSet(sce), lwd = 2)
dev.off()

# 9. UMAP colored by pseudotime (PNG via ragg)
ragg::agg_png("trajectory_umap_pseudotime.png", width = 800, height = 600, units = "px")
print(
  plotUMAP(sce, colour_by = "pseudotime") +
    ggtitle("UMAP Colored by Pseudotime")
)
dev.off()

# 10. Pseudotime distribution histogram (PNG via ragg)
ragg::agg_png("pseudotime_histogram.png", width = 800, height = 600, units = "px")
hist(
  pseudotime_values,
  breaks = 30,
  main   = "Pseudotime Distribution",
  xlab   = "Pseudotime",
  col    = "lightblue"
)
dev.off()


