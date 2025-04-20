# 1. Load required packages for differential expression and enrichment analysis
library(DESeq2)       
library(ggplot2)       
library(ggrepel)      
library(pheatmap)      
library(clusterProfiler)  
library(org.Mm.eg.db)    
library(AnnotationDbi)   

# 2. Set working directory to data directory 
setwd("C:/Users/MuktharM/OneDrive - AGR-AGR/Documents/DEG_project/")

# 3. Load count data and sample metadata
count_data <- read.csv("count_data.csv", row.names = 1, stringsAsFactors = FALSE)
col_data   <- read.csv("col_data.csv", row.names = 1, stringsAsFactors = FALSE)

# count_data columns and col_data rows
col_data$group <- factor(col_data$group, levels = c("group1", "group2"))
# count matrix columns to match col_data row order 
count_data <- count_data[, rownames(col_data)]

# 4.  DESeq2 dataset and run differential expression analysis
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = col_data, design = ~ group)
dds$group <- relevel(dds$group, ref = "group1")  
dds <- DESeq(dds)  
# 5.  results for the group comparison (group1 vs group2)
res <- results(dds, contrast = c("group", "group2", "group1"))
res_df <- as.data.frame(res)

# results by adjusted p-value for easier inspection
res_df <- res_df[order(res_df$padj, na.last=TRUE), ]
# gene identifier column
res_df$EnsemblID <- rownames(res_df)
res_df <- res_df[, c(ncol(res_df), 1:(ncol(res_df)-1))]

# 6. DEGs by significance criteria
res_sig <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, ]

if(!"EnsemblID" %in% colnames(res_sig) && nrow(res_sig)>0){
  res_sig$EnsemblID <- rownames(res_sig)
}
res_sig <- res_sig[, c("EnsemblID", setdiff(colnames(res_sig), "EnsemblID"))]

# 7. Save results to CSV files
write.csv(res_df, "DESeq2_results_unfiltered.csv", row.names = FALSE)
write.csv(res_sig, "DESeq2_results_filtered.csv", row.names = FALSE)


# 8. Convert Ensembl IDs to Entrez IDs and perform GO enrichment 
if(nrow(res_sig) > 0) {
  
  ens_ids <- sub("\\..*$", "", res_sig$EnsemblID)
  ens_ids <- unique(ens_ids)
  
  entrez_ids <- mapIds(org.Mm.eg.db, keys = ens_ids, column = "ENTREZID",
                       keytype = "ENSEMBL", multiVals = "first")
 
  entrez_ids <- unique(entrez_ids[!is.na(entrez_ids)])
  if(length(entrez_ids) > 0) {
    # GO enrichment for Biological Process ontology
    ego <- enrichGO(gene = entrez_ids,
                    OrgDb = org.Mm.eg.db,
                    keyType = "ENTREZID",
                    ont = "BP",              
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05)
    ego_res <- as.data.frame(ego)
    # Save GO enrichment results 
    if(nrow(ego_res) > 0) {
      write.csv(ego_res, "GO_BP_enrichment_results.csv", row.names = FALSE)
    }
  }
}



# 9a. Volcano Plot of differential expression

volcano_df <- res_df
volcano_df$padj[is.na(volcano_df$padj)] <- 1  

volcano_df$Status <- "NotSig"
volcano_df$Status[volcano_df$padj < 0.05 & volcano_df$log2FoldChange > 1] <- "Up"
volcano_df$Status[volcano_df$padj < 0.05 & volcano_df$log2FoldChange < -1] <- "Down"
volcano_df$Status <- factor(volcano_df$Status, levels = c("Up","Down","NotSig"))

# top 5 up and top 5 down by log2FC (among significant DEGs)
label_genes <- c()
if(sum(volcano_df$Status=="Up") > 0) {
  up_genes <- head(volcano_df[volcano_df$Status=="Up", ][order(volcano_df[volcano_df$Status=="Up", "log2FoldChange"], decreasing=TRUE), ], 5)$EnsemblID
  label_genes <- c(label_genes, as.character(up_genes))
}
if(sum(volcano_df$Status=="Down") > 0) {
  down_genes <- head(volcano_df[volcano_df$Status=="Down", ][order(volcano_df[volcano_df$Status=="Down", "log2FoldChange"], decreasing=FALSE), ], 5)$EnsemblID
  label_genes <- c(label_genes, as.character(down_genes))
}
volcano_df$Label <- ifelse(volcano_df$EnsemblID %in% label_genes, volcano_df$EnsemblID, NA)

# Create volcano plot with ggplot2
volcano_plot <- ggplot(volcano_df, aes(x = log2FoldChange, y = -log10(padj), color = Status, label = Label)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Up"="red", "Down"="blue", "NotSig"="grey")) +
  geom_text_repel(min.segment.length = 0, box.padding = 0.3, max.overlaps = Inf) +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 adjusted p-value") +
  theme_bw(base_size = 14)
# Save volcano plot to file
ggsave("volcano_plot.png", plot = volcano_plot, width = 6, height = 6)

# 9b. PCA Plot of samples 
vsd <- vst(dds, blind = FALSE)  
# Prepare PCA data with sample labels
pcaData <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

if(!"name" %in% colnames(pcaData)) {
  pcaData$name <- rownames(pcaData)  
}
# Plot PCA with samples labeled
pca_plot <- ggplot(pcaData, aes(x = PC1, y = PC2, color = group, label = name)) +
  geom_point(size = 3) +
  geom_text_repel(min.segment.length = 0) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +  
  labs(title = "PCA Plot of Samples") +
  theme_bw(base_size = 14)
ggsave("PCA_plot.png", plot = pca_plot, width = 6, height = 5.5)

# 9c. Heatmap of sample-to-sample distances

sample_dists <- dist(t(assay(vsd)))
sample_dist_matrix <- as.matrix(sample_dists)
rownames(sample_dist_matrix) <- colnames(sample_dist_matrix) <- colnames(vsd)

dist_colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

pheatmap(sample_dist_matrix,
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists,
         col = dist_colors, border_color = NA,
         main = "Sample Distance Heatmap",
         filename = "sample_distance_heatmap.png", width = 5.5, height = 5.5)

