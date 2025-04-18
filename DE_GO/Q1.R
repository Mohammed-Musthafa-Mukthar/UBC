.libPaths("/gpfs/fs7/aafc/phenocart/PhenomicsProjects/UFPSGPSCProject/4_Assets/R/custom_lib")

library(edgeR)
library(clusterProfiler)
library(org.Mm.eg.db)

# Load and clean counts matrix
counts <- read.delim(
  "/gpfs/fs7/aafc/phenocart/PhenomicsProjects/RNASeq/DE_GO/Data/counts.txt",
  comment.char = "#",
  row.names = 1,
  check.names = FALSE
)

# Automatically select sample columns (those with SAMN IDs)
count_data <- counts[, grep("SAMN", colnames(counts))]

# Define group for the 4 samples
group <- factor(c("KO", "KO", "WT", "WT"))

# Validate
stopifnot(length(group) == ncol(count_data))

# DGE processing
dge <- DGEList(counts = count_data, group = group)
dge <- dge[filterByExpr(dge), , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)

# Differential expression
design <- model.matrix(~ group)
dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)
qlf <- glmQLFTest(fit, coef = 2)
deg <- topTags(qlf, n = Inf)$table

# Filter significant DEGs
sig_genes <- deg[deg$FDR < 0.05 & abs(deg$logFC) > 1, ]
ensembl_ids <- rownames(sig_genes)

# Strip Ensembl version (e.g., ENSMUSG00000102693.2 → ENSMUSG00000102693)
ensembl_ids_clean <- sub("\\..*", "", ensembl_ids)

# Map Ensembl to Entrez
gene_entrez <- bitr(
  ensembl_ids_clean,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

# GO enrichment
ego <- enrichGO(
  gene = gene_entrez$ENTREZID,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE
)

# Plot top GO terms
barplot(ego, showCategory = 20)
dotplot(ego, showCategory = 20)
