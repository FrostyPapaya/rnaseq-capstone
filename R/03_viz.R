# R/03_viz.R
# Visualisations: PCA, Volcano, Heatmap
library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)
library(pheatmap)
library(RColorBrewer)
# Load
load("data/processed/deseq2.RData")
meta <- readRDS("data/processed/meta.rds")
dir.create("results/figures", showWarnings = FALSE)
# ── 1. PCA ──────────────────────────────────────────────────────────────────
pca_data <- plotPCA(rld, intgroup = "tissue", returnData = TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"))
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2,
                                  color = tissue, label = name)) +
  geom_point(size = 4) +
  geom_text(vjust = -0.8, size = 3) +
  xlab(paste0("PC1: ", pct_var[1], "% variance")) +
  ylab(paste0("PC2: ", pct_var[2], "% variance")) +
  scale_color_brewer(palette = "Set1") +
  theme_bw(base_size = 14) +
  ggtitle("PCA: Cerebellum vs Spinal Cord")
ggsave("results/figures/pca.pdf", pca_plot, width = 7, height = 5)
ggsave("results/figures/pca.png", pca_plot, width = 7, height = 5, dpi = 300)
cat("PCA saved\n")
# ── 2. Volcano plot ──────────────────────────────────────────────────────────
res_df <- read.csv("results/tables/de_results.csv", row.names = 1)
# Label top 10 significant genes by padj
top_genes <- rownames(res_df[order(res_df$padj), ])[1:10]
pdf("results/figures/volcano.pdf", width = 8, height = 6)
EnhancedVolcano(res_df,
                lab            = rownames(res_df),
                x              = "log2FoldChange",
                y              = "padj",
                selectLab      = top_genes,
                pCutoff        = 0.05,
                FCcutoff       = 1.0,
                pointSize      = 2,
                labSize        = 3,
                title          = "Spinal Cord vs Cerebellum",
                subtitle       = "Influenza infection mouse study (GSE96870)",
                col            = c("grey30", "forestgreen", "royalblue", "red2"),
                legendPosition = "right")
dev.off()
cat("Volcano saved\n")
# ── 3. Heatmap ───────────────────────────────────────────────────────────────
# Top N DE genes by padj (only genes that exist in rld)
sig_genes <- rownames(res)[!is.na(res$padj) & res$padj < 0.05]
sig_genes <- sig_genes[sig_genes %in% rownames(rld)]  # Only keep genes in rld
topN      <- sig_genes[1:min(20, length(sig_genes))]

cat("Significant genes in heatmap:", length(topN), "\n")

if (length(topN) > 0) {
  mat <- assay(rld)[topN, ]
mat <- mat - rowMeans(mat)  # centre rows
# Annotation bar
ann_col <- data.frame(
  Tissue    = meta$tissue,
  Infection = meta$infection,
  row.names = meta$sample
)
ann_colors <- list(
  Tissue    = c(Cerebellum = "#E41A1C", SpinalCord = "#377EB8"),
  Infection = c(Influenza  = "#FF7F00", NonInfected = "#4DAF4A")
)
pdf("results/figures/heatmap.pdf", width = 8, height = 10)
pheatmap(mat,
         annotation_col  = ann_col,
         annotation_colors = ann_colors,
         color           = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
         cluster_cols    = TRUE,
         cluster_rows    = TRUE,
         show_rownames   = TRUE,
         fontsize_row    = 7,
         fontsize_col    = 10,
         main            = "Top 20 DE Genes (rlog scaled)")
dev.off()
# PNG version
png("results/figures/heatmap.png", width = 800, height = 1000, res = 120)
pheatmap(mat,
         annotation_col  = ann_col,
         annotation_colors = ann_colors,
         color           = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
         cluster_cols    = TRUE,
         cluster_rows    = TRUE,
         show_rownames   = TRUE,
         fontsize_row    = 7,
         fontsize_col    = 10,
         main            = "Top 20 DE Genes (rlog scaled)")
dev.off()
} else {
  cat("No significant genes to plot in heatmap\n")
}
cat("Heatmap saved\n")
cat("\nAll figures saved to results/figures/\n")
