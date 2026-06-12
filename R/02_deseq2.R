# R/02_deseq2.R
# Differential expression: Cerebellum vs SpinalCord
library(DESeq2)
library(apeglm)
# Load objects from Week 6
txi  <- readRDS("data/processed/txi.rds")
meta <- readRDS("data/processed/meta.rds")
# Create DESeq2 object
# Design: compare by tissue (Cerebellum is reference level)
dds <- DESeqDataSetFromTximport(txi,
                                colData = meta,
                                design  = ~ tissue)
# Pre-filter: remove genes with very low counts
# Keep genes with >= 10 counts in at least half the samples
keep <- rowSums(counts(dds) >= 10) >= (ncol(dds) / 2)
dds  <- dds[keep, ]
cat("Genes after filtering:", nrow(dds), "\n")
# Run DESeq2 (fits NB model, estimates dispersions, tests)
dds <- DESeq(dds)
# Extract results: SpinalCord vs Cerebellum
res <- results(dds,
               contrast = c("tissue", "SpinalCord", "Cerebellum"),
               alpha    = 0.05)
# Shrink LFC estimates (better for ranking and visualization)
res <- lfcShrink(dds,
                 coef = "tissue_SpinalCord_vs_Cerebellum",
                 type = "apeglm")
# Summary
cat("\nDESeq2 results summary:\n")
summary(res)
# rlog transform for visualization (blind=FALSE uses the design)
rld <- rlog(dds, blind = FALSE)
# Save everything
save(dds, res, rld, file = "data/processed/deseq2.RData")
# Export results table
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df <- res_df[order(res_df$padj, na.last = TRUE), ]
write.csv(res_df, "results/tables/de_results.csv", row.names = FALSE)
# Export normalised counts
norm_counts <- counts(dds, normalized = TRUE)
write.csv(norm_counts, "results/tables/normalised_counts.csv")
cat("\nSignificant genes (padj < 0.05):", sum(res$padj < 0.05, na.rm = TRUE), "\n")
cat("Saved: deseq2.RData, de_results.csv, normalised_counts.csv\n")
