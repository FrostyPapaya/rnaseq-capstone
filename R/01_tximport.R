# R/01_tximport.R
# Import Salmon quantification results into R
library(tximport)
# Load metadata
meta <- read.csv("data/metadata.csv", stringsAsFactors = FALSE)
meta$tissue <- factor(meta$tissue, levels = c("Cerebellum", "SpinalCord"))
rownames(meta) <- meta$sample
# Point to quant.sf files
files <- file.path("data/processed/salmon", meta$SRRid, "quant.sf")
names(files) <- meta$sample
# Check all files exist
stopifnot(all(file.exists(files)))
cat("All quant.sf files found\n")
# Load tx2gene
tx2gene <- read.table("data/tx2gene.tsv",
                      header = FALSE,
                      col.names = c("tx_id", "gene_id"))
# Strip version numbers from IDs (e.g. ENSMUST00000177564.1 -> ENSMUST00000177564)
tx2gene$tx_id   <- sub("\\..*", "", tx2gene$tx_id)
tx2gene$gene_id <- sub("\\..*", "", tx2gene$gene_id)
# Import
txi <- tximport(files,
                type       = "salmon",
                tx2gene    = tx2gene,
                ignoreTxVersion = TRUE)
# Sanity checks
cat("Count matrix dimensions:", dim(txi$counts), "\n")
cat("Sample names:", colnames(txi$counts), "\n")
cat("First few genes:\n")
print(head(txi$counts))
# Save
saveRDS(txi, "data/processed/txi.rds")
saveRDS(meta, "data/processed/meta.rds")
cat("Saved txi.rds and meta.rds\n")
