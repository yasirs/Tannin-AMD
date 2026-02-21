# RefSeq Mapping via biomaRt (Alternative Approach)
# Since org.Hs.eg.db doesn't recognize the RefSeq IDs, try biomaRt

library(dplyr)

# Check if biomaRt is available
if (!requireNamespace("biomaRt", quietly = TRUE)) {
  message("⚠️  biomaRt is not installed.")
  message("To install: BiocManager::install('biomaRt')")
  message("\nCannot proceed with RefSeq mapping via biomaRt.")
  quit(status = 1)
}

library(biomaRt)

message("=== RefSeq Mapping via biomaRt ===\n")

# Load annotation
probe_file <- "results/cohort-GSE29801/dry_amd_de/probe_annotation_improved.csv"
all_probes <- read.csv(probe_file, stringsAsFactors = FALSE)

n_before <- sum(!is.na(all_probes$gene_symbol) & all_probes$gene_symbol != "", na.rm = TRUE)
message(sprintf("Starting with %d mapped probes (%.1f%%)\n",
                n_before, 100 * n_before / nrow(all_probes)))

# Connect to Ensembl
message("Connecting to Ensembl biomaRt...")
tryCatch({
  ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  message("  ✅ Connected\n")
}, error = function(e) {
  message(sprintf("  ⚠️  Connection failed: %s", e$message))
  message("  Trying alternative host...")
  ensembl <<- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl",
                         host = "useast.ensembl.org")
  message("  ✅ Connected to alternative host\n")
})

# Get unmapped probes with RefSeq IDs
unmapped_idx <- which(is.na(all_probes$gene_symbol) | all_probes$gene_symbol == "")
refseq_ids <- all_probes$REFSEQ[unmapped_idx]

# Clean RefSeq IDs
refseq_clean <- gsub("\"", "", refseq_ids)
refseq_clean <- gsub("\\..*$", "", refseq_clean)
refseq_clean <- trimws(refseq_clean)
refseq_clean <- refseq_clean[!is.na(refseq_clean) & refseq_clean != ""]

unique_refseq <- unique(refseq_clean)
message(sprintf("Mapping %d unique RefSeq IDs via biomaRt...\n", length(unique_refseq)))

# Query biomaRt
refseq_mapping <- getBM(
  attributes = c('refseq_mrna', 'hgnc_symbol'),
  filters = 'refseq_mrna',
  values = unique_refseq,
  mart = ensembl
)

# Remove empty symbols
refseq_mapping <- refseq_mapping[refseq_mapping$hgnc_symbol != "", ]
message(sprintf("  Retrieved %d mappings from biomaRt\n", nrow(refseq_mapping)))

# Map back to probes
for (i in unmapped_idx) {
  refseq_id <- all_probes$REFSEQ[i]
  if (!is.na(refseq_id) && refseq_id != "") {
    # Clean ID
    refseq_clean_id <- gsub("\"", "", refseq_id)
    refseq_clean_id <- gsub("\\..*$", "", refseq_clean_id)
    refseq_clean_id <- trimws(refseq_clean_id)
    
    # Find mapping
    match_idx <- which(refseq_mapping$refseq_mrna == refseq_clean_id)
    if (length(match_idx) > 0) {
      all_probes$gene_symbol[i] <- refseq_mapping$hgnc_symbol[match_idx[1]]
    }
  }
}

n_after <- sum(!is.na(all_probes$gene_symbol) & all_probes$gene_symbol != "", na.rm = TRUE)
improvement <- n_after - n_before

message(sprintf("\n=== RefSeq Mapping via biomaRt Complete ==="))
message(sprintf("Before: %d mapped (%.1f%%)", n_before, 100 * n_before / nrow(all_probes)))
message(sprintf("After:  %d mapped (%.1f%%)", n_after, 100 * n_after / nrow(all_probes)))
message(sprintf("Improvement: +%d probes (%.1f%% → %.1f%%)\n",
                improvement,
                100 * n_before / nrow(all_probes),
                100 * n_after / nrow(all_probes)))

if (improvement > 0) {
  # Save
  write.csv(all_probes, probe_file, row.names = FALSE)
  message(sprintf("✅ Saved: %s", probe_file))
  
  complete_file <- "results/cohort-GSE29801/dry_amd_de/probe_annotation_complete.csv"
  write.csv(all_probes, complete_file, row.names = FALSE)
  message(sprintf("✅ Saved: %s", complete_file))
} else {
  message("⚠️  No additional mappings found via biomaRt")
}
