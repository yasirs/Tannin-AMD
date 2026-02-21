# RefSeq ID Mapping Fix for GSE29801
# Resolves the RefSeq mapping issue by properly formatting IDs

library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)

message("=== RefSeq ID Mapping (Fixed Approach) ===\n")

# Load current annotation  
probe_file <- "results/cohort-GSE29801/dry_amd_de/probe_annotation_improved.csv"
all_probes <- read.csv(probe_file, stringsAsFactors = FALSE)

n_before <- sum(!is.na(all_probes$gene_symbol) & all_probes$gene_symbol != "", na.rm = TRUE)
message(sprintf("Starting with %d mapped probes (%.1f%%)\n",
                n_before, 100 * n_before / nrow(all_probes)))

# Get unmapped probes with RefSeq IDs
unmapped_idx <- which(is.na(all_probes$gene_symbol) | all_probes$gene_symbol == "")
refseq_ids <- all_probes$REFSEQ[unmapped_idx]

# Clean RefSeq IDs: remove quotes, version numbers, whitespace
refseq_clean <- refseq_ids
refseq_clean <- gsub("\"", "", refseq_clean)  # Remove quotes
refseq_clean <- gsub("\\..*$", "", refseq_clean)  # Remove version (.1, .2, etc)
refseq_clean <- trimws(refseq_clean)  # Remove whitespace
refseq_clean <- refseq_clean[!is.na(refseq_clean) & refseq_clean != ""]

# Get unique RefSeq IDs to map
unique_refseq <- unique(refseq_clean)
message(sprintf("Found %d unique RefSeq IDs to map\n", length(unique_refseq)))

# Get all valid REFSEQ keys from org.Hs.eg.db
valid_refseq_keys <- keys(org.Hs.eg.db, keytype = "REFSEQ")
message(sprintf("org.Hs.eg.db has %d valid REFSEQ keys\n", length(valid_refseq_keys)))

# Find which of our RefSeq IDs are valid
valid_to_map <- intersect(unique_refseq, valid_refseq_keys)
message(sprintf("  %d of our RefSeq IDs are valid in org.Hs.eg.db (%.1f%%)\n",
                length(valid_to_map),
                100 * length(valid_to_map) / length(unique_refseq)))

if (length(valid_to_map) > 0) {
  message(sprintf("Mapping %d valid RefSeq IDs...\n", length(valid_to_map)))
  
  # Map valid RefSeq to Symbols
  refseq_symbols <- suppressMessages(mapIds(org.Hs.eg.db,
                                           keys = valid_to_map,
                                           column = "SYMBOL",
                                           keytype = "REFSEQ",
                                           multiVals = "first"))
  
  message(sprintf("  Retrieved %d gene symbols\n", sum(!is.na(refseq_symbols))))
  
  # Map back to probes
  for (i in unmapped_idx) {
    refseq_id <- all_probes$REFSEQ[i]
    if (!is.na(refseq_id) && refseq_id != "") {
      # Clean the RefSeq ID
      refseq_clean_id <- gsub("\"", "", refseq_id)
      refseq_clean_id <- gsub("\\..*$", "", refseq_clean_id)
      refseq_clean_id <- trimws(refseq_clean_id)
      
      # Look up symbol
      if (refseq_clean_id %in% names(refseq_symbols)) {
        symbol <- refseq_symbols[refseq_clean_id]
        if (!is.na(symbol)) {
          all_probes$gene_symbol[i] <- symbol
        }
      }
    }
  }
  
  n_after <- sum(!is.na(all_probes$gene_symbol) & all_probes$gene_symbol != "", na.rm = TRUE)
  improvement <- n_after - n_before
  
  message(sprintf("\n=== RefSeq Mapping Complete ==="))
  message(sprintf("Before: %d mapped (%.1f%%)", n_before, 100 * n_before / nrow(all_probes)))
  message(sprintf("After:  %d mapped (%.1f%%)", n_after, 100 * n_after / nrow(all_probes)))
  message(sprintf("Improvement: +%d probes (%.1f%% → %.1f%%)\n",
                  improvement,
                  100 * n_before / nrow(all_probes),
                  100 * n_after / nrow(all_probes)))
  
  # Save updated annotation
  write.csv(all_probes, probe_file, row.names = FALSE)
  message(sprintf("✅ Saved: %s", probe_file))
  
  # Also save to complete file
  complete_file <- "results/cohort-GSE29801/dry_amd_de/probe_annotation_complete.csv"
  write.csv(all_probes, complete_file, row.names = FALSE)
  message(sprintf("✅ Saved: %s", complete_file))
} else {
  message("⚠️  No valid RefSeq IDs found that match org.Hs.eg.db")
}
