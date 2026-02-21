# Additional Gene Bank Mapping in Batches
# Process GenBank accessions in smaller batches to avoid timeout

library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)

message("=== GenBank Batch Mapping (Strategy 5) ===\n")

# Load current annotation
probe_file <- "results/cohort-GSE29801/dry_amd_de/probe_annotation_complete.csv"
all_probes <- read.csv(probe_file, stringsAsFactors = FALSE)

n_before <- sum(!is.na(all_probes$gene_symbol) & all_probes$gene_symbol != "", na.rm = TRUE)
message(sprintf("Starting with %d mapped probes (%.1f%%)\n",
                n_before, 100 * n_before / nrow(all_probes)))

# Get unmapped probes with GenBank IDs
unmapped_idx <- which(is.na(all_probes$gene_symbol) | all_probes$gene_symbol == "")
genbank_ids <- all_probes$GB_ACC[unmapped_idx]
genbank_ids <- genbank_ids[!is.na(genbank_ids) & genbank_ids != ""]

# Create mapping from index to GenBank ID
idx_to_gb <- setNames(unmapped_idx[genbank_ids != "" & !is.na(genbank_ids)],
                      genbank_ids[genbank_ids != "" & !is.na(genbank_ids)])

message(sprintf("Found %d unique GenBank IDs to map\n", length(unique(names(idx_to_gb)))))

# Process in batches of 500
batch_size <- 500
unique_gb <- unique(names(idx_to_gb))
n_batches <- ceiling(length(unique_gb) / batch_size)
total_mapped <- 0

message(sprintf("Processing in %d batches of %d...\n", n_batches, batch_size))

for (batch_num in 1:n_batches) {
  start_idx <- (batch_num - 1) * batch_size + 1
  end_idx <- min(batch_num * batch_size, length(unique_gb))
  batch_gb <- unique_gb[start_idx:end_idx]
  
  message(sprintf("  Batch %d/%d: mapping %d GenBank IDs...", 
                  batch_num, n_batches, length(batch_gb)))
  
  tryCatch({
    gb_symbols <- suppressMessages(mapIds(org.Hs.eg.db,
                         keys = batch_gb,
                         column = "SYMBOL",
                         keytype = "ACCNUM",
                         multiVals = "first"))
    
    # Map back to probes
    batch_mapped <- 0
    for (gb_id in names(gb_symbols)) {
      symbol <- gb_symbols[gb_id]
      if (!is.na(symbol)) {
        # Find all probes with this GenBank ID
        probe_indices <- idx_to_gb[names(idx_to_gb) == gb_id]
        for (idx in probe_indices) {
          if (is.na(all_probes$gene_symbol[idx]) || all_probes$gene_symbol[idx] == "") {
            all_probes$gene_symbol[idx] <- symbol
            batch_mapped <- batch_mapped + 1
          }
        }
      }
    }
    
    total_mapped <- total_mapped + batch_mapped
    message(sprintf("    Mapped %d probes in this batch", batch_mapped))
    
  }, error = function(e) {
    message(sprintf("    ⚠️  Batch %d failed: %s", batch_num, e$message))
  })
}

n_after <- sum(!is.na(all_probes$gene_symbol) & all_probes$gene_symbol != "", na.rm = TRUE)
improvement <- n_after - n_before

message(sprintf("\n=== GenBank Batch Mapping Complete ==="))
message(sprintf("Before: %d mapped (%.1f%%)", n_before, 100 * n_before / nrow(all_probes)))
message(sprintf("After:  %d mapped (%.1f%%)", n_after, 100 * n_after / nrow(all_probes)))
message(sprintf("Improvement: +%d probes (%.1f%% → %.1f%%)\n",
                improvement,
                100 * n_before / nrow(all_probes),
                100 * n_after / nrow(all_probes)))

# Save updated annotation
backup_file <- gsub("\\.csv$", "_backup2.csv", probe_file)
file.copy(probe_file, backup_file, overwrite = TRUE)
write.csv(all_probes, probe_file, row.names = FALSE)

improved_file <- "results/cohort-GSE29801/dry_amd_de/probe_annotation_improved.csv"
write.csv(all_probes, improved_file, row.names = FALSE)

message(sprintf("✅ Saved: %s", probe_file))
message(sprintf("✅ Saved: %s", improved_file))
