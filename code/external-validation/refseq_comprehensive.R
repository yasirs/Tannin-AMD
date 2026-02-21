# Comprehensive RefSeq Mapping - Multiple Strategies
# Tries different RefSeq types and multiple data sources

library(dplyr)

message("=== Comprehensive RefSeq Mapping Strategy ===\n")

# Load annotation
probe_file <- "results/cohort-GSE29801/dry_amd_de/probe_annotation_improved.csv"
all_probes <- read.csv(probe_file, stringsAsFactors = FALSE)

n_start <- sum(!is.na(all_probes$gene_symbol) & all_probes$gene_symbol != "", na.rm = TRUE)
message(sprintf("Starting: %d mapped (%.1f%%)\n", n_start, 100 * n_start / nrow(all_probes)))

# Get unmapped with RefSeq
unmapped_idx <- which(is.na(all_probes$gene_symbol) | all_probes$gene_symbol == "")
refseq_raw <- all_probes$REFSEQ[unmapped_idx]

# Clean RefSeq IDs
refseq_clean <- gsub("\"", "", refseq_raw)
refseq_clean <- trimws(refseq_clean)
refseq_clean <- refseq_clean[!is.na(refseq_clean) & refseq_clean != ""]

# Analyze RefSeq types
refseq_types <- data.frame(
  type = c("NM_", "NR_", "XM_", "XR_", "Other"),
  count = c(
    sum(grepl("^NM_", refseq_clean)),
    sum(grepl("^NR_", refseq_clean)),
    sum(grepl("^XM_", refseq_clean)),
    sum(grepl("^XR_", refseq_clean)),
    sum(!grepl("^(NM_|NR_|XM_|XR_)", refseq_clean))
  )
)

message("RefSeq ID types found:")
print(refseq_types)
message("")

unique_refseq <- unique(refseq_clean)
message(sprintf("Total unique RefSeq IDs: %d\n", length(unique_refseq)))

#%% Strategy 1: Try org.Hs.eg.db with different keytypes
if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  
  message("=== Strategy 1: org.Hs.eg.db (all RefSeq types) ===")
  
  # Try with version numbers removed
  refseq_no_version <- gsub("\\..*$", "", unique_refseq)
  
  # Check overlap
  valid_keys <- keys(org.Hs.eg.db, keytype = "REFSEQ")
  overlap <- intersect(refseq_no_version, valid_keys)
  message(sprintf("  Match with no version: %d / %d (%.1f%%)",
                  length(overlap), length(unique_refseq),
                  100 * length(overlap) / length(unique_refseq)))
  
  if (length(overlap) > 0) {
    symbols <- suppressMessages(mapIds(org.Hs.eg.db,
                                      keys = overlap,
                                      column = "SYMBOL", 
                                      keytype = "REFSEQ",
                                      multiVals = "first"))
    
    # Map back
    for (i in unmapped_idx) {
      refseq_id <- gsub("\"", "", all_probes$REFSEQ[i])
      refseq_id <- gsub("\\..*$", "", refseq_id)
      if (refseq_id %in% names(symbols) && !is.na(symbols[refseq_id])) {
        all_probes$gene_symbol[i] <- symbols[refseq_id]
      }
    }
    
    n_now <- sum(!is.na(all_probes$gene_symbol) & all_probes$gene_symbol != "")
    message(sprintf("  Mapped: +%d probes\n", n_now - n_start))
    n_start <- n_now
  } else {
    message("  No matches found\n")
  }
}

#%% Strategy 2: biomaRt with current Ensembl
if (requireNamespace("biomaRt", quietly = TRUE)) {
  library(biomaRt)
  
  message("=== Strategy 2: biomaRt (current Ensembl) ===")
  
  tryCatch({
    ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
    message("  Connected to current Ensembl")
    
    # Try mapping
    refseq_to_map <- unique(gsub("\\..*$", "", refseq_clean))
    
    mapping <- getBM(
      attributes = c('refseq_mrna', 'hgnc_symbol'),
      filters = 'refseq_mrna',
      values = refseq_to_map,
      mart = ensembl
    )
    
    mapping <- mapping[mapping$hgnc_symbol != "", ]
    message(sprintf("  Retrieved %d mappings", nrow(mapping)))
    
    # Also try predicted RefSeq
    mapping_pred <- getBM(
      attributes = c('refseq_mrna_predicted', 'hgnc_symbol'),
      filters = 'refseq_mrna_predicted',
      values = refseq_to_map,
      mart = ensembl
    )
    
    mapping_pred <- mapping_pred[mapping_pred$hgnc_symbol != "", ]
    if (nrow(mapping_pred) > 0) {
      colnames(mapping_pred)[1] <- "refseq_mrna"
      mapping <- rbind(mapping, mapping_pred)
      mapping <- unique(mapping)
    }
    message(sprintf("  Total with predicted: %d mappings", nrow(mapping)))
    
    # Map back
    unmapped_idx_now <- which(is.na(all_probes$gene_symbol) | all_probes$gene_symbol == "")
    for (i in unmapped_idx_now) {
      refseq_id <- gsub("\"", "", all_probes$REFSEQ[i])
      refseq_id <- gsub("\\..*$", "", refseq_id)
      
      matches <- mapping[mapping$refseq_mrna == refseq_id, ]
      if (nrow(matches) > 0) {
        all_probes$gene_symbol[i] <- matches$hgnc_symbol[1]
      }
    }
    
    n_now <- sum(!is.na(all_probes$gene_symbol) & all_probes$gene_symbol != "")
    improvement <- n_now - n_start
    message(sprintf("  Mapped: +%d probes\n", improvement))
    n_start <- n_now
    
  }, error = function(e) {
    message(sprintf("  Error: %s\n", e$message))
  })
  
  #%% Strategy 3: Try Ensembl archive (older version from ~2010)
  message("=== Strategy 3: biomaRt (Ensembl archive GRCh37) ===")
  
  tryCatch({
    # Try GRCh37 (older assembly, may match older GPL platform)
    ensembl_old <- useEnsembl(biomart = "genes", 
                              dataset = "hsapiens_gene_ensembl",
                              GRCh = 37)
    message("  Connected to Ensembl GRCh37")
    
    unmapped_idx_now <- which(is.na(all_probes$gene_symbol) | all_probes$gene_symbol == "")
    refseq_unmapped <- all_probes$REFSEQ[unmapped_idx_now]
    refseq_unmapped <- gsub("\"", "", refseq_unmapped)
    refseq_unmapped <- gsub("\\..*$", "", refseq_unmapped)
    refseq_unmapped <- unique(refseq_unmapped[refseq_unmapped != "" & !is.na(refseq_unmapped)])
    
    if (length(refseq_unmapped) > 0) {
      mapping_old <- getBM(
        attributes = c('refseq_mrna', 'hgnc_symbol'),
        filters = 'refseq_mrna',
        values = refseq_unmapped,
        mart = ensembl_old
      )
      
      mapping_old <- mapping_old[mapping_old$hgnc_symbol != "", ]
      message(sprintf("  Retrieved %d mappings from GRCh37", nrow(mapping_old)))
      
      # Map back
      for (i in unmapped_idx_now) {
        refseq_id <- gsub("\"", "", all_probes$REFSEQ[i])
        refseq_id <- gsub("\\..*$", "", refseq_id)
        
        matches <- mapping_old[mapping_old$refseq_mrna == refseq_id, ]
        if (nrow(matches) > 0) {
          all_probes$gene_symbol[i] <- matches$hgnc_symbol[1]
        }
      }
      
      n_now <- sum(!is.na(all_probes$gene_symbol) & all_probes$gene_symbol != "")
      improvement <- n_now - n_start
      message(sprintf("  Mapped: +%d probes\n", improvement))
      n_start <- n_now
    }
    
  }, error = function(e) {
    message(sprintf("  Error with GRCh37: %s\n", e$message))
  })
}

# Final results
n_final <- sum(!is.na(all_probes$gene_symbol) & all_probes$gene_symbol != "")
total_improvement <- n_final - (n_start - (n_final - n_start))

message("\n=== Final Results ===")
message(sprintf("Total improvement from RefSeq: +%d probes", total_improvement))
message(sprintf("Final: %d mapped (%.1f%%)\n", n_final, 100 * n_final / nrow(all_probes)))

if (total_improvement > 0) {
  write.csv(all_probes, probe_file, row.names = FALSE)
  message(sprintf("✅ Saved: %s", probe_file))
  
  complete_file <- "results/cohort-GSE29801/dry_amd_de/probe_annotation_complete.csv"
  write.csv(all_probes, complete_file, row.names = FALSE)
  message(sprintf("✅ Saved: %s", complete_file))
}
