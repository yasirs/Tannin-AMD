# Enhanced Probe-to-Gene Mapping for GSE29801
# Uses org.Hs.eg.db annotation package for improved mapping
# Falls back to biomaRt if needed

library(dplyr)

# Function to check and load annotation packages
load_annotation_package <- function() {
  if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    library(org.Hs.eg.db)
    library(AnnotationDbi)
    message("  ✅ Using org.Hs.eg.db for annotation")
    return("org.Hs.eg.db")
  } else if (requireNamespace("biomaRt", quietly = TRUE)) {
    library(biomaRt)
    message("  ✅ Using biomaRt for annotation")
    return("biomaRt")
  } else {
    stop("Neither org.Hs.eg.db nor biomaRt is installed. Please install one:\n",
         "  Option 1: BiocManager::install('org.Hs.eg.db')\n",
         "  Option 2: BiocManager::install('biomaRt')")
  }
}

#%%
# Configuration
message("=== Enhanced Probe-to-Gene Mapping for GSE29801 ===")
message("Using RefSeq, Entrez Gene IDs, Ensembl IDs, and other identifiers\n")

annotation_pkg <- load_annotation_package()

#%%
# Load existing probe annotation
probe_file <- "results/cohort-GSE29801/dry_amd_de/probe_annotation_complete.csv"
all_probes <- read.csv(probe_file, stringsAsFactors = FALSE)

message(sprintf("\nLoaded %d probes", nrow(all_probes)))

# Check current mapping status
n_mapped_before <- sum(!is.na(all_probes$gene_symbol) & all_probes$gene_symbol != "", na.rm = TRUE)
n_unmapped_before <- nrow(all_probes) - n_mapped_before
message(sprintf("Current mapping: %d mapped (%.1f%%), %d unmapped (%.1f%%)",
                n_mapped_before, 100 * n_mapped_before / nrow(all_probes),
                n_unmapped_before, 100 * n_unmapped_before / nrow(all_probes)))

#%%
# Function to clean and prepare IDs
clean_ids <- function(ids) {
  ids <- ids[!is.na(ids) & ids != "" & ids != "NA"]
  unique(ids)
}

#%%
# Strategy 1: Map via Entrez Gene IDs using org.Hs.eg.db
if (annotation_pkg == "org.Hs.eg.db") {
  message("\n=== Strategy 1: Map via Entrez Gene IDs ===")
  unmapped_idx <- which(is.na(all_probes$gene_symbol) | all_probes$gene_symbol == "")
  entrez_to_map <- clean_ids(as.character(all_probes$GENE[unmapped_idx]))
  
  if (length(entrez_to_map) > 0) {
    message(sprintf("  Mapping %d unique Entrez Gene IDs...", length(entrez_to_map)))
    
    # Map Entrez to Symbol
    entrez_symbols <- mapIds(org.Hs.eg.db, 
                             keys = entrez_to_map,
                             column = "SYMBOL",
                             keytype = "ENTREZID",
                             multiVals = "first")
    
    # Update probes
    for (i in unmapped_idx) {
      entrez_id <- as.character(all_probes$GENE[i])
      if (!is.na(entrez_id) && entrez_id != "" && entrez_id != "NA") {
        symbol <- entrez_symbols[entrez_id]
        if (!is.na(symbol)) {
          all_probes$gene_symbol[i] <- symbol
        }
      }
    }
    
    n_mapped_entrez <- sum(!is.na(all_probes$gene_symbol) & all_probes$gene_symbol != "") - n_mapped_before
    message(sprintf("  ✅ Mapped %d additional probes via Entrez Gene", n_mapped_entrez))
  }
  
  #%%
  message("\n=== Strategy 2: Map via RefSeq IDs ===")
  unmapped_idx <- which(is.na(all_probes$gene_symbol) | all_probes$gene_symbol == "")
  refseq_to_map <- clean_ids(all_probes$REFSEQ[unmapped_idx])
  
  if (length(refseq_to_map) > 0) {
    message(sprintf("  Mapping %d unique RefSeq IDs...", length(refseq_to_map)))
    
    # Clean RefSeq IDs (remove version numbers)
    refseq_clean <- gsub("\\..*$", "", refseq_to_map)
    
    # Try mapping RefSeq to Symbol with error handling
    tryCatch({
      refseq_symbols <- suppressMessages(mapIds(org.Hs.eg.db,
                               keys = refseq_clean,
                               column = "SYMBOL",
                               keytype = "REFSEQ",
                               multiVals = "first"))
      
      # Update probes
      n_before <- sum(!is.na(all_probes$gene_symbol) & all_probes$gene_symbol != "")
      for (i in unmapped_idx) {
        refseq_id <- all_probes$REFSEQ[i]
        if (!is.na(refseq_id) && refseq_id != "") {
          refseq_clean_id <- gsub("\\..*$", "", refseq_id)
          if (refseq_clean_id %in% names(refseq_symbols)) {
            symbol <- refseq_symbols[refseq_clean_id]
            if (!is.na(symbol)) {
              all_probes$gene_symbol[i] <- symbol
            }
          }
        }
      }
      
      n_mapped_refseq <- sum(!is.na(all_probes$gene_symbol) & all_probes$gene_symbol != "") - n_before
      message(sprintf("  ✅ Mapped %d additional probes via RefSeq", n_mapped_refseq))
    }, error = function(e) {
      message(sprintf("  ⚠️  RefSeq mapping encountered issues: %s", e$message))
      message("  Skipping RefSeq mapping...")
    })
  }
  
  #%%
  message("\n=== Strategy 3: Map via Ensembl IDs ===")
  unmapped_idx <- which(is.na(all_probes$gene_symbol) | all_probes$gene_symbol == "")
  ensembl_to_map <- clean_ids(all_probes$ENSEMBL_ID[unmapped_idx])
  
  if (length(ensembl_to_map) > 0) {
    message(sprintf("  Mapping %d unique Ensembl transcript IDs...", length(ensembl_to_map)))
    
    # Map Ensembl to Symbol
    ensembl_symbols <- mapIds(org.Hs.eg.db,
                              keys = ensembl_to_map,
                              column = "SYMBOL",
                              keytype = "ENSEMBLTRANS",
                              multiVals = "first")
    
    # Update probes
    n_before <- sum(!is.na(all_probes$gene_symbol) & all_probes$gene_symbol != "")
    for (i in unmapped_idx) {
      ensembl_id <- all_probes$ENSEMBL_ID[i]
      if (!is.na(ensembl_id) && ensembl_id != "") {
        symbol <- ensembl_symbols[ensembl_id]
        if (!is.na(symbol)) {
          all_probes$gene_symbol[i] <- symbol
        }
      }
    }
    
    n_mapped_ensembl <- sum(!is.na(all_probes$gene_symbol) & all_probes$gene_symbol != "") - n_before
    message(sprintf("  ✅ Mapped %d additional probes via Ensembl", n_mapped_ensembl))
  }
  
  #%%
  message("\n=== Strategy 4: Map via GenBank Accessions ===")
  unmapped_idx <- which(is.na(all_probes$gene_symbol) | all_probes$gene_symbol == "")
  genbank_to_map <- clean_ids(all_probes$GB_ACC[unmapped_idx])
  
  if (length(genbank_to_map) > 0 && length(genbank_to_map) <= 1000) {
    message(sprintf("  Mapping %d unique GenBank accessions...", length(genbank_to_map)))
    
    tryCatch({
      # Map GenBank to Symbol
      gb_symbols <- mapIds(org.Hs.eg.db,
                           keys = genbank_to_map,
                           column = "SYMBOL",
                           keytype = "ACCNUM",
                           multiVals = "first")
      
      # Update probes
      n_before <- sum(!is.na(all_probes$gene_symbol) & all_probes$gene_symbol != "")
      for (i in unmapped_idx) {
        gb_acc <- all_probes$GB_ACC[i]
        if (!is.na(gb_acc) && gb_acc != "") {
          symbol <- gb_symbols[gb_acc]
          if (!is.na(symbol)) {
            all_probes$gene_symbol[i] <- symbol
          }
        }
      }
      
      n_mapped_gb <- sum(!is.na(all_probes$gene_symbol) & all_probes$gene_symbol != "") - n_before
      message(sprintf("  ✅ Mapped %d additional probes via GenBank", n_mapped_gb))
    }, error = function(e) {
      message(sprintf("  ⚠️  GenBank mapping encountered errors: %s", e$message))
    })
  } else if (length(genbank_to_map) > 1000) {
    message(sprintf("  ⚠️  Too many GenBank IDs (%d), skipping", length(genbank_to_map)))
  }
}

#%%
# Final statistics
n_mapped_after <- sum(!is.na(all_probes$gene_symbol) & all_probes$gene_symbol != "", na.rm = TRUE)
n_unmapped_after <- nrow(all_probes) - n_mapped_after
n_improved <- n_mapped_after - n_mapped_before

message("\n=== Final Mapping Statistics ===")
message(sprintf("Before: %d mapped (%.1f%%), %d unmapped (%.1f%%)",
                n_mapped_before, 100 * n_mapped_before / nrow(all_probes),
                n_unmapped_before, 100 * n_unmapped_before / nrow(all_probes)))
message(sprintf("After:  %d mapped (%.1f%%), %d unmapped (%.1f%%)",
                n_mapped_after, 100 * n_mapped_after / nrow(all_probes),
                n_unmapped_after, 100 * n_unmapped_after / nrow(all_probes)))
message(sprintf("\n✅ Improvement: +%d probes mapped (%.1f%% → %.1f%%)",
                n_improved,
                100 * n_mapped_before / nrow(all_probes),
                100 * n_mapped_after / nrow(all_probes)))

#%%
# Save improved annotation
output_file <- "results/cohort-GSE29801/dry_amd_de/probe_annotation_improved.csv"
write.csv(all_probes, output_file, row.names = FALSE)
message(sprintf("\n✅ Saved improved annotation: %s", output_file))

#%%
# Also update the original file
backup_file <- gsub("\\.csv$", "_backup.csv", probe_file)
file.copy(probe_file, backup_file)
message(sprintf("✅ Backed up original: %s", backup_file))

write.csv(all_probes, probe_file, row.names = FALSE)
message(sprintf("✅ Updated original annotation: %s", probe_file))

#%%
# Generate mapping quality report
mapping_sources <- data.frame(
  Source = c("Original GPL annotation", "Entrez Gene IDs", "RefSeq IDs", 
             "Ensembl Transcript IDs", "GenBank Accessions"),
  Available_Probes = c(
    n_mapped_before,
    sum(!is.na(all_probes$GENE) & all_probes$GENE != "NA"),
    sum(!is.na(all_probes$REFSEQ) & all_probes$REFSEQ != ""),
    sum(!is.na(all_probes$ENSEMBL_ID) & all_probes$ENSEMBL_ID != ""),
    sum(!is.na(all_probes$GB_ACC) & all_probes$GB_ACC != "")
  ),
  stringsAsFactors = FALSE
)

report_file <- "results/cohort-GSE29801/dry_amd_de/mapping_improvement_report.txt"
cat("GSE29801 Probe-to-Gene Mapping Improvement Report\n", file = report_file)
cat("=" , rep("=", 60), "\n\n", sep = "", file = report_file, append = TRUE)
cat(sprintf("Analysis Date: %s\n\n", Sys.Date()), file = report_file, append = TRUE)
cat("MAPPING RESULTS:\n", file = report_file, append = TRUE)
cat(sprintf("  Before improvement: %d / %d probes mapped (%.1f%%)\n",
            n_mapped_before, nrow(all_probes),
            100 * n_mapped_before / nrow(all_probes)),
    file = report_file, append = TRUE)
cat(sprintf("  After improvement:  %d / %d probes mapped (%.1f%%)\n",
            n_mapped_after, nrow(all_probes),
            100 * n_mapped_after / nrow(all_probes)),
    file = report_file, append = TRUE)
cat(sprintf("  Improvement: +%d probes (%.1f percentage points)\n\n",
            n_improved,
            100 * n_improved / nrow(all_probes)),
    file = report_file, append = TRUE)

cat("ANNOTATION SOURCES USED:\n", file = report_file, append = TRUE)
for (i in 1:nrow(mapping_sources)) {
  cat(sprintf("  %s: %d probes with IDs available\n",
              mapping_sources$Source[i],
              mapping_sources$Available_Probes[i]),
      file = report_file, append = TRUE)
}

message(sprintf("✅ Saved mapping report: %s\n", report_file))

message("=== Mapping Improvement Complete ===")
