# Final solution: Extract probe-to-gene mapping from raw Agilent data files
# Create complete annotation for GSE29801 results

library(dplyr)

#%%
# Step 1: Extract all probe IDs from a raw data file in correct order
message("=== Step 1: Extracting Probe IDs ===")

raw_file <- "data/external/geo/GSE29801/raw_data/GSM738525_RPEChoroid_93.txt.gz"
raw_data <- read.delim(gzfile(raw_file), stringsAsFactors = FALSE, header = TRUE)

# Extract probe IDs
all_probes <- data.frame(
  row_number = 1:nrow(raw_data),
  agilent_id = raw_data$ID,
  stringsAsFactors = FALSE
)

message(sprintf("Extracted %d probe IDs from raw data", nrow(all_probes)))

# Show mapping for row numbers used in results (12-100)
message("\nExample row number -> probe ID mapping:")
print(all_probes[c(12, 50, 100), ])

#%%
# Step 2: Try to download GPL platform annotation for gene symbols
message("\n=== Step 2: Downloading GPL4133 Annotation ===")

# Try the newer format
gpl_url2 <- "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL4133&targ=self&form=text&view=data"
gpl_file2 <- "data/external/geo/GSE29801/GPL4133_annotation.txt"

if (!file.exists(gpl_file2)) {
  tryCatch({
    download.file(gpl_url2, gpl_file2, method = "auto", quiet = FALSE)
    message("  Downloaded GPL4133 platform data")
  }, error = function(e) {
    message(sprintf("  Could not download: %s", e$message))
  })
}

#%%
# Step 3: Parse GPL annotation if available
if (file.exists(gpl_file2)) {
  message("\n===Step 3: Parsing GPL Annotation ===")
  
  # Read GPL file
  gpl_lines <- readLines(gpl_file2)
  
  # Find table section
  table_start <- which(grepl("^!platform_table_begin", gpl_lines))
  table_end <- which(grepl("^!platform_table_end", gpl_lines))
  
  if (length(table_start) > 0 && length(table_end) > 0) {
    # Extract table data
    table_lines <- gpl_lines[(table_start + 1):(table_end - 1)]
    
    # Parse as TSV
    gpl_table <- read.table(text = table_lines, sep = "\t", header = TRUE, 
                            quote = "", fill = TRUE, stringsAsFactors = FALSE)
    
    message(sprintf("  Parsed GPL table: %d rows, columns: %s",
                    nrow(gpl_table), paste(colnames(gpl_table), collapse=", ")))
    
    # Merge with probe list
    all_probes <- all_probes %>%
      left_join(gpl_table, by = c("agilent_id" = "ID"))
    
    # Check if we got gene symbols
    if ("GENE_SYMBOL" %in% colnames(all_probes)) {
      all_probes <- all_probes %>%
        rename(gene_symbol = GENE_SYMBOL)
    } else if ("Gene.Symbol" %in% colnames(all_probes)) {
      all_probes <- all_probes %>%
        rename(gene_symbol = Gene.Symbol)
    }
    
    n_mapped <- sum(!is.na(all_probes$gene_symbol) & all_probes$gene_symbol != "", na.rm = TRUE)
    message(sprintf("  Mapped %d / %d probes to gene symbols (%.1f%%)",
                    n_mapped, nrow(all_probes),
                    100 * n_mapped / nrow(all_probes)))
  }
} else {
  message("\n⚠️  GPL annotation file not available")
  message("   Will save probe IDs only - you can add gene symbols later")
  all_probes$gene_symbol <- NA
}

#%%
# Save complete probe annotation
write.csv(all_probes,
          "results/cohort-GSE29801/dry_amd_de/probe_annotation_complete.csv",
          row.names = FALSE)

message("\n✅ Saved complete probe annotation")

#%%
# Function to annotate result files
annotate_results <- function(results_file, output_file = NULL) {
  
  results <- read.csv(results_file)
  
  if (!"gene" %in% colnames(results)) {
    stop("No 'gene' column found")
  }
  
  results$gene <- as.integer(results$gene)
  
  # Merge with probe annotation
  cols_to_add <- c("row_number", "agilent_id")
  if ("gene_symbol" %in% colnames(all_probes)) {
    cols_to_add <- c(cols_to_add, "gene_symbol")
  }
  
  results_annotated <- results %>%
    left_join(all_probes[, cols_to_add], by = c("gene" = "row_number"))  %>%
    select(gene, agilent_id, everything())
  
  # Move gene_symbol after agilent_id if it exists
  if ("gene_symbol" %in% colnames(results_annotated)) {
    results_annotated <- results_annotated %>%
      select(gene, agilent_id, gene_symbol, everything())
  }
  
  if (!is.null(output_file)) {
    write.csv(results_annotated, output_file, row.names = FALSE)
    message(sprintf("  ✅ Saved: %s", basename(output_file)))
  }
  
  return(results_annotated)
}

#%%
# Annotate all key result files
message("\n=== Annotating Result Files ===")

key_files <- c(
  "meta_combined_intermediate.csv",
  "meta_combined_early.csv",
  "meta_all_dry_amd.csv",
  "macular_variance_analysis.csv",
  "extramacular_variance_analysis.csv",
  "extramacular_age_main.csv",
  "macular_age_main.csv"
)

for (file in key_files) {
  input_path <- file.path("results/cohort-GSE29801/dry_amd_de", file)
  output_path <- sub("\\.csv$", "_annotated.csv", input_path)
  
  if (file.exists(input_path)) {
    message(sprintf("\n📝 %s", file))
    tryCatch({
      result <- annotate_results(input_path, output_path)
      
      # Show top results
      if ("meta_fdr" %in% colnames(result)) {
        top <- result %>% arrange(meta_p) %>%
          select(any_of(c("agilent_id", "gene_symbol", "meta_logFC", "meta_p", "meta_fdr"))) %>%
          head(10)
      } else if ("f_fdr" %in% colnames(result)) {
        top <- result %>% arrange(f_pvalue) %>%
          select(any_of(c("agilent_id", "gene_symbol", "var_fc", "f_pvalue", "f_fdr"))) %>%
          head(10)
      } else if ("P.Value" %in% colnames(result)) {
        top <- result %>% arrange(P.Value) %>%
          select(any_of(c("agilent_id", "gene_symbol", "logFC", "P.Value", "adj.P.Val"))) %>%
          head(10)
      }
      
      if (exists("top") && ncol(top) > 0) {
        print(top, row.names = FALSE)
      }
    }, error = function(e) {
      message(sprintf("  ❌ Error: %s", e$message))
    })
  }
}

message("\n=== Annotation Complete ===")
message("All _annotated.csv files now have:")
message("  - agilent_id: Agilent probe ID (e.g., A_24_P66027)")
message("  - gene_symbol: Gene symbol (if GPL annotation available)")
