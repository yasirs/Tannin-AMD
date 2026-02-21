# GSE29801 Dry AMD Data Preparation Script
# Filters to RPE-choroid tissues and dry AMD subtypes only

#%%
# Load libraries
library(dplyr)
library(tidyr)

#%%
# Configuration
BASE_DIR <- "/home/ysuhail/work/Tannin-AMD"
DATA_DIR <- file.path(BASE_DIR, "data/external/geo/GSE29801")
OUTPUT_DIR <- file.path(BASE_DIR, "results/cohort-GSE29801/dry_amd_de")

# Create output directory
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

#%%
# Parse series matrix to extract metadata
parse_series_matrix <- function(path) {
  con <- gzfile(path, "rt")
  lines <- readLines(con)
  close(con)
  
  # Find sample characteristics lines
  char_lines <- lines[grepl("^!Sample_characteristics_ch1", lines)]
  geo_line <- lines[grepl("^!Sample_geo_accession", lines)]
  
  # Parse sample IDs
  sample_ids <- strsplit(geo_line, "\t")[[1]][-1]
  sample_ids <- gsub("\"", "", sample_ids)
  
  # Initialize metadata data frame
  metadata <- data.frame(sample_id = sample_ids, stringsAsFactors = FALSE)
  
  # Parse each characteristic line
  for (line in char_lines) {
    parts <- strsplit(line, "\t")[[1]][-1]
    parts <- gsub("\"", "", parts)
    
    # Extract key-value pairs
    if (length(parts) > 0 && grepl(":", parts[1])) {
      key <- trimws(sub(":.*", "", parts[1]))
      values <- sapply(parts, function(x) trimws(sub("^[^:]+:\\s*", "", x)))
      metadata[[key]] <- values
    }
  }
  
  return(metadata)
}

#%%
# Load metadata
message("Loading metadata from series matrix...")
metadata <- parse_series_matrix(file.path(DATA_DIR, "GSE29801_series_matrix.txt.gz"))

# Clean up column names
colnames(metadata) <- gsub(" ", "_", colnames(metadata))
colnames(metadata) <- gsub("[()]", "", colnames(metadata))

message(sprintf("Loaded metadata for %d samples", nrow(metadata)))

#%%
# Filter to RPE-choroid tissues only
message("\nFiltering to RPE-choroid tissues...")
rpe_choroid_samples <- metadata %>%
  filter(grepl("RPE-choroid", tissue))

message(sprintf("  Retained %d RPE-choroid samples (excluded %d retina samples)", 
                nrow(rpe_choroid_samples), 
                nrow(metadata) - nrow(rpe_choroid_samples)))

#%%
# Filter to dry AMD and controls only (exclude wet AMD and unspecified)
message("\nFiltering to dry AMD and controls...")

# Define dry AMD subtypes and controls
dry_amd_types <- c("normal", "MD1", "MD2", "dry AMD", "GA")

dry_samples <- rpe_choroid_samples %>%
  filter(amd_classification %in% dry_amd_types)

excluded <- rpe_choroid_samples %>%
  filter(!amd_classification %in% dry_amd_types)

message(sprintf("  Retained %d samples:", nrow(dry_samples)))
table(dry_samples$amd_classification) %>% print()

message(sprintf("\n  Excluded %d samples:", nrow(excluded)))
if (nrow(excluded) > 0) {
  table(excluded$amd_classification) %>% print()
}

#%%
# Create disease status and stage variables
dry_samples <- dry_samples %>%
  mutate(
    disease_status = ifelse(amd_classification == "normal", "Control", "AMD"),
    dry_stage = case_when(
      amd_classification == "normal" ~ "Control",
      amd_classification == "MD1" ~ "MD1",
      amd_classification == "MD2" ~ "MD2",
      amd_classification == "dry AMD" ~ "dry_AMD",
      amd_classification == "GA" ~ "GA"
    ),
    dry_stage = factor(dry_stage, levels = c("Control", "MD1", "MD2", "dry_AMD", "GA"))
  )

#%%
# Parse age and sex
dry_samples <- dry_samples %>%
  mutate(
    age = as.numeric(age_years),
    sex = factor(gender, levels = c("male", "female"))
  )

# Check for RIN (RNA integrity number)
if ("rna_integrity_number_rin" %in% colnames(dry_samples)) {
  dry_samples <- dry_samples %>%
    mutate(rin = as.numeric(rna_integrity_number_rin))
  message("\n✅ RIN values available")
} else {
  dry_samples$rin <- NA
  message("\n⚠️  RIN values not available")
}

#%%
# Separate by tissue location
macular_samples <- dry_samples %>%
  filter(tissue == "macular RPE-choroid")

extramacular_samples <- dry_samples %>%
  filter(tissue == "extramacular RPE-choroid")

message(sprintf("\nSample distribution:"))
message(sprintf("  Macular RPE-choroid: %d samples", nrow(macular_samples)))
table(macular_samples$dry_stage) %>% print()

message(sprintf("\n  Extramacular RPE-choroid: %d samples", nrow(extramacular_samples)))
table(extramacular_samples$dry_stage) %>% print()

#%%
# Load expression matrix
message("\nLoading expression matrix...")
expr_raw <- read.csv(gzfile(file.path(DATA_DIR, "GSE29801_expression_matrix.csv.gz")), 
                     row.names = 1, check.names = FALSE)

# Separate expression data from Group column
if ("Group" %in% colnames(expr_raw)) {
  expr_data <- expr_raw[, colnames(expr_raw) != "Group"]
} else {
  expr_data <- expr_raw
}

message(sprintf("Expression matrix: %d samples × %d probes", nrow(expr_data), ncol(expr_data)))

#%%
# Map probes to gene symbols
message("\nMapping probes to gene symbols...")

# Load GPL annotation
gpl_path <- file.path(DATA_DIR, "GSE29801_extracted/GPL4133_old_annotations.txt.gz")

if (file.exists(gpl_path)) {
  # Parse GPL
  id_to_sym <- list()
  con <- gzfile(gpl_path, 'rt')
  lines <- readLines(con)
  close(con)
  
  header_found <- FALSE
  for (line in lines) {
    line <- trimws(line)
    if (!header_found) {
      if (grepl("^ID", line) && grepl("GENE_SYMBOL", line)) {
        header_found <- TRUE
        headers <- strsplit(line, "\t")[[1]]
      }
      next
    }
    
    parts <- strsplit(line, "\t")[[1]]
    if (length(parts) < length(headers)) next
    
    tryCatch({
      id_idx <- which(headers == "ID")
      sym_idx <- which(headers == "GENE_SYMBOL")
      
      probe_id <- parts[id_idx]
      sym <- parts[sym_idx]
      
      if (length(sym) > 0 && sym != "" && sym != "nan") {
        id_to_sym[[probe_id]] <- sym
      }
    }, error = function(e) {})
  }
  
  message(sprintf("  Mapped %d probes to gene symbols", length(id_to_sym)))
  
  # Map column names
  probe_cols <- colnames(expr_data)
  gene_symbols <- sapply(probe_cols, function(p) {
    if (p %in% names(id_to_sym)) id_to_sym[[p]] else p
  })
  
  colnames(expr_data) <- gene_symbols
  
  # Collapse duplicate gene symbols by averaging
  expr_t <- as.data.frame(t(expr_data))
  expr_t$gene <- rownames(expr_t)
  
  expr_collapsed <- expr_t %>%
    group_by(gene) %>%
    summarise(across(everything(), mean, na.rm = TRUE)) %>%
    as.data.frame()
  
  rownames(expr_collapsed) <- expr_collapsed$gene
  expr_collapsed$gene <- NULL
  
  expr_data <- as.data.frame(t(expr_collapsed))
  
  message(sprintf("  After collapsing: %d unique genes", ncol(expr_data)))
} else {
  message("  ⚠️  GPL annotation file not found, using probe IDs")
}

#%%
# Filter expression data to dry AMD samples
macular_expr <- expr_data[macular_samples$sample_id, ]
extramacular_expr <- expr_data[extramacular_samples$sample_id, ]

# Ensure metadata order matches expression data
macular_samples <- macular_samples[match(rownames(macular_expr), macular_samples$sample_id), ]
extramacular_samples <- extramacular_samples[match(rownames(extramacular_expr), extramacular_samples$sample_id), ]

#%%
# Save prepared data
message("\nSaving prepared datasets...")

# Macular
saveRDS(list(
  expr = macular_expr,
  metadata = macular_samples
), file.path(OUTPUT_DIR, "macular_rpe_choroid_data.rds"))

# Extramacular
saveRDS(list(
  expr = extramacular_expr,
  metadata = extramacular_samples
), file.path(OUTPUT_DIR, "extramacular_rpe_choroid_data.rds"))

# Combined (for reference)
combined_samples <- rbind(
  macular_samples %>% mutate(tissue_location = "macular"),
  extramacular_samples %>% mutate(tissue_location = "extramacular")
)
combined_expr <- rbind(macular_expr, extramacular_expr)

saveRDS(list(
  expr = combined_expr,
  metadata = combined_samples
), file.path(OUTPUT_DIR, "combined_rpe_choroid_data.rds"))

#%%
# Summary statistics
message("\n=== Data Preparation Complete ===")
message(sprintf("Output directory: %s", OUTPUT_DIR))
message("\nSaved files:")
message("  - macular_rpe_choroid_data.rds")
message("  - extramacular_rpe_choroid_data.rds")
message("  - combined_rpe_choroid_data.rds")

message("\nSample summary:")
message(sprintf("  Macular: %d samples, %d genes", nrow(macular_expr), ncol(macular_expr)))
message(sprintf("  Extramacular: %d samples, %d genes", nrow(extramacular_expr), ncol(extramacular_expr)))

message("\nCovariate summary:")
message(sprintf("  Age range: %.0f - %.0f years", 
                min(combined_samples$age, na.rm=TRUE), 
                max(combined_samples$age, na.rm=TRUE)))
message(sprintf("  Sex distribution: %d male, %d female", 
                sum(combined_samples$sex == "male"), 
                sum(combined_samples$sex == "female")))
if (any(!is.na(combined_samples$rin))) {
  message(sprintf("  RIN range: %.1f - %.1f", 
                  min(combined_samples$rin, na.rm=TRUE), 
                  max(combined_samples$rin, na.rm=TRUE)))
}
