#!/usr/bin/env Rscript
# ============================================================
# Intron confidence scoring from a protein alignment
# ------------------------------------------------------------
# What this script does:
#   1) Reads an intron table that contains:
#        - Gene
#        - CodonIndexBefore   (codon index immediately BEFORE the intron)
#        - Intron             (optional but recommended)
#   2) Reads a protein multiple sequence alignment (FASTA)
#   3) Maps CodonIndexBefore -> alignment column (aln_col)
#   4) Computes, for each intron, alignment-quality metrics in a ±WINDOW window:
#        - gap_frac: fraction of gap characters ("-") across all sequences in the window
#        - entropy: mean Shannon entropy across columns in the window (ignoring gaps)
#   5) Applies simple confidence cutoffs for gap_frac and entropy, then a combined label:
#        - conf_gap: high / medium / low
#        - conf_entropy: high / medium / low
#        - conf_final: high / medium / low
#   6) Writes a single output table with the original columns + mapping + metrics + labels

# ============================================================

suppressPackageStartupMessages({
  library(Biostrings)
  library(readr)
  library(dplyr)
  library(tibble)
})

# ----------------------------
# 0) USER SETTINGS
# ----------------------------

# Intron table (TSV). Must contain:
#   - Gene
#   - CodonIndexBefore
#   - Intron (optional but helpful)
IN_INTRONS_TSV <- "results/introns.tsv"

# Protein alignment (FASTA)
IN_ALN_FASTA <- "results/protein_alignment.aln.fasta"

# Output table (TSV)
OUT_CONF_TSV <- "results/introns_confidence.tsv"

# Window radius for metrics: ±WINDOW alignment columns
WINDOW <- 10

# Cutoffs for gap fraction (lower is better)
# Example interpretation:
#   <= 0.10  -> high confidence
#   <= 0.30  -> medium confidence
#   >  0.30  -> low confidence
GAP_CUTOFFS <- c(high = 0.10, medium = 0.30)

# Cutoffs for entropy (lower is better)
# Example interpretation:
#   <= 0.4   -> high confidence
#   <= 0.8   -> medium confidence
#   >  0.8   -> low confidence
ENTROPY_CUTOFFS <- c(high = 0.4, medium = 0.8)

dir.create(dirname(OUT_CONF_TSV), showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# 1) LOAD INPUTS
# ----------------------------

intr <- read_tsv(IN_INTRONS_TSV, col_types = cols(.default = "c"))
stopifnot(all(c("Gene", "CodonIndexBefore") %in% names(intr)))

# Convert CodonIndexBefore to integer and drop unusable rows
intr <- intr %>%
  mutate(CodonIndexBefore = suppressWarnings(as.integer(CodonIndexBefore))) %>%
  filter(!is.na(Gene), !is.na(CodonIndexBefore), CodonIndexBefore >= 1)

aln <- readAAStringSet(IN_ALN_FASTA)

# Use only the first token of each FASTA header as the gene name
names(aln) <- sub("\\s.*$", "", names(aln))

# Keep only introns for genes present in the alignment
missing_genes <- sort(setdiff(unique(intr$Gene), names(aln)))
if (length(missing_genes) > 0) {
  message(
    "NOTE: Dropping ", length(missing_genes),
    " gene(s) from intron table because they are missing from the alignment."
  )
}

intr <- intr %>% filter(Gene %in% names(aln))
stopifnot(nrow(intr) > 0)

# ----------------------------
# 2) MAP CodonIndexBefore -> ALIGNMENT COLUMN
# ----------------------------
# Given an aligned protein sequence (with gaps) and a residue index (1-based, ungapped),
# return the alignment column where that residue occurs.

aa_col_of_residue <- function(aln_seq_string, residue_index_1based) {
  if (is.na(residue_index_1based) || residue_index_1based < 1) return(NA_integer_)
  
  chars <- strsplit(aln_seq_string, "", fixed = TRUE)[[1]]
  count <- 0L
  
  for (i in seq_along(chars)) {
    if (chars[i] != "-") {
      count <- count + 1L
      if (count == residue_index_1based) return(i)
    }
  }
  
  NA_integer_
}

gene_aln_df <- tibble(
  Gene       = names(aln),
  aln_string = as.character(aln)
)

intr <- intr %>%
  left_join(gene_aln_df, by = "Gene") %>%
  mutate(
    aln_col = purrr::map2_int(aln_string, CodonIndexBefore, aa_col_of_residue)
  )

# If anything fails to map, show and stop (it usually means CodonIndexBefore is too large)
if (any(is.na(intr$aln_col))) {
  bad <- intr %>%
    filter(is.na(aln_col)) %>%
    select(Gene, Intron, CodonIndexBefore) %>%
    distinct()
  print(bad)
  stop("Some introns could not be mapped into the alignment (NA aln_col). Fix those rows and rerun.")
}

# ----------------------------
# 3) PREP ALIGNMENT MATRIX FOR WINDOW METRICS
# ----------------------------
# Convert the alignment into a character matrix: [sequence x column]
aln_mat <- do.call(rbind, strsplit(as.character(aln), "", fixed = TRUE))
L <- ncol(aln_mat)

# ----------------------------
# 4) GAP FRACTION + SHANNON ENTROPY IN ±WINDOW
# ----------------------------

# Gap fraction over a window around a focal alignment column:
# fraction of characters that are "-" across ALL sequences and ALL columns in the window
gap_frac_window <- function(col) {
  lo <- max(1L, col - WINDOW)
  hi <- min(L,  col + WINDOW)
  mean(aln_mat[, lo:hi, drop = FALSE] == "-")
}

# Shannon entropy of a single alignment column (ignoring gaps)
shannon_entropy_col <- function(x) {
  keep <- x[x != "-"]
  if (length(keep) == 0) return(NA_real_)
  p <- table(keep) / length(keep)
  -sum(p * log2(p))
}

# Mean entropy across all columns in the window (ignoring gaps per column)
entropy_window <- function(col) {
  lo <- max(1L, col - WINDOW)
  hi <- min(L,  col + WINDOW)
  cols <- aln_mat[, lo:hi, drop = FALSE]
  mean(apply(cols, 2, shannon_entropy_col), na.rm = TRUE)
}

intr <- intr %>%
  mutate(
    gap_frac = vapply(aln_col, gap_frac_window, numeric(1)),
    entropy  = vapply(aln_col, entropy_window,  numeric(1))
  )

# ----------------------------
# 5) APPLY CUTOFFS -> CONFIDENCE LABELS
# ----------------------------

# Gap confidence
intr <- intr %>%
  mutate(
    conf_gap = case_when(
      gap_frac <= GAP_CUTOFFS[["high"]]   ~ "high",
      gap_frac <= GAP_CUTOFFS[["medium"]] ~ "medium",
      TRUE                                ~ "low"
    ),
    conf_entropy = case_when(
      entropy <= ENTROPY_CUTOFFS[["high"]]   ~ "high",
      entropy <= ENTROPY_CUTOFFS[["medium"]] ~ "medium",
      TRUE                                    ~ "low"
    ),
    # Final confidence rule:
    # - high only if BOTH metrics are high
    # - low if EITHER metric is low
    # - otherwise medium
    conf_final = case_when(
      conf_gap == "high" & conf_entropy == "high" ~ "high",
      conf_gap == "low"  | conf_entropy == "low"  ~ "low",
      TRUE                                        ~ "medium"
    )
  )

# Quick summary prints
message("conf_gap counts:")
print(table(intr$conf_gap, useNA = "ifany"))
message("conf_entropy counts:")
print(table(intr$conf_entropy, useNA = "ifany"))
message("conf_final counts:")
print(table(intr$conf_final, useNA = "ifany"))

# ----------------------------
# 6) WRITE OUTPUT TABLE
# ----------------------------
# Keep original columns and add:
#   aln_col, gap_frac, entropy, conf_gap, conf_entropy, conf_final
# Drop the alignment string column before writing.

out <- intr %>%
  select(-aln_string) %>%
  arrange(Gene, aln_col)

write_tsv(out, OUT_CONF_TSV)
message("Saved intron confidence table to: ", OUT_CONF_TSV)
