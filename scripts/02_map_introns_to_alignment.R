#!/usr/bin/env Rscript
# ============================================================
# Map intron positions onto a protein multiple-sequence alignment
# ------------------------------------------------------------
# What this script does (high level):
#   1) Reads a table of introns that includes a codon index
#      indicating where each intron occurs in the coding sequence
#   2) Reads a protein multiple sequence alignment (FASTA)
#   3) Converts each intron’s codon position into an alignment column
#      by counting non-gap residues in the aligned sequence
#   4) Builds:
#        - a tidy table of introns mapped to alignment columns
#        - a wide binary presence/absence matrix (genes × sites)
#
#
# Requirements:
#   Biostrings, dplyr, readr, tidyr, purrr, tibble
# ============================================================

suppressPackageStartupMessages({
  library(Biostrings)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(purrr)
  library(tibble)
})

# ----------------------------
# 0) USER SETTINGS
# ----------------------------

# Intron table (TSV). Must contain at least:
#   - Gene
#   - CodonIndexBefore  (codon index immediately before the intron)
IN_INTRONS <- "results/introns.tsv"

# Protein multiple-sequence alignment (FASTA)
IN_ALN <- "results/protein_alignment.aln.fasta"

# Output directory
OUT_DIR <- "results/intron_summaries"

# Output files
OUT_MAT  <- file.path(OUT_DIR, "intron_presence_matrix.csv")
OUT_TIDY <- file.path(OUT_DIR, "intron_positions_tidy.csv")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# 1) LOAD INTRON TABLE
# ----------------------------

# Read intron table as character columns to avoid parsing issues
introns <- read_tsv(IN_INTRONS, col_types = cols(.default = "c"))

# Basic sanity check for required columns
stopifnot(all(c("Gene", "CodonIndexBefore") %in% names(introns)))

# Convert codon index to integer and drop unusable rows
introns <- introns %>%
  mutate(
    CodonIndexBefore = suppressWarnings(as.integer(CodonIndexBefore))
  ) %>%
  filter(
    !is.na(Gene),
    !is.na(CodonIndexBefore)
  )

# ----------------------------
# 2) LOAD PROTEIN ALIGNMENT
# ----------------------------

# Read alignment as an AAStringSet
aln <- readAAStringSet(IN_ALN)

# Use only the first token of each FASTA header as the gene name
names(aln) <- sub("\\s.*$", "", names(aln))

# Store alignment sequences in a data frame
gene_aln_df <- tibble(
  Gene       = names(aln),
  aln_string = as.character(aln)
)

# Keep only introns for genes present in the alignment
introns <- introns %>%
  filter(Gene %in% gene_aln_df$Gene)

# ----------------------------
# 3) MAP CODON INDEX → ALIGNMENT COLUMN
# ----------------------------

# Helper function:
# Given an aligned protein sequence (with gaps) and a residue index
# (1-based, ungapped), return the alignment column where that residue occurs.
aa_col_of_residue <- function(aln_seq_string, residue_index_1based) {
  if (is.na(residue_index_1based) || residue_index_1based < 1) {
    return(NA_integer_)
  }
  
  chars <- strsplit(aln_seq_string, "", fixed = TRUE)[[1]]
  count <- 0L
  
  for (i in seq_along(chars)) {
    if (chars[i] != "-") {
      count <- count + 1L
      if (count == residue_index_1based) {
        return(i)
      }
    }
  }
  
  NA_integer_
}

# Join introns to alignment sequences and compute alignment columns
introns_mapped <- introns %>%
  left_join(gene_aln_df, by = "Gene") %>%
  mutate(
    # Clamp codon index just in case
    idx_before = pmax(1L, CodonIndexBefore),
    aln_col    = map2_int(aln_string, idx_before, aa_col_of_residue)
  ) %>%
  filter(!is.na(aln_col)) %>%
  select(Gene, Intron, CodonIndexBefore, aln_col)

# ----------------------------
# 4) BUILD BINARY PRESENCE / ABSENCE MATRIX
# ----------------------------

# Each unique alignment column is treated as a site
sites <- sort(unique(introns_mapped$aln_col))
genes <- sort(unique(introns_mapped$Gene))

# Tidy long-format table: one row per gene × site
intron_tidy <- introns_mapped %>%
  distinct(Gene, aln_col) %>%   # collapse duplicate mappings if present
  mutate(present = 1L) %>%
  complete(
    Gene    = genes,
    aln_col = sites,
    fill    = list(present = 0L)
  )

# Convert to wide matrix (genes × sites)
mat <- intron_tidy %>%
  mutate(col_name = sprintf("pos_%03d", aln_col)) %>%
  select(Gene, col_name, present) %>%
  pivot_wider(
    names_from  = col_name,
    values_from = present,
    values_fill = 0L
  ) %>%
  arrange(Gene)

# ----------------------------
# 5) SAVE OUTPUTS
# ----------------------------

# Binary presence/absence matrix
write_csv(mat, OUT_MAT)

# Tidy intron position table (useful for plotting/debugging)
intron_tidy_out <- introns_mapped %>%
  arrange(aln_col, Gene) %>%
  mutate(col_name = sprintf("pos_%03d", aln_col)) %>%
  select(Gene, Intron, CodonIndexBefore, aln_col, col_name)

write_csv(intron_tidy_out, OUT_TIDY)

message("Saved intron presence/absence matrix to: ", OUT_MAT)
message("Saved tidy intron position table to: ", OUT_TIDY)
