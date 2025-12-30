suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})

# ============================================================
# Intron positional conservation statistics (permutation tests)
# ------------------------------------------------------------
# Input:
#   - A binary intron presence/absence matrix (CSV)
#   - One row per gene
#   - One column per intron “site” (0/1), plus a "Gene" column
#
# Output:
#   - A CSV summary table with observed statistics, null expectations,
#     effect sizes, and one-sided p-values.
#
# What gets tested:
#   (1) Mean pairwise Jaccard similarity across genes (higher = more similar)
#   (2) Number of sites shared by >= MIN_OCC genes (higher = more shared sites)
#
# Null model:
#   - Keeps each gene’s intron count fixed (row sums)
#   - Randomly redistributes introns across the available sites
# ============================================================

# ----------------------------
# 0) USER SETTINGS
# ----------------------------

# Input intron matrix (CSV)
IN_MAT  <- "results/intron_summaries/intron_matrix.csv"

# Output results (CSV)
OUT_CSV <- "results/intron_summaries/intron_conservation_permtest.csv"

# Number of permutations for the null distribution
NPERM <- 10000

# Random seed for reproducibility
SEED <- 1

# Define a “shared” intron site as present in at least this many genes
MIN_OCC <- 2

dir.create(dirname(OUT_CSV), showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# 1) LOAD MATRIX
# ----------------------------

df <- read_csv(IN_MAT, show_col_types = FALSE)
stopifnot("Gene" %in% names(df))

# Move Gene into rownames so the remaining columns are the matrix
mat_df <- df |>
  as.data.frame() |>
  column_to_rownames("Gene")

# Detect site columns:
# - Accept "pos_###" style columns
# - Also accept purely numeric columns like "906"
cn <- colnames(mat_df)
is_pos <- grepl("^pos_[0-9]+$", cn)
is_num <- grepl("^[0-9]+$", cn)
site_cols <- cn[is_pos | is_num]
stopifnot(length(site_cols) > 0)

# Sort site columns numerically so output is stable
site_index <- ifelse(
  grepl("^pos_", site_cols),
  as.integer(sub("^pos_", "", site_cols)),
  as.integer(site_cols)
)
site_cols <- site_cols[order(site_index)]

# Convert to a numeric 0/1 matrix
mat <- as.matrix(mat_df[, site_cols, drop = FALSE])
mode(mat) <- "numeric"
mat[is.na(mat)] <- 0
mat[mat != 0] <- 1

# Basic sanity checks
G <- nrow(mat)  # number of genes
S <- ncol(mat)  # number of sites
stopifnot(G >= 2, S >= 2)

row_sums <- rowSums(mat)  # introns per gene
stopifnot(sum(row_sums) > 0)

set.seed(SEED)

# ----------------------------
# 2) DEFINE TEST STATISTICS
# ----------------------------

# (A) Mean pairwise Jaccard similarity across genes
# Jaccard(A,B) = |A ∩ B| / |A ∪ B|
mean_pairwise_jaccard <- function(M) {
  inter <- M %*% t(M)                 # intersection sizes for all gene pairs
  rs <- rowSums(M)
  union <- outer(rs, rs, "+") - inter # union sizes for all gene pairs
  
  J <- inter / union
  diag(J) <- NA_real_
  
  # If two genes both have 0 introns, union=0 (undefined); ignore those pairs
  J[union == 0] <- NA_real_
  
  mean(J, na.rm = TRUE)
}

# (B) Count how many sites are shared by at least MIN_OCC genes
shared_sites_count <- function(M, min_occ = 2L) {
  sum(colSums(M) >= min_occ)
}

# Observed statistics from the real matrix
obs_jacc   <- mean_pairwise_jaccard(mat)
obs_shared <- shared_sites_count(mat, min_occ = MIN_OCC)

# ----------------------------
# 3) PERMUTATION NULL MODEL
# ----------------------------
# Build a random matrix that:
#   - keeps each row sum fixed (same intron count per gene)
#   - places those introns uniformly across the S sites

permute_matrix_preserve_row_sums <- function(row_sums, S) {
  G <- length(row_sums)
  M <- matrix(0L, nrow = G, ncol = S)
  
  for (i in seq_len(G)) {
    k <- row_sums[i]
    if (k > 0) {
      idx <- sample.int(S, size = k, replace = FALSE)
      M[i, idx] <- 1L
    }
  }
  
  M
}

# Store permutation results
perm_jacc   <- numeric(NPERM)
perm_shared <- integer(NPERM)

for (b in seq_len(NPERM)) {
  Mp <- permute_matrix_preserve_row_sums(row_sums, S)
  perm_jacc[b] <- mean_pairwise_jaccard(Mp)
  perm_shared[b] <- shared_sites_count(Mp, min_occ = MIN_OCC)
}

# ----------------------------
# 4) P-VALUES + EFFECT SIZES
# ----------------------------

# One-sided p-values: "greater" = more conservation than expected by chance
p_jacc   <- (sum(perm_jacc >= obs_jacc) + 1) / (NPERM + 1)
p_shared <- (sum(perm_shared >= obs_shared) + 1) / (NPERM + 1)

# Simple z-like effect sizes relative to the null
z_jacc   <- (obs_jacc - mean(perm_jacc)) / sd(perm_jacc)
z_shared <- (obs_shared - mean(perm_shared)) / sd(as.numeric(perm_shared))

# ----------------------------
# 5) SAVE SUMMARY TABLE
# ----------------------------

summ <- tibble(
  matrix_file   = IN_MAT,
  n_genes       = G,
  n_sites       = S,
  total_introns = sum(row_sums),
  nperm         = NPERM,
  seed          = SEED,
  min_occ       = MIN_OCC,
  
  stat = c(
    "mean_pairwise_jaccard",
    paste0("n_sites_with_occupancy_ge_", MIN_OCC)
  ),
  
  observed  = c(obs_jacc, obs_shared),
  null_mean = c(mean(perm_jacc), mean(perm_shared)),
  null_sd   = c(sd(perm_jacc), sd(as.numeric(perm_shared))),
  effect_z  = c(z_jacc, z_shared),
  
  p_value_greater = c(p_jacc, p_shared)
)

write_csv(summ, OUT_CSV)

print(summ)
message("Saved permutation test summary to: ", OUT_CSV)
