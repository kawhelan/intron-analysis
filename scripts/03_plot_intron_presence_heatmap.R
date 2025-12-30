#!/usr/bin/env Rscript
# ============================================================
# Intron presence/absence heatmap (ComplexHeatmap)
# ------------------------------------------------------------
# What this script does:
#   1) Reads a CSV intron presence/absence matrix
#        - First column must be "Gene"
#        - Remaining columns are sites (0/1 values)
#   2) Detects which columns are "site columns" (by default: columns that contain digits)
#   3) Converts values to 0/1 and drops site columns that are all zeros (optional)
#   4) Plots a heatmap showing intron presence (1) vs absence (0)
#   5) Saves the plot as a PNG
#
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# ----------------------------
# 0) USER SETTINGS
# ----------------------------

# Input matrix (CSV):
#   - must contain a column named "Gene"
#   - all other columns should be intron site columns with 0/1 (or numeric) values
IN_CSV  <- "results/intron_summaries/intron_presence_matrix.csv"

# Output figure (PNG)
OUT_PNG <- "results/intron_summaries/intron_presence_heatmap.png"

# If TRUE, cluster rows by similarity
CLUSTER_ROWS <- FALSE

# If TRUE, remove site columns where every value is 0
DROP_EMPTY_COLUMNS <- TRUE

# Optional: cap how many column labels to print (to avoid unreadable axes)
MAX_COLUMN_LABELS <- 40

# ----------------------------
# 1) LOAD MATRIX
# ----------------------------

df <- read_csv(IN_CSV, show_col_types = FALSE)

# Make sure the Gene column exists
stopifnot("Gene" %in% names(df))

# Convert to a matrix with Gene as row names
mat0 <- df %>%
  as.data.frame() %>%
  column_to_rownames("Gene")

# ----------------------------
# 2) IDENTIFY SITE COLUMNS
# ----------------------------
# This tries to be robust to site column names like:
#   pos_027, X906, 906.1, site105, etc.
# Rule: a "site column" is any column name that contains at least one digit.

all_cols <- colnames(mat0)

extract_pos <- function(x) {
  m <- regexpr("[0-9]+", x)
  ifelse(m == -1, NA_integer_, as.integer(regmatches(x, m)))
}

pos <- extract_pos(all_cols)

site_cols <- all_cols[!is.na(pos)]
site_pos  <- pos[!is.na(pos)]

if (length(site_cols) == 0) {
  stop(
    "No site columns detected. Your columns must contain digits (e.g., pos_001).\n",
    "Columns found:\n  ",
    paste(head(all_cols, 30), collapse = ", ")
  )
}

# Sort site columns by their numeric position so the heatmap runs left-to-right
ord <- order(site_pos)
site_cols <- site_cols[ord]
site_pos  <- site_pos[ord]

# ----------------------------
# 3) CLEAN VALUES TO 0/1
# ----------------------------

# Pull out only site columns and force numeric
mat <- as.matrix(mat0[, site_cols, drop = FALSE])
mat <- suppressWarnings(matrix(
  as.numeric(mat),
  nrow = nrow(mat),
  ncol = ncol(mat),
  dimnames = dimnames(mat)
))

# Replace NA with 0, then force everything to 0/1
mat[is.na(mat)] <- 0
mat <- ifelse(mat != 0, 1, 0)

# Optionally drop columns that have no introns in any gene
if (DROP_EMPTY_COLUMNS) {
  keep_cols <- apply(mat, 2, function(x) any(x == 1))
  mat <- mat[, keep_cols, drop = FALSE]
  site_pos <- site_pos[keep_cols]
}

if (ncol(mat) == 0) {
  stop("No site columns left after filtering. Check that your matrix contains 1s.")
}

# ----------------------------
# 4) PREP LABELS (OPTIONAL)
# ----------------------------
# If you have hundreds of site columns, labeling every column becomes unreadable.
# This prints only ~MAX_COLUMN_LABELS labels spread across the heatmap.

lab <- as.character(site_pos)
step <- if (length(lab) > MAX_COLUMN_LABELS) ceiling(length(lab) / MAX_COLUMN_LABELS) else 1
lab[((seq_along(lab) - 1) %% step) != 0] <- ""

# ----------------------------
# 5) BUILD A SIMPLE HEATMAP
# ----------------------------
# White = absent (0)
# Dark = present (1)

mat_chr <- ifelse(mat == 1, "1", "0")

col_fun <- c(
  "0" = "white",
  "1" = "#1a9850"  # one color for "present" (edit if you want)
)

# Optional: add a small bar showing the number of introns per gene
row_counts <- rowSums(mat, na.rm = TRUE)

row_anno <- rowAnnotation(
  `Intron count` = row_counts,
  col = list(`Intron count` = colorRamp2(
    c(min(row_counts), max(row_counts)),
    c("#f0f0f0", "#636363")
  )),
  gp = gpar(col = NA),
  width = unit(6, "mm"),
  show_annotation_name = FALSE
)

ht <- Heatmap(
  mat_chr,
  name = "Intron",
  col = col_fun,
  
  cluster_rows = CLUSTER_ROWS,
  cluster_columns = FALSE,
  
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 10),
  
  show_column_names = TRUE,
  column_labels = lab,
  column_names_gp = gpar(fontsize = 8),
  column_names_rot = 90,
  
  left_annotation = row_anno,
  
  # Light grid so cell boundaries are visible
  rect_gp = gpar(col = "#808080", lwd = 0.4),
  
  row_title = NULL,
  use_raster = TRUE
)

# ----------------------------
# 6) SAVE FIGURE
# ----------------------------

dir.create(dirname(OUT_PNG), showWarnings = FALSE, recursive = TRUE)

png(OUT_PNG, width = 2200, height = 1400, res = 200)
draw(ht)
dev.off()

message("Saved heatmap to: ", OUT_PNG)
