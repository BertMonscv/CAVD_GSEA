#!/usr/bin/env Rscript
# =============================================================================
# verify_gene_set_overlap.R
#
# Pre-analysis plan Section 4.1 #3 compliance check.
# Verifies pairwise Jaccard overlap among curated gene sets is <= 0.30.
#
# Scope:
#   PRIMARY (plan literal): 4 custom sets listed in plan Section 4 table
#     - CUSTOM_LEPTIN_MTOR_SIGNALING
#     - CUSTOM_MTOR_AUTOPHAGY_LIPOPHAGY_AXIS
#     - CUSTOM_LIPOPHAGY_CORE
#     - CUSTOM_VIC_OSTEOGENIC_DIFFERENTIATION
#     -> 6 pairs
#
#   SUPPLEMENTARY (completeness): also include CUSTOM_S100A9_RAGE_CALCIFICATION
#     (5-set pool member but listed as supplementary in plan)
#     -> 4 additional pairs (10 total)
#
# Input:
#   Supp_Table_S1_v2.tsv  (master gene-set table, one row per gene)
#     required columns: gene_symbol, gene_set
#
# Output:
#   gene_set_overlap_log.tsv  (one row per pair, both scopes)
#
# Exit status:
#   0 if every pair <= 0.30
#   1 if any pair fails (downstream fgsea call must not proceed)
#
# Usage:
#   Rscript verify_gene_set_overlap.R \
#       --input  Supp_Table_S1_v2.tsv \
#       --output gene_set_overlap_log.tsv
#
# Dependencies: base R only (no external packages required).
# =============================================================================

# Dependencies: base R only -- no library() calls needed.

# ---- args --------------------------------------------------------------------
parse_args <- function(argv) {
  defaults <- list(
    input  = "Supp_Table_S1_v2.tsv",
    output = "gene_set_overlap_log.tsv",
    threshold = 0.30
  )
  i <- 1
  while (i <= length(argv)) {
    a <- argv[i]
    if (a == "--input")     { defaults$input  <- argv[i + 1]; i <- i + 2; next }
    if (a == "--output")    { defaults$output <- argv[i + 1]; i <- i + 2; next }
    if (a == "--threshold") { defaults$threshold <- as.numeric(argv[i + 1]); i <- i + 2; next }
    if (a %in% c("-h", "--help")) {
      cat("Usage: Rscript verify_gene_set_overlap.R [--input FILE] [--output FILE] [--threshold N]\n")
      quit(save = "no", status = 0)
    }
    stop("Unknown argument: ", a)
  }
  defaults
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

cat("verify_gene_set_overlap.R\n")
cat("  input:     ", args$input, "\n")
cat("  output:    ", args$output, "\n")
cat("  threshold: ", args$threshold, "\n\n")

# ---- read --------------------------------------------------------------------
if (!file.exists(args$input)) {
  stop("Input file not found: ", args$input)
}

df <- read.table(
  args$input,
  sep = "\t", header = TRUE,
  stringsAsFactors = FALSE,
  quote = "", comment.char = "",
  encoding = "UTF-8"
)

required_cols <- c("gene_symbol", "gene_set")
missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) {
  stop("Input missing required columns: ", paste(missing_cols, collapse = ", "))
}

# ---- scope definitions -------------------------------------------------------
PLAN_FOUR <- c(
  "CUSTOM_LEPTIN_MTOR_SIGNALING",
  "CUSTOM_MTOR_AUTOPHAGY_LIPOPHAGY_AXIS",
  "CUSTOM_LIPOPHAGY_CORE",
  "CUSTOM_VIC_OSTEOGENIC_DIFFERENTIATION"
)
SUPPLEMENTARY_FIFTH <- "CUSTOM_S100A9_RAGE_CALCIFICATION"
ALL_FIVE <- c(PLAN_FOUR, SUPPLEMENTARY_FIFTH)

# sanity: every expected set must be present
present_sets <- unique(df$gene_set)
missing_sets <- setdiff(ALL_FIVE, present_sets)
if (length(missing_sets) > 0) {
  stop("Input is missing gene_set values: ", paste(missing_sets, collapse = ", "))
}

# build gene_set -> unique symbol set
gene_sets <- lapply(ALL_FIVE, function(s) {
  unique(df$gene_symbol[df$gene_set == s])
})
names(gene_sets) <- ALL_FIVE

cat("Set sizes:\n")
for (s in ALL_FIVE) cat("  ", s, ": ", length(gene_sets[[s]]), "\n", sep = "")
cat("\n")

# ---- jaccard helper ----------------------------------------------------------
jaccard_row <- function(a_name, b_name, scope_label, threshold) {
  a <- gene_sets[[a_name]]
  b <- gene_sets[[b_name]]
  inter <- intersect(a, b)
  uni   <- union(a, b)
  jac   <- if (length(uni) == 0) 0 else length(inter) / length(uni)
  data.frame(
    scope        = scope_label,
    set_A        = a_name,
    set_B        = b_name,
    n_A          = length(a),
    n_B          = length(b),
    n_intersect  = length(inter),
    n_union      = length(uni),
    jaccard      = sprintf("%.4f", jac),
    threshold    = sprintf("%.2f", threshold),
    passes       = if (jac <= threshold) "PASS" else "FAIL",
    shared_genes = paste(sort(inter), collapse = ","),
    stringsAsFactors = FALSE
  )
}

# ---- compute pairs -----------------------------------------------------------
compute_pairs <- function(names, scope_label, threshold) {
  combs <- utils::combn(names, 2, simplify = FALSE)
  do.call(rbind, lapply(combs, function(p) {
    jaccard_row(p[1], p[2], scope_label, threshold)
  }))
}

primary <- compute_pairs(PLAN_FOUR, "plan_section_4.1_four_custom", args$threshold)
all_pairs <- compute_pairs(ALL_FIVE, "tmp", args$threshold)

# keep only the supplementary-extra pairs (those not already in primary)
primary_key <- paste(pmin(primary$set_A, primary$set_B),
                     pmax(primary$set_A, primary$set_B), sep = "||")
all_key     <- paste(pmin(all_pairs$set_A, all_pairs$set_B),
                     pmax(all_pairs$set_A, all_pairs$set_B), sep = "||")
supp <- all_pairs[!(all_key %in% primary_key), , drop = FALSE]
supp$scope <- "all_five_custom_including_S100A9"

# ---- report ------------------------------------------------------------------
print_table <- function(rows, label) {
  cat(strrep("=", 88), "\n", sep = "")
  cat(label, "\n", sep = "")
  cat(strrep("=", 88), "\n", sep = "")
  for (i in seq_len(nrow(rows))) {
    r <- rows[i, ]
    cat(sprintf("  %-40s  vs  %-40s  jac=%s  %s\n",
                r$set_A, r$set_B, r$jaccard, r$passes))
  }
  cat("\n")
}

print_table(primary,
            sprintf("PRIMARY: plan Section 4.1 - %d pairs, threshold %.2f",
                    nrow(primary), args$threshold))
print_table(supp,
            sprintf("SUPPLEMENTARY: %d additional S100A9 pairs", nrow(supp)))

all_results <- rbind(primary, supp)
n_fail <- sum(all_results$passes == "FAIL")

primary_pass <- all(primary$passes == "PASS")
all_pass     <- all(all_results$passes == "PASS")

cat(strrep("=", 88), "\n", sep = "")
cat("VERDICT (primary 4-custom scope, 6 pairs): ",
    if (primary_pass) "PASS" else "FAIL", "\n", sep = "")
cat("VERDICT (full 5-custom scope, 10 pairs):   ",
    if (all_pass) "PASS" else "FAIL", "\n", sep = "")
cat("Max Jaccard observed: ",
    sprintf("%.4f", max(as.numeric(all_results$jaccard))), "\n", sep = "")
cat(strrep("=", 88), "\n\n", sep = "")

# ---- write -------------------------------------------------------------------
write.table(
  all_results,
  file = args$output,
  sep = "\t", quote = FALSE,
  row.names = FALSE, col.names = TRUE,
  fileEncoding = "UTF-8",
  eol = "\n"
)
cat("Wrote: ", args$output, "  (", nrow(all_results), " pair rows)\n", sep = "")

# ---- exit --------------------------------------------------------------------
# Primary scope (4 custom) determines compliance with the pre-analysis plan literal.
# If primary fails, downstream fgsea must not proceed (per plan Section 4.1 #3).
if (!primary_pass) {
  cat("\nFAIL: at least one primary-scope pair exceeds threshold. ",
      "fgsea must NOT proceed until gene sets are revised (plan Section 4.1 #3).\n",
      sep = "")
  quit(save = "no", status = 1)
}
quit(save = "no", status = 0)
