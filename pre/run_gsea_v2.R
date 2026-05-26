# =============================================================================
# run_gsea_v2.R
#
# Frozen analysis script for CAVD GSEA Pre-Analysis Plan v2.
# (Plan: pre/CAVD_GSEA_PreAnalysisPlan_v2.md, frozen 2026-05-26,
#  git commit c76f64489e94180e77741d8cbf0aca9a37134e9c, tag prereg-frozen-v2)
#
# This script is committed BEFORE any fgsea call is executed against any
# cohort (per plan Section 10.1 #4). No edits permitted after the freeze
# commit except via a dated amendment file.
#
# -----------------------------------------------------------------------------
# Cohort coverage and sample selection (per plan Section 2.1):
#
# Cohort 1: GSE51472 (discovery, primary ranking; 5 / 5)
#   GSE51472 contains 3 groups of 5 samples each:
#     Normal:    GSM1246204-08
#     Sclerotic: GSM1246209-13   <-- EXCLUDED from Normal-vs-Calcified contrast
#     Calcified: GSM1246214-18   (a.k.a. "stenotic" in source metadata)
#   Plan §2.1 table specifies "5 / 5 Normal/Calcified"; the Sclerotic group
#   is excluded by this sample count. No plan amendment required.
#
# Cohort 2: GSE51472 + GSE12644 integrated (sensitivity; 15 / 15 post-ComBat)
#   GSE12644 contributes 10 Normal + 10 Calcified (all tricuspid aortic
#   valves, GPL570, RMA-processed). Combined with GSE51472's 5 + 5 gives
#   15 / 15. ComBat is applied to the combined expression matrix with
#   calcification status preserved as the biological variable (plan §2.2).
#
# Cohort 3: GSE83453 (cross-platform validation; 8 / 9)
#   GSE83453 contains 3 groups:
#      8 TAVn  (tricuspid, no stenosis)        --> Normal
#      9 TAVc  (tricuspid, with stenosis)      --> Calcified
#     10 BAV   (bicuspid, calcified)           <-- EXCLUDED
#   Rationale for BAV exclusion: BAV is an anatomical variant whose
#   calcification mechanism differs from tricuspid calcification
#   (Mathieu et al.). Including BAV would confound calcification effect
#   with anatomical type, and would make GSE83453 incomparable to GSE51472
#   and GSE12644 (both tricuspid-only). Plan §2.1's "8/9" reflects this
#   design; no plan amendment required.
#   GSE83453 platform = GPL10558 (Illumina HumanHT-12 v4), uses
#   illuminaHumanv4.db annotation, data already log2-transformed (lumi).
# =============================================================================


# ---- Stage 1: environment + packages ----------------------------------------

.required_packages <- list(
  fgsea       = "1.28.0",
  limma       = "3.62.0",
  sva         = "3.50.0",
  msigdbr     = "7.5.1",
  HGNChelper  = NA,
  digest      = NA
)

for (pkg in names(.required_packages)) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Required package '%s' is not installed.", pkg))
  }
  required_ver <- .required_packages[[pkg]]
  if (!is.na(required_ver)) {
    actual_ver <- as.character(utils::packageVersion(pkg))
    if (utils::compareVersion(actual_ver, required_ver) < 0) {
      stop(sprintf("Package '%s' version %s is below required %s.",
                   pkg, actual_ver, required_ver))
    }
  }
}

# Strict pin on msigdbr at exactly 7.5.1 (plan §5.2)
msigdbr_ver <- as.character(utils::packageVersion("msigdbr"))
if (msigdbr_ver != "7.5.1") {
  stop(sprintf("msigdbr must be exactly 7.5.1 per plan §5.2 (got %s)", msigdbr_ver))
}

# Annotation packages (Bioconductor)
for (annot_pkg in c("BiocManager", "AnnotationDbi",
                    "hgu133plus2.db", "illuminaHumanv4.db",
                    "BiocParallel")) {
  if (!requireNamespace(annot_pkg, quietly = TRUE)) {
    stop(sprintf("Required annotation/Bioconductor package '%s' is not installed.",
                 annot_pkg))
  }
}

suppressPackageStartupMessages({
  library(fgsea)
  library(limma)
  library(sva)
  library(msigdbr)
  library(HGNChelper)
  library(AnnotationDbi)
  library(hgu133plus2.db)
  library(illuminaHumanv4.db)
})

.PLAN_SEED <- 42L

# Toggle: TRUE = run 100-bootstrap NES_SD per pathway per cohort (final, plan §7.1)
#         FALSE = skip bootstrap, NES_SD = NA (dry-run)
RUN_BOOTSTRAP <- TRUE

# Paths (working dir = repo root)
.PATH_GSE51472    <- "csv/GSE51472_series_matrix.txt"
.PATH_GSE12644    <- "csv/GSE12644_series_matrix.txt"
.PATH_GSE83453    <- "csv/GSE83453_series_matrix.txt"
.PATH_GENESETS    <- "pre/Supp_Table_S1_v2.tsv"
.PATH_OUT_DIR     <- "results"
.PATH_SESSIONINFO <- "pre/sessionInfo_v2.txt"

dir.create(.PATH_OUT_DIR, showWarnings = FALSE, recursive = TRUE)
writeLines(capture.output(sessionInfo()), .PATH_SESSIONINFO)
cat(sprintf("[stage 1] sessionInfo written to %s\n", .PATH_SESSIONINFO))
cat(sprintf("[stage 1] RUN_BOOTSTRAP = %s\n", RUN_BOOTSTRAP))


# ---- Stage 2: series matrix loader (generic) --------------------------------

# Reads a GEO series matrix file and returns:
#   $expr     : numeric matrix, rows=probes, cols=GSM IDs
#   $samples  : data.frame with sample_id, all !Sample_characteristics_ch1
#                rows, and !Sample_source_name_ch1
read_series_matrix <- function(path) {
  lines <- readLines(path)
  begin <- grep("^!series_matrix_table_begin", lines)
  end   <- grep("^!series_matrix_table_end",   lines)
  stopifnot(length(begin) == 1L, length(end) == 1L, end > begin)

  geo_line <- lines[grep("^!Sample_geo_accession", lines)]
  stopifnot(length(geo_line) == 1L)
  geo_ids <- gsub("\"", "", strsplit(geo_line, "\t")[[1]][-1])

  src_line <- lines[grep("^!Sample_source_name_ch1", lines)]
  if (length(src_line) == 1L) {
    src <- gsub("\"", "", strsplit(src_line, "\t")[[1]][-1])
  } else {
    src <- rep(NA_character_, length(geo_ids))
  }

  # All !Sample_characteristics_ch1 rows; stack them
  char_lines <- lines[grep("^!Sample_characteristics_ch1", lines)]
  char_mat <- if (length(char_lines) > 0L) {
    do.call(rbind, lapply(char_lines, function(L) {
      vals <- gsub("\"", "", strsplit(L, "\t")[[1]][-1])
      vals
    }))
  } else matrix("", nrow = 0, ncol = length(geo_ids))

  title_line <- lines[grep("^!Sample_title", lines)]
  titles <- if (length(title_line) == 1L) {
    gsub("\"", "", strsplit(title_line, "\t")[[1]][-1])
  } else rep(NA_character_, length(geo_ids))

  samples <- data.frame(
    sample_id = geo_ids,
    title = titles,
    source_name = src,
    stringsAsFactors = FALSE
  )
  if (nrow(char_mat) > 0L) {
    char_df <- as.data.frame(t(char_mat), stringsAsFactors = FALSE)
    colnames(char_df) <- paste0("char_", seq_len(ncol(char_df)))
    samples <- cbind(samples, char_df)
  }

  expr_lines <- lines[(begin + 1L):(end - 1L)]
  con <- textConnection(expr_lines)
  on.exit(close(con), add = TRUE)
  expr <- read.table(con, header = TRUE, sep = "\t", row.names = 1,
                     check.names = FALSE, quote = "\"",
                     stringsAsFactors = FALSE, comment.char = "")
  expr <- as.matrix(expr)

  stopifnot(identical(colnames(expr), geo_ids))
  list(expr = expr, samples = samples)
}

# Log2-transform if data appear linear (plan §2.2)
log2_if_needed <- function(expr, label) {
  qs <- quantile(as.vector(expr), probs = c(0, 0.5, 0.99, 1.0), na.rm = TRUE)
  needs <- unname(qs["99%"]) >= 100
  cat(sprintf("  [%s] quantiles 0/50/99/100%%: %.2f / %.2f / %.2f / %.2f\n",
              label, qs[1], qs[2], qs[3], qs[4]))
  if (needs) {
    cat(sprintf("  [%s] LINEAR scale detected --> log2(x + 1) applied\n", label))
    expr <- log2(expr + 1)
  } else {
    cat(sprintf("  [%s] already log2 scale, no transform\n", label))
  }
  expr
}


# ---- Stage 2a: load GSE51472 -----------------------------------------------

cat("\n[stage 2a] loading GSE51472...\n")
gse51472 <- read_series_matrix(.PATH_GSE51472)

# Find the "disease state" row in characteristics
state_col <- which(apply(gse51472$samples, 2, function(v)
  any(grepl("disease state:", v, fixed = TRUE))))
stopifnot(length(state_col) == 1L)
gse51472$samples$disease_state <- sub("^disease state:\\s*", "",
                                       gse51472$samples[[state_col]])

keep_normal    <- gse51472$samples$disease_state == "Normal aortic valve"
keep_calcified <- gse51472$samples$disease_state == "Calcified aortic valve"
keep_idx <- keep_normal | keep_calcified
stopifnot(sum(keep_normal) == 5L, sum(keep_calcified) == 5L)
cat(sprintf("  excluded (Sclerotic): %s\n",
            paste(gse51472$samples$sample_id[!keep_idx], collapse = ", ")))

expr_51472 <- gse51472$expr[, keep_idx, drop = FALSE]
samp_51472 <- gse51472$samples[keep_idx, , drop = FALSE]
samp_51472$group <- factor(
  ifelse(samp_51472$disease_state == "Normal aortic valve", "Normal", "Calcified"),
  levels = c("Normal", "Calcified"))
samp_51472$batch <- "GSE51472"
expr_51472 <- log2_if_needed(expr_51472, "GSE51472")


# ---- Stage 2b: load GSE12644 ------------------------------------------------

cat("\n[stage 2b] loading GSE12644...\n")
gse12644 <- read_series_matrix(.PATH_GSE12644)

# GSE12644 uses source_name_ch1 = "Normal aortic valve" / "Calcified aortic valve"
gse12644$samples$disease_state <- gse12644$samples$source_name

keep_normal_12 <- gse12644$samples$disease_state == "Normal aortic valve"
keep_calc_12   <- gse12644$samples$disease_state == "Calcified aortic valve"
stopifnot(sum(keep_normal_12) == 10L, sum(keep_calc_12) == 10L)
keep_idx_12 <- keep_normal_12 | keep_calc_12

expr_12644 <- gse12644$expr[, keep_idx_12, drop = FALSE]
samp_12644 <- gse12644$samples[keep_idx_12, , drop = FALSE]
samp_12644$group <- factor(
  ifelse(samp_12644$disease_state == "Normal aortic valve", "Normal", "Calcified"),
  levels = c("Normal", "Calcified"))
samp_12644$batch <- "GSE12644"
expr_12644 <- log2_if_needed(expr_12644, "GSE12644")


# ---- Stage 2c: load GSE83453 ------------------------------------------------

cat("\n[stage 2c] loading GSE83453...\n")
gse83453 <- read_series_matrix(.PATH_GSE83453)

# GSE83453: classify via source_name_ch1
gse83453$samples$disease_state <- gse83453$samples$source_name

# source_name values: "human normal aortic valve", "human stenotic aortic valve",
# "human bicuspid aortic valve" (BAV; excluded)
is_normal <- grepl("normal", gse83453$samples$source_name, ignore.case = TRUE)
is_TAVc   <- grepl("stenotic", gse83453$samples$source_name, ignore.case = TRUE)
is_BAV    <- grepl("bicuspid", gse83453$samples$source_name, ignore.case = TRUE)

stopifnot(sum(is_normal) == 8L)
stopifnot(sum(is_TAVc)   == 9L)
stopifnot(sum(is_BAV)    == 10L)

cat(sprintf("  excluded (BAV): %d samples (%s)\n",
            sum(is_BAV),
            paste(head(gse83453$samples$sample_id[is_BAV], 3), "...",
                  collapse = ", ")))

keep_idx_83 <- is_normal | is_TAVc
expr_83453 <- gse83453$expr[, keep_idx_83, drop = FALSE]
samp_83453 <- gse83453$samples[keep_idx_83, , drop = FALSE]
samp_83453$group <- factor(
  ifelse(is_normal[keep_idx_83], "Normal", "Calcified"),
  levels = c("Normal", "Calcified"))
samp_83453$batch <- "GSE83453"
expr_83453 <- log2_if_needed(expr_83453, "GSE83453")


# ---- Stage 3: probe -> gene with IQR collapse (per platform) ----------------

# Generic IQR collapse function
collapse_probes_by_iqr <- function(expr, probe2sym, label) {
  # probe2sym: data.frame with columns PROBEID, SYMBOL
  cat(sprintf("  [%s] probes in input: %d\n", label, nrow(expr)))

  p2s <- probe2sym[!is.na(probe2sym$SYMBOL), ]
  multi <- unique(p2s$PROBEID[duplicated(p2s$PROBEID)])
  p2s <- p2s[!p2s$PROBEID %in% multi, ]
  p2s <- p2s[p2s$PROBEID %in% rownames(expr), ]

  cat(sprintf("  [%s] probes with unique SYMBOL annotation: %d\n",
              label, nrow(p2s)))
  cat(sprintf("  [%s] probes dropped (multi-SYMBOL):         %d\n",
              label, length(multi)))

  expr_annot <- expr[p2s$PROBEID, , drop = FALSE]
  probe_iqr <- apply(expr_annot, 1, IQR)

  tbl <- data.frame(
    probe  = p2s$PROBEID,
    symbol = p2s$SYMBOL,
    iqr    = probe_iqr,
    stringsAsFactors = FALSE
  )
  tbl <- tbl[order(tbl$symbol, -tbl$iqr), ]
  collapsed <- tbl[!duplicated(tbl$symbol), ]

  expr_gene <- expr_annot[collapsed$probe, , drop = FALSE]
  rownames(expr_gene) <- collapsed$symbol

  cat(sprintf("  [%s] unique gene symbols after IQR collapse: %d\n",
              label, nrow(expr_gene)))
  expr_gene
}

# Affymetrix GPL570 probe -> gene (GSE51472, GSE12644)
cat("\n[stage 3] probe -> gene collapse\n")
probes_affy <- unique(c(rownames(expr_51472), rownames(expr_12644)))
p2s_affy <- AnnotationDbi::select(
  hgu133plus2.db,
  keys = probes_affy,
  columns = c("SYMBOL"),
  keytype = "PROBEID"
)

expr_gene_51472 <- collapse_probes_by_iqr(expr_51472, p2s_affy, "GSE51472")
expr_gene_12644 <- collapse_probes_by_iqr(expr_12644, p2s_affy, "GSE12644")

# Illumina GPL10558 probe -> gene (GSE83453)
probes_illumina <- rownames(expr_83453)
p2s_illumina <- AnnotationDbi::select(
  illuminaHumanv4.db,
  keys = probes_illumina,
  columns = c("SYMBOL"),
  keytype = "PROBEID"
)
expr_gene_83453 <- collapse_probes_by_iqr(expr_83453, p2s_illumina, "GSE83453")


# ---- Stage 4a: limma DE per cohort -----------------------------------------

run_limma_DE <- function(expr_gene, samples, label, seed) {
  cat(sprintf("\n[stage 4 - %s] limma DE...\n", label))
  design <- model.matrix(~ 0 + samples$group)
  colnames(design) <- levels(samples$group)
  fit <- lmFit(expr_gene, design)
  cont <- makeContrasts(Calcified_vs_Normal = Calcified - Normal, levels = design)
  fit2 <- contrasts.fit(fit, cont)
  fit2 <- eBayes(fit2)
  tt <- topTable(fit2, coef = "Calcified_vs_Normal",
                 number = Inf, sort.by = "none")
  cat(sprintf("  [%s] %d genes, logFC range [%.3f, %.3f], median %.3f\n",
              label, nrow(tt), min(tt$logFC), max(tt$logFC), median(tt$logFC)))

  # Stable jitter for tie-breaking (plan §3.2)
  set.seed(seed)
  jitter_vec <- runif(nrow(tt), min = -1e-10, max = 1e-10)
  ranks <- tt$logFC + jitter_vec
  names(ranks) <- rownames(tt)
  ranks <- sort(ranks, decreasing = TRUE)
  cat(sprintf("  [%s] ranked vector head: %s = %.3f, tail: %s = %.3f\n",
              label, names(ranks)[1], unname(ranks[1]),
              names(ranks)[length(ranks)], unname(ranks[length(ranks)])))
  ranks
}

ranks_GSE51472 <- run_limma_DE(expr_gene_51472, samp_51472, "GSE51472",
                                seed = .PLAN_SEED)


# ---- Stage 4b: integrated cohort GSE51472 + GSE12644 with ComBat -----------

cat("\n[stage 4b] integrating GSE51472 + GSE12644 with ComBat...\n")
common_genes <- intersect(rownames(expr_gene_51472), rownames(expr_gene_12644))
cat(sprintf("  common genes: %d\n", length(common_genes)))
expr_comb <- cbind(expr_gene_51472[common_genes, , drop = FALSE],
                   expr_gene_12644[common_genes, , drop = FALSE])
samp_comb <- rbind(
  data.frame(sample_id = samp_51472$sample_id, group = samp_51472$group,
             batch = samp_51472$batch, stringsAsFactors = FALSE),
  data.frame(sample_id = samp_12644$sample_id, group = samp_12644$group,
             batch = samp_12644$batch, stringsAsFactors = FALSE)
)
stopifnot(identical(colnames(expr_comb), samp_comb$sample_id))
cat(sprintf("  combined matrix: %d genes x %d samples (Normal=%d, Calcified=%d)\n",
            nrow(expr_comb), ncol(expr_comb),
            sum(samp_comb$group == "Normal"),
            sum(samp_comb$group == "Calcified")))

# PCA before ComBat
pca_pre <- prcomp(t(expr_comb), scale. = FALSE)

# Apply ComBat: batch = study, mod = ~ group (preserve biological signal)
mod_combat <- model.matrix(~ samp_comb$group)
expr_combat <- ComBat(dat = expr_comb,
                      batch = samp_comb$batch,
                      mod = mod_combat,
                      par.prior = TRUE, prior.plots = FALSE)
cat("  ComBat applied (preserving group as biological variable)\n")

# PCA after ComBat
pca_post <- prcomp(t(expr_combat), scale. = FALSE)

# Save PCA check plot (plan §10.1 item 7)
.PATH_COMBAT_PCA <- file.path(.PATH_OUT_DIR, "combat_pca_check.pdf")
pdf(.PATH_COMBAT_PCA, width = 10, height = 5)
par(mfrow = c(1, 2))
plot_pca <- function(pca, title) {
  v <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  cols <- ifelse(samp_comb$batch == "GSE51472", "red", "blue")
  shapes <- ifelse(samp_comb$group == "Normal", 1, 19)
  plot(pca$x[, 1], pca$x[, 2],
       col = cols, pch = shapes, cex = 1.5,
       xlab = sprintf("PC1 (%.1f%%)", v[1]),
       ylab = sprintf("PC2 (%.1f%%)", v[2]),
       main = title)
  legend("topright", legend = c("GSE51472", "GSE12644", "Normal", "Calcified"),
         col = c("red", "blue", "black", "black"),
         pch = c(15, 15, 1, 19), cex = 0.7)
}
plot_pca(pca_pre,  "Before ComBat")
plot_pca(pca_post, "After ComBat")
dev.off()
cat(sprintf("  PCA check saved: %s\n", .PATH_COMBAT_PCA))

# limma DE on the integrated cohort (use ranks_INT name)
ranks_INTEGRATED <- run_limma_DE(expr_combat, samp_comb, "INTEGRATED",
                                  seed = .PLAN_SEED)


# ---- Stage 4c: validation cohort GSE83453 ----------------------------------

ranks_GSE83453 <- run_limma_DE(expr_gene_83453, samp_83453, "GSE83453",
                                seed = .PLAN_SEED)


# ---- Stage 5: gene sets + coverage check + fgsea ---------------------------

cat("\n[stage 5] loading gene sets, coverage check, fgsea...\n")

# 5a: Hallmark from msigdbr 7.5.1, with SHA256 anchor
hallmark_df <- msigdbr(species = "Homo sapiens", category = "H")
stopifnot(length(unique(hallmark_df$gs_name)) == 50L)
hallmark_sets <- split(hallmark_df$human_gene_symbol, hallmark_df$gs_name)
hallmark_sets <- lapply(hallmark_sets, function(x) sort(unique(x)))

.MTORC1_SHA256_EXPECTED <- "af61682e7fa9e484cfb294f25742a0c8afe0d8f1cdc6d03bfe9f9da45aba7e70"
mtorc1_genes <- hallmark_sets[["HALLMARK_MTORC1_SIGNALING"]]
stopifnot(length(mtorc1_genes) == 200L)
mtorc1_sha <- digest::digest(mtorc1_genes, algo = "sha256")
if (mtorc1_sha != .MTORC1_SHA256_EXPECTED) {
  stop(sprintf("MTORC1 SHA256 mismatch: expected %s, got %s",
               .MTORC1_SHA256_EXPECTED, mtorc1_sha))
}
cat(sprintf("  Hallmark MTORC1 SHA256 anchor: OK\n"))

# 5b: custom sets
gs_df <- read.table(.PATH_GENESETS, sep = "\t", header = TRUE,
                    stringsAsFactors = FALSE, quote = "",
                    comment.char = "", encoding = "UTF-8")
custom_sets <- split(gs_df$gene_symbol, gs_df$gene_set)
custom_sets <- lapply(custom_sets, function(x) sort(unique(trimws(x))))

.EXPECTED_SIZES <- list(
  CUSTOM_LEPTIN_MTOR_SIGNALING            = 22L,
  CUSTOM_MTOR_AUTOPHAGY_LIPOPHAGY_AXIS    = 20L,
  CUSTOM_LIPOPHAGY_CORE                   = 22L,
  CUSTOM_VIC_OSTEOGENIC_DIFFERENTIATION   = 27L,
  CUSTOM_S100A9_RAGE_CALCIFICATION        = 17L
)
for (nm in names(.EXPECTED_SIZES)) {
  if (length(custom_sets[[nm]]) != .EXPECTED_SIZES[[nm]]) {
    stop(sprintf("Set '%s' size mismatch: expected %d got %d",
                 nm, .EXPECTED_SIZES[[nm]], length(custom_sets[[nm]])))
  }
}

# 5c: pools
hypothesis_pool_names <- c(
  "HALLMARK_MTORC1_SIGNALING",
  "CUSTOM_LEPTIN_MTOR_SIGNALING",
  "CUSTOM_MTOR_AUTOPHAGY_LIPOPHAGY_AXIS",
  "CUSTOM_LIPOPHAGY_CORE",
  "CUSTOM_VIC_OSTEOGENIC_DIFFERENTIATION"
)
hypothesis_sets <- c(
  list(HALLMARK_MTORC1_SIGNALING = hallmark_sets[["HALLMARK_MTORC1_SIGNALING"]]),
  custom_sets[hypothesis_pool_names[-1]]
)
supplementary_sets <- list(
  CUSTOM_S100A9_RAGE_CALCIFICATION =
    custom_sets[["CUSTOM_S100A9_RAGE_CALCIFICATION"]]
)
hallmark_pool <- hallmark_sets

# 5d: coverage check across all 3 cohorts
compute_coverage <- function(sets, ranks, pool_label, cohort_label) {
  do.call(rbind, lapply(names(sets), function(nm) {
    set_genes <- sets[[nm]]
    n_input  <- length(set_genes)
    n_mapped <- sum(set_genes %in% names(ranks))
    decision <- if (n_mapped >= 8L) "proceed" else "insufficient_coverage"
    data.frame(
      cohort = cohort_label, pool = pool_label, gene_set = nm,
      N_input = n_input, N_mapped = n_mapped,
      coverage_decision = decision, stringsAsFactors = FALSE
    )
  }))
}

cohort_ranks <- list(
  GSE51472   = ranks_GSE51472,
  INTEGRATED = ranks_INTEGRATED,
  GSE83453   = ranks_GSE83453
)

coverage_all <- do.call(rbind, lapply(names(cohort_ranks), function(cn) {
  rbind(
    compute_coverage(hypothesis_sets,    cohort_ranks[[cn]], "hypothesis_5pool", cn),
    compute_coverage(supplementary_sets, cohort_ranks[[cn]], "supplementary",    cn),
    compute_coverage(hallmark_pool,      cohort_ranks[[cn]], "hallmark_50pool",  cn)
  )
}))

cat("\n  coverage check (hypothesis pool + supplementary, all 3 cohorts):\n")
print(subset(coverage_all, pool %in% c("hypothesis_5pool", "supplementary")))

.PATH_COVERAGE <- file.path(.PATH_OUT_DIR, "coverage_check_v2.tsv")
write.table(coverage_all, file = .PATH_COVERAGE,
            sep = "\t", quote = FALSE, row.names = FALSE,
            fileEncoding = "UTF-8", eol = "\n")
cat(sprintf("  coverage table written: %s\n", .PATH_COVERAGE))


# 5e: fgseaMultilevel runner (per plan §5.1)
run_fgsea_pool <- function(pathways, ranks, label) {
  cat(sprintf("\n  fgsea call [%s] (n_pathways = %d)...\n",
              label, length(pathways)))
  set.seed(.PLAN_SEED)
  res <- fgseaMultilevel(
    pathways    = pathways,
    stats       = ranks,
    minSize     = 8L,
    maxSize     = 500L,
    eps         = 0,
    sampleSize  = 101L,
    BPPARAM     = BiocParallel::SerialParam()
  )
  res$cohort <- label
  cat(sprintf("    pathways tested: %d\n", nrow(res)))
  res
}

# Call 1-3: GSE51472 (hypothesis_5pool + supplementary + hallmark)
res_h_GSE51472 <- run_fgsea_pool(hypothesis_sets,    ranks_GSE51472,  "GSE51472_hypothesis_5pool")
res_s_GSE51472 <- run_fgsea_pool(supplementary_sets, ranks_GSE51472,  "GSE51472_supplementary")
res_H_GSE51472 <- run_fgsea_pool(hallmark_pool,      ranks_GSE51472,  "GSE51472_hallmark_50pool")

# Call 4-5: INTEGRATED (hypothesis + supplementary)
res_h_INT      <- run_fgsea_pool(hypothesis_sets,    ranks_INTEGRATED, "INTEGRATED_hypothesis_5pool")
res_s_INT      <- run_fgsea_pool(supplementary_sets, ranks_INTEGRATED, "INTEGRATED_supplementary")

# Call 6-7: GSE83453 (hypothesis + supplementary)
res_h_GSE83453 <- run_fgsea_pool(hypothesis_sets,    ranks_GSE83453,   "GSE83453_hypothesis_5pool")
res_s_GSE83453 <- run_fgsea_pool(supplementary_sets, ranks_GSE83453,   "GSE83453_supplementary")

# BH within pool (already per-call in fgsea, re-verify explicitly)
bh_within_pool <- function(res, pool_label) {
  res$padj_within_pool <- p.adjust(res$pval, method = "BH")
  res$bh_pool <- pool_label
  res
}

res_h_GSE51472 <- bh_within_pool(res_h_GSE51472, "hypothesis_5pool")
res_s_GSE51472 <- bh_within_pool(res_s_GSE51472, "supplementary_solo")
res_H_GSE51472 <- bh_within_pool(res_H_GSE51472, "hallmark_50pool")
res_h_INT      <- bh_within_pool(res_h_INT,      "hypothesis_5pool")
res_s_INT      <- bh_within_pool(res_s_INT,      "supplementary_solo")
res_h_GSE83453 <- bh_within_pool(res_h_GSE83453, "hypothesis_5pool")
res_s_GSE83453 <- bh_within_pool(res_s_GSE83453, "supplementary_solo")


# ---- Stage 6: NES_SD bootstrap + output ------------------------------------

cat("\n[stage 6] writing outputs...\n")

# Combine all results (3 cohorts hypothesis + supplementary + 1 Hallmark)
all_results <- rbind(
  res_h_GSE51472, res_s_GSE51472, res_H_GSE51472,
  res_h_INT,      res_s_INT,
  res_h_GSE83453, res_s_GSE83453
)
all_results$N_input  <- coverage_all$N_input[match(
  paste(all_results$cohort, all_results$pathway),
  paste(
    ifelse(grepl("GSE51472",   coverage_all$cohort), "GSE51472",
    ifelse(grepl("INTEGRATED", coverage_all$cohort), "INTEGRATED",
                                                     "GSE83453")),
    coverage_all$gene_set))]
# Simpler: rebuild N_input/N_mapped from coverage_all by cohort name in result
result_cohort_short <- ifelse(grepl("^GSE51472",   all_results$cohort), "GSE51472",
                       ifelse(grepl("^INTEGRATED", all_results$cohort), "INTEGRATED",
                                                                        "GSE83453"))
all_results$cohort_short <- result_cohort_short
all_results$N_input <- mapply(function(ch, pw) {
  m <- coverage_all$N_input[coverage_all$cohort == ch & coverage_all$gene_set == pw]
  if (length(m) >= 1) m[1] else NA_integer_
}, result_cohort_short, all_results$pathway, USE.NAMES = FALSE)
all_results$N_mapped <- mapply(function(ch, pw) {
  m <- coverage_all$N_mapped[coverage_all$cohort == ch & coverage_all$gene_set == pw]
  if (length(m) >= 1) m[1] else NA_integer_
}, result_cohort_short, all_results$pathway, USE.NAMES = FALSE)


# 6a: Bootstrap NES_SD (plan §7.1 fallback)
bootstrap_nes_sd <- function(pathway_genes, ranks, n_boot = 100L, seed_offset) {
  nes_boot <- rep(NA_real_, n_boot)
  for (i in seq_len(n_boot)) {
    seed_i <- .PLAN_SEED + seed_offset + i
    set.seed(seed_i)
    boot_ranks <- sample(ranks, length(ranks), replace = TRUE)
    names(boot_ranks) <- names(ranks)
    boot_ranks <- sort(boot_ranks, decreasing = TRUE)
    set.seed(seed_i)
    res_b <- tryCatch(
      fgseaMultilevel(pathways = list(p = pathway_genes), stats = boot_ranks,
                      minSize = 8L, maxSize = 500L, eps = 0, sampleSize = 101L,
                      BPPARAM = BiocParallel::SerialParam()),
      error = function(e) NULL
    )
    if (!is.null(res_b) && nrow(res_b) == 1L) nes_boot[i] <- res_b$NES
  }
  sd(nes_boot, na.rm = TRUE)
}

all_results$NES_SD <- NA_real_
if (RUN_BOOTSTRAP) {
  cat("  bootstrap NES_SD (100 resamples per pathway per cohort)...\n")
  cat(sprintf("  total: %d pathway-cohort combinations\n", nrow(all_results)))
  t0 <- Sys.time()
  for (i in seq_len(nrow(all_results))) {
    pw <- all_results$pathway[i]
    ch <- result_cohort_short[i]
    ranks_use <- cohort_ranks[[ch]]
    pw_genes <- if (pw %in% names(hypothesis_sets))    hypothesis_sets[[pw]]
              else if (pw %in% names(supplementary_sets)) supplementary_sets[[pw]]
              else if (pw %in% names(hallmark_pool))      hallmark_pool[[pw]]
              else NULL
    if (is.null(pw_genes) || is.null(ranks_use)) next
    all_results$NES_SD[i] <- bootstrap_nes_sd(
      pw_genes, ranks_use, n_boot = 100L, seed_offset = i * 1000L
    )
    if (i %% 20 == 0) {
      elapsed <- as.numeric(Sys.time() - t0, units = "secs")
      eta <- elapsed / i * (nrow(all_results) - i)
      cat(sprintf("    %d / %d done, elapsed %.1fs, ETA %.1fs\n",
                  i, nrow(all_results), elapsed, eta))
    }
  }
  cat(sprintf("  bootstrap finished in %.1fs\n",
              as.numeric(Sys.time() - t0, units = "secs")))
} else {
  cat("  bootstrap SKIPPED (RUN_BOOTSTRAP = FALSE); NES_SD = NA\n")
}

all_results$NES_lower_95CI <- all_results$NES - 1.96 * all_results$NES_SD
all_results$NES_upper_95CI <- all_results$NES + 1.96 * all_results$NES_SD

# Leading edge formatting
all_results$leading_edge_size <- vapply(all_results$leadingEdge, length, integer(1))
all_results$leading_edge <- vapply(all_results$leadingEdge,
                                    function(x) paste(x, collapse = ","),
                                    character(1))
all_results$ranking_metric <- "limma_logFC_with_stable_jitter_seed42"

out_cols <- c("pathway", "cohort", "cohort_short", "bh_pool",
              "N_input", "N_mapped",
              "NES", "NES_SD", "NES_lower_95CI", "NES_upper_95CI",
              "pval", "padj", "padj_within_pool",
              "leading_edge", "leading_edge_size",
              "ranking_metric")
all_results_out <- as.data.frame(all_results)[, out_cols, drop = FALSE]

.PATH_S2 <- file.path(.PATH_OUT_DIR, "Supp_Table_S2_v2.tsv")
write.table(all_results_out, file = .PATH_S2,
            sep = "\t", quote = FALSE, row.names = FALSE,
            fileEncoding = "UTF-8", eol = "\n")
cat(sprintf("  Supp_Table_S2 written: %s (%d rows)\n", .PATH_S2, nrow(all_results_out)))

# Forest plot data (plan §7.2): hypothesis pool across 3 cohorts
forestplot_data <- subset(all_results_out,
                          pathway %in% hypothesis_pool_names &
                          cohort_short %in% c("GSE51472", "INTEGRATED", "GSE83453") &
                          bh_pool == "hypothesis_5pool")
.PATH_FOREST <- file.path(.PATH_OUT_DIR, "forestplot_data.tsv")
write.table(forestplot_data, file = .PATH_FOREST,
            sep = "\t", quote = FALSE, row.names = FALSE,
            fileEncoding = "UTF-8", eol = "\n")
cat(sprintf("  forestplot_data written: %s (%d rows)\n",
            .PATH_FOREST, nrow(forestplot_data)))


# ---- Stage 7: cross-cohort consistency classification (plan §8) ------------

cat("\n[stage 7] cross-cohort consistency classification (plan §8)...\n")

classify_consistency <- function(nes_vec, padj_vec) {
  # nes_vec and padj_vec are length-3 vectors (GSE51472, INTEGRATED, GSE83453)
  if (any(is.na(nes_vec))) return("incomplete_data")
  signs <- sign(nes_vec)
  same_sign <- length(unique(signs[signs != 0])) <= 1L
  if (!same_sign) return("Inconsistent")
  n_sig <- sum(padj_vec < 0.05, na.rm = TRUE)
  if (n_sig >= 2L) "Strong consistency"
  else "Directional consistency"
}

consistency_tbl <- do.call(rbind, lapply(hypothesis_pool_names, function(pw) {
  rows <- subset(all_results_out,
                 pathway == pw &
                 cohort_short %in% c("GSE51472", "INTEGRATED", "GSE83453") &
                 bh_pool == "hypothesis_5pool")
  rows <- rows[match(c("GSE51472", "INTEGRATED", "GSE83453"),
                     rows$cohort_short), ]
  nes_vec  <- rows$NES
  padj_vec <- rows$padj_within_pool
  data.frame(
    pathway = pw,
    NES_GSE51472   = nes_vec[1],
    NES_INTEGRATED = nes_vec[2],
    NES_GSE83453   = nes_vec[3],
    padj_GSE51472   = padj_vec[1],
    padj_INTEGRATED = padj_vec[2],
    padj_GSE83453   = padj_vec[3],
    classification = classify_consistency(nes_vec, padj_vec),
    stringsAsFactors = FALSE
  )
}))

.PATH_CONSIST <- file.path(.PATH_OUT_DIR, "cross_cohort_consistency_v2.tsv")
write.table(consistency_tbl, file = .PATH_CONSIST,
            sep = "\t", quote = FALSE, row.names = FALSE,
            fileEncoding = "UTF-8", eol = "\n")
cat(sprintf("  consistency table written: %s\n", .PATH_CONSIST))

# Re-write sessionInfo at end
writeLines(capture.output(sessionInfo()), .PATH_SESSIONINFO)

cat("\n[stage 1-7 OK] all 3 cohorts complete.\n\n")

# Console summary: hypothesis pool across 3 cohorts
cat("=== Cross-cohort hypothesis-driven pool results ===\n")
print(consistency_tbl)

cat("\n=== Supplementary set (CUSTOM_S100A9_RAGE_CALCIFICATION) ===\n")
supp_summary <- subset(all_results_out,
                       pathway == "CUSTOM_S100A9_RAGE_CALCIFICATION" &
                       cohort_short %in% c("GSE51472", "INTEGRATED", "GSE83453"))
print(supp_summary[, c("cohort_short", "NES", "pval", "padj_within_pool",
                       "N_mapped")])
