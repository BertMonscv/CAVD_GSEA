# =============================================================================
# draw_figure1.R
#
# Regenerates Figure 1 (baseline transcriptomic landscape, panels A-I) using
# preprocessing IDENTICAL to the frozen analysis script pre/run_gsea_v2.R
# (commit 6d6a7cd / tag prereg-frozen-v2 pipeline): same series-matrix parsing,
# same log2(x+1) linear-scale detection, same probe->gene IQR collapse, same
# limma model (~ 0 + group, Calcified - Normal).
#
# This guarantees that Figure 1 (DEG, PCA, markers) and Figure 2E (GSEA forest)
# share one preprocessing pipeline; no divergence between the DEG-level and
# pathway-level analyses.
#
# OUTPUT: results/Figure1_panels/ (individual panels) + a console summary of
# all numbers (PCA variance, DEG counts, marker logFCs, top GO/KEGG terms) so
# the manuscript Section 3.1 text can be updated to match.
#
# Panels:
#   1A  PCA of each cohort separately (Normal vs Calcified)
#   1B  PCA of integrated cohort before/after ComBat (group separation)
#   1C  Volcano: discovery (GSE51472)
#   1D  Volcano: integrated (GSE51472+GSE12644, post-ComBat)
#   1E  Osteogenic markers boxplot (discovery)
#   1F  Osteogenic markers boxplot (validation GSE83453)
#   1G  Cross-cohort logFC of key transcripts (forest)
#   1H  GO enrichment (discovery DEGs)
#   1I  KEGG enrichment (discovery DEGs)
#
# DEG threshold: |log2FC| > 1 AND adj.P.Val < 0.05  (standard; matches v22 text)
# =============================================================================

suppressPackageStartupMessages({
  library(limma)
  library(sva)
  library(ggplot2)
  library(AnnotationDbi)
  library(hgu133plus2.db)
  library(illuminaHumanv4.db)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
  library(patchwork)
})

set.seed(42)

# Global theme: center & bold all plot titles (and subtitles) for a clean look
theme_set(theme_bw(base_size = 10) +
          theme(plot.title    = element_text(hjust = 0.5, face = "bold"),
                plot.subtitle = element_text(hjust = 0.5)))

.PATH_GSE51472 <- "csv/GSE51472_series_matrix.txt"
.PATH_GSE12644 <- "csv/GSE12644_series_matrix.txt"
.PATH_GSE83453 <- "csv/GSE83453_series_matrix.txt"
.OUT <- "results/Figure1_panels"
dir.create(.OUT, showWarnings = FALSE, recursive = TRUE)

DEG_LFC <- 1.0
DEG_PADJ <- 0.05


# ---- preprocessing functions (identical to run_gsea_v2.R) -------------------

read_series_matrix <- function(path) {
  lines <- readLines(path)
  begin <- grep("^!series_matrix_table_begin", lines)
  end   <- grep("^!series_matrix_table_end",   lines)
  geo_line <- lines[grep("^!Sample_geo_accession", lines)]
  geo_ids <- gsub("\"", "", strsplit(geo_line, "\t")[[1]][-1])
  src_line <- lines[grep("^!Sample_source_name_ch1", lines)]
  src <- if (length(src_line) == 1L) gsub("\"", "", strsplit(src_line, "\t")[[1]][-1]) else rep(NA, length(geo_ids))
  char_lines <- lines[grep("^!Sample_characteristics_ch1", lines)]
  char_mat <- if (length(char_lines) > 0L) do.call(rbind, lapply(char_lines, function(L) gsub("\"", "", strsplit(L, "\t")[[1]][-1]))) else matrix("", 0, length(geo_ids))
  samples <- data.frame(sample_id = geo_ids, source_name = src, stringsAsFactors = FALSE)
  if (nrow(char_mat) > 0L) {
    cd <- as.data.frame(t(char_mat), stringsAsFactors = FALSE)
    colnames(cd) <- paste0("char_", seq_len(ncol(cd)))
    samples <- cbind(samples, cd)
  }
  expr_lines <- lines[(begin + 1L):(end - 1L)]
  con <- textConnection(expr_lines); on.exit(close(con), add = TRUE)
  expr <- as.matrix(read.table(con, header = TRUE, sep = "\t", row.names = 1,
                               check.names = FALSE, quote = "\"", comment.char = ""))
  list(expr = expr, samples = samples)
}

log2_if_needed <- function(expr) {
  qs <- quantile(as.vector(expr), 0.99, na.rm = TRUE)
  if (qs >= 100) log2(expr + 1) else expr
}

collapse_iqr <- function(expr, annot_db) {
  p2s <- AnnotationDbi::select(annot_db, keys = rownames(expr),
                               columns = "SYMBOL", keytype = "PROBEID")
  p2s <- p2s[!is.na(p2s$SYMBOL), ]
  multi <- unique(p2s$PROBEID[duplicated(p2s$PROBEID)])
  p2s <- p2s[!p2s$PROBEID %in% multi, ]
  p2s <- p2s[p2s$PROBEID %in% rownames(expr), ]
  ex <- expr[p2s$PROBEID, , drop = FALSE]
  iqr <- apply(ex, 1, IQR)
  tbl <- data.frame(probe = p2s$PROBEID, symbol = p2s$SYMBOL, iqr = iqr,
                    stringsAsFactors = FALSE)
  tbl <- tbl[order(tbl$symbol, -tbl$iqr), ]
  keep <- tbl[!duplicated(tbl$symbol), ]
  out <- ex[keep$probe, , drop = FALSE]
  rownames(out) <- keep$symbol
  out
}

run_limma <- function(expr_gene, group) {
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  fit <- lmFit(expr_gene, design)
  cont <- makeContrasts(Calcified_vs_Normal = Calcified - Normal, levels = design)
  fit2 <- eBayes(contrasts.fit(fit, cont))
  topTable(fit2, coef = "Calcified_vs_Normal", number = Inf, sort.by = "none")
}


# ---- load 3 cohorts (identical exclusions to run_gsea_v2.R) -----------------

cat("[load] GSE51472...\n")
g1 <- read_series_matrix(.PATH_GSE51472)
state_col <- which(apply(g1$samples, 2, function(v) any(grepl("disease state:", v, fixed = TRUE))))
g1$samples$disease <- sub("^disease state:\\s*", "", g1$samples[[state_col]])
keep1 <- g1$samples$disease %in% c("Normal aortic valve", "Calcified aortic valve")
e1 <- log2_if_needed(g1$expr[, keep1])
grp1 <- factor(ifelse(g1$samples$disease[keep1] == "Normal aortic valve", "Normal", "Calcified"),
               levels = c("Normal", "Calcified"))
eg1 <- collapse_iqr(e1, hgu133plus2.db)

cat("[load] GSE12644...\n")
g2 <- read_series_matrix(.PATH_GSE12644)
g2$samples$disease <- g2$samples$source_name
keep2 <- g2$samples$disease %in% c("Normal aortic valve", "Calcified aortic valve")
e2 <- log2_if_needed(g2$expr[, keep2])
grp2 <- factor(ifelse(g2$samples$disease[keep2] == "Normal aortic valve", "Normal", "Calcified"),
               levels = c("Normal", "Calcified"))
eg2 <- collapse_iqr(e2, hgu133plus2.db)

cat("[load] GSE83453...\n")
g3 <- read_series_matrix(.PATH_GSE83453)
g3$samples$disease <- g3$samples$source_name
is_n <- grepl("normal", g3$samples$source_name, ignore.case = TRUE)
is_c <- grepl("stenotic", g3$samples$source_name, ignore.case = TRUE)
keep3 <- is_n | is_c
e3 <- log2_if_needed(g3$expr[, keep3])
grp3 <- factor(ifelse(is_n[keep3], "Normal", "Calcified"), levels = c("Normal", "Calcified"))
eg3 <- collapse_iqr(e3, illuminaHumanv4.db)


# ---- integrated cohort + ComBat (identical to run_gsea_v2.R) ----------------

cat("[combat] integrating GSE51472 + GSE12644...\n")
common <- intersect(rownames(eg1), rownames(eg2))
ecomb <- cbind(eg1[common, ], eg2[common, ])
batch <- c(rep("GSE51472", ncol(eg1)), rep("GSE12644", ncol(eg2)))
grpC <- factor(c(as.character(grp1), as.character(grp2)), levels = c("Normal", "Calcified"))
mod <- model.matrix(~ grpC)
ecombat <- ComBat(dat = ecomb, batch = batch, mod = mod, par.prior = TRUE, prior.plots = FALSE)


# ---- limma DE per cohort ----------------------------------------------------

cat("[limma] DE...\n")
tt1 <- run_limma(eg1, grp1)         # discovery
ttI <- run_limma(ecombat, grpC)     # integrated
tt3 <- run_limma(eg3, grp3)         # validation

deg_count <- function(tt) {
  up <- sum(tt$logFC > DEG_LFC & tt$adj.P.Val < DEG_PADJ)
  dn <- sum(tt$logFC < -DEG_LFC & tt$adj.P.Val < DEG_PADJ)
  c(up = up, down = dn, total = up + dn)
}
cat("\n=== DEG COUNTS (|log2FC|>1 & adj.P<0.05) ===\n")
cat(sprintf("  Discovery  (GSE51472):       up=%d down=%d total=%d\n", deg_count(tt1)[1], deg_count(tt1)[2], deg_count(tt1)[3]))
cat(sprintf("  Integrated (51472+12644):    up=%d down=%d total=%d\n", deg_count(ttI)[1], deg_count(ttI)[2], deg_count(ttI)[3]))
cat(sprintf("  Validation (GSE83453):       up=%d down=%d total=%d\n", deg_count(tt3)[1], deg_count(tt3)[2], deg_count(tt3)[3]))


# ---- PCA helper -------------------------------------------------------------

pca_panel <- function(expr_gene, group, title) {
  p <- prcomp(t(expr_gene), scale. = FALSE)
  v <- 100 * p$sdev^2 / sum(p$sdev^2)
  df <- data.frame(PC1 = p$x[,1], PC2 = p$x[,2], group = group)
  ggplot(df, aes(PC1, PC2, color = group)) +
    geom_point(size = 2.5, alpha = 0.85) +
    scale_color_manual(values = c(Normal = "#1F77B4", Calcified = "#D62728")) +
    labs(title = title,
         x = sprintf("PC1 (%.1f%%)", v[1]),
         y = sprintf("PC2 (%.1f%%)", v[2])) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold"))
}

cat("\n=== PCA VARIANCE (PC1%) ===\n")
for (nm in c("GSE51472", "GSE12644", "GSE83453")) {
  eg <- switch(nm, GSE51472 = eg1, GSE12644 = eg2, GSE83453 = eg3)
  p <- prcomp(t(eg), scale. = FALSE)
  v <- 100 * p$sdev^2 / sum(p$sdev^2)
  cat(sprintf("  %s: PC1=%.1f%% PC2=%.1f%%\n", nm, v[1], v[2]))
}


# ---- Panel 1A: per-cohort PCA -----------------------------------------------

p1A <- (pca_panel(eg1, grp1, "GSE51472 (discovery)") |
        pca_panel(eg2, grp2, "GSE12644") |
        pca_panel(eg3, grp3, "GSE83453 (validation)")) +
       plot_layout(guides = "collect") & theme(legend.position = "bottom")
ggsave(file.path(.OUT, "Figure1A_PCA_percohort.pdf"), p1A, width = 12, height = 4)


# ---- Panel 1B: integrated PCA before/after ComBat ---------------------------

pca_combat <- function(mat, title) {
  p <- prcomp(t(mat), scale. = FALSE)
  v <- 100 * p$sdev^2 / sum(p$sdev^2)
  df <- data.frame(PC1 = p$x[,1], PC2 = p$x[,2], group = grpC, batch = batch)
  ggplot(df, aes(PC1, PC2, color = group, shape = batch)) +
    geom_point(size = 2.8, alpha = 0.85) +
    scale_color_manual(values = c(Normal = "#1F77B4", Calcified = "#D62728")) +
    labs(title = title, x = sprintf("PC1 (%.1f%%)", v[1]), y = sprintf("PC2 (%.1f%%)", v[2])) +
    theme_bw(base_size = 10) + theme(legend.position = "bottom", panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold"))
}
p1B <- (pca_combat(ecomb, "Integrated: before ComBat") | pca_combat(ecombat, "Integrated: after ComBat")) +
       plot_layout(guides = "collect") & theme(legend.position = "bottom")
ggsave(file.path(.OUT, "Figure1B_PCA_combat.pdf"), p1B, width = 9, height = 4.2)


# ---- Panel 1C/1D: volcano ---------------------------------------------------

# osteogenic markers we ALWAYS want labelled on the volcano (the study's focus),
# so they are not crowded out by the more-significant immune genes
OSTEO_LABELS <- c("SPP1", "IBSP", "RUNX2", "BMP2", "ALPL", "SOX9", "COL1A1", "ENPP1")

volcano <- function(tt, title) {
  tt$sig <- "ns"
  tt$sig[tt$logFC >  DEG_LFC & tt$adj.P.Val < DEG_PADJ] <- "up"
  tt$sig[tt$logFC < -DEG_LFC & tt$adj.P.Val < DEG_PADJ] <- "down"
  tt$gene <- rownames(tt)
  # label = top 10 by adj.P  UNION  osteogenic markers that are present & significant
  topg  <- head(tt[order(tt$adj.P.Val), ], 10)
  osteo <- tt[rownames(tt) %in% OSTEO_LABELS & tt$adj.P.Val < DEG_PADJ & abs(tt$logFC) > DEG_LFC, ]
  lab   <- unique(rbind(topg, osteo))
  # color osteogenic labels distinctly (orange) vs the rest (black)
  lab$labcol <- ifelse(lab$gene %in% OSTEO_LABELS, "#E69F00", "black")
  ggplot(tt, aes(logFC, -log10(adj.P.Val), color = sig)) +
    geom_point(size = 1, alpha = 0.5) +
    scale_color_manual(values = c(ns = "grey75", up = "#D62728", down = "#1F77B4")) +
    geom_vline(xintercept = c(-DEG_LFC, DEG_LFC), linetype = "dashed", color = "grey50", linewidth = 0.3) +
    geom_hline(yintercept = -log10(DEG_PADJ), linetype = "dashed", color = "grey50", linewidth = 0.3) +
    ggrepel::geom_text_repel(data = lab, aes(label = gene), color = lab$labcol,
                             size = 2.6, fontface = ifelse(lab$gene %in% OSTEO_LABELS, "bold", "plain"),
                             max.overlaps = 30, min.segment.length = 0, seed = 42) +
    labs(title = title, x = "log2 fold-change (Calcified vs Normal)", y = "-log10 adj.P") +
    theme_bw(base_size = 10) +
    theme(legend.position = "none", panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold"))
}
if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")
library(ggrepel)
ggsave(file.path(.OUT, "Figure1C_volcano_discovery.pdf"), volcano(tt1, "Discovery (GSE51472)"), width = 5.5, height = 5)
ggsave(file.path(.OUT, "Figure1D_volcano_validation.pdf"), volcano(tt3, "Validation (GSE83453)"), width = 5.5, height = 5)


# ---- Panel 1E/1F: osteogenic markers ----------------------------------------

markers <- c("SPP1", "IBSP", "RUNX2", "BMP2", "ALPL", "SOX9")
marker_box <- function(expr_gene, group, title) {
  present <- markers[markers %in% rownames(expr_gene)]
  df <- do.call(rbind, lapply(present, function(g)
    data.frame(gene = g, expr = expr_gene[g, ], group = group)))
  df$gene <- factor(df$gene, levels = markers)
  ggplot(df, aes(group, expr, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    geom_jitter(width = 0.15, size = 1, alpha = 0.7) +
    facet_wrap(~ gene, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = c(Normal = "#1F77B4", Calcified = "#D62728")) +
    labs(title = title, x = NULL, y = "log2 expression") +
    theme_bw(base_size = 9) + theme(legend.position = "none",
                                    axis.text.x = element_text(angle = 30, hjust = 1),
                                    plot.title = element_text(hjust = 0.5, face = "bold"))
}
ggsave(file.path(.OUT, "Figure1E_markers_discovery.pdf"), marker_box(eg1, grp1, "Discovery (GSE51472)"), width = 10, height = 3)
ggsave(file.path(.OUT, "Figure1F_markers_validation.pdf"), marker_box(eg3, grp3, "Validation (GSE83453)"), width = 10, height = 3)

cat("\n=== KEY MARKER logFC (Calcified vs Normal) ===\n")
for (g in markers) {
  l1 <- if (g %in% rownames(tt1)) tt1[g, "logFC"] else NA
  lI <- if (g %in% rownames(ttI)) ttI[g, "logFC"] else NA
  l3 <- if (g %in% rownames(tt3)) tt3[g, "logFC"] else NA
  cat(sprintf("  %-8s discovery=%+.2f integrated=%+.2f validation=%+.2f\n", g, l1, lI, l3))
}


# ---- Panel 1G: cross-cohort logFC forest of key transcripts -----------------

key_tx <- c("SPP1", "IBSP", "RUNX2", "ALPL", "COL1A1", "S100A8", "S100A9", "MMP9", "SOX9")
forest_df <- do.call(rbind, lapply(key_tx, function(g) {
  rbind(
    data.frame(gene = g, cohort = "Discovery",  lfc = if (g %in% rownames(tt1)) tt1[g,"logFC"] else NA),
    data.frame(gene = g, cohort = "Integrated", lfc = if (g %in% rownames(ttI)) ttI[g,"logFC"] else NA),
    data.frame(gene = g, cohort = "Validation", lfc = if (g %in% rownames(tt3)) tt3[g,"logFC"] else NA)
  )
}))
forest_df$gene <- factor(forest_df$gene, levels = rev(key_tx))
forest_df$cohort <- factor(forest_df$cohort, levels = c("Discovery", "Integrated", "Validation"))
p1G <- ggplot(forest_df, aes(lfc, gene, color = cohort)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_point(size = 2.5, position = position_dodge(width = 0.6)) +
  scale_color_manual(values = c(Discovery = "#D62728", Integrated = "#1F77B4", Validation = "#2CA02C")) +
  labs(title = "Cross-cohort log2FC of key transcripts", x = "log2 fold-change", y = NULL) +
  theme_bw(base_size = 10) + theme(legend.position = "bottom", panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave(file.path(.OUT, "Figure1G_crosscohort_forest.pdf"), p1G, width = 6, height = 5)


# ---- Panel 1H/1I: GO + KEGG enrichment (discovery DEGs) ---------------------

cat("\n[enrich] GO + KEGG on discovery DEGs...\n")
deg_up <- rownames(tt1)[tt1$logFC > DEG_LFC & tt1$adj.P.Val < DEG_PADJ]
deg_all <- rownames(tt1)[abs(tt1$logFC) > DEG_LFC & tt1$adj.P.Val < DEG_PADJ]
entrez <- bitr(deg_all, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID

ego <- enrichGO(entrez, OrgDb = org.Hs.eg.db, ont = "BP",
                pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.1,
                readable = TRUE)
ekegg <- enrichKEGG(entrez, organism = "hsa", pvalueCutoff = 0.05)

if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
  go_df <- head(as.data.frame(ego), 12)
  go_df$Description <- factor(go_df$Description, levels = rev(go_df$Description))
  p1H <- ggplot(go_df, aes(-log10(p.adjust), Description, size = Count, color = -log10(p.adjust))) +
    geom_point() + scale_color_viridis_c() +
    labs(title = "GO biological process (discovery DEGs)", x = "-log10 adj.P", y = NULL) +
    theme_bw(base_size = 9) + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  ggsave(file.path(.OUT, "Figure1H_GO.pdf"), p1H, width = 7, height = 5)
  cat("\n=== TOP 5 GO BP ===\n")
  print(head(as.data.frame(ego)[, c("Description", "p.adjust", "Count")], 5))
}

if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
  kegg_df <- head(as.data.frame(ekegg), 12)
  kegg_df$Description <- factor(kegg_df$Description, levels = rev(kegg_df$Description))
  p1I <- ggplot(kegg_df, aes(-log10(p.adjust), Description, size = Count, color = -log10(p.adjust))) +
    geom_point() + scale_color_viridis_c() +
    labs(title = "KEGG pathway (discovery DEGs)", x = "-log10 adj.P", y = NULL) +
    theme_bw(base_size = 9) + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  ggsave(file.path(.OUT, "Figure1I_KEGG.pdf"), p1I, width = 7, height = 5)
  cat("\n=== TOP 5 KEGG ===\n")
  print(head(as.data.frame(ekegg)[, c("Description", "p.adjust", "Count")], 5))
}

cat("\n[done] All panels written to results/Figure1_panels/\n")
cat("Review the console numbers above and update manuscript Section 3.1 to match.\n")
