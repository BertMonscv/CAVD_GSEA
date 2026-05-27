# =============================================================================
# draw_figure2E_forestplot.R   (v3)
#
# Forest plot of GSEA NES with 95% CI for 5 hypothesis-driven pathways
# across 3 cohorts. Reads results/forestplot_data.tsv (frozen-run output
# from git tag gsea-run-v2, commit 18f0f7f).
#
# v3 changes:
#   - Classification labels moved OUT of the plot panel into a right-side
#     annotation strip (via patchwork), eliminating overlap with CI bars
#   - Title and subtitle centered above the plot
#   - Subtitle reworded: explain filled vs hollow, and what 95% CI is
#   - x-axis label rewritten without Unicode arrows
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    install.packages("patchwork")
  }
  library(patchwork)
})

# ---- read frozen-run data ---------------------------------------------------

fp <- read.table("results/forestplot_data.tsv",
                 sep = "\t", header = TRUE, stringsAsFactors = FALSE,
                 quote = "", comment.char = "")
stopifnot(nrow(fp) == 15L)

consist <- read.table("results/cross_cohort_consistency_v2.tsv",
                      sep = "\t", header = TRUE, stringsAsFactors = FALSE,
                      quote = "", comment.char = "")


# ---- ordering and labels ---------------------------------------------------

pathway_levels <- c(
  "HALLMARK_MTORC1_SIGNALING",
  "CUSTOM_LEPTIN_MTOR_SIGNALING",
  "CUSTOM_MTOR_AUTOPHAGY_LIPOPHAGY_AXIS",
  "CUSTOM_LIPOPHAGY_CORE",
  "CUSTOM_VIC_OSTEOGENIC_DIFFERENTIATION"
)
pathway_labels <- c(
  "Hallmark\nMTORC1",
  "CUSTOM\nLEPTIN_MTOR",
  "CUSTOM\nMTOR_AUTOPHAGY\n_LIPOPHAGY",
  "CUSTOM\nLIPOPHAGY",
  "CUSTOM\nVIC_OSTEOGENIC"
)
fp$pathway_label <- factor(
  pathway_labels[match(fp$pathway, pathway_levels)],
  levels = rev(pathway_labels)
)

cohort_levels <- c("GSE51472", "INTEGRATED", "GSE83453")
cohort_pretty <- c("Discovery (GSE51472)",
                   "Sensitivity (Integrated)",
                   "Validation (GSE83453)")
fp$cohort_label <- factor(
  cohort_pretty[match(fp$cohort_short, cohort_levels)],
  levels = cohort_pretty
)

fp$significant <- fp$padj_within_pool < 0.05

classification_map <- setNames(consist$classification, consist$pathway)
fp$classification <- classification_map[fp$pathway]


# ---- colors ----------------------------------------------------------------

cohort_colors <- c(
  "Discovery (GSE51472)"     = "#D62728",
  "Sensitivity (Integrated)" = "#1F77B4",
  "Validation (GSE83453)"    = "#2CA02C"
)
class_colors <- c(
  "Strong consistency"      = "#1A9850",
  "Directional consistency" = "#FDAE61",
  "Inconsistent"            = "#D73027"
)


# ---- x-axis range ----------------------------------------------------------

xmin_data <- min(fp$NES_lower_95CI, na.rm = TRUE)
xmax_data <- max(fp$NES_upper_95CI, na.rm = TRUE)
xlim_left  <- floor(xmin_data) - 0.3
xlim_right <- ceiling(xmax_data) + 0.3


# ---- main forest plot ------------------------------------------------------

dodge <- position_dodge(width = 0.65)

fp$point_fill <- ifelse(fp$significant,
                        cohort_colors[as.character(fp$cohort_label)],
                        NA_character_)

p_main <- ggplot(fp, aes(x = NES, y = pathway_label,
                         color = cohort_label, group = cohort_label)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "grey50", linewidth = 0.4) +
  geom_errorbar(aes(xmin = NES_lower_95CI, xmax = NES_upper_95CI),
                width = 0, linewidth = 0.7,
                position = dodge, orientation = "y") +
  geom_point(aes(fill = point_fill),
             size = 3.5, shape = 21, stroke = 1.0,
             position = dodge) +
  scale_color_manual(values = cohort_colors, name = "Cohort") +
  scale_fill_identity(guide = "none", na.value = NA) +
  scale_x_continuous(limits = c(xlim_left, xlim_right),
                     breaks = seq(-3, 3, by = 1),
                     expand = c(0, 0)) +
  labs(
    x = "Normalized Enrichment Score (NES)  (positive = up in calcified valves)",
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.y = element_text(size = 9, lineheight = 0.85),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 9),
    plot.margin = margin(t = 4, r = 4, b = 4, l = 8)
  ) +
  guides(
    color = guide_legend(
      override.aes = list(fill = unname(cohort_colors), shape = 21, size = 4)
    )
  )


# ---- right-side classification strip ---------------------------------------

# A separate "plot" that only renders text labels at the right pathway position.
# Same y-axis as p_main so rows align vertically when patchworked together.

class_df <- unique(fp[, c("pathway_label", "classification")])
class_df$class_color <- class_colors[class_df$classification]

p_class <- ggplot(class_df, aes(x = 1, y = pathway_label)) +
  geom_text(aes(label = classification, color = classification),
            hjust = 0.5, size = 3.5, fontface = "bold") +
  scale_color_manual(values = class_colors, guide = "none") +
  scale_x_continuous(limits = c(0.5, 1.5), expand = c(0, 0)) +
  theme_void(base_size = 11) +
  theme(
    # Add bottom space so the strip aligns with main panel (which has x-axis text)
    plot.margin = margin(t = 4, r = 8, b = 25, l = 4)
  )


# ---- combine with patchwork ------------------------------------------------

combined <- p_main + p_class +
  plot_layout(widths = c(5, 1.3)) +
  plot_annotation(
    title = "Figure 2E. GSEA NES across three cohorts (hypothesis-driven pool)",
    subtitle = paste0(
      "Five pre-specified pathways (plan section 4). ",
      "Filled circles: padj < 0.05; hollow: padj >= 0.05 (BH within 5-pool).\n",
      "Horizontal bars: 95% CI from 100-bootstrap NES_SD (plan section 7.1). ",
      "Right column: consistency class (plan section 8)."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      plot.subtitle = element_text(size = 9.5, hjust = 0.5,
                                    color = "grey25",
                                    lineheight = 1.2,
                                    margin = margin(b = 8))
    )
  )


# ---- write file ------------------------------------------------------------

out_pdf <- "results/Figure2E_forestplot.pdf"
ggsave(out_pdf, plot = combined, width = 11, height = 6, units = "in")
cat(sprintf("written: %s\n", out_pdf))

out_png <- "results/Figure2E_forestplot.png"
ggsave(out_png, plot = combined, width = 11, height = 6, units = "in", dpi = 300)
cat(sprintf("written: %s\n", out_png))


# ---- console summary -------------------------------------------------------

cat("\n--- classification per pathway ---\n")
print(consist[, c("pathway", "classification")])

cat("\n--- significance flags ---\n")
print(table(fp$cohort_short, fp$significant))
