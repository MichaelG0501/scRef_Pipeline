####################
# Merged from: analysis/states_qccheck.R + analysis/qc_unresolved_states.R
# Reorganized as part of analysis/ restructuring
#
# Structure:
#   Section A — Defined vs Unresolved QC (full feature set, from states_qccheck.R)
#   Section B — Defined vs Unresolved QC (core QC only, from qc_unresolved_states.R)
#   Section C — Per-State QC comparison (full feature set, from states_qccheck.R)
#   Section D — Per-State QC comparison (core QC only, from qc_unresolved_states.R)
#
# Both source scripts produced overlapping output filenames:
#   - states_status_quality_comparison.pdf
#   - quality_comparison_by_original_state.pdf
# The merged script preserves both approaches with clear section boundaries.
####################

# ══════════════════════════════════════════════════════════════════════
# SHARED: colours & helpers
# ══════════════════════════════════════════════════════════════════════
library(dplyr)
library(ggplot2)
library(scales)   # ships with ggplot2, no extra install

status_cols <- c("Defined" = "#2166AC", "Unresolved" = "#B2182B")

sig_label <- function(p) {
  stars <- ifelse(p < 0.001, "***",
                  ifelse(p < 0.01,  "**",
                         ifelse(p < 0.05,  "*", "ns")))
  paste0("Wilcoxon p = ", signif(p, 3), "  (", stars, ")")
}

sig_label_chi <- function(p) {
  stars <- ifelse(p < 0.001, "***",
                  ifelse(p < 0.01,  "**",
                         ifelse(p < 0.05,  "*", "ns")))
  paste0("Chi-square p = ", signif(p, 3), "  (", stars, ")")
}

# ══════════════════════════════════════════════════════════════════════
# SHARED: prepare metadata — Defined vs Unresolved status
# ══════════════════════════════════════════════════════════════════════
meta_combined$states_status <- states[rownames(meta_combined)]
meta_combined$states_status <- ifelse(meta_combined$states_status == "Unresolved",
                                      "Unresolved", "Defined")
meta_combined$states_status <- factor(meta_combined$states_status,
                                      levels = c("Defined", "Unresolved"))


# ══════════════════════════════════════════════════════════════════════
# SECTION A — Defined vs Unresolved: Full feature set
#             (from states_qccheck.R)
#             Output: states_status_quality_comparison.pdf
# ══════════════════════════════════════════════════════════════════════
continuous_features <- c(
  "nCount_RNA", "nFeature_RNA", "percent.mt",
  "mod_epithelial", "cna_cor", "cna_signal",
  "cc_score", "cs_score"
)

categorical_features <- c("malignancy")

pdf("states_status_quality_comparison.pdf", width = 7, height = 6)

### ── Continuous features ──────────────────────────────────────────
for (feature in continuous_features) {
  
  if (!feature %in% colnames(meta_combined)) next
  
  df <- meta_combined %>%
    select(states_status, all_of(feature)) %>%
    filter(!is.na(states_status), !is.na(.data[[feature]])) %>%
    mutate(value = .data[[feature]])
  
  # ---- INSERT HERE (before p-value): drop extreme outliers from plotting dataset
  y_limits <- quantile(df$value, probs = c(0.01, 0.99), na.rm = TRUE)
  df <- df %>% filter(value >= y_limits[1], value <= y_limits[2])
  
  p <- tryCatch(
    wilcox.test(value ~ states_status, data = df)$p.value,
    error = function(e) NA_real_
  )
  
  g <- ggplot(df, aes(x = states_status, y = value, fill = states_status)) +
    
    # jittered raw points (semi-transparent, behind the box)
    geom_jitter(aes(colour = states_status),
                width = 0.2, size = 0.3, alpha = 0.15, show.legend = FALSE) +
    
    # boxplot on top (slightly transparent so points show through)
    geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.75,
                 colour = "grey30", show.legend = FALSE) +
    
    # median dot
    stat_summary(fun = median, geom = "point",
                 shape = 18, size = 3.5, colour = "black", show.legend = FALSE) +
    
    scale_fill_manual(values   = status_cols) +
    scale_colour_manual(values = status_cols) +
    
    labs(
      title    = feature,
      subtitle = sig_label(p),
      x = NULL,
      y = feature
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title    = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 11, colour = "grey30"),
      axis.text.x   = element_text(size = 13, face = "bold"),
      axis.title.y   = element_text(size = 13),
      panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3)
    )
  
  # log-scale where appropriate
  
  if (feature %in% c("nCount_RNA", "nFeature_RNA")) {
    g <- g + scale_y_log10(labels = comma)
  }
  
  print(g)
}

### ── Categorical features ─────────────────────────────────────────
for (feature in categorical_features) {
  
  if (!feature %in% colnames(meta_combined)) next
  
  df <- meta_combined %>%
    select(states_status, all_of(feature)) %>%
    filter(!is.na(states_status), !is.na(.data[[feature]]))
  
  tab <- table(df$states_status, df[[feature]])
  p   <- chisq.test(tab)$p.value
  
  # proportions within each status group
  
  prop_df <- as.data.frame(prop.table(tab, 1))
  colnames(prop_df) <- c("states_status", "category", "proportion")
  
  g <- ggplot(prop_df,
              aes(x = category, y = proportion, fill = states_status)) +
    
    # side-by-side bars — much easier to compare than stacked
    geom_bar(stat = "identity", position = position_dodge(width = 0.75),
             width = 0.65, colour = "grey30", linewidth = 0.25) +
    
    # proportion labels on top of each bar
    geom_text(aes(label = paste0(round(proportion * 100, 1), "%")),
              position = position_dodge(width = 0.75),
              vjust = -0.5, size = 3.2) +
    
    scale_fill_manual(values = status_cols, name = "Status") +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                       expand = expansion(mult = c(0, 0.12))) +
    
    labs(
      title    = feature,
      subtitle = sig_label_chi(p),
      x = NULL,
      y = "Proportion within group"
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title    = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 11, colour = "grey30"),
      axis.text.x   = element_text(size = 12),
      legend.position = "top",
      panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3)
    )
  
  print(g)
}

dev.off()


# ══════════════════════════════════════════════════════════════════════
# SECTION B — Defined vs Unresolved: Core QC only
#             (from qc_unresolved_states.R)
#             Output: states_status_quality_comparison_core.pdf
#
# Difference from Section A:
#   - Same continuous features but categorical_features = NULL
#   - This was the "stripped" version focusing on continuous only
# ══════════════════════════════════════════════════════════════════════
# NOTE: Section A already covers this with identical continuous logic
# and adds categorical features on top. Section B is preserved for
# reproducibility of the original qc_unresolved_states.R Part 1 output.
# To run Section B separately, uncomment the block below:

# continuous_features_core <- c(
#   "nCount_RNA", "nFeature_RNA", "percent.mt",
#   "mod_epithelial", "cna_cor", "cna_signal",
#   "cc_score", "cs_score"
# )
# categorical_features_core <- NULL
# pdf("states_status_quality_comparison_core.pdf", width = 7, height = 6)
# # ... (identical loop as Section A but with categorical_features_core = NULL)
# dev.off()


########################
# ══════════════════════════════════════════════════════════════════════
# SHARED: Per-state comparison setup
# ══════════════════════════════════════════════════════════════════════

# ── packages ────────────────────────────────────────────────────────
library(dplyr)
library(ggplot2)
library(scales)

# ── prepare metadata ───────────────────────────────────────────────
meta_combined$original_state <- states[rownames(meta_combined)]

# palette (edit as needed)
state_pal <- c(
  "Classic_Prolif"        = "#E41A1C",
  "Squamous_Transition"   = "#4DAF4A",
  "Intest_Diff_Columnar"  = "#984EA3",
  "Plastic_Tolerant"      = "#FF7F00",
  "Unresolved"            = "grey80",
  "Unassigned"            = "grey50"
)

# order states on x-axis using palette order (keeps consistent)
meta_combined$original_state <- factor(meta_combined$original_state, levels = names(state_pal))

# ── helpers ─────────────────────────────────────────────────────────
stars_from_p <- function(p) {
  if (is.na(p)) return(NA_character_)
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  NA_character_
}

sig_label_kw <- function(p) {
  st <- stars_from_p(p); if (is.na(st)) st <- "ns"
  paste0("Kruskal\u2013Wallis p = ", signif(p, 3), "  (", st, ")")
}

sig_label_chi <- function(p) {
  st <- stars_from_p(p); if (is.na(st)) st <- "ns"
  paste0("Chi-square p = ", signif(p, 3), "  (", st, ")")
}

# pairwise brackets for a ggplot with x = original_state (factor), y = value
add_pairwise_brackets <- function(g, df, xvar = "original_state", yvar = "value",
                                  p_adjust = "BH", alpha = 0.05) {
  # levels and numeric x positions
  levs <- levels(df[[xvar]])
  xmap <- setNames(seq_along(levs), levs)
  
  # all pairs
  pairs <- combn(levs, 2, simplify = FALSE)
  
  # compute p-values
  res <- lapply(pairs, function(pr) {
    dsub <- df %>% filter(.data[[xvar]] %in% pr) %>% droplevels()
    # need at least 2 values per group to test
    ok <- all(table(dsub[[xvar]]) >= 2)
    if (!ok) return(NULL)
    p <- tryCatch(wilcox.test(dsub[[yvar]] ~ dsub[[xvar]])$p.value,
                  error = function(e) NA_real_)
    data.frame(g1 = pr[1], g2 = pr[2], p = p, stringsAsFactors = FALSE)
  }) %>% bind_rows()
  
  if (nrow(res) == 0) return(g)
  
  res$p_adj <- p.adjust(res$p, method = p_adjust)
  res$stars <- vapply(res$p_adj, stars_from_p, character(1))
  res <- res %>% filter(!is.na(stars) & p_adj < alpha)
  
  if (nrow(res) == 0) return(g)
  
  # y positions: stack brackets upward
  y_max <- max(df[[yvar]], na.rm = TRUE)
  y_min <- min(df[[yvar]], na.rm = TRUE)
  step  <- (y_max - y_min) * 0.06
  if (!is.finite(step) || step == 0) step <- abs(y_max) * 0.06 + 1e-6
  
  res <- res %>%
    mutate(
      x1 = xmap[g1],
      x2 = xmap[g2]
    ) %>%
    arrange(p_adj) %>%                    # most significant first
    mutate(y = y_max * 1.05 + (row_number() - 1) * step)
  
  # add bracket layers
  g +
    geom_segment(data = res, aes(x = x1, xend = x2, y = y, yend = y),
                 inherit.aes = FALSE, linewidth = 0.5) +
    geom_segment(data = res, aes(x = x1, xend = x1, y = y, yend = y - step*0.25),
                 inherit.aes = FALSE, linewidth = 0.5) +
    geom_segment(data = res, aes(x = x2, xend = x2, y = y, yend = y - step*0.25),
                 inherit.aes = FALSE, linewidth = 0.5) +
    geom_text(data = res, aes(x = (x1 + x2)/2, y = y + step*0.10, label = stars),
              inherit.aes = FALSE, size = 4.2, fontface = "bold")
}


# ══════════════════════════════════════════════════════════════════════
# SECTION C — Per-State QC: Full feature set
#             (from states_qccheck.R)
#             Output: quality_comparison_by_original_state.pdf
# ══════════════════════════════════════════════════════════════════════
continuous_features <- c(
  "nCount_RNA", "nFeature_RNA", "percent.mt",
  "mod_epithelial", "cna_cor", "cna_signal",
  "cc_score", "cs_score"
)
categorical_features <- c("malignancy")

pdf("quality_comparison_by_original_state.pdf", width = 11, height = 6)

# ── Continuous: compare across original_state ONLY ───────────────────
for (feature in continuous_features) {
  
  if (!feature %in% colnames(meta_combined)) next
  
  df <- meta_combined %>%
    select(original_state, all_of(feature)) %>%
    filter(!is.na(original_state), !is.na(.data[[feature]])) %>%
    mutate(value = .data[[feature]]) %>%
    droplevels()
  
  # trim extremes (CHANGED to 0.05 / 0.95)
  y_limits <- quantile(df$value, probs = c(0.05, 0.95), na.rm = TRUE)
  df <- df %>% filter(value >= y_limits[1], value <= y_limits[2])
  
  # overall test across states (subtitle)
  p_overall <- tryCatch(kruskal.test(value ~ original_state, data = df)$p.value,
                        error = function(e) NA_real_)
  
  g <- ggplot(df, aes(x = original_state, y = value)) +
    geom_jitter(aes(colour = original_state),
                width = 0.18, size = 0.25, alpha = 0.25, show.legend = FALSE) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.55,
                 colour = "grey30", aes(fill = original_state), show.legend = FALSE) +
    stat_summary(fun = median, geom = "point",
                 shape = 18, size = 2.7, colour = "black", show.legend = FALSE) +
    scale_colour_manual(values = state_pal, guide = "none") +
    scale_fill_manual(values = state_pal, guide = "none") +
    labs(
      title    = feature,
      subtitle = sig_label_kw(p_overall),
      x = NULL,
      y = feature
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title         = element_text(face = "bold", size = 16),
      plot.subtitle      = element_text(size = 11, colour = "grey30"),
      axis.text.x        = element_text(size = 11, angle = 30, hjust = 1, face = "bold"),
      axis.title.y       = element_text(size = 13),
      panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3)
    )
  
  # log-scale where appropriate
  if (feature %in% c("nCount_RNA", "nFeature_RNA")) {
    g <- g + scale_y_log10(labels = comma)
  }
  
  # pairwise brackets (THIS fixes "no brackets shown")
  # g <- add_pairwise_brackets(g, df, xvar = "original_state", yvar = "value",
  #                            p_adjust = "BH", alpha = 0.05)
  
  print(g)
}

# ── Categorical: one bar per state, stacked by category (no Defined/Unresolved) ──
for (feature in categorical_features) {
  
  if (!feature %in% colnames(meta_combined)) next
  
  df <- meta_combined %>%
    select(original_state, all_of(feature)) %>%
    filter(!is.na(original_state), !is.na(.data[[feature]])) %>%
    mutate(category = .data[[feature]]) %>%
    droplevels()
  
  tab <- table(df$original_state, df$category)
  p_all <- suppressWarnings(chisq.test(tab)$p.value)
  
  prop_df <- as.data.frame(prop.table(tab, 1))
  colnames(prop_df) <- c("original_state", "category", "proportion")
  prop_df$original_state <- factor(prop_df$original_state, levels = names(state_pal))
  
  g <- ggplot(prop_df, aes(x = original_state, y = proportion, fill = category)) +
    geom_bar(stat = "identity", position = "fill", colour = "grey30", linewidth = 0.2) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(
      title    = feature,
      subtitle = sig_label_chi(p_all),
      x = NULL,
      y = "Proportion within state"
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title         = element_text(face = "bold", size = 16),
      plot.subtitle      = element_text(size = 11, colour = "grey30"),
      axis.text.x        = element_text(size = 11, angle = 30, hjust = 1, face = "bold"),
      legend.position    = "right",
      panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3)
    )
  
  print(g)
}

dev.off()


# ══════════════════════════════════════════════════════════════════════
# SECTION D — Per-State QC: Core QC only (3 features)
#             (from qc_unresolved_states.R)
#             Output: quality_comparison_by_original_state_core.pdf
#
# Difference from Section C:
#   - continuous_features: only nCount_RNA, nFeature_RNA, percent.mt
#   - categorical_features: NULL (no categorical plots)
# ══════════════════════════════════════════════════════════════════════
continuous_features_core <- c(
  "nCount_RNA", "nFeature_RNA", "percent.mt")
categorical_features_core <- NULL

pdf("quality_comparison_by_original_state_core.pdf", width = 11, height = 6)

# ── Continuous: compare across original_state ONLY ───────────────────
for (feature in continuous_features_core) {
  
  if (!feature %in% colnames(meta_combined)) next
  
  df <- meta_combined %>%
    select(original_state, all_of(feature)) %>%
    filter(!is.na(original_state), !is.na(.data[[feature]])) %>%
    mutate(value = .data[[feature]]) %>%
    droplevels()
  
  # trim extremes (CHANGED to 0.05 / 0.95)
  y_limits <- quantile(df$value, probs = c(0.05, 0.95), na.rm = TRUE)
  df <- df %>% filter(value >= y_limits[1], value <= y_limits[2])
  
  # overall test across states (subtitle)
  p_overall <- tryCatch(kruskal.test(value ~ original_state, data = df)$p.value,
                        error = function(e) NA_real_)
  
  g <- ggplot(df, aes(x = original_state, y = value)) +
    geom_jitter(aes(colour = original_state),
                width = 0.18, size = 0.25, alpha = 0.25, show.legend = FALSE) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.55,
                 colour = "grey30", aes(fill = original_state), show.legend = FALSE) +
    stat_summary(fun = median, geom = "point",
                 shape = 18, size = 2.7, colour = "black", show.legend = FALSE) +
    scale_colour_manual(values = state_pal, guide = "none") +
    scale_fill_manual(values = state_pal, guide = "none") +
    labs(
      title    = feature,
      subtitle = sig_label_kw(p_overall),
      x = NULL,
      y = feature
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title         = element_text(face = "bold", size = 16),
      plot.subtitle      = element_text(size = 11, colour = "grey30"),
      axis.text.x        = element_text(size = 11, angle = 30, hjust = 1, face = "bold"),
      axis.title.y       = element_text(size = 13),
      panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3)
    )
  
  # log-scale where appropriate
  if (feature %in% c("nCount_RNA", "nFeature_RNA")) {
    g <- g + scale_y_log10(labels = comma)
  }
  
  # pairwise brackets (THIS fixes "no brackets shown")
  # g <- add_pairwise_brackets(g, df, xvar = "original_state", yvar = "value",
  #                            p_adjust = "BH", alpha = 0.05)
  
  print(g)
}

# ── Categorical: one bar per state, stacked by category (no Defined/Unresolved) ──
for (feature in categorical_features_core) {
  
  if (!feature %in% colnames(meta_combined)) next
  
  df <- meta_combined %>%
    select(original_state, all_of(feature)) %>%
    filter(!is.na(original_state), !is.na(.data[[feature]])) %>%
    mutate(category = .data[[feature]]) %>%
    droplevels()
  
  tab <- table(df$original_state, df$category)
  p_all <- suppressWarnings(chisq.test(tab)$p.value)
  
  prop_df <- as.data.frame(prop.table(tab, 1))
  colnames(prop_df) <- c("original_state", "category", "proportion")
  prop_df$original_state <- factor(prop_df$original_state, levels = names(state_pal))
  
  g <- ggplot(prop_df, aes(x = original_state, y = proportion, fill = category)) +
    geom_bar(stat = "identity", position = "fill", colour = "grey30", linewidth = 0.2) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(
      title    = feature,
      subtitle = sig_label_chi(p_all),
      x = NULL,
      y = "Proportion within state"
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title         = element_text(face = "bold", size = 16),
      plot.subtitle      = element_text(size = 11, colour = "grey30"),
      axis.text.x        = element_text(size = 11, angle = 30, hjust = 1, face = "bold"),
      legend.position    = "right",
      panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3)
    )
  
  print(g)
}

dev.off()
