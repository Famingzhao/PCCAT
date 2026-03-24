suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(IOBR)
  library(survival)
})


# Read a CIBERSORTx signature matrix exported from the web server.
read_cibersortx_signature <- function(path) {
  read.table(
    file = path,
    header = TRUE,
    sep = "\t",
    row.names = 1,
    check.names = FALSE
  )
}


# Run bulk deconvolution with an existing CN signature matrix.
run_bulk_cn_deconvolution <- function(
  expr_matrix,
  signature_matrix,
  abs_method = "sig.score",
  absolute = FALSE,
  perm = 0
) {
  IOBR::CIBERSORT(
    mixture_file = expr_matrix,
    abs_method = abs_method,
    sig_matrix = signature_matrix,
    absolute = absolute,
    perm = perm
  )
}


# Run deconvolution for a named list of bulk expression matrices.
run_bulk_cn_deconvolution_list <- function(
  expr_list,
  signature_matrix,
  abs_method = "sig.score",
  absolute = FALSE,
  perm = 0
) {
  lapply(expr_list, function(expr_matrix) {
    run_bulk_cn_deconvolution(
      expr_matrix = expr_matrix,
      signature_matrix = signature_matrix,
      abs_method = abs_method,
      absolute = absolute,
      perm = perm
    )
  })
}


# Get CN columns from a deconvolution result table.
get_cn_columns <- function(df, prefix = "^CN") {
  grep(prefix, colnames(df), value = TRUE)
}


# Bind deconvolution results to matching clinical metadata.
attach_cn_to_clinical <- function(
  clinical_df,
  deconvolution_df,
  sample_id_col = NULL
) {
  cn_cols <- get_cn_columns(deconvolution_df)
  deconv_sub <- deconvolution_df[, cn_cols, drop = FALSE]

  if (is.null(sample_id_col)) {
    common_ids <- intersect(rownames(clinical_df), rownames(deconv_sub))
    clinical_df <- clinical_df[common_ids, , drop = FALSE]
    deconv_sub <- deconv_sub[common_ids, , drop = FALSE]
    return(cbind(clinical_df, deconv_sub))
  }

  clinical_df <- clinical_df %>%
    tibble::rownames_to_column(".sample_id_internal") %>%
    left_join(
      deconv_sub %>% tibble::rownames_to_column(sample_id_col),
      by = sample_id_col
    )
  rownames(clinical_df) <- clinical_df$.sample_id_internal
  clinical_df$.sample_id_internal <- NULL
  clinical_df
}


# Add a simple CN contrast score.
add_cn_ratio_score <- function(
  df,
  numerator = "CN7",
  denominator = "CN8",
  score_name = "CN7_CN8",
  use_scaled_difference = FALSE
) {
  if (!all(c(numerator, denominator) %in% colnames(df))) {
    stop("Missing CN columns required for score: ", numerator, ", ", denominator)
  }

  if (isTRUE(use_scaled_difference)) {
    df[[score_name]] <- as.numeric(scale(df[[numerator]])) - as.numeric(scale(df[[denominator]]))
  } else {
    df[[score_name]] <- df[[numerator]] - df[[denominator]]
  }
  df
}


# Build a boxplot for one CN feature across a clinical grouping variable.
plot_cn_by_group <- function(
  df,
  feature,
  group_col,
  fill_values = NULL,
  title = NULL,
  ylab = NULL
) {
  ggplot(df, aes(x = .data[[group_col]], y = .data[[feature]], fill = .data[[group_col]])) +
    geom_boxplot(outlier.shape = NA, width = 0.8) +
    geom_jitter(width = 0.1, height = 0, size = 0.8, alpha = 0.8) +
    labs(
      x = NULL,
      y = ylab %||% feature,
      title = title %||% feature
    ) +
    theme_bw() +
    theme(legend.position = "none")
}


# Convert a list of deconvolution outputs into mean CN proportions by dataset.
summarize_cn_proportions <- function(
  deconv_list,
  dataset_order = names(deconv_list),
  cn_levels = NULL
) {
  mean_df <- lapply(names(deconv_list), function(dataset_name) {
    tmp <- deconv_list[[dataset_name]]
    cn_cols <- get_cn_columns(tmp)
    data.frame(
      Dataset = dataset_name,
      t(colMeans(tmp[, cn_cols, drop = FALSE], na.rm = TRUE)),
      check.names = FALSE
    )
  }) %>%
    bind_rows()

  mean_long <- mean_df %>%
    mutate(across(starts_with("CN"), as.numeric)) %>%
    pivot_longer(
      cols = starts_with("CN"),
      names_to = "CN",
      values_to = "Proportion"
    ) %>%
    mutate(
      Proportion = Proportion * 100,
      Dataset = factor(Dataset, levels = dataset_order)
    )

  if (!is.null(cn_levels)) {
    mean_long$CN <- factor(mean_long$CN, levels = cn_levels)
  }

  mean_long
}


# Plot a stacked CN proportion barplot across bulk cohorts.
plot_cn_proportion_barplot <- function(
  mean_long,
  cn_colors = NULL
) {
  ggplot(mean_long, aes(x = Dataset, y = Proportion, fill = CN)) +
    geom_bar(stat = "identity", width = 0.85) +
    labs(y = "Average community proportion (%)", x = "Datasets") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    {
      if (!is.null(cn_colors)) scale_fill_manual(values = cn_colors) else NULL
    }
}


# Build a median-split survival data.frame for one CN feature.
prepare_survival_groups <- function(
  df,
  feature,
  time_col,
  status_col,
  group_col = "group"
) {
  keep <- df[!is.na(df[[time_col]]) & !is.na(df[[status_col]]) & !is.na(df[[feature]]), , drop = FALSE]
  cutoff <- stats::median(keep[[feature]], na.rm = TRUE)
  keep[[group_col]] <- ifelse(keep[[feature]] > cutoff, "High", "Low")
  keep[[group_col]] <- factor(keep[[group_col]], levels = c("Low", "High"))
  keep
}


# Run a univariable Cox model for one CN feature.
run_univariable_cox <- function(
  df,
  feature,
  time_col,
  status_col
) {
  keep <- df[!is.na(df[[time_col]]) & !is.na(df[[status_col]]) & !is.na(df[[feature]]), , drop = FALSE]
  fit <- survival::coxph(
    stats::as.formula(paste0("Surv(", time_col, ", ", status_col, ") ~ ", feature)),
    data = keep
  )
  stat_df <- summary(fit)$coefficients
  conf_df <- summary(fit)$conf.int
  data.frame(
    feature = feature,
    hazard_ratio = conf_df[1, "exp(coef)"],
    conf_low = conf_df[1, "lower .95"],
    conf_high = conf_df[1, "upper .95"],
    p_value = stat_df[1, "Pr(>|z|)"],
    n = nrow(keep),
    row.names = NULL,
    check.names = FALSE
  )
}


# Run univariable Cox models for all available CN features in one cohort.
run_univariable_cox_for_cn_panel <- function(
  df,
  time_col,
  status_col,
  features = NULL
) {
  if (is.null(features)) {
    features <- get_cn_columns(df)
  }
  bind_rows(lapply(features, function(feature) {
    run_univariable_cox(df, feature = feature, time_col = time_col, status_col = status_col)
  }))
}


# Run univariable Cox models across multiple clinical cohorts.
run_univariable_cox_for_cohort_list <- function(
  clinical_list,
  time_col,
  status_col,
  features = NULL
) {
  bind_rows(lapply(names(clinical_list), function(dataset_name) {
    out <- run_univariable_cox_for_cn_panel(
      df = clinical_list[[dataset_name]],
      time_col = time_col,
      status_col = status_col,
      features = features
    )
    out$Dataset <- dataset_name
    out
  }))
}


# Null-coalescing helper for plotting labels.
`%||%` <- function(a, b) if (!is.null(a)) a else b


# Example workflow:
# signature_matrix <- read_cibersortx_signature(
#   "./Step4.CibersortX/CIBERSORTx_signature_matrix.txt"
# )
#
# deconv_list <- run_bulk_cn_deconvolution_list(
#   expr_list = bulk_expr_list,
#   signature_matrix = signature_matrix
# )
#
# clinical_with_cn <- lapply(names(clinical_list), function(dataset_name) {
#   attach_cn_to_clinical(
#     clinical_df = clinical_list[[dataset_name]],
#     deconvolution_df = deconv_list[[dataset_name]]
#   )
# }) %>% setNames(names(clinical_list))
#
# clinical_with_cn$TCGA <- add_cn_ratio_score(
#   clinical_with_cn$TCGA,
#   numerator = "CN7",
#   denominator = "CN8",
#   score_name = "CN7_CN8"
# )
#
# p1 <- plot_cn_by_group(
#   df = clinical_with_cn$TCGA,
#   feature = "CN7_CN8",
#   group_col = "GS_group",
#   title = "TCGA",
#   ylab = "CN7 - CN8"
# )
#
# mean_long <- summarize_cn_proportions(
#   deconv_list = deconv_list,
#   dataset_order = names(deconv_list),
#   cn_levels = paste0("CN", 1:11)
# )
# p2 <- plot_cn_proportion_barplot(mean_long)
#
# cox_df <- run_univariable_cox_for_cohort_list(
#   clinical_list = clinical_with_cn,
#   time_col = "BCR.time",
#   status_col = "BCR",
#   features = paste0("CN", 1:11)
# )
