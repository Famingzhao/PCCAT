suppressPackageStartupMessages({
  library(Seurat)
  library(qs)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})


# Read a Seurat object from .qs or .rds.
read_seurat_object <- function(path) {
  if (grepl("\\.qs$", path, ignore.case = TRUE)) {
    return(qs::qread(path))
  }
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    return(readRDS(path))
  }
  stop("Unsupported input format: ", path)
}


# Save a Seurat object to .qs or .rds.
save_seurat_object <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (grepl("\\.qs$", path, ignore.case = TRUE)) {
    qs::qsave(object, file = path)
    return(invisible(path))
  }
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    saveRDS(object, file = path)
    return(invisible(path))
  }
  stop("Unsupported output format: ", path)
}


# Normalize user-specified RCTD cell-state names when needed.
standardize_rctd_labels <- function(
  rctd_matrix,
  rename_map = c("myoCAFs" = "vCAFs"),
  merge_map = list("mCAFs" = c("mCAFs", "dCAFs"))
) {
  stopifnot(is.matrix(rctd_matrix) || is.data.frame(rctd_matrix))
  rctd_matrix <- as.matrix(rctd_matrix)

  rownames(rctd_matrix) <- dplyr::recode(rownames(rctd_matrix), !!!rename_map)
  if (any(duplicated(rownames(rctd_matrix)))) {
    rctd_matrix <- rowsum(rctd_matrix, group = rownames(rctd_matrix), reorder = FALSE)
  }

  for (target_name in names(merge_map)) {
    source_names <- intersect(merge_map[[target_name]], rownames(rctd_matrix))
    if (length(source_names) == 0) {
      next
    }
    merged_signal <- colSums(rctd_matrix[source_names, , drop = FALSE])
    rctd_matrix <- rctd_matrix[setdiff(rownames(rctd_matrix), source_names), , drop = FALSE]
    rctd_matrix <- rbind(rctd_matrix, merged_signal)
    rownames(rctd_matrix)[nrow(rctd_matrix)] <- target_name
  }

  rctd_matrix
}


# Replace or create the RCTD assay inside a Seurat object.
set_rctd_assay <- function(
  object,
  rctd_matrix,
  assay_name = "RCTD"
) {
  rctd_matrix <- as.matrix(rctd_matrix)
  object[[assay_name]] <- SeuratObject::CreateAssayObject(counts = rctd_matrix)
  object <- SeuratObject::SetAssayData(
    object,
    assay = assay_name,
    slot = "data",
    new.data = rctd_matrix
  )
  object
}


# Summarize selected epithelial programs from the RCTD assay.
add_cn_program_scores <- function(
  object,
  assay_name = "RCTD",
  ne_features = c("LPCs-HMMR", "NE-EZH2", "NE-CHG+", "NE-FOXA2"),
  lum_features = c("Lum-DPP4", "Lum-ERG", "Lum-ETV1", "Lum-FAM110B", "Lum-KLK4", "Lum-PCA3")
) {
  assay_data <- GetAssayData(object, assay = assay_name, slot = "data")

  present_ne <- intersect(ne_features, rownames(assay_data))
  present_lum <- intersect(lum_features, rownames(assay_data))

  object$NE_cell <- if (length(present_ne) > 0) {
    colSums(assay_data[present_ne, , drop = FALSE])
  } else {
    0
  }
  object$Lum_cell <- if (length(present_lum) > 0) {
    colSums(assay_data[present_lum, , drop = FALSE])
  } else {
    0
  }

  object
}


# Run CN clustering directly from an RCTD assay.
run_cn_clustering_from_rctd <- function(
  object,
  assay_name = "RCTD",
  dims = 1:10,
  npcs = 30,
  resolutions = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0),
  cluster_graph_name = "RCTD_snn"
) {
  DefaultAssay(object) <- assay_name

  object <- object %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(npcs = npcs, verbose = FALSE, features = rownames(object)) %>%
    RunUMAP(reduction = "pca", dims = dims, verbose = FALSE) %>%
    FindNeighbors(reduction = "pca", dims = dims, graph.name = cluster_graph_name, verbose = FALSE)

  for (res in resolutions) {
    object <- FindClusters(
      object,
      resolution = res,
      algorithm = 1,
      graph.name = cluster_graph_name,
      verbose = FALSE
    )
  }

  object
}


# Optionally run BBKNN if the package is installed.
run_bbknn_if_available <- function(
  object,
  batch_key = "sampleID",
  reduction = "pca",
  umap_name = "umap",
  n_pcs = 10
) {
  if (!requireNamespace("bbknnR", quietly = TRUE)) {
    message("Package bbknnR is not installed. Skipping BBKNN.")
    return(object)
  }

  bbknnR::RunBBKNN(
    object,
    batch_key = batch_key,
    reduction = reduction,
    UMAP_name = umap_name,
    run_TSNE = FALSE,
    n_pcs = n_pcs
  )
}


# Convert a clustering column into a CN label.
assign_cn_labels <- function(
  object,
  cluster_col,
  cn_col = "CN",
  prefix = "CN"
) {
  if (!cluster_col %in% colnames(object@meta.data)) {
    stop("Cluster column not found in metadata: ", cluster_col)
  }

  cluster_values <- as.character(object@meta.data[[cluster_col]])
  unique_clusters <- sort(unique(cluster_values))
  cn_levels <- paste0(prefix, seq_along(unique_clusters))
  names(cn_levels) <- unique_clusters

  object[[cn_col]] <- factor(
    cn_levels[cluster_values],
    levels = cn_levels
  )
  object
}


# Build a multi-resolution UMAP panel.
plot_cn_resolution_panel <- function(
  object,
  cluster_cols,
  ncol = 4
) {
  plot_list <- lapply(cluster_cols, function(cluster_col) {
    DimPlot(object, reduction = "umap", group.by = cluster_col, label = TRUE) +
      NoAxes() +
      ggtitle(cluster_col)
  })
  patchwork::wrap_plots(plotlist = plot_list, ncol = ncol)
}


# Plot the final CN labels on the UMAP.
plot_cn_umap <- function(
  object,
  cn_col = "CN",
  cols = NULL
) {
  DimPlot(
    object,
    reduction = "umap",
    group.by = cn_col,
    label = TRUE,
    cols = cols,
    label.box = TRUE
  ) + NoAxes()
}


# Plot CN labels in tissue space for selected images.
plot_cn_spatial <- function(
  object,
  images = NULL,
  cn_col = "CN",
  cols = NULL,
  pt_size_factor = 2.5
) {
  SpatialDimPlot(
    object,
    group.by = cn_col,
    images = images,
    label = FALSE,
    cols = cols,
    pt.size.factor = pt_size_factor
  )
}


# Export compact metadata for downstream plotting or statistics.
export_cn_metadata <- function(
  object,
  outfile,
  keep_cols = c("sampleID", "orig.ident", "NE_cell", "Lum_cell", "CN")
) {
  present_cols <- intersect(keep_cols, colnames(object@meta.data))
  cn_res <- object@meta.data[, present_cols, drop = FALSE]
  dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
  write.csv(cn_res, file = outfile, row.names = TRUE)
  invisible(cn_res)
}


# Example workflow:
# st.data <- read_seurat_object("./Outdata/StepF.For_PCCAT_stRNAseq_merge_n43.qs")
# rctd_matrix <- GetAssayData(st.data, assay = "RCTD", slot = "counts")
# rctd_matrix <- standardize_rctd_labels(rctd_matrix)
# st.data <- set_rctd_assay(st.data, rctd_matrix, assay_name = "RCTD")
# st.data <- run_cn_clustering_from_rctd(st.data, assay_name = "RCTD", dims = 1:10)
# st.data <- add_cn_program_scores(st.data, assay_name = "RCTD")
# st.data <- assign_cn_labels(st.data, cluster_col = "RCTD_snn_res.0.4", cn_col = "CN")
# p_res <- plot_cn_resolution_panel(
#   st.data,
#   cluster_cols = c(
#     "RCTD_snn_res.0.1", "RCTD_snn_res.0.2", "RCTD_snn_res.0.3",
#     "RCTD_snn_res.0.4", "RCTD_snn_res.0.5", "RCTD_snn_res.0.6",
#     "RCTD_snn_res.0.8", "RCTD_snn_res.1"
#   )
# )
# p_cn <- plot_cn_umap(st.data, cn_col = "CN")
# ggsave("./Outplot/Step1.CN_resolution_panel.png", p_res, width = 16, height = 10)
# ggsave("./Outplot/Step1.CN_umap.png", p_cn, width = 7, height = 5)
# export_cn_metadata(st.data, "./Outdata/Step1.ST_cn_res.csv")
# save_seurat_object(st.data, "./Outdata/Step1.Spatial_CN_from_RCTD.qs")
