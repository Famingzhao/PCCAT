suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(purrr)
  library(qs)
  library(future)
  library(ggplot2)
  library(patchwork)
  library(spacexr)
})


# Set parallel workers for functions that support future or multicore execution.
set_spatial_future_plan <- function(workers = 4) {
  future::plan("multisession", workers = workers)
  invisible(workers)
}


# Run a standard Seurat preprocessing workflow for Visium objects.
preprocess_visium_samples <- function(
  st_data_list,
  dims = 1:30,
  resolutions = c(0.05, 0.2, 0.3, 0.5, 0.8, 1.0)
) {
  lapply(st_data_list, function(st_data) {
    st_data <- SCTransform(st_data, assay = "Spatial", verbose = FALSE)
    st_data <- st_data %>%
      RunPCA(assay = "SCT", verbose = FALSE) %>%
      FindNeighbors(reduction = "pca", dims = dims) %>%
      RunUMAP(reduction = "pca", dims = dims)

    for (res in resolutions) {
      st_data <- FindClusters(
        st_data,
        resolution = res,
        algorithm = 1,
        verbose = FALSE
      )
    }

    st_data
  })
}


# Convert a Seurat Visium object into a spacexr SpatialRNA object.
build_spatial_rna_from_visium <- function(st_data, flip_x = TRUE, image_width = 540) {
  counts <- GetAssayData(st_data, slot = "counts")
  coords <- GetTissueCoordinates(st_data)[, 1:2]

  if (isTRUE(flip_x)) {
    coords$x <- image_width - coords$x
  }

  colnames(coords) <- c("ycoord", "xcoord")
  coords <- coords[, c("xcoord", "ycoord")]
  nUMI <- colSums(counts)

  spacexr::SpatialRNA(coords, counts, nUMI)
}


# Add normalized RCTD weights back into a Seurat object.
add_rctd_assay_to_seurat <- function(st_data, rctd_result, assay_name = "RCTD") {
  signature_exp <- as.data.frame(
    t(as.data.frame(spacexr::normalize_weights(rctd_result@results$weights)))
  )

  common_barcodes <- intersect(colnames(st_data), colnames(signature_exp))
  st_data <- subset(st_data, cells = common_barcodes)
  signature_exp <- signature_exp[, common_barcodes, drop = FALSE]

  st_data[[assay_name]] <- SeuratObject::CreateAssayObject(counts = signature_exp)
  st_data <- SeuratObject::SetAssayData(
    st_data,
    slot = "scale.data",
    new.data = as.matrix(signature_exp),
    assay = assay_name
  )

  st_data
}


# Plot all inferred RCTD cell-type weights for one Visium sample.
plot_rctd_spatial_features <- function(
  st_data,
  sample_id,
  outdir,
  assay_name = "RCTD",
  width = 20,
  height = 16
) {
  features_to_plot <- rownames(st_data[[assay_name]])
  p <- SpatialFeaturePlot(st_data, features = features_to_plot, ncol = 5)
  ggplot2::ggsave(
    filename = file.path(outdir, paste0("Step3.", sample_id, "_AllCellTypes_RCTD.png")),
    plot = p,
    width = width,
    height = height
  )
  invisible(p)
}


# Run RCTD for a single Visium sample and optionally save intermediate objects.
run_rctd_for_visium_sample <- function(
  st_data,
  sample_id,
  sc_reference,
  outdir = "Step3.RCTD_Res",
  celltype_id = "ct.L2",
  max_cores = 4,
  doublet_mode = "full",
  save_intermediate = TRUE,
  make_plot = TRUE
) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  reference <- if (inherits(sc_reference, "character")) {
    qs::qread(sc_reference)
  } else {
    sc_reference
  }

  puck <- build_spatial_rna_from_visium(st_data)

  if (isTRUE(save_intermediate)) {
    qs::qsave(puck, file = file.path(outdir, paste0("Step3.", sample_id, "_Stpuck.qs")))
  }

  if (!is.null(celltype_id)) {
    message("Using reference annotation label: ", celltype_id)
  }

  my_rctd <- spacexr::create.RCTD(
    puck = puck,
    reference = reference,
    max_cores = max_cores,
    CELL_MIN_INSTANCE = 1,
    test_mode = FALSE
  )
  my_rctd <- spacexr::run.RCTD(my_rctd, doublet_mode = doublet_mode)

  if (isTRUE(save_intermediate)) {
    qs::qsave(
      my_rctd,
      file = file.path(outdir, paste0("Step3.", sample_id, "_results_", doublet_mode, ".qs"))
    )
  }

  st_data <- add_rctd_assay_to_seurat(st_data, my_rctd, assay_name = "RCTD")

  if (isTRUE(make_plot)) {
    plot_rctd_spatial_features(st_data, sample_id = sample_id, outdir = outdir)
  }

  message("Finished RCTD: ", sample_id)
  st_data
}


# Run the full RCTD workflow for a named list of Visium objects.
run_rctd_for_visium_list <- function(
  st_data_list,
  sc_reference,
  outdir = "Step3.RCTD_Res",
  celltype_id = "ct.L2",
  max_cores = 4,
  doublet_mode = "full",
  save_intermediate = TRUE,
  make_plot = TRUE
) {
  purrr::imap(st_data_list, function(st_data, sample_id) {
    run_rctd_for_visium_sample(
      st_data = st_data,
      sample_id = sample_id,
      sc_reference = sc_reference,
      outdir = outdir,
      celltype_id = celltype_id,
      max_cores = max_cores,
      doublet_mode = doublet_mode,
      save_intermediate = save_intermediate,
      make_plot = make_plot
    )
  })
}


# Filter CosMx/Nanostring Seurat objects using simple QC thresholds.
filter_cosmx_samples <- function(
  nano_list,
  min_features = 25,
  min_counts = 50,
  min_area = 10,
  max_area = 500,
  assay_prefix = "Nanostring"
) {
  feature_col <- paste0("nFeature_", assay_prefix)
  count_col <- paste0("nCount_", assay_prefix)

  lapply(nano_list, function(x) {
    meta <- x@meta.data
    keep_cells <- rownames(meta)[
      meta[[feature_col]] >= min_features &
        meta[[count_col]] >= min_counts &
        meta[["Area.um2"]] >= min_area &
        meta[["Area.um2"]] <= max_area
    ]
    subset(x, cells = keep_cells)
  })
}


# Plot CosMx QC features after filtering.
plot_cosmx_qc <- function(
  nano_list,
  outdir = "./Outplot",
  features = c("nCount_Nanostring", "nFeature_Nanostring"),
  max_cutoff = "q95",
  width = 12,
  height = 5
) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  invisible(lapply(names(nano_list), function(sample_id) {
    obj <- nano_list[[sample_id]]
    p <- ImageFeaturePlot(
      obj,
      fov = sample_id,
      features = features,
      max.cutoff = max_cutoff
    )
    ggplot2::ggsave(
      plot = p,
      filename = file.path(
        outdir,
        paste0("Step3.", sample_id, "_After_QC_nCount_nFeature_ImageFeaturePlot.png")
      ),
      height = height,
      width = width
    )
  }))
}


# Example usage:
# set_spatial_future_plan(workers = 4)
#
# st.data.list <- readr::read_rds("./Outdata/Step1.stRNAseq_data.rds")
# st.data.list <- preprocess_visium_samples(st.data.list)
# rctd.res.list <- run_rctd_for_visium_list(
#   st_data_list = st.data.list,
#   sc_reference = "./Step3.SCRef_downsample.qs",
#   outdir = "Step3.RCTD_Res",
#   celltype_id = "ct.L2",
#   max_cores = 4
# )
#
# nano.list <- qs::qread("./Outdata/Step2.Nano_Seurat_Obj_list_sample_phe.qs")
# nano.list.filter <- filter_cosmx_samples(
#   nano_list = nano.list,
#   min_features = 25,
#   min_counts = 50,
#   min_area = 10,
#   max_area = 500
# )
# qs::qsave(
#   nano.list.filter,
#   file = "./Outdata/Step3.Nano_Seurat_Obj_list_sample_phe_after_QC.qs"
# )
# plot_cosmx_qc(nano.list.filter, outdir = "./Outplot")
