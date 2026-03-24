suppressPackageStartupMessages({
  library(Seurat)
  library(qs)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(future)
  library(CellChat)
})


# Read a Seurat object from .qs or .rds.
read_seurat_spatial_object <- function(path) {
  if (grepl("\\.qs$", path, ignore.case = TRUE)) {
    return(qs::qread(path))
  }
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    return(readRDS(path))
  }
  stop("Unsupported file format: ", path)
}


# Merge or simplify detailed cell-type labels before CellChat.
collapse_cellchat_labels <- function(
  object,
  input_col,
  output_col = "cellchat_group",
  recode_map = NULL
) {
  if (!input_col %in% colnames(object@meta.data)) {
    stop("Metadata column not found: ", input_col)
  }

  object[[output_col]] <- as.character(object@meta.data[[input_col]])
  if (!is.null(recode_map)) {
    object[[output_col]] <- dplyr::recode(object[[output_col]], !!!recode_map)
  }
  object[[output_col]] <- as.factor(object[[output_col]])
  object
}


# Prepare the inputs required by CellChat for spatial transcriptomics.
prepare_cellchat_spatial_inputs <- function(
  object,
  assay_name = "Spatial",
  group_col = "cellchat_group",
  scale_factors = NULL
) {
  if (!group_col %in% colnames(object@meta.data)) {
    stop("Grouping column not found: ", group_col)
  }

  data_input <- Seurat::GetAssayData(object, slot = "data", assay = assay_name)
  meta <- data.frame(
    labels = object@meta.data[[group_col]],
    row.names = colnames(object),
    stringsAsFactors = FALSE
  )
  spatial_locs <- Seurat::GetTissueCoordinates(
    object,
    scale = NULL,
    cols = c("imagerow", "imagecol")
  )

  if (is.null(scale_factors)) {
    scale_factors <- list(
      spot.diameter = 65,
      spot = 65,
      fiducial = 65,
      hires = 1,
      lowres = 1
    )
  }

  list(
    data_input = data_input,
    meta = meta,
    spatial_locs = spatial_locs,
    scale_factors = scale_factors
  )
}


# Create a spatial CellChat object using the human or mouse database.
create_spatial_cellchat_object <- function(
  data_input,
  meta,
  spatial_locs,
  scale_factors,
  species = c("human", "mouse")
) {
  species <- match.arg(species)

  cellchat <- createCellChat(
    object = data_input,
    meta = meta,
    group.by = "labels",
    datatype = "spatial",
    coordinates = spatial_locs,
    scale.factors = scale_factors
  )

  cellchat@DB <- if (species == "human") CellChatDB.human else CellChatDB.mouse
  cellchat
}


# Run the standard CellChat workflow for spatial transcriptomics.
run_spatial_cellchat_pipeline <- function(
  cellchat,
  workers = 4,
  future_max_size = 100 * 1024^3,
  min_cells = 10,
  compute_type = "truncatedMean",
  trim = 0.1,
  distance_use = TRUE,
  scale_distance = 0.1
) {
  future::plan("multisession", workers = workers)
  options(future.globals.maxSize = future_max_size)

  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(
    cellchat,
    type = compute_type,
    trim = trim,
    distance.use = distance_use,
    scale.distance = scale_distance
  )
  cellchat <- filterCommunication(cellchat, min.cells = min_cells)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  cellchat
}


# Plot global CellChat circle plots for interaction counts and weights.
plot_cellchat_circle_summary <- function(
  cellchat,
  color_use = NULL,
  file = NULL,
  width = 10,
  height = 5
) {
  if (!is.null(file)) {
    pdf(file, width = width, height = height)
    on.exit(dev.off(), add = TRUE)
  }

  par(mfrow = c(1, 2), xpd = TRUE)
  netVisual_circle(
    cellchat@net$count,
    vertex.weight = rowSums(cellchat@net$count),
    color.use = color_use,
    weight.scale = TRUE,
    label.edge = FALSE,
    title.name = "Number of interactions"
  )
  netVisual_circle(
    cellchat@net$weight,
    vertex.weight = rowSums(cellchat@net$weight),
    color.use = color_use,
    weight.scale = TRUE,
    label.edge = FALSE,
    title.name = "Interaction strength"
  )
}


# Plot a heatmap summary of CellChat interactions.
plot_cellchat_heatmap_summary <- function(
  cellchat,
  color_use = NULL
) {
  p1 <- netVisual_heatmap(
    cellchat,
    measure = "count",
    color.heatmap = "Reds",
    color.use = color_use
  )
  p2 <- netVisual_heatmap(
    cellchat,
    measure = "weight",
    color.heatmap = "Reds",
    color.use = color_use
  )
  patchwork::wrap_plots(p1, p2, ncol = 2)
}


# Plot one signaling pathway using circle, hierarchy, or spatial layout.
plot_cellchat_pathway <- function(
  cellchat,
  signaling,
  layout = c("circle", "hierarchy", "spatial"),
  vertex_receiver = NULL,
  color_use = NULL,
  sources_use = NULL,
  targets_use = NULL,
  alpha_image = 0.3
) {
  layout <- match.arg(layout)

  if (layout == "hierarchy") {
    return(netVisual_aggregate(
      cellchat,
      signaling = signaling,
      vertex.receiver = vertex_receiver,
      layout = "hierarchy"
    ))
  }

  if (layout == "spatial") {
    return(netVisual_aggregate(
      cellchat,
      signaling = signaling,
      layout = "spatial",
      color.use = color_use,
      sources.use = sources_use,
      targets.use = targets_use,
      alpha.image = alpha_image,
      edge.width.max = 2,
      vertex.size.max = 1
    ))
  }

  netVisual_aggregate(cellchat, signaling = signaling, layout = "circle")
}


# Plot a selected ligand-receptor pair on tissue coordinates.
plot_cellchat_lr_spatial <- function(
  cellchat,
  pair_lr,
  point_size = 1,
  cutoff = 0.05,
  enriched_only = FALSE
) {
  spatialFeaturePlot(
    cellchat,
    pairLR.use = pair_lr,
    point.size = point_size,
    do.binary = TRUE,
    cutoff = cutoff,
    enriched.only = enriched_only,
    color.heatmap = "Reds",
    direction = 1
  )
}


# Save a CellChat object.
save_cellchat_object <- function(cellchat, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (grepl("\\.qs$", path, ignore.case = TRUE)) {
    qs::qsave(cellchat, file = path)
  } else if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    saveRDS(cellchat, file = path)
  } else {
    stop("Unsupported output format: ", path)
  }
  invisible(path)
}


# Example workflow:
# st.data <- read_seurat_spatial_object("./Outdata/spatial_seurat_for_cellchat.qs")
# st.data <- collapse_cellchat_labels(
#   object = st.data,
#   input_col = "ct.l2",
#   output_col = "cellchat_group"
# )
#
# inputs <- prepare_cellchat_spatial_inputs(
#   object = st.data,
#   assay_name = "Spatial",
#   group_col = "cellchat_group"
# )
#
# cellchat <- create_spatial_cellchat_object(
#   data_input = inputs$data_input,
#   meta = inputs$meta,
#   spatial_locs = inputs$spatial_locs,
#   scale_factors = inputs$scale_factors,
#   species = "human"
# )
#
# cellchat <- run_spatial_cellchat_pipeline(cellchat, workers = 4)
# save_cellchat_object(cellchat, "./Outdata/cellchat_spatial_results.qs")
#
# plot_cellchat_circle_summary(
#   cellchat,
#   file = "./Outplot/cellchat_circle_summary.pdf"
# )
#
# p_heat <- plot_cellchat_heatmap_summary(cellchat)
# ggsave("./Outplot/cellchat_heatmap_summary.pdf", p_heat, width = 12, height = 6)
