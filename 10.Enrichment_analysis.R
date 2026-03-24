suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(clusterProfiler)
  library(enrichplot)
})

# Read a GMT file and remove empty genes.
read_term2gene_clean <- function(gmt_file) {
  term2gene <- clusterProfiler::read.gmt(gmt_file)
  term2gene$gene <- dplyr::na_if(term2gene$gene, "")
  stats::na.omit(term2gene)
}


# Run GSEA from a named ranked vector.
run_gsea_enrichment <- function(
  ranked_gene_list,
  gmt_file,
  pvalue_cutoff = 0.05
) {
  term2gene <- read_term2gene_clean(gmt_file)
  clusterProfiler::GSEA(
    geneList = ranked_gene_list,
    TERM2GENE = term2gene,
    pvalueCutoff = pvalue_cutoff
  )
}


# Plot a clusterProfiler GSEA result.
plot_gsea_dotplot <- function(
  gsea_result,
  title = "GSEA analysis",
  font_size = 8,
  label_format = 80
) {
  enrichplot::dotplot(
    gsea_result,
    x = "NES",
    font.size = font_size,
    label_format = label_format
  ) +
    ggplot2::ggtitle(title) +
    ggplot2::scale_size(range = c(1.5, 3.5)) +
    ggplot2::theme(legend.key.size = grid::unit(0.5, "cm"))
}


# Run GO BP enrichment for a DEG list with a user-supplied enrichment helper.
# The enrich_fun should accept the same arguments as the original Myenrich helper.
run_go_bp_enrichment <- function(
  deg_list,
  enrich_fun,
  category = "go",
  geneid = "SYMBOL"
) {
  lapply(deg_list, enrich_fun, category = category, geneid = geneid)
}


# Plot non-empty GO enrichment results.
plot_go_bp_results <- function(
  go_results,
  show_category = 10,
  font_size = 12,
  ncol = 2
) {
  plots <- lapply(names(go_results), function(name) {
    result <- go_results[[name]]
    if (is.null(result) || nrow(as.data.frame(result)) == 0) {
      return(NULL)
    }
    clusterProfiler::barplot(
      result,
      font.size = font_size,
      showCategory = show_category,
      label_format = 80
    ) +
      ggplot2::ggtitle(name)
  })

  plots <- Filter(Negate(is.null), plots)
  if (length(plots) == 0) {
    return(NULL)
  }

  cowplot::plot_grid(plotlist = plots, ncol = ncol)
}


# Quantify pathway activity and store the result as a Seurat assay.
sc_pathway_seurat <- function(
  obj,
  method = "VISION",
  imputation = FALSE,
  ncores = 2,
  spatial = FALSE,
  gene_list = "KEGG",
  assay_name = "pathway",
  set_default_assay = FALSE
) {
  signatures_kegg <- system.file(
    "data",
    "KEGG_metabolism_nc.gmt",
    package = "scMetabolism"
  )
  signatures_reactome <- system.file(
    "data",
    "REACTOME_metabolism.gmt",
    package = "scMetabolism"
  )

  if (gene_list == "KEGG") {
    gmt_file <- signatures_kegg
    cat("Using KEGG gene sets\n")
  } else if (gene_list == "REACTOME") {
    gmt_file <- signatures_reactome
    cat("Using REACTOME gene sets\n")
  } else {
    if (!file.exists(gene_list)) {
      stop("Custom GMT file does not exist: ", gene_list)
    }
    gmt_file <- gene_list
    cat("Using custom gene sets\n")
  }

  if (spatial) {
    count_exp <- obj@assays$Spatial@counts
  } else {
    count_exp <- obj@assays$RNA@counts
  }
  count_exp <- as.data.frame(as.matrix(count_exp))

  if (isTRUE(imputation)) {
    library(rsvd)
    cat("Running ALRA imputation\n")
    result_completed <- alra(as.matrix(t(count_exp)))
    count_exp2 <- result_completed[[3]]
    rownames(count_exp2) <- colnames(count_exp)
  } else {
    count_exp2 <- count_exp
  }

  cat("Calculating pathway activity with ", method, "\n", sep = "")

  if (method == "VISION") {
    library(VISION)
    n_umi <- colSums(count_exp2)
    scaled_counts <- t(t(count_exp2) / n_umi) * stats::median(n_umi)
    vis <- Vision(scaled_counts, signatures = gmt_file)
    options(mc.cores = ncores)
    vis <- analyze(vis)
    signature_exp <- data.frame(t(vis@SigScores))
  } else if (method == "AUCell") {
    library(AUCell)
    library(GSEABase)
    cells_rankings <- AUCell_buildRankings(
      as.matrix(count_exp2),
      nCores = ncores,
      plotStats = FALSE
    )
    gene_sets <- getGmt(gmt_file)
    cells_auc <- AUCell_calcAUC(gene_sets, cells_rankings)
    signature_exp <- data.frame(getAUC(cells_auc))
  } else if (method == "ssGSEA") {
    library(GSVA)
    library(GSEABase)
    gene_sets <- getGmt(gmt_file)
    gsva_es <- gsva(
      as.matrix(count_exp2),
      gene_sets,
      method = "ssgsea",
      kcdf = "Poisson",
      parallel.sz = ncores
    )
    signature_exp <- data.frame(gsva_es)
  } else if (method == "gsva") {
    library(GSVA)
    library(GSEABase)
    gene_sets <- getGmt(gmt_file)
    gsva_es <- gsva(
      as.matrix(count_exp2),
      gene_sets,
      method = "gsva",
      kcdf = "Poisson",
      parallel.sz = ncores
    )
    signature_exp <- data.frame(gsva_es)
  } else {
    stop("Unsupported method: ", method)
  }

  colnames(signature_exp) <- rownames(obj@meta.data)
  obj[[assay_name]] <- SeuratObject::CreateAssayObject(counts = signature_exp)
  obj <- SeuratObject::SetAssayData(
    obj,
    slot = "scale.data",
    new.data = as.matrix(signature_exp),
    assay = assay_name
  )

  if (isTRUE(set_default_assay)) {
    DefaultAssay(obj) <- assay_name
  }

  obj
}


# Backward-compatible alias for older notebooks/scripts.
sc.Pathway.Seurat <- sc_pathway_seurat


# Plot pathway scores from a Seurat object.
plot_pathway_violin <- function(
  obj,
  target_pathways,
  assay_name = "pathway",
  ncol = 4,
  pt_size = 0,
  cols = NULL
) {
  Seurat::VlnPlot(
    obj,
    features = target_pathways,
    pt.size = pt_size,
    stack = FALSE,
    ncol = ncol,
    cols = cols,
    assay = assay_name
  ) &
    ggplot2::labs(x = NULL, y = "Pathway activity score") &
    ggplot2::stat_summary(
      fun = "mean",
      geom = "point",
      position = ggplot2::position_dodge(0.9)
    ) &
    ggplot2::stat_summary(
      fun.data = "mean_sd",
      geom = "errorbar",
      width = 0.15,
      size = 0.3,
      position = ggplot2::position_dodge(0.9)
    ) &
    ggplot2::theme(legend.position = "none")
}


# Example usage:
# ranked_gene_list <- sce.list
# gsea_result <- run_gsea_enrichment(
#   ranked_gene_list = ranked_gene_list,
#   gmt_file = "./h.all.v7.2.symbols_PCa_gene.list.gmt"
# )
# p_gsea <- plot_gsea_dotplot(gsea_result, title = "GSEA analysis of High vs. Low risk")
#
# sig_gobp <- run_go_bp_enrichment(mye.DEG.list, enrich_fun = Myenrich)
# p_gobp <- plot_go_bp_results(sig_gobp)
#
# m.seurat <- sc_pathway_seurat(
#   obj = m.seurat,
#   method = "AUCell",
#   imputation = FALSE,
#   ncores = 10,
#   assay_name = "pathway",
#   gene_list = "./h.all.v7.2.symbols_PCa_gene.list.gmt"
# )
# p_violin <- plot_pathway_violin(
#   obj = m.seurat,
#   target_pathways = target.pathways,
#   assay_name = "pathway",
#   cols = g.colSet$cluster.name.main
# )
