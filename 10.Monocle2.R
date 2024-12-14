library(Seurat)
library(dplyr)
library(readr)
library(tidyverse)
library(patchwork)
library(monocle)
library(qs)
library(ggplot2)
library(ggpubr)
library(ggridges)
library(ggsci)
library(patchwork)
library(openxlsx)

#### 1. inoput data
seurat.data = qread("./Outdata/Step7.CD8T_annotation_AUCell.qs")

#### 2. Run monocle2
if(F){
  #2.1 Extract data, phenotype data, and feature data from the SeuratObject
  data <- as(as.matrix(seurat.data.subset@assays$RNA@counts), 'sparseMatrix')
  pd <- new('AnnotatedDataFrame', data = seurat.data.subset@meta.data)
  fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
  fd <- new('AnnotatedDataFrame', data = fData)
  
  HSMM <- newCellDataSet(data,
                         phenoData = pd,
                         featureData = fd,
                         lowerDetectionLimit = 0.5,
                         expressionFamily = negbinomial.size())
 
  
  ### 2.2 Estimate size factors and dispersions
  HSMM <- HSMM %>% 
    estimateSizeFactors() %>% 
    estimateDispersions()
  
  ### 2.3 QC
  HSMM <- detectGenes(HSMM, min_expr = 1)
  print(head(fData(HSMM)))
  expressed_genes <- row.names(subset(fData(HSMM),
                                      num_cells_expressed >= 10))
  
  head(pData(HSMM))
  
  ### 2.4 select variable genes
  HSMM <- setOrderingFilter(HSMM, VariableFeatures(seurat.data.subset))
  plot_ordering_genes(HSMM)
  
  ### 2.5 reduceDimension
  HSMM <- reduceDimension(HSMM, max_components = 2,
                          # num_dim = 10,
                          reduction_method = 'DDRTree',
                          residualModelFormulaStr = "~SampleID",verbose = F) 
 
  HSMM <- orderCells(HSMM)
  pData(HSMM) %>% head()
  table(HSMM$State, as.character(HSMM$ct.L3))
  
  HSMM <- orderCells(HSMM, root_state = c(6,7), num_paths = NULL, reverse = T)
  
  #### 2.6 Vis
  colour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
           "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
           "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
           "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
  
  a1 <- plot_cell_trajectory(HSMM, color_by = "ct.L3",
                             cell_size = 0.1,) + 
    scale_color_manual(values = colour);a1
  head(a1$data)
  
  plotdf2=as.data.frame(t(HSMM@reducedDimS))
  colnames(plotdf2)=c("component1","component2")
  plotdf2$ct.L3 = HSMM$ct.L3
  plotdf2$Pseudotime = HSMM$Pseudotime
  
  a2 <- plot_cell_trajectory(HSMM, color_by = "State", cell_size = 0.01) + scale_color_manual(values = colour)
  a2
  
  a3 <- plot_cell_trajectory(HSMM, color_by = "Pseudotime", cell_size = 0.01)  + 
    ggsci::scale_color_gsea()
  a3
  
  a1 + a3
  
  a4 <- plot_cell_trajectory(HSMM, color_by = "Group3", cell_size = 0.01)+ 
    scale_color_manual(values = g.colSet$group3)
  a4
  (a1 + a2) / (a3 + a4)
  
  HSMM$ct.L3 = factor(HSMM$ct.L3, levels = c("CD8Tn_CCR7", "CD8Tm_IL7R", 
                                             "CD8Tem_GZMK_early", "CD8Tex_GZMK",
                                             "CD8Tex_terminal"))
  input.data = data.frame(barcode = HSMM$barcode,
                          ct.L3 = HSMM$ct.L3,
                          Pseudotime = HSMM$Pseudotime,
                          Branch = HSMM$State)
  
  ggplot(data = input.data, aes(x = ct.L3,  y = Pseudotime, fill = ct.L3))+ 
    geom_boxplot(outlier.alpha = 0)+ 
    scale_fill_manual(values = g.colSet$T.ct.L3)+ 
    theme_classic()+
    mytheme +  labs(x=NULL) + 
    theme( axis.text.x = element_text(angle = text.angle, hjust = text.hjust))
 }


#### 3.Run BEAM analysis
BEAM_res=BEAM(HSMM,
              branch_point = 2,
              cores = 10, 
              progenitor_method = "duplicate")
BEAM_res_CD8T=BEAM_res[,c("gene_short_name","pval","qval")]

# Vis
tmp1=plot_genes_branched_heatmap(HSMM[row.names(subset(BEAM_res_CD8T,qval<1e-4)),],
                                 branch_point = 2,
                                 num_clusters = 4, 
                                 cores = 10,
                                 branch_labels = c("Cell fate 1", "Cell fate 2"),
                                 hmcols = colorRampPalette(rev(brewer.pal(9, "PRGn")))(62),
                                 branch_colors = c("#979797", "#F05662", "#7990C8"),
                                 use_gene_short_name = T,
                                 show_rownames = F,
                                 return_heatmap = T 
)
tmp1$ph_res

gene_group=tmp1$annotation_row
gene_group$gene=rownames(gene_group)

### Using ClusterGVis for Vis
library(ClusterGVis)
library(RColorBrewer)

df <- plot_genes_branched_heatmap2(HSMM[row.names(subset(BEAM_res_CD8T,qval < 1e-4)),],
                                   branch_point = 2,
                                   num_clusters = 4,
                                   cores = 5,
                                   use_gene_short_name = T,
                                   show_rownames = T)
visCluster(object = df,plot.type = "heatmap",
           pseudotime_col = c("#ff0000","#cdcdcd","#0000ff"))

# enrich for each cluster
all_degs = data.frame(gene = row.names(BEAM_res_CD8T),
                      q.value = BEAM_res_CD8T$qval)
all_degs$log10q = -log10(all_degs$q.value)
all_degs = all_degs[order(all_degs$log10q, decreasing = T),]

if(T){
  test.gene = df$long.res
  test.gene = test.gene[!duplicated(test.gene$gene),]
  test.gene = dplyr::select(test.gene,cluster,gene,cluster_name)
  test.gene = dplyr::left_join(test.gene, all_degs)
  test.gene = split(test.gene, test.gene$cluster)
  str(test.gene)
  test.gene = lapply(test.gene, function(x){
    x = x[order(x$log10q, decreasing = T),]
    x$rank = 1:length(x$log10q)
    return(x)
  })
  str(test.gene)
  names(test.gene) = c("C1 (125)", "C2 (116)", "C3 (249)", "C4 (621)")
  
  gene.data = df$long.res
  gene.data$cluster = paste0("C",gene.data$cluster)
  table(gene.data$cluster)
  
  gene.data = split(gene.data, gene.data$cluster)
  
  # GOBP 
  gobp.CD8T = lapply(gene.data, function(x){
    input.gene = unique(x$gene)
    print(length(input.gene))
    Myenrich(genes =input.gene, category = "gobp")
  })
  
  Siggobpplot <- lapply(names(gobp.CD8T), function(z)barplot(gobp.CD8T[[z]],
                                                             font.size = 12,
                                                             showCategory = 10) +
                          ggtitle(z) + scale_y_discrete(labels=function(x) str_wrap(x, width=30)))
  
  p.gobp = cowplot::plot_grid(plotlist = Siggobpplot, ncol = 4) 
  p.gobp
  
  # KEGG
  kegg.endo = lapply(gene.data, function(x){
    input.gene = x$gene
    Myenrich(genes =input.gene, category = "kegg")
  })
  Sigkeggplot <- lapply(names(kegg.endo), function(z)barplot(kegg.endo[[z]],
                                                             font.size = 12,
                                                             showCategory = 10) +
                          ggtitle(z) + scale_y_discrete(labels=function(x) str_wrap(x, width=30)))
  
  keppcol <- which(sapply(Sigkeggplot, function(z)nrow(z$data)) != 0)
  if (length(keppcol) != 0){
    Sigkeggplot <- Sigkeggplot[keppcol]
    p.kepp = cowplot::plot_grid(plotlist = Sigkeggplot, ncol = 4) 
  }
  
  # ReactomePA
  library(ReactomePA)
  Reac.list <- lapply(gene.data, function(x){
    input.gene = x$gene
    sig_ID <- bitr(input.gene,fromType = 'SYMBOL',
                   toType = 'ENTREZID',
                   OrgDb = "org.Hs.eg.db")
    Reac <- enrichPathway(gene = sig_ID$ENTREZID,
                          organism = 'human',
                          pvalueCutoff = 0.05)
    Reac <- setReadable(Reac, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
    return(Reac)
  })
  reactome.plot <- lapply(names(Reac.list), function(z)barplot(Reac.list[[z]],
                                                               font.size = 12,
                                                               showCategory = 10) +
                            ggtitle(z) + scale_y_discrete(labels=function(x) str_wrap(x, width=30)))
  
  p.react = wrap_plots(reactome.plot, ncol = 4)
  p.C4 = wrap_plots(p.gobp, p.kepp,p.react, ncol = 1)
}