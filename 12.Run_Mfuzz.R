library(Seurat)
library(dplyr)
library(patchwork)
library(qs)
library(openxlsx)
library(ClusterGVis)
library(org.Hs.eg.db)
library(ggpubr)
library(readr)
library(ggplot2)

############# 1.load data
endo.data <- qread(file = "./Step9.Mfuzz/Step9.Mfuzz.Endo.qs")

############# 2. Mfuzz
endo.data.list= SplitObject(endo.data, split.by = "Group")
expr.dat = endo.data@assays$RNA@data

# filter low expression genes
pct.exp <- apply(X = t(expr.dat), MARGIN = 2, FUN = PercentAbove, 
                 threshold = 0)
select.gene = names(pct.exp)[pct.exp > 0.1]

input.data = AverageExpression(object = endo.data, assays = "RNA", group.by = "Group")
input.data = input.data$RNA[select.gene,] %>% as.data.frame()

### 2.1 Run Mfuzz
head(input.data)
exps = input.data

if(F){
  # check optimal cluster numbers
  getClusters(exp = exps)
  
  # using mfuzz for clustering
  cm <- clusterData(exp = exps,
                    cluster.method = "mfuzz",
                    cluster.num = 6)
  
  # plot line only
  visCluster(object = cm,
             plot.type = "line",
             ms.col = c("#08519C", "white", "#A50F15"))
  
  # plot heatmap only
  visCluster(object = cm,
             column_names_rot = 45,
             plot.type = "both",
             sample.col = g.colSet$Group[levels(cm$long.res$cell_type)],
             ctAnno.col = g.colSet$color_DataSets)
  
  # enrich for clusters
  enrich <- enrichCluster(object = cm,
                          OrgDb = org.Hs.eg.db,
                          type = "BP",
                          pvalueCutoff = 0.05,
                          topn = 50,
                          seed = 5201314)
  
  check
  head(enrich,3)
  
  # annotation for specific clusters
  pdf('./Step9.Mfuzz/Endo_Mfuzz_gene_cluster_subc.pdf',height = 10,width = 10,onefile = F)
  visCluster(object = cm,
             plot.type = "both",
             column_names_rot = 45,
             show_row_dend = F,
             # markGenes = markGenes,
             markGenes.side = "left",
             # genes.gp = c('italic',fontsize = 12,col = "black"),
             annoTerm.data = enrich,
             line.side = "left",
             # go.col = rep(ggsci::pal_d3()(11),each = 5),
             go.size = "pval")
  dev.off()
}
