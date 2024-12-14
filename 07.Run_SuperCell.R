rm(list=ls())
library(data.table)
library(ggpubr)
library(Seurat)
library(dplyr)
library(patchwork)
library(readr)
library(ggplot2)
library(RColorBrewer)
library(future)
library(clustree)
library(cowplot)
library(stringr)
library(harmony)
library(rliger)
library(reshape2)
library(scales)
library(NMF)
library(ggsci)
library(pheatmap)
library(AUCell)
library(tidyverse)
library(corrplot)
library(qs)

# change the current plan to access parallelization
plan("multisession", workers = 2)
# plan()=
options(future.globals.maxSize = 300 * 1024^3)

#### 1.load data
seurat.data = qread("Outdata/StepF.All.Cells.qs")

#### 2.SuperCell
if(T){
  library(SuperCell)
  gamma <- 20 # Graining level
  k.knn <- 5
  n.pc  <- 10 # number of PCs
  
  ### 2.1 Run SuperCell
  SC <- SCimplify(
    GetAssayData(seurat.data),  # gene expression matrix 
    gamma = gamma, # graining level
    cell.annotation = as.character(seurat.data$ct.L3),
    cell.split.condition = seurat.data$SampleID, # metacell do not mix cells from different cell lines
    n.var.genes = 1000, # number of the top variable genes to use for dimentionality reduction 
    n.pc = n.pc) # number of proncipal components to use
  
  genes.use <- SC$genes.use
  gc()
  
  SC$ct.L2 <- supercell_assign(seurat.data$ct.L2,
                               supercell_membership = SC$membership)
  head(SC$ct.L2)
  
  SC$ct.L3 <- supercell_assign(seurat.data$ct.L3,
                               supercell_membership = SC$membership)
  head(SC$ct.L3)
  
  SC$SampleID <- supercell_assign(seurat.data$SampleID,
                                  supercell_membership = SC$membership)
  head(SC$SampleID)
  
  SC$celltype.sub.article <- supercell_assign(seurat.data$celltype.sub.article,
                                              supercell_membership = SC$membership)
  head(SC$celltype.sub.article)
  
  SC.GE <- supercell_GE(GetAssayData(seurat.data), 
                        groups = SC$membership)
  dim(SC.GE)
  
  genes.use <- SC$genes.use
  
  SC$SC_PCA <- supercell_prcomp(
    Matrix::t(SC.GE),
    supercell_size = SC$supercell_size, 
    genes.use = genes.use)
  
  SC$SC_UMAP <- supercell_UMAP(
    SC, 
    n_neighbors = 10)
  
  supercell_plot_UMAP(
    SC,
    group = "ct.L2",
    title = paste0("Combined construction of metacells")
  )
  
  supercell_plot_UMAP(
    SC,
    group = "ct.L3",
    title = paste0("Combined construction of metacells")
  )
  
  # plot network of metacells colored by cell line assignment
  SC$ct.L3 = m.seurat$ct.L3
  supercell_plot(SC$graph.supercells,
                 group = SC$ct.L3,
                 color.use = g.colSet$Epi.level3,
                 seed = 1,
                 main = "Metacells colored by cell line assignment")
  
  
  SC$ct.L2 = m.seurat$ct.L2
  supercell_plot(SC$graph.supercells,
                 group = SC$ct.L2,
                 color.use = g.colSet$Epi.main,
                 seed = 1,
                 main = "Metacells colored by cell line assignment")
  
  supercell_plot(SC$graph.supercells,
                 group = SC$Group3,
                 color.use = g.colSet$group3,
                 seed = 1,
                 main = "Metacells colored by Sample Types")
  
  ### 2.2 SuperCell to Seurat
  m.seurat <- supercell_2_Seurat(SC.GE = as.matrix(SC.GE),
                                 SC = SC,
                                 fields = c("SampleID","ct.L2","ct.L3"))
  table(m.seurat$SampleID)
  table(m.seurat$ct.L2)
  head(m.seurat)
  
  meta.data = dplyr::select(seurat.data@meta.data,SampleID,Group:Core_or_extension)
  meta.data = meta.data[!duplicated(meta.data$SampleID),]
  
  m.seurat.meta = m.seurat@meta.data
  m.seurat.meta$barcode = row.names(m.seurat.meta)
  head(m.seurat.meta)
  m.seurat@meta.data= m.seurat.meta
  m.seurat.meta = left_join(m.seurat.meta,meta.data)
  row.names(m.seurat.meta) = m.seurat.meta$barcode
  
  m.seurat = CreateSeuratObject(counts = m.seurat@assays$RNA@counts,
                                meta.data = m.seurat@meta.data,
                                min.cells = 0,
                                min.features = 0)
  
  m.seurat <- m.seurat %>% NormalizeData(verbose = F) %>%
    FindVariableFeatures(selection.method = "vst", 
                         nfeatures = 1000, verbose = F) %>%
    ScaleData(verbose = F) %>%
    RunPCA(npcs = 50, verbose = F)
}