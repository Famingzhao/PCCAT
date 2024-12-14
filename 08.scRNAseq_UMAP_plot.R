library(Seurat)
library(dplyr)
library(patchwork)
library(readr)
library(ggplot2)
library(RColorBrewer)
library(future)
library(cowplot)
library(stringr)
library(ggpubr)
library(qs)
plan("multisession", workers = 4)
plan()
options(future.globals.maxSize = 10 * 1024^3)
if(T){
  text.size = 8
  text.angle = 45
  text.hjust = 1
  legend.position = "right"
  mytheme <- theme(plot.title = element_text(size = text.size+2,color="black",hjust = 0.5),
                   axis.ticks = element_line(color = "black"),
                   axis.title = element_text(size = text.size,color ="black"), 
                   axis.text = element_text(size=text.size,color = "black"),
                   panel.grid=element_blank(), 
                   legend.position = legend.position,
                   legend.text = element_text(size= text.size),
                   legend.title= element_text(size= text.size)
  )
  
  DotPlot_2 <- function(object,
                        features,
                        assay = NULL, 
                        scale = T,
                        dot.range.min = 0,
                        dot.range.max = 3.5,
                        Combine=F,
                        legend.position = "right",
                        label.size = 4,
                        label_widths = 0.1,
                        x.lab = NULL,
                        y.lab = NULL,
                        title = NULL,
                        text.size = 8,
                        text.angle = 90,
                        text.vjust = 0.5,
                        text.hjust = 1,
                        group.by = NULL,
                        color.use = NULL,
                        cols = c("lightgrey", 
                                 "blue"),
                        legend.key.size = 0.5,
                        col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6, 
                        idents = NULL, split.by = NULL, cluster.idents = FALSE, 
                        scale.by = "radius", scale.min = NA, scale.max = NA,
                        ...
  ){
    library(dplyr)
    library(ggplot2)
    library(aplot)
    mytheme <- theme(plot.title = element_text(size = text.size+2,color="black",hjust = 0.5),
                     axis.ticks = element_line(color = "black"),
                     axis.title = element_text(size = text.size,color ="black"), 
                     axis.text = element_text(size=text.size,color = "black"),
                     axis.text.x = element_text(angle = text.angle, hjust = text.hjust, 
                                                vjust = text.vjust),
                     panel.grid=element_blank(), 
                     legend.position = legend.position,
                     legend.text = element_text(size= text.size),
                     legend.title= element_text(size= text.size)
    )
    
    if(!is.null(group.by)){
      Idents(object) = group.by
    }
    
    if(Combine){
      if(is.null(color.use)){
        color.use <- alpha(paletteer::paletteer_d('ggsci::category20c_d3'), 1)
      }
      df <- data.frame(x = 0, y = levels(object), stringsAsFactors = F )
      df$y <- factor(df$y, levels = df$y )
      p1 <- ggplot(df, aes(x, y, color = factor(y))) +
        geom_point(size = label.size, show.legend = F) +
        scale_color_manual(values = color.use) +
        theme_classic() +
        scale_x_continuous(expand = c(0,0)) + mytheme + 
        theme(
          plot.margin = margin(r=0),
          axis.title = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 0,color="white"),
          axis.ticks = element_blank(),
          axis.line = element_blank()
        )
    }
    
    p2 = DotPlot(object = object,cols = cols,
                 assay = assay,
                 col.min = col.min, col.max = col.max, dot.min = dot.min, dot.scale = dot.scale, 
                 idents = idents,split.by = split.by, cluster.idents = cluster.idents, 
                 scale.by = scale.by, scale.min = scale.min, scale.max = scale.max,
                 features = features,scale = scale,...)+theme_bw()+
      mytheme + labs(x = x.lab,y = y.lab,title = title)  +
      scale_size(range = c(dot.range.min,dot.range.max))+
      theme(legend.key.size = unit(legend.key.size, "cm"))+ 
      guides(size = guide_legend(title = "Per.Exp"),
             color = guide_colorbar(title = paste("Aver.Exp.","Scaled",sep = "\n")))
    if(Combine){
      p3 = p2 %>% insert_left(p1, width=label_widths)
      return(p3)
    }else{
      return(p2) 
    }
  }
  
  plot.clusters.group = function (data = seurat_data,
                                  clusters = seurat_clusters,
                                  legend.position = "top",
                                  group = orig.ident,widths = c(3,1),
                                  log =TRUE,
                                  order=T,
                                  width = 0.9,
                                  text.size = 12,
                                  legend.title = "Group",
                                  color = 1,
                                  xlab = "",cell.counts = T,
                                  legend.key.size = 0.5,...){ 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    library(paletteer)
    mytheme = theme(plot.title = element_text(size = text.size,color="black",hjust = 0.5),
                    axis.title = element_text(size = c(text.size-2),color ="black"), 
                    axis.text = element_text(size=c(text.size-2),color = "black"),
                    panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank(),
                    legend.text = element_text(size= c(text.size-4)),
                    legend.title= element_text(size= c(text.size-4)),
                    legend.key.size = unit(legend.key.size, "cm") 
    )
    
    count_table <- table(data@meta.data[,clusters], data@meta.data[,group])
    count_mtx <- as.data.frame.matrix(count_table)
    count_mtx$cluster <- rownames(count_mtx)
    melt_mtx <- melt(count_mtx)
    melt_mtx$cluster <- melt_mtx$cluster
    
    cluster_size <- aggregate(value ~ cluster, data = melt_mtx, FUN = sum)
    
    if(!is.factor(data@meta.data[,clusters])){
      data@meta.data[,clusters] = as.factor(data@meta.data[,clusters])
      cluster_size$cluster <- factor(cluster_size$cluster,levels = levels(data@meta.data[,clusters]))
      melt_mtx$cluster <- factor(melt_mtx$cluster,levels = levels(data@meta.data[,clusters]))
    }
    if("0" %in% cluster_size$cluster){
      sorted_labels <- paste(sort(as.integer(levels(cluster_size$cluster)),decreasing = T))
    }else{
      sorted_labels <- paste(cluster_size$cluster[order(cluster_size$value)])
    }
    if(order){
      cluster_size$cluster <- factor(cluster_size$cluster,levels = sorted_labels)
      melt_mtx$cluster <- factor(melt_mtx$cluster,levels = sorted_labels)
    }else{
      cluster_size$cluster <- factor(cluster_size$cluster,levels = levels(data@meta.data[,clusters]))
      melt_mtx$cluster <- factor(melt_mtx$cluster,levels = levels(data@meta.data[,clusters]))
    }
    colnames(melt_mtx)[2] <- "dataset"
    
    if(log){
      p1 <- ggplot(cluster_size, aes(y= cluster,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
        theme_bw() + scale_x_log10() + xlab("Cells per cluster") + ylab("") + mytheme
    }else{
      p1 <- ggplot(cluster_size, aes(y= cluster,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
        theme_bw() + xlab("Cells per cluster") + ylab("") + mytheme
    }
    
    if(color[1]==1){
      if(length(unique(melt_mtx$dataset)) < 21){
        p2 <- ggplot(melt_mtx,aes(x=cluster,y=value,fill=dataset)) + 
          geom_bar(position="fill", stat="identity",width = width) + theme_bw()  + 
          scale_fill_manual(values = alpha(paletteer::paletteer_d('ggsci::category20c_d3'), 0.9))+
          ylab(paste0("Fraction of cells in each ",tolower(legend.title))) + xlab(xlab) + 
          theme(legend.position=legend.position) + guides(fill = guide_legend(title = legend.title)) +
          scale_y_continuous(labels = function(x) x * 100, expand = c(0.01, 0.01))+ mytheme
      }else{
        warning("The color limit is <21")
        p2 <- ggplot(melt_mtx,aes(x=cluster,y=value,fill=dataset)) + 
          geom_bar(position="fill", stat="identity",width = width) + theme_bw()  + 
          ylab(paste0("Fraction of cells in each ",tolower(legend.title))) + xlab(xlab) +
          theme(legend.position=legend.position) + guides(fill = guide_legend(title = legend.title)) +
          scale_y_continuous(labels =function(x) x * 100, expand = c(0.01, 0.01))+ mytheme
      }
    }
    ###########################
    if(color[1]==2){
      if(length(unique(melt_mtx$dataset)) < 9){
        p2 <- ggplot(melt_mtx,aes(x=cluster,y=value,fill=dataset)) + 
          geom_bar(position="fill", stat="identity",width = width) + theme_bw() + 
          scale_fill_brewer(palette = "Set2")+
          ylab(paste0("Fraction of cells in each ",tolower(legend.title))) + xlab(xlab) + 
          theme(legend.position=legend.position) + guides(fill = guide_legend(title = legend.title)) +
          scale_y_continuous(labels =function(x) x * 100, expand = c(0.01, 0.01))+ mytheme
      }else{
        warning("The color limit is <9")
        p2 <- ggplot(melt_mtx,aes(x=cluster,y=value,fill=dataset)) + 
          geom_bar(position="fill", stat="identity",width = width) + theme_bw()  + 
          ylab(paste0("Fraction of cells in each ",tolower(legend.title))) + xlab(xlab) +
          theme(legend.position=legend.position) + guides(fill = guide_legend(title = legend.title)) +
          scale_y_continuous(labels =function(x) x * 100, expand = c(0.01, 0.01))+ mytheme
      }
    }
    
    if (!is.numeric(color)) {
      p2 <- ggplot(melt_mtx,aes(x=cluster,y=value,fill=dataset)) + 
        geom_bar(position="fill", stat="identity",width = width) + theme_bw()  + 
        ylab(paste0("Fraction of cells in each ",tolower(legend.title))) + xlab(xlab) +
        theme(legend.position=legend.position) + guides(fill = guide_legend(title = legend.title)) +
        scale_y_continuous(labels =function(x) x * 100, expand = c(0.01, 0.01))+ mytheme
      p2 = p2 + scale_fill_manual(values = color)
    }
    
    if(cell.counts){
      p2 = wrap_plots(ncol = 2,p2 + coord_flip(),p1,widths = widths)
    }else{p2}
    return(p2)
  } 
}

########## All cells---------------
seurat.data.all = qread("../Outdata/StepF.All.Cells.qs")
DimPlot(seurat.data.all,raster=T,
        label = T,
        group.by = "ct.L1")
DimPlot(seurat.data.all,raster=T,
        label = T,
        group.by = "ct.L2")
DimPlot(seurat.data.all,raster=T,
        label = T,
        group.by = "ct.L3")

markers = c('EPCAM',"KRT8","KRT18","KRT19","KLK3", #Epi
            "COL1A1","DCN","LUM",'ACTA2','FBLN1',
            "PLVAP",'VWF',"RAMP2",'CLDN5',"PECAM1",
            "LYZ","C1QA","C1QB","CD163","CSF1R",
            'TPSAB1',"KIT",'MS4A2','TPSB2','GATA2',
            "CD3D","CD3E","CD3G","CCL5","IL7R",
            "CD79A","MS4A1", "IGHM","CD19", "CD22"
)
DotPlot_2(seurat.data.all,
          features = check_genes,
          Combine = T, dot.range.min = 0,
          color.use = g.colSet$cluster.name.main)

#### 2. Cell proportion
### 2.1 total proportion bar plot
plot.clusters.group(data = seurat.data.all,
                    clusters = "group",
                    legend.position = "right",
                    group = "ct.L2",
                    widths = c(2,1),
                    log = F,
                    text.size = 12,
                    legend.title = "Group",
                    color = 1,
                    xlab = "",
                    cell.counts = F,
                    order = F)+
  scale_fill_manual(values = g.colSet$cluster.name.main)+
  mytheme

## 2.2 OR Ro/e
meta.tb <- seurat.data.all@meta.data
meta.tb$sampleID = as.character(meta.tb$sampleID)
table( meta.tb$sampleID)
res_all = gg.tissueDist(cellInfo.tb = meta.tb,
                        meta.cluster = meta.tb$ct.L2,
                        colname.patient = "SampleID",
                        loc = meta.tb$group2,
                        cuts =c(0, 0.1, 1, 1.5, 2, 5, Inf),
                        bin.label = c("-", "+/-", "+", "++", "+++", "++++"),
                        verbose = 1,
                        z.hi.OR = 4,
                        z.hi.ROE = 2,
                        text_color = g.colSet$cluster.name.main,
                        legend.key.size = 0.3
                        
)

## 2.3 Propotion box plot
## boxplot
library(data.table)
dat.plot = data.table(meta.tb)
dat.plot$Celltype = dat.plot$ct.L2
dat.plot$Group = dat.plot$group2
dat.plot$SampleID = dat.plot$sampleID

dat.plot <- dat.plot[,.(N=.N),by=c("SampleID","Celltype","Group")]
sum.data <- dat.plot[,{.(NTotal = sum(.SD$N))},by=c("SampleID","Group")]

dat.plot = left_join(dat.plot,sum.data)
dat.plot$freq = dat.plot$N / dat.plot$NTotal *100

out.prefix.loc.cmp = "./Outplot/Vis/"
p.list.perMcls <- llply(unique(sort(dat.plot$Celltype)),function(mcls){
  text.size = 8
  text.angle = 45
  text.hjust = 1
  legend.position = "none"
  p <- ggboxplot(dat.plot[Celltype== mcls ,],x="Group",y="freq",
                 color = "Group", legend="none",title=mcls,
                 #fill = "loc",alpha=0.8,
                 #add = "none",outlier.shape=NA) +
                 xlab="",ylab="frequency",font.label = list(size = text.size, color = "black"),
                 add = "jitter",outlier.shape=NA,add.params = list(size = 1)) +
    scale_color_manual(values=g.colSet$cluster.name.main) +
    stat_compare_means(label="p.format",size = 2.8,method = "kruskal.test") +
    coord_cartesian(clip="off") +
    theme(plot.title = element_text(size = text.size+2,color="black",hjust = 0.5),
          axis.ticks = element_line(color = "black"),
          axis.title = element_text(size = text.size,color ="black"), 
          axis.text = element_text(size=text.size,color = "black"),
          axis.text.x = element_text(angle = text.angle, hjust = text.hjust ), #,vjust = 0.5
          legend.position = legend.position,
          legend.text = element_text(size= text.size),
          legend.title= element_text(size= text.size),
          strip.background = element_rect(color="black",size= 1, linetype="solid") 
    )
  return(p)
},.parallel=T)
p.prop.2 = wrap_plots(p.list.perMcls,ncol = 6);p.prop.2

#### 3.AUCell score
sc.Pathway.Seurat = function (obj, method = "VISION", imputation = F, ncores = 2, Spatial=F,
                                geneList = "KEGG",assay.names = "pathway",DefaultAssay=F) 
{
  signatures_KEGG_metab <- system.file("data", "KEGG_metabolism_nc.gmt", 
                                       package = "scMetabolism")
  signatures_REACTOME_metab <- system.file("data", "REACTOME_metabolism.gmt", 
                                           package = "scMetabolism")
  
  if (geneList == "KEGG") {
    gmtFile <- signatures_KEGG_metab
    cat("Your choice is: KEGG\n")
  }
  if (geneList == "REACTOME") {
    gmtFile <- signatures_REACTOME_metab
    cat("Your choice is: REACTOME\n")
  }
  ##########    
  if(!geneList %in% c("KEGG","REACTOME")){
    library(GSEABase)
    gmtFile <- geneList
    cat("Custom gene sets by users\n")
    if(! file.exists(gmtFile)) stop('File ',gmtFile,' not available\n')
    geneSets <- getGmt(gmtFile)
  }
  ########## 
  # library(scMetabolism)
  if(Spatial){
    countexp <- obj@assays$Spatial@counts
    countexp <- data.frame(as.matrix(countexp))
  }else{
    countexp <- obj@assays$RNA@counts
    countexp <- data.frame(as.matrix(countexp))
  }
  ##########    
  if (imputation == F) {
    countexp2 <- countexp
  }
  if (imputation == T) {
    library(rsvd)
    cat("Start imputation...\n")
    cat("Citation: George C. Linderman, Jun Zhao, Yuval Kluger. Zero-preserving imputation of scRNA-seq data using low-rank approximation. bioRxiv. doi: https://doi.org/10.1101/397588 \n")
    result.completed <- alra(as.matrix(t(countexp)))
    countexp2 <- result.completed[[3]]
    row.names(countexp2) <- colnames(countexp)
  }
  cat("Start quantify the pathway activity...\n")
  if (method == "VISION") {
    library(VISION)
    n.umi <- colSums(countexp2)
    scaled_counts <- t(t(countexp2)/n.umi) * median(n.umi)
    vis <- Vision(scaled_counts, signatures = gmtFile)
    options(mc.cores = ncores)
    vis <- analyze(vis)
    signature_exp <- data.frame(t(vis@SigScores))
  }
  if (method == "AUCell") {
    library(AUCell)
    library(GSEABase)
    cells_rankings <- AUCell_buildRankings(as.matrix(countexp2), 
                                           nCores = ncores, plotStats = F)
    geneSets <- getGmt(gmtFile)
    cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
    signature_exp <- data.frame(getAUC(cells_AUC))
  }
  if (method == "ssGSEA") {
    library(GSVA)
    library(GSEABase)
    geneSets <- getGmt(gmtFile)
    gsva_es <- gsva(as.matrix(countexp2), geneSets, method = c("ssgsea"), 
                    kcdf = c("Poisson"), parallel.sz = ncores)
    signature_exp <- data.frame(gsva_es)
  }
  if (method == "gsva") {
    library(GSVA)
    library(GSEABase)
    geneSets <- getGmt(gmtFile)
    gsva_es <- gsva(as.matrix(countexp2), geneSets, method = c("gsva"), 
                    kcdf = c("Poisson"), parallel.sz = ncores)
    signature_exp <- data.frame(gsva_es)
  }
  colnames(signature_exp) = row.names(obj@meta.data)
  obj[[assay.names]] <- SeuratObject::CreateAssayObject(counts = signature_exp)
  obj <- SeuratObject::SetAssayData(obj, slot = "scale.data",
                                    new.data = as.matrix(signature_exp), assay = assay.names)
  if(DefaultAssay == T){
    DefaultAssay(obj) <- assay.names
  }
  obj
}
m.seurat <- sc.Pathway.Seurat(obj = m.seurat, 
                                method = "AUCell", 
                                imputation = F,
                                ncores = 10,
                                assay.names = "pathway",
                                geneList = "./h.all.v7.2.symbols_PCa_gene.list.gmt")

VlnPlot(seurat.data.all,features = target.pathways,
        pt.size = 0,stack = F,ncol = 4, 
        cols = g.colSet$cluster.name.main,assay = "pathway")&
  labs(x=NULL,y="AUCell scores")&mytheme&
  stat_summary(fun = "mean", geom = "point",position = position_dodge(0.9)) &
  stat_summary(fun.data = "mean_sd", geom = "errorbar", width = .15,
               size = 0.3,position = position_dodge(0.9))&
  theme(legend.position = "none")

