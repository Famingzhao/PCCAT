rm(list=ls())
library(Seurat)
library(RColorBrewer)
library(tidyr)
library(purrr)
library(ggplot2)
library(patchwork)
library(dplyr)
library(readr)
library(future)
library(cowplot)
library(stringr)
library(ggpubr)
library(tidyverse)
library(ggsci)
library(spacexr)
# change the current plan to access parallelization
plan("multisession", workers = 4)
plan()

### 1.input
st.data.list <- read_rds("./Outdata/Step1.stRNAseq_data.rds")

### 2.SCT
st.data.list = lapply(st.data.list, function(st.data){
  st.data <- SCTransform(st.data, assay = "Spatial", verbose = FALSE)
  st.data <- st.data %>%  RunPCA( assay = "SCT", verbose = FALSE) %>% 
    FindNeighbors(reduction = "pca", dims = 1:30) %>% 
    RunUMAP(reduction = "pca", dims = 1:30)
  
  for (res in c(0.05, 0.2, 0.3, 0.5, 0.8, 1)) {
    # res=0.01
    print(res)
    st.data <- FindClusters(st.data, resolution = res, algorithm = 1, verbose = FALSE)
  }
  
  return(st.data)
})

#### 3. Run RCTD
RCTD.res.list <- imap(st.data.list, function(st.data, name){
  SAMPLE_ID = name
  CELLTYPE_ID = "ct.L2"
  SC_Reference = "./Step3.SCRef_downsample.qs"
  
  ### 3.1 RDCT list Run--------------------
  ## 3.1.1 Create the Reference object
  reference <- qread(SC_Reference)
  
  ## 3.1.2 Create the STpuck object
  counts = GetAssayData(st.data, slot = "counts")
  img <- GetTissueCoordinates(st.data)[,1:2]
  img$x = 540-img$x
  # head(img)
  colnames(img) = c("ycoord", "xcoord")
  coords = img[,c("xcoord",  "ycoord")]
  nUMI <- colSums(counts) # In this case, total counts per pixel is nUMI
  
  puck <- SpatialRNA(coords, counts, nUMI)
  ## Examine SpatialRNA object (optional)
  # print(dim(puck@counts)) # observe Digital Gene Expression matrix
  # hist(log(puck@nUMI,2)) # histogram of log_2 nUMI
  
  # print(head(puck@coords)) # start of coordinate data.frame
  barcodes <- colnames(puck@counts) # pixels to be used (a list of barcode names).
  
  # This list can be restricted if you want to crop the puck e.g.
  # puck <- restrict_puck(puck, barcodes) provides a basic plot of the nUMI of each pixel
  # on the plot:
  # plot_puck_continuous(puck, barcodes, puck@nUMI, ylimit = c(0,round(quantile(puck@nUMI,0.9))),
  #                      title ='plot of nUMI')
  qsave(puck,file = paste0("Step3.RCTD_Res/Step3.",SAMPLE_ID,"_Stpuck.qs"))
  
  ### 3.1.3 Creating RCTD Object and Run RCTD
  myRCTD <- create.RCTD(puck, reference, 
                        max_cores = 4, test_mode = FALSE) # here puck is the SpatialRNA object, and reference is the Reference object.
  
  myRCTD <- run.RCTD(myRCTD, doublet_mode = "full")
  qsave(myRCTD, file = paste0("Step3.RCTD_Res/Step3.",SAMPLE_ID,"_results_full.qs"))
  
  ### 3.1.4 results
  signature_exp <- as.data.frame(normalize_weights(myRCTD@results$weights)) %>% t() %>% as.data.frame()
  st.data = subset(st.data, barcode %in% colnames(signature_exp))
  
  st.data[["RCTD"]] <- SeuratObject::CreateAssayObject(counts = signature_exp)
  st.data <- SeuratObject::SetAssayData(st.data, slot = "scale.data",
                                        new.data = as.matrix(signature_exp), assay = "RCTD")
  
  #Vis
  p1 = SpatialFeaturePlot(st.data, features = row.names(st.data[["RCTD"]]), ncol = 5)
  ggplot2::ggsave(p1, filename = paste0("Step3.RCTD_Res/Step3.",SAMPLE_ID,"_AllCellTypes_RCTD.png"), width = 20, height = 16)
  print(paste0("End: ", SAMPLE_ID))
  return(st.data)
})
