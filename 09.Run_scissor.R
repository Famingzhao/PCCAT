rm(list=ls())
library(Scissor)
library(Seurat)
library(qs)
library(dplyr)
library(ggsci)
library(future)
library(ggplot2)
library(patchwork)
library(readr)
# change the current plan to access parallelization
plan("multisession", workers = 2)
options(future.globals.maxSize = 300 * 1024^3)
Scissor_fmz <- function (bulk_dataset, sc_dataset, phenotype, tag = NULL, alpha = NULL, 
                         cutoff = 0.2, family = c("gaussian", "binomial", "cox"), 
                         Save_file = "Scissor_inputs.RData", Load_file = NULL,netw = NULL) 
{
  library(Seurat)
  library(preprocessCore)
  library(qs)
  if (is.null(Load_file)) {
    common <- intersect(rownames(bulk_dataset), rownames(sc_dataset))
    if (length(common) == 0) {
      stop("There is no common genes between the given single-cell and bulk samples.")
    }
    if (class(sc_dataset) == "Seurat") {
      sc_exprs <- as.matrix(sc_dataset@assays$RNA@data)
      
      if(is.null(netw)){
        network <- as.matrix(sc_dataset@graphs$RNA_snn)
      }
    }
    else {
      sc_exprs <- as.matrix(sc_dataset)
      Seurat_tmp <- CreateSeuratObject(sc_dataset)
      Seurat_tmp <- FindVariableFeatures(Seurat_tmp, selection.method = "vst", 
                                         verbose = F)
      Seurat_tmp <- ScaleData(Seurat_tmp, verbose = F)
      Seurat_tmp <- RunPCA(Seurat_tmp, features = VariableFeatures(Seurat_tmp), 
                           verbose = F)
      Seurat_tmp <- FindNeighbors(Seurat_tmp, dims = 1:10, 
                                  verbose = F)
      network <- as.matrix(Seurat_tmp@graphs$RNA_snn)
    }
    if(is.null(netw)){
      diag(network) <- 0
      network[which(network != 0)] <- 1
    }
    dataset0 <- cbind(bulk_dataset[common, ], sc_exprs[common, 
    ])
    
    dataset1 <- normalize.quantiles(as.matrix(dataset0))
    rownames(dataset1) <- rownames(dataset0)
    colnames(dataset1) <- colnames(dataset0)
    Expression_bulk <- dataset1[, 1:ncol(bulk_dataset)]
    Expression_cell <- dataset1[, (ncol(bulk_dataset) + 
                                     1):ncol(dataset1)]
    X <- cor(Expression_bulk, Expression_cell)
    quality_check <- quantile(X)
    print("|**************************************************|")
    print("Performing quality-check for the correlations")
    print("The five-number summary of correlations:")
    print(quality_check)
    print("|**************************************************|")
    if (quality_check[3] < 0.01) {
      warning("The median correlation between the single-cell and bulk samples is relatively low.")
    }
    if (family == "binomial") {
      Y <- as.numeric(phenotype)
      z <- table(Y)
      if (length(z) != length(tag)) {
        stop("The length differs between tags and phenotypes. Please check Scissor inputs and selected regression type.")
      }
      else {
        print(sprintf("Current phenotype contains %d %s and %d %s samples.", 
                      z[1], tag[1], z[2], tag[2]))
        print("Perform logistic regression on the given phenotypes:")
      }
    }
    if (family == "gaussian") {
      Y <- as.numeric(phenotype)
      z <- table(Y)
      if (length(z) != length(tag)) {
        stop("The length differs between tags and phenotypes. Please check Scissor inputs and selected regression type.")
      }
      else {
        tmp <- paste(z, tag)
        print(paste0("Current phenotype contains ", 
                     paste(tmp[1:(length(z) - 1)], collapse = ", "), 
                     ", and ", tmp[length(z)], " samples."))
        print("Perform linear regression on the given phenotypes:")
      }
    }
    if (family == "cox") {
      Y <- as.matrix(phenotype)
      if (ncol(Y) != 2) {
        stop("The size of survival data is wrong. Please check Scissor inputs and selected regression type.")
      }
      else {
        print("Perform cox regression on the given clinical outcomes:")
      }
    }
    if(is.null(netw)){
      save(X, Y,Expression_bulk, Expression_cell, 
           file = Save_file)
      qsave(network, file = "network.qs")
    }else{
      save(X, Y, Expression_bulk, Expression_cell, 
           file = Save_file)
    }
  }
  else {
    load(Load_file)
  }
  
  if(!is.null(netw)){
    network <- qread(file = paste0(netw,"network.qs"))
  }
  
  if (is.null(alpha)) {
    alpha <- c(0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 
               0.6, 0.7, 0.8, 0.9)
  }
  for (i in 1:length(alpha)) {
    set.seed(123)
    fit0 <- APML1(X, Y, family = family, penalty = "Net", 
                  alpha = alpha[i], Omega = network, nlambda = 100, 
                  nfolds = min(10, nrow(X)))
    fit1 <- APML1(X, Y, family = family, penalty = "Net", 
                  alpha = alpha[i], Omega = network, lambda = fit0$lambda.min)
    if (family == "binomial") {
      Coefs <- as.numeric(fit1$Beta[2:(ncol(X) + 1)])
    }
    else {
      Coefs <- as.numeric(fit1$Beta)
    }
    Cell1 <- colnames(X)[which(Coefs > 0)]
    Cell2 <- colnames(X)[which(Coefs < 0)]
    percentage <- (length(Cell1) + length(Cell2))/ncol(X)
    print(sprintf("alpha = %s", alpha[i]))
    print(sprintf("Scissor identified %d Scissor+ cells and %d Scissor- cells.", 
                  length(Cell1), length(Cell2)))
    print(sprintf("The percentage of selected cell is: %s%%", 
                  formatC(percentage * 100, format = "f", digits = 3)))
    if (percentage < cutoff) {
      break
    }
    cat("\n")
  }
  print("|**************************************************|")
  return(list(para = list(alpha = alpha[i], lambda = fit0$lambda.min, 
                          family = family), Coefs = Coefs, Scissor_pos = Cell1, 
              Scissor_neg = Cell2))
}

### 1. load scRNA-seq data
m.seurat <- qread(file = "./Step4.SuperCell/Step4.m.seurat.BBKNN.qs")

### 2. load Bulk data
load("./Rawdata/TCGA_PRAD_Expression.Rdata")
tcga_tpm[1:10,1:10]
tcga_tpm = 2^tcga_tpm-1
colnames(tcga_tpm) == tcga_phe$barcode
 
#### 3. Run Scissor
### 3.1 binomial
phenotype <- ifelse(tcga_phe$group == "Adj",0,1)
tag <- c('Adj','Tumor')
binomial_res <- Scissor_fmz(bulk_dataset = tcga_tpm,
                            sc_dataset =  m.seurat, 
                            phenotype =  phenotype,  
                            tag = tag, 
                            cutoff = 0.2,
                            alpha = c(0.005, 0.01, 0.02, 0.03, 0.04,0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 
                                      0.6, 0.7, 0.8, 0.9),
                            family = "binomial", netw = "./Step7.Scissor_Res/",
                            Save_file = './Step7.Scissor_Res/TCGA_Scissor_PRAD_Tumor_vs_Adj.RData')

### 3.2 cox Survival
phenotype_PFI <- select(tcga_phe[tcga_phe$group == "Tumor",],PFI.time,PFI)
colnames(phenotype_PFI) <- c("time", "status")
tcga_tpm_tumor = tcga_tpm[,row.names(phenotype_PFI)]
infos_PFI <- Scissor_fmz(bulk_dataset = tcga_tpm_tumor,
                         sc_dataset =  m.seurat, 
                         phenotype =  phenotype_PFI,  
                         cutoff = 0.2,
                         alpha = c(0.005, 0.01, 0.02, 0.03, 0.04,0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 
                                   0.6, 0.7, 0.8, 0.9),
                         netw = "./Step7.Scissor_Res/", 
                         family = "cox",
                         Save_file = './Step7.Scissor_Res/TCGA_Scissor_PFI_survival.RData')

### 3.3 gaussian
(load("./PMID34857732_PCaProfiler.Rdata"))
pcatlas.data[1:10,1:10]
table(colnames(pcatlas.data) == pcatlas.phe.data$ID)
row.names(pcatlas.phe.data) = pcatlas.phe.data$ID
pcatlas.data = 2^pcatlas.data-1
table(pcatlas.phe.data$Sample.Type)
pcatlas.phe.data$GS_level2 = recode(pcatlas.phe.data$Sample.Type,
                                    "NORMAL"=0,
                                    "Adjacent"=1,
                                    "GS<7"=2,
                                    "GS=7"=3,
                                    "GS>7"=4,
                                    "CRPC"=5,
                                    "NEPC"=6)

pcatlas.phe = pcatlas.phe.data[!is.na(pcatlas.phe.data$GS_level),]
table(pcatlas.phe$GS_level2)
pcatlas.data_tumor = pcatlas.data[,row.names(pcatlas.phe)]

phenotype <- pcatlas.phe$GS_level2
tag <- c("NORMAL"=0,
         "Adjacent"=1,
         "GS<7"=2,
         "GS=7"=3,
         "GS>7"=4,
         "CRPC"=5,
         "NEPC"=6)
infos_pcatlas <- Scissor_fmz(bulk_dataset = pcatlas.data_tumor,
                             sc_dataset =  m.seurat, 
                             phenotype =  phenotype,  
                             tag = tag, 
                             cutoff = 0.2,
                             alpha = c(0.005, 0.01, 0.02, 0.03, 0.04,0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 
                                      0.6, 0.7, 0.8, 0.9),
                             family = "gaussian", 
                             netw = "./Step7.Scissor_Res/",
                             Save_file = './Step7.Scissor_Res/PCaProfiler_Scissor_Gleason_level.RData')