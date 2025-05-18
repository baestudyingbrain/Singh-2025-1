######################################################
## scRNA-seq Preprocessing Pipeline (Seurat + Harmony)
## Project: ASPM_expression_scRNA
######################################################

suppressMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(sctransform)
  library(harmony)
})

dir <- "~/ASPM_project/"

######################################################
## Load and Filter Raw Seurat Objects
######################################################

obj_list <- readRDS(file.path(dir, "rds/obj_all_list.rds"))
obj_list <- lapply(obj_list, function(obj) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
  subset(obj, subset = nCount_RNA >= 5000 & percent.mt <= 5)
})
obj.srt <- merge(x = obj_list[[1]], y = obj_list[-1], add.cell.ids = names(obj_list))

######################################################
## Load and Merge Singlet (Doublet-Removed) Seurat Objects
######################################################

singlet_files <- list.files(file.path(dir, "rds"), pattern = "singlet", full.names = TRUE)
singlet_list <- lapply(singlet_files, readRDS)
obj.srt <- merge(x = singlet_list[[1]], y = singlet_list[-1], add.cell.ids = names(singlet_list))

######################################################
## Normalization, Feature Selection, Scaling (on all genes)
######################################################

obj.srt <- NormalizeData(obj.srt)
obj.srt <- FindVariableFeatures(obj.srt, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(obj.srt)
obj.srt <- ScaleData(obj.srt, features = all.genes)

######################################################
## PCA and Harmony Batch Correction
######################################################

obj.srt <- RunPCA(obj.srt, features = VariableFeatures(obj.srt), npcs = 10)
obj.srt <- RunHarmony(obj.srt, group.by.vars = "orig.ident")

######################################################
## Nearest Neighbors, Clustering, UMAP
######################################################

obj.srt <- FindNeighbors(obj.srt, reduction = "harmony", dims = 1:10)
obj.srt <- FindClusters(obj.srt, resolution = c(0.1, 0.2, 0.3, 0.4))
obj.srt <- RunUMAP(obj.srt, reduction = "harmony", dims = 1:10)

######################################################
## Save Final Processed Seurat Object
######################################################

saveRDS(obj.srt, file = file.path(dir, "rds", "obj_singlet_harmony_umap.rds"))


