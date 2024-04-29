# liger integration (non-negative matrix factorization)
#install.packages('rliger')
library(rliger)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)

set.seed(2023)

wong_cohort <- readRDS("~/datasets/Wong/wong_cohort.rds")

# rna assay (as shown in tutorial)
DefaultAssay(wong_cohort) <- "RNA"
DietSeurat(wong_cohort, assays="RNA")
wong_cohort <- NormalizeData(wong_cohort)
wong_cohort <- FindVariableFeatures(wong_cohort)
wong_cohort <- ScaleData(wong_cohort, split.by="sample_id", do.center=FALSE)
wong_cohort <- RunOptimizeALS(wong_cohort, k=20, lambda=5, split.by="sample_id")
wong_cohort <- RunQuantileNorm(wong_cohort, split.by="sample_id")
wong_cohort <- FindNeighbors(wong_cohort, reduction="iNMF", dims=1:20)
wong_cohort <- FindClusters(wong_cohort, resolution=0.5)
wong_cohort <- RunUMAP(wong_cohort, dims=1:ncol(wong_cohort[["iNMF"]]), reduction="iNMF")
DimPlot(wong_cohort, group.by = c("celltype_manual", "sample_id"))
saveRDS(wong_cohort, "wong_inmf_integrated.rds")