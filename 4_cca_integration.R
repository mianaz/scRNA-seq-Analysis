# read in dataframe
library(Seurat)
set.seed(2023)

# CCA integration
wong_cohort <- readRDS("~/datasets/Wong/wong_cohort.rds")
# convert to v5 assay
wong_cohort[["RNA"]] <- as(wong_cohort[["RNA"]], "Assay5")
wong_cohort[["RNA"]] <- split(wong_cohort[["RNA"]], f=wong_cohort$sample_id)
wong_cohort <- SCTransform(wong_cohort)
wong_cohort <- RunPCA(wong_cohort)
wong_cohort <- RunUMAP(wong_cohort, dims = 1:30)
DimPlot(wong_cohort, reduction = "umap", group.by = c("celltype_manual", "sample_id"))
ggsave("wong_no_integration_sct.png", width = 12, height=5)
wong_cohort <- IntegrateLayers(object = wong_cohort, method = CCAIntegration, normalization.method = "SCT", 
                               orig.reduction = "pca", new.reduction = "integrated.cca")

wong_cohort <- FindNeighbors(wong_cohort, reduction = "integrated.cca", dims = 1:30)
wong_cohort <- FindClusters(wong_cohort, resolution = 1)
wong_cohort <- RunUMAP(wong_cohort, dims = 1:30, reduction = "integrated.cca")
DimPlot(wong_cohort, reduction="umap", group.by=c("celltype_manual", "sample_id"))
ggsave("cca_integration_sct.png", width = 12, height=5)
saveRDS(wong_cohort, "wong_cca_sct_integrated.rds")