# load required libraries
library(Seurat)
library(ggplot2)
library(scDblFinder)
library(BiocParallel)
library(dplyr)
library(cowplot)
wd="/Users/miana/Desktop/Fall 2023/scRNA course/Final Project"
setwd(wd)
set.seed(2023)

# ICC/IDC: Wong cohort (7 patients, paired normal/tumor)
wong <- readRDS("GSE185344_PH_scRNA.final.rds")
wong_cohort <- wong[['rawobj']]
rm(wong)

# examine raw counts
VlnPlot(wong_cohort, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, alpha=0.05, pt.size = 0.2)
ggsave("wong_preQC.png", width=15, height=5)
table(wong_cohort$orig.ident)
wong_cohort$sample_id <- wong_cohort$orig.ident
wong_cohort$group <- sub(".*_.*_","",wong_cohort$sample_id)
wong_cohort$patient_id <- gsub("_Tumor$|_Benign$","\\1",wong_cohort$sample_id)

# we only need RNA assay for now
wong_cohort <- DietSeurat(wong_cohort, assay="RNA")
wong_cohort[["percent.mt"]] <- PercentageFeatureSet(wong_cohort, pattern = "^MT-")
wong_cohort[["log1p_count"]] <- log1p(wong_cohort$nCount_RNA)
wong_cohort[["log1p_feature"]] <- log1p(wong_cohort$nFeature_RNA)

# split samples into separate objects
wong_list <- SplitObject(wong_cohort, split.by="sample_id")
rm(wong_cohort)

# filter low quality cells
is.outlier <- function(obj, metric, nmads){
  M = obj@meta.data[,metric]
  outlier = (M < median(M) - nmads * mad(M)) | (
    median(M) + nmads *mad(M) < M
  )
  return(outlier)
}

wong_list <- lapply(wong_list, function(x){
  x$outlier <- (is.outlier(x, "log1p_count", 5)|is.outlier(x, "log1p_feature", 5)|is.outlier(x, "percent.mt", 3))
  x <- subset(x, outlier==F)
  })

# additionally filter out cells with more than 20% mt genes
wong_list <- lapply(wong_list, function(x) subset(x, percent.mt<20))

# doublet removal
wong_list_sce <- lapply(wong_list, as.SingleCellExperiment)
wong_list_sce <- lapply(wong_list_sce, scDblFinder)
for(i in names(wong_list_sce)){
  wong_list[[i]] <- AddMetaData(wong_list[[i]], wong_list_sce[[i]]$scDblFinder.class, col.name="scDblFinder.class")
}
wong_list <- lapply(wong_list, function(x) subset(x, scDblFinder.class=="singlet"))

rm(wong_list_sce)

# merge list
wong_merged <- merge(wong_list[[1]], wong_list[2:length(wong_list)])
VlnPlot(wong_merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, alpha=0.05, pt.size = 0.2)
ggsave("wong_postQC.png", width=15, height=5)

# save qc'ed list and merged obj
saveRDS(wong_merged, "wong_filtered.rds")






