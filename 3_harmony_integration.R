# load required libraries
library(Seurat)
library(ggplot2)
library(scDblFinder)
library(BiocParallel)
library(dplyr)
library(cowplot)
library(harmony)
wd="/Users/miana/Desktop/Fall 2023/scRNA course/Final Project"
setwd(wd)
set.seed(2023)

wong_cohort <- readRDS("wong_cohort_annotated.rds")

t1 <- Sys.time()
# Harmony
wong_cohort <- RunHarmony(wong_cohort, group.by.vars = "sample_id") %>% 
  FindNeighbors(dims = 1:30, reduction="harmony") %>%
  RunUMAP(dims = 1:30, reduction="harmony") %>%
  FindClusters()

t2 <- Sys.time()

t2-t1

DimPlot(wong_cohort, group.by=c("celltype_manual", "sample_id")) 
ggsave("wong_harmony.png", width=12, height=5)
saveRDS(wong_cohort, "wong_harmony.rds")