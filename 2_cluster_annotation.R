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

wong_cohort <- readRDS("wong_merged_filtered.rds")

# unsupervised clustering
wong_cohort <- SCTransform(wong_cohort, vars.to.regress = c("percent.mt")) %>%
  RunPCA(dims=1:50) %>%
  FindNeighbors(dims = 1:30) %>%
  RunUMAP(dims = 1:30, reduction.key = "umap.unintegrated") %>%
  FindClusters(cluster.name="unintegrated_clusters")

# manual annotation
wong.markers <- FindAllMarkers(wong_cohort)
wong.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> wong.top10
#DoHeatmap(wong_cohort, features = top10$gene) + NoLegend()
#save_plot(filename="heatmap.pdf", plot=last_plot(), base_height=20, base_width=40)
wong.top10 %>% print(n=+Inf)
DotPlot(wong_cohort, features=c("EPCAM", "KRT8", "KRT15", "VWF","SELE", "ACTA2","MYH11","FBLN1", "PTPRC","CD3E", "CD4", "CD8A", "CD1E", "S100A9", "CSF3R","CD19", "CD79A", "PLP1"))
ggsave("marker_exp.png", width=15, height=5)
Idents(wong_cohort) <- wong_cohort$seurat_clusters

# annotate clusters based on top markers
wong_cohort <- RenameIdents(wong_cohort, '0'="Endo", '1'="CD8 T" , '2'="CD4 T", '3'="MF/DC", '4'="Club", 
                            '5'="CD4 T", '6'="LE", '7'="CD8 T", '8'="CD8 T", '9'="CD8 T", '10'="LE", '11'="BE", 
                            '12'="SM", '13'="LE", '14'="B", '15'="LE", '16'="B", '17'="B", 
                            '18'="LE", '19'="NF", '20'="CD8 T", '21'="Mast", '22'="LE", '23'="FB", '24'="Endo", 
                            '25'="CD4 T", '26'="Neuron", '27'="Other")

wong_cohort$celltype_manual <- Idents(wong_cohort)

DimPlot(wong_cohort, group.by=c("celltype_manual", "sample_id")) 
ggsave("wong_unintegrated.png", width=12, height=5)
saveRDS(wong_cohort, "wong_cohort_annotated.rds")
