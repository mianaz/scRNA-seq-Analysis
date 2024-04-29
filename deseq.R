# for gsea
library(org.Hs.eg.db)
library(clusterProfiler)

# for DE
library(DESeq2)

# read counts data
counts <- read.csv("~/Downloads/lncap_lascpc_rawcounts.csv", header = T, row.names = 1)
head(counts)

# generate coldata
coldata <- data.frame(cell_line=c(rep("LNCaP", 3), rep("LASCPC", 3)))
# convert to factor
coldata$cell_line <- as.factor(coldata$cell_line)
rownames(coldata) <- colnames(counts) 
# double check rownames of coldata matches colnames of counts                     
all(rownames(coldata) == colnames(counts))

# create deseq dataset
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~cell_line)

# pre-filtering(optional)
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds_filter <- dds[keep,]

# perform DE
#dds <- DESeq(dds)
#res <- results(dds, name="lncap_vs_lascpc")
#res
dds_filter <- DESeq(dds_filter)
res <- results(dds_filter)
res
resultsNames(dds_filter)
# LFC shrinkage
resLFC <- lfcShrink(dds_filter, coef="cell_line_LNCaP_vs_LASCPC", type="apeglm")
# order by p-value
resOrdered <- resLFC[order(resLFC$pvalue),]
write.csv(resOrdered, "deseq_res.csv")

# perform enrichment analysis

# convert ensembl to entrez
entrez <- bitr(rownames(resOrdered), fromType="ENSEMBL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")

# prepare genelist
genelist <- resOrdered$log2FoldChange
names(genelist) <- rownames(resOrdered)
# filter unmapped genes
genelist <- genelist[names(genelist) %in% entrez[,1]]
# conver mapped genes to to entrezid
names(genelist) <- entrez[match(names(genelist), entrez[,1]), 2]
genelist <- sort(genelist, decreasing=T)

# kegg gsea
kegg_gse <- gseKEGG(genelist, organism="hsa", eps=0, seed=2023)
kegg_gse_res <- kegg_gse@result # note that we consider pathways with NES>0 to be overexpressed/activated; NES<0 to be repressed.
saveRDS(kegg_res, "kegg_gse_res.rds")

colnames(kegg_res)

# kegg enrichment
# select those with fold change > 2
de <- names(genelist)[abs(genelist) > 2]
kegg_en <- enrichKEGG(de, organism = "hsa")
kegg_en_res <- kegg_en@result
saveRDS(kegg_en@result, "kegg_en_res.rds")

