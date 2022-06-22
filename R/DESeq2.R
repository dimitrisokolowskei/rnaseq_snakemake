# Import Libraries ----
library(tidyverse)
library(tximport)
library(EnsDb.Hsapiens.v86)
library(DESeq2)
library(rhdf5)
library(apeglm)
library(gplots)
library(RColorBrewer)
library(limma)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi) 
library(DOSE)
library(enrichplot)
library(ggnewscale)

# Import expression matrices and metadata ----
targets <- read_tsv("metadata.txt") # Study metadata
path_transcripts <- file.path(targets$sample, "abundance.tsv") # Path to samples

import_data <- function (reference, path) {
  Tx <- transcripts(reference, columns=c("tx_id", "gene_name"))
  Tx <- as_tibble(Tx)
  Tx <- dplyr::rename(Tx, target_id = tx_id)
  Tx <- dplyr::select(Tx, "target_id", "gene_name")
  Txi_gene <- tximport(path, 
                       type = "kallisto", 
                       tx2gene = Tx, 
                       txOut = FALSE, #determines whether your data represented at transcript or gene level
                       countsFromAbundance = "lengthScaledTPM",
                       ignoreTxVersion = TRUE)
}

txi <- import_data(EnsDb.Hsapiens.v86, path_transcripts)

# Count Matrix Input ----
sampleLabels <- targets$sample
colnames(txi$counts) <- sampleLabels

dds <- DESeqDataSetFromTximport(txi,
                                colData = targets,
                                design = ~ group)

# Pre-filtering & Note factor levels ----
keep <- rowSums(counts(dds)) >= 10 # Put minium read lenght
dds <- dds[keep, ]

dds$group <- relevel(dds$group, ref = "healthy") # Put the name of your metadata label control

# Differential Expression Analysis ----
dds <- DESeq(dds)
res <- results(dds, contrast = c("group", "healthy", "treated"))

# Log Fold Change shrinkage (LFC) ----
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="group_treated_vs_healthy", type="apeglm")
resLFC

# Pvalue Ordering ----
resOrdered <- res[order(res$pvalue),]
summary(resOrdered)

# DGE Normalization
summary(results(dds, alpha=0.05))
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)

# MA plot ----
plotMA(dds, ylim=c(-2,2))
plotMA(resLFC, ylim=c(-2,2))
# Volcano Plot ----

diff_genes <- res
diff_genes <- data.frame(diff_genes)
diff_genes$diffexpressed <- "Down"
diff_genes$diffexpressed[diff_genes$log2FoldChange>0.58 & diff_genes$padj<0.01] <- "Up"
diff_genes$diffexpressed[diff_genes$log2FoldChange<0.58 & diff_genes$padj<0.01] <- "Down"

diff_genes$deLabel <- NA
p <- ggplot(diff_genes, aes(x=log2FoldChange, y=-log10(pvalue)), col=diff_genes, label=deLabel) +
  geom_point() +
  theme_minimal() 

mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
p2 <- p + scale_colour_manual(values = mycolors)

 scale_color_manual(values = c("blue","black","red")) +
   theme(text=element_text(size=20))



# Principal Component Analysis ----
group <- factor(targets$group)

pca.res <- prcomp(t(normalized_counts), scale.=F, retx=T)
pc.var <- pca.res$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
pc.per <- round(pc.var/sum(pc.var)*100, 1) 
pca.res.df <- as_tibble(pca.res$x)
pca.plot <- ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, label=sampleLabels, color = group) +
  geom_point(size=4) +
  stat_ellipse() +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="PCA plot") +
  coord_fixed() +
  theme_bw()

# Heatmap ----
myheatcolors <- rev(brewer.pal(name="RdBu", n=11))
clustRows <- hclust(as.dist(1-cor(t(normalized_counts), method="pearson")), method="complete") #cluster rows by pearson correlation
clustColumns <- hclust(as.dist(1-cor(normalized_counts, method="spearman")), method="complete")
module.assign <- cutree(clustRows, k=2)
module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(module.assign)] 
heatmap.2(normalized_counts, 
          Rowv=as.dendrogram(clustRows), 
          Colv=as.dendrogram(clustColumns),
          RowSideColors=module.color,
          col=myheatcolors, scale='row', keysize = 1, labRow=rownames(normalized_counts),
          density.info="none", trace="none",  
          cexRow=1, cexCol=1, margins=c(8,20)) 




# Gene Enrichment Analysis ----
res.enrich <- res

GO_results <- enrichGO(gene = res.enrich, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP")
as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 20)) 
dotplot(GO_results, showCategory = 20)


edo <- pairwise_termsim(GO_results)
p1 <- emapplot(edo, layout="kk")
