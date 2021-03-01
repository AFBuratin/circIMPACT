## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval = FALSE-------------------------------------------------------------
#  devtools::install_github("AFBuratin/circIMPACT")

## ----setup, warning = FALSE, message = FALSE----------------------------------
library(circIMPACT)
library(DESeq2)
library(data.table)
library(plyr)
library(Rtsne)
library(ggplot2)
library(ggrepel)
library(plotly)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(knitr)
library(kableExtra)
library(formattable)
library(htmltools)
library(sparkline)
library(tidyverse)
library(RColorBrewer)
library(purrr)
library(magrittr)
library(webshot)

## -----------------------------------------------------------------------------
data("circularData")
data("count.matrix")
data("meta")
data("coldata.df")
## define sample/condition colors
intgroup.dt <- meta[, .(.N), by = .(sample, 
                                      condition)][order(sample), 
                                                  .(sample), by = condition]
samples.per.condition <- data.frame(intgroup.dt[, .N, by = .(condition)], row.names = "condition")
  
n.conditions <- nrow(samples.per.condition)
hues <- circIMPACT::gg_color_hue(n.conditions)
for(i in 1:n.conditions){
      n.hues <- samples.per.condition[i, "N"] + 2
      col.hues <- colorRampPalette(colors = c(hues[i], 
                                              "white"))(n.hues)[1:(n.hues-2)]
      
      intgroup.dt[condition == rownames(samples.per.condition)[i], `:=`(color = col.hues,
                                                                        hue = hues[i])]
    }
## create a deseq dataset object 
dds.circular <- suppressMessages(DESeqDataSetFromMatrix(countData = ceiling(circularData[, coldata.df$sample[order(coldata.df$condition)]]),
                                   colData = coldata.df[order(coldata.df$condition),],
                                   design = ~ condition))
dds.circular <- suppressMessages(estimateSizeFactors(dds.circular))
sf <- sizeFactors(dds.circular)

## -----------------------------------------------------------------------------
data.filt <- circularData[rowSums(circularData >= 5) >= 2,]
dds.filt.expr <- suppressMessages(DESeqDataSetFromMatrix(countData = ceiling(data.filt[,coldata.df$sample[order(coldata.df$condition)]]),
                                   colData = coldata.df[order(coldata.df$condition),],
                                   design = ~ condition))
dds.filt.expr <- suppressMessages(estimateSizeFactors(dds.filt.expr))
sf.filt <- sizeFactors(dds.filt.expr)
circNormDeseq <- counts(dds.filt.expr, normalized = T)

## ----message=FALSE, warning=FALSE, include=TRUE-------------------------------

circIMPACT <- marker.selection(dat = circNormDeseq, dds = dds.filt.expr, sf = sf.filt, p.cutoff = 0.1, lfc.cutoff = 1, 
                                 method.d = "spearman", method.c = "complete", k = 2, median = TRUE)

## -----------------------------------------------------------------------------
circMark <- circIMPACT$circ.targetIDS[5]
circMark_group.df <- circIMPACT$group.df[circIMPACT$group.df$circ_id==circMark,]
circMark_group.df$counts <- merge(circMark_group.df, reshape2::melt(circNormDeseq[circMark,]), by.x = "sample_id", by.y = "row.names")[,"value"]
mu <- ddply(circMark_group.df, "group", summarise, Mean=mean(counts), Median=median(counts), Variance=var(counts))

p <- ggplot(circMark_group.df, aes(x=counts, color=group, fill=group)) +
  geom_density(alpha=0.3) + 
  geom_vline(data=mu, aes(xintercept=Median, color=group),
             linetype="dashed") +
  geom_text(data=mu, aes(x=Median[group=="g1"] - 0.55, 
                         label=paste0("Median:", round(Median[group=="g1"], 3), " Variance:", round(Variance[group=="g1"]), 3), y=0.15),
            colour="black", angle=90, text=element_text(size=9)) +
  geom_text(data=mu, aes(x=Median[group=="g2"] - 0.55, 
                       label=paste0("Median:", round(Median[group=="g2"], 3), " Variance:", round(Variance[group=="g2"]), 3 ), y=0.15), 
          colour="black", angle=90, text=element_text(size=11)) +  scale_fill_brewer(palette="Dark2") + 
  scale_color_brewer(palette="Dark2") + 
  labs(title=paste0("circMarker (", circMark, ")", " counts density curve"), x = "Normalized read counts", y = "Density") + 
  theme_classic()
p


## ----message=FALSE------------------------------------------------------------
markers.circrnas <- circIMPACT$circ.targetIDS
mat.filt.mark <- circularData[markers.circrnas, ]

dds.filt.mark <- DESeqDataSetFromMatrix(countData = ceiling(mat.filt.mark[,coldata.df$sample]),
                                   colData = coldata.df,
                                   design = ~ 1)
dds.filt.vst <- varianceStabilizingTransformation(dds.filt.mark, fitType = "local", blind = F)
norm.counts.filt <- assay(dds.filt.vst)

## -----------------------------------------------------------------------------
dt <- norm.counts.filt


## -----------------------------------------------------------------------------
## Rtsne function may take some minutes to complete...
set.seed(9)
mydist <- dist(t(norm.counts.filt))
## t-SNE representation
# set a perplexity parameter consistent with the number of samples
tsne_data <- Rtsne(mydist, pca = F, perplexity=1, max_iter=5000)

## getting the two dimension matrix
d_tsne_1 = as.data.frame(tsne_data$Y)
rownames(d_tsne_1) <- colnames(norm.counts.filt)


## ----clustering, fig.cap = "t-SNE dimensionality reduction representation. K-means and hierarchical clustering are compared."----

## keeping original data
d_tsne_1_original=d_tsne_1

## Creating k-means clustering model, and assigning the result to the data used to create the tsne
fit_cluster_kmeans=kmeans(scale(d_tsne_1), 2)
d_tsne_1_original$cl_kmeans = factor(fit_cluster_kmeans$cluster)

## Creating hierarchical cluster model, and assigning the result to the data used to create the tsne
fit_cluster_hierarchical=hclust(dist(scale(d_tsne_1)))

## setting 2 clusters as output
d_tsne_1_original$cl_hierarchical = factor(cutree(fit_cluster_hierarchical, k=2))

# Plotting the cluster models onto t-SNE output

plot_cluster=function(data, var_cluster, palette)
{
  ggplot(data, aes_string(x="V1", y="V2", color=var_cluster)) +
  geom_point(size=3) +
  guides(colour=guide_legend(override.aes=list(size=3))) +
  geom_text_repel(aes(label = rownames(data)), 
                  hjust = 0.5, vjust = -1) +
  xlab("") + ylab("") +
  ggtitle("") +
  theme_light(base_size=11) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        legend.direction = "horizontal", 
        legend.position = "bottom",
        legend.box = "horizontal") + 
    scale_colour_brewer(palette = palette) 
}


plot_k=plot_cluster(d_tsne_1_original, "cl_kmeans", "Dark2")
plot_h=plot_cluster(d_tsne_1_original, "cl_hierarchical", "Set1")

## and finally: putting the plots side by side with gridExtra lib...
library(gridExtra)
grid.arrange(plot_k, plot_h,  ncol=2)


## -----------------------------------------------------------------------------
pca <- prcomp(x = t(norm.counts.filt), center = T)
d <- data.frame(pca$x[rownames(coldata.df), c("PC1", "PC2")], coldata.df)
PC1.var <- summary(pca)$importance["Proportion of Variance", 1]
PC2.var <- summary(pca)$importance["Proportion of Variance", 2]
g1 <- ggplot(data = d, 
       mapping = aes(x = PC1, y = PC2)) +
    geom_point(size = 4) +
    coord_fixed(ratio = 1) +
    xlab(paste0("PC1: ", percent(PC1.var))) +
    ylab(paste0("PC2: ", percent(PC2.var))) +
    theme_classic() + 
    theme(legend.position = "bottom", 
          plot.title = element_text(hjust = .5))
library(factoextra)
#### compute contribution 
contrib <- function(ind.coord, comp.sdev, n.ind){
  100*(1/n.ind)*ind.coord^2/comp.sdev^2
}
ind.contrib <- t(apply(pca$x, 1, contrib, 
                       pca$sdev, nrow(pca$x)))
g3 <- fviz_contrib(pca, choice="var", axes = 1:3, top = 8)
library(ggpubr)

ggpubr::ggarrange(ggarrange(g1, g3, 
                            ncol = 2, 
                            labels = c("A", "B")),
                  nrow = 1)  

## -----------------------------------------------------------------------------

set.seed(201)
dds.filt.mark <- estimateSizeFactors(dds.filt.mark)
circNormDeseq <- counts(dds.filt.mark, normalized = T)

base_mean = log2(rowMeans(circNormDeseq)+0.001)
mat_scaled = t(apply(dt, 1, function(x) scale(x = x, center = T, scale = T)))
colnames(mat_scaled) <- colnames(dt)
cond = colData(dds.filt.expr)$condition
## choice of kmeans results as cluster of samples
clus = d_tsne_1_original$cl_kmeans
cond.colors <- unique(intgroup.dt$hue)
names(cond.colors) <- unique(intgroup.dt$condition)
ha = HeatmapAnnotation(df = data.frame(condition = cond, cluster = clus),
                       col = list(condition = cond.colors),
                       show_annotation_name = F,
                       annotation_legend_param = list(condition = list(nrow = 2, direction = "horizontal")))

mat.dend <- as.dendrogram(fit_cluster_hierarchical)
fit_cluster_kmeans$cluster  
ht <- Heatmap(mat_scaled, name = "expression", 
        # km = 2,
        # column_km = 2,
        column_order = names(fit_cluster_kmeans$cluster[order(fit_cluster_kmeans$cluster)]),
        col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
        top_annotation = ha, 
        # top_annotation_height = unit(4, "mm"),
        clustering_distance_columns = "euclidean",
        clustering_method_column = "complete",
        cluster_columns = F,
        clustering_distance_rows = "spearman",#"minkowski",
        clustering_method_rows = "ward.D2",
        cluster_rows = T,
        # row_dend_side = "right",
        # row_names_side = "left",
        show_row_names = T, 
        show_column_names = F, 
        width = unit(9, "cm"),
        show_row_dend = T,
        show_column_dend = T,
        # row_dend_reorder = TRUE,
        row_names_gp = gpar(fontsize = 5),
        heatmap_legend_param = list(direction = "horizontal")) +
Heatmap(base_mean, name = "log2(base mean)", show_row_names = F, width = unit(2, "mm"), col = inferno(255), show_column_names = F, row_names_gp = gpar(fontsize = 5), heatmap_legend_param = list(direction = "horizontal"))

draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")


## ----Filterlinear_data--------------------------------------------------------
## filter out genes low expressed 
min.count <- 20
min.col <- 5

filt.mat <- count.matrix[rowSums(count.matrix >= min.count) >= min.col, ]


## ----message=FALSE, warning=FALSE---------------------------------------------
#normalized circRNAs data 
circNormDeseq <- counts(dds.filt.expr, normalized = T) %>% as.data.frame()
circNormDeseq$circ_id <- rownames(circNormDeseq)

library(doParallel)
no_cores <- detectCores() - 1  
registerDoParallel(cores=no_cores)  

# gene_mark <- foreach::foreach(i=1:5, .combine = rbind) %dopar% {
# 
#   results.temp <- data.frame(geneexpression(circ_idofinterest = markers.circrnas[i], circRNAs = circNormDeseq, 
#                                        linearRNAs = filt.mat, colData = coldata.df, padj = 0.1, 
#                                        group = circIMPACT$group.df[circIMPACT$group.df$circ_id%in%markers.circrnas[i],],
#                                        covariates = NULL), circIMPACT = markers.circrnas[i])
# }

gene_mark_hipk3 <- data.frame(geneexpression(circ_idofinterest = "11:33286412-33287511", circRNAs = circNormDeseq, 
                                       linearRNAs = filt.mat, colData = coldata.df, padj = 0.1, 
                                       group = circIMPACT$group.df[circIMPACT$group.df$circ_id%in%"11:33286412-33287511",],
                                       covariates = NULL), circIMPACT = "11:33286412-33287511")


## -----------------------------------------------------------------------------

gene_mark <- as.data.table(gene_mark_hipk3)
gene_mark %>% dplyr::rename("Gene" = "gene_id", "logFC" = "log2FoldChange") %>% 
  arrange(padj) %>% 
  select(circIMPACT, Gene, logFC) %>% head(20) %>% 
  formattable::formattable(., align = c("c","c","c"), list(
          gene_id = formattable::formatter("span", style = ~ formattable::style(color = "grey", font.weight = "bold")),
          circIMPACT = formattable::formatter("span", style = ~ formattable::style(color = "grey", font.weight = "bold")),
          logFC = circIMPACT::color_tile3(digits = 3, n = 18, fun = "comma", palette = "PiYG")))

# knitr::kable(gene_mark %>% dplyr::group_by(circRNA_markers, n.degs) %>% 
# dplyr::summarise(DEGs = paste(sort(gene_id),collapse=", ")),
#       escape = F, align = "c", row.names = T, caption = "circRNA-DEGs assosiation") %>% kable_styling(c("striped"), full_width = T)
gene_mark[gene_mark$gene_id=="HPSE",]

