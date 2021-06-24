## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval = FALSE-------------------------------------------------------------
#  devtools::install_github("AFBuratin/circIMPACT", force = TRUE)

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
library(randomForest)
library(DaMiRseq)
library(webshot)
library(ggpubr)
library(tidyr)
library(factoextra)
library(NbClust)

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
data.filt <- circularData[rowSums(circularData >= 10) >= 2,]
dds.filt.expr <- suppressMessages(DESeqDataSetFromMatrix(countData = ceiling(data.filt[,coldata.df$sample[order(coldata.df$condition)]]),
                                   colData = coldata.df[order(coldata.df$condition),],
                                   design = ~ condition))
dds.filt.expr <- suppressMessages(estimateSizeFactors(dds.filt.expr))
sf.filt <- sizeFactors(dds.filt.expr)
circNormDeseq <- counts(dds.filt.expr, normalized = T)

## ----message=FALSE, warning=FALSE, include=TRUE-------------------------------

circIMPACT <- circIMPACT::marker.selection(dat = circNormDeseq, 
                               dds = dds.filt.expr, 
                               sf = sf.filt, p.cutoff = 0.05, lfc.cutoff = 1)

circIMPACT_Kmeans <- marker.selection(dat = circNormDeseq,
                                      dds = dds.filt.expr,
                                      sf = sf.filt, p.cutoff = 0.05, lfc.cutoff = NULL, 
                                      method.d = "euclidean", method.c = "ward.D2",
                                      median = FALSE, choose.k = TRUE, index.m = "silhouette")

## -----------------------------------------------------------------------------
circMark <-"11:33286412-33287511"
circIMPACT_Kmeans$group.df[circIMPACT_Kmeans$group.df$circ_id==circMark,]
circMark_group.df <- circIMPACT$group.df[circIMPACT$group.df$circ_id==circMark,]
circMark_group.df$counts <- merge(circMark_group.df, reshape2::melt(circNormDeseq[circMark,]), by.x = "sample_id", by.y = "row.names")[,"value"]
mu <- ddply(circMark_group.df, "group", summarise, Mean=mean(counts), Median=median(counts), Variance=sd(counts))

p <- ggplot(circMark_group.df, aes(x=counts, color=group, fill=group)) +
  geom_density(alpha=0.3) + 
  ylim(c(0,0.02)) +
  geom_vline(data=mu, aes(xintercept=Median, color=group),
             linetype="dashed") +
  geom_text(data=mu, aes(x=Median[group=="g1"] - 1, 
                         label=paste0("Median:", round(Median[group=="g1"], 2), "/ Variance:", round(Variance[group=="g1"], 2)), y=0.012),
            colour="black", angle=90, text=element_text(size=12)) +
  geom_text(data=mu, aes(x=Median[group=="g2"] - 1, 
                       label=paste0("Median:", round(Median[group=="g2"], 3), "/ Variance:", round(as.numeric(Variance[group=="g2"]), 2)), y=0.012), 
          colour="black", angle=90, text=element_text(size=11)) +  scale_fill_brewer(palette="Dark2") + 
  scale_color_brewer(palette="Dark2") + 
  labs(title=paste0("circHIPK3 (", circMark, ")", " counts density curve"), x = "Normalized read counts", y = "Density") + 
  theme_classic() +
  theme(axis.text.x = element_text(size=18), 
        axis.text.y = element_text(size=18), 
        axis.title = element_text(size=20), 
        legend.text = element_text(size=15), 
        legend.title = element_text(size=15),
        plot.title = element_text(size = 18),
        plot.subtitle = element_text(size = 13, color = "darkslategrey"))
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
       mapping = aes(x = PC1, y = PC2, shape = condition)) +
    geom_point(size = 5) +
    coord_fixed(ratio = 1) +
    xlab(paste0("PC1: ", percent(PC1.var))) +
    ylab(paste0("PC2: ", percent(PC2.var))) +
    scale_shape_manual("", values = c(8,16)) +
    theme_classic() + 
    theme(legend.position = "bottom",
          legend.direction = "vertical",
          plot.title = element_text(hjust = .5),
          text = element_text(size=20),
          axis.text.x = element_text(size=20),
          axis.text.y = element_text(size=20),
          aspect.ratio = 1)


circRNArank = circIMPACT::circ_contrib(matrix_pc = pca, K = 25)


ggpubr::ggarrange(ggarrange(g1, circRNArank$plot,
                            ncol = 2,
                            labels = c("PCA", "")),
                  nrow = 1)


## -----------------------------------------------------------------------------
library(factoextra)
res = fviz_contrib(pca, choice ="var", axes = 1, top = 25)
top25 = res$data$name[order(res$data$contrib, decreasing=T)][1:25]
circAnnotation = read.csv("/media/Data/Li/circRNA_expression_per_sample.csv", header = T)

circAnnotation$circCHR =  sub(":.*", "", circAnnotation$circ_id)
circAnnotation$circstart = gsub(".*:(.*)\\-.*", "\\1", circAnnotation$circ_id)
circAnnotation$circend = sub(".*-", "\\1", circAnnotation$circ_id)
circAnnotation$circ_id <- paste0(circAnnotation$circCHR, ":", as.numeric(circAnnotation$circstart)-1, "-", circAnnotation$circend)

rownames(pca$rotation) = paste0("circ", circAnnotation$gene_names[match( rownames(norm.counts.filt), circAnnotation$circ_id)], "_", rownames(norm.counts.filt))

p = fviz_pca_biplot(pca, #select.ind = list(contrib = 5), 
               select.var = list(contrib = 25),
               ggtheme = theme_minimal(),
               title = "CircRNAs PCA loadings",
               #habillage=coldata.df$condition,
               # Individuals
                geom.ind = "point",
                #fill.ind = coldata.df$condition, 
               col.ind = "black",
              #pointshape = coldata.df$condition, 
               pointsize = 3,
                palette = "jco",
                addEllipses = FALSE,
               repel = TRUE,
                # Variables
               # alpha.var ="contrib", 
               col.var = "contrib",
                # gradient.cols = c("#C51B7D", "#E9A3C9", "#A1D76A", "#4D9221"), 
                gradient.cols = c("darkviolet", "deeppink2", "deeppink4"),
                legend.title = list(shape = "Condition", color = "Contrib")) +
                                    #alpha = "Contrib")) +
  labs(subtitle = "Top 25 variables with highest contribution of the event to PC1 and PC2",
       caption = "") +
  theme(axis.text.x = element_text(size=18), 
        axis.text.y = element_text(size=18), 
        axis.title = element_text(size=20), 
        legend.text = element_text(size=15), 
        legend.title = element_text(size=15),
        plot.title = element_text(size = 18),
        plot.subtitle = element_text(size = 13, color = "gray28"))
p

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

# library(doParallel)
# no_cores <- detectCores() - 5
# registerDoParallel(cores=no_cores)
# gene_mark <- foreach::foreach(i=1:25, .combine = rbind) %dopar% {
# 
#   results.temp <- data.frame(gene.expression(circ_idofinterest = top25[i], circRNAs = circNormDeseq,
#                                        linearRNAs = filt.mat,
#                                        colData = coldata.df, padj = 0.05,
#                                        group = circIMPACT$group.df[circIMPACT$group.df$circ_id%in%top25[i],],
#                                        covariates = NULL), circIMPACT = top25[i])
# }

gene_mark_hipk3 <- data.frame(gene.expression(circ_idofinterest = "11:33286412-33287511",
                                             circRNAs = circNormDeseq,
                                             linearRNAs = filt.mat, 
                                             colData = coldata.df, 
                                             padj = 0.1, 
                                             group = circIMPACT$group.df[circIMPACT$group.df$circ_id%in%"11:33286412-33287511",],
                                             covariates = NULL), 
                              circIMPACT = "11:33286412-33287511")


## -----------------------------------------------------------------------------

gene_mark_hip <- as.data.table(gene_mark_hipk3)
gene_mark_tab = gene_mark_hip %>% dplyr::group_by(circIMPACT, n.degs) %>% dplyr::mutate(UP=sum(log2FoldChange>0), DW=sum(log2FoldChange<0)) %>% dplyr::select(circIMPACT, n.degs, UP, DW)

gene_mark_hip %>% dplyr::rename("Gene" = "gene_id", "logFC" = "log2FoldChange") %>% 
  arrange(padj) %>% 
  select(circIMPACT, Gene, logFC) %>% head(20) %>% 
  formattable::formattable(., align = c("c","c","c"), list(
          gene_id = formattable::formatter("span", style = ~ formattable::style(color = "grey", font.weight = "bold")),
          circIMPACT = formattable::formatter("span", style = ~ formattable::style(color = "grey", font.weight = "bold")),
          logFC = circIMPACT::color_tile3(digits = 3, n = 18, fun = "comma", palette = "PiYG")))

# knitr::kable(gene_mark %>% dplyr::group_by(circRNA_markers, n.degs) %>% 
# dplyr::summarise(DEGs = paste(sort(gene_id),collapse=", ")),
#       escape = F, align = "c", row.names = T, caption = "circRNA-DEGs assosiation") %>% kable_styling(c("striped"), full_width = T)
gene_mark_hip[gene_mark_hip$gene_id=="HPSE",]
# head(gene_mark_hip)
# Make a basic volcano plot
gene_mark_hipk3$expression = ifelse(gene_mark_hipk3$padj <= 0.01 & abs(gene_mark_hipk3$log2FoldChange) >= 1, 
                     ifelse(gene_mark_hipk3$log2FoldChange >= 1 ,'Up','Down'),
                     'Stable')
p <- ggplot(data = gene_mark_hipk3, 
            aes(x = log2FoldChange, 
                y = -log10(padj), 
                colour=expression,
                label = gene_id)) +
  geom_point(alpha=0.4, size=3.5) +
  geom_text_repel(aes(label=ifelse((abs(log2FoldChange) >= 5)&(padj<=0.01), gene_id, "")), color = "black", size = 6) +
  scale_color_manual("",
                     guide = guide_legend(title.position = "top",
                                          ncol = 1,
                                          title.hjust = .5),
                     values = setNames(c("blue", "grey","red"),
                                       nm = c("Down", "Stable","Up"))) +
  # scale_color_manual(values=c("blue", "grey","red"))+
  xlim(c(-6.5, 6.5)) +
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = 2,lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",
       y="-log10 (adj.p-value)",
       title="")  +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank(),
        legend.text = element_text(size=20),
        text = element_text(size=20),
        axis.text.x = element_text(size=25),
        axis.text.y = element_text(size=25))

p
#ggplotly(p)


## -----------------------------------------------------------------------------
my_data = data.frame(circHIPK3 = as.vector(t(circNormDeseq["11:33286412-33287511",-7])), 
                     HPSE = filt.mat["HPSE",])
ggscatter(my_data, x = "circHIPK3", y = "HPSE", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "circHIPK3 expression", ylab = "QKI gene expression")

## -----------------------------------------------------------------------------
my_data = data.frame(circHIPK3 = as.vector(t(circNormDeseq["11:33286412-33287511",-7])), 
                     HPSE = filt.mat["HPSE",])
my_data$sample_id = rownames(my_data)
my_data$condition = meta$condition[match(my_data$sample_id, meta$sample)]
ggplot(reshape2::melt(my_data, id.var = c("sample_id","condition")),
       aes(x=condition, y=value, group=condition)) + 
  geom_boxplot(aes(fill=condition)) + 
  facet_grid( variable ~ ., scales='free') +
  scale_fill_manual(values = c("#1B9E77","#D95F02")) +
theme_bw() +
  xlab("") + 
  ylab("Normalized Expression") +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="none", 
        legend.title = element_blank(),
        legend.text = element_text(size=20),
        text = element_text(size=20),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(20),
        strip.text.x = element_text(size = 20))


## -----------------------------------------------------------------------------
library(doParallel)
library(caret)
library(dplyr)
library(tidyr)
library(data.table)
# no_cores <- detectCores() - 6
# registerDoParallel(cores=no_cores)
# gene.class <- foreach::foreach(i=1:25, .combine = list) %dopar% {
# 
#   results.temp <- gene_class(circ_idofinterest = top25[i], circRNAs = circNormDeseq,
#                                        linearRNAs = filt.mat, colData = coldata.df,
#                                        group = circIMPACT$group.df[circIMPACT$group.df$circ_id%in%top25[i],],
#                                        covariates = NULL, th.corr = 0.3)
# }

gene_class_hipk3 <- circIMPACT::gene.class(circ_idofinterest = "11:33286412-33287511", 
                                           circRNAs = circNormDeseq,
                                       linearRNAs = filt.mat,
                                       colData = coldata.df,
                                       group =
              circIMPACT$group.df[circIMPACT$group.df$circ_id%in%"11:33286412-33287511",],
                                       covariates = NULL,
              th.corr = 0.3)
#p
#ggplotly(p)

## -----------------------------------------------------------------------------
VI <- importance(gene_class_hipk3$RF)
VI.mat <- as.data.frame(VI)
r <- rownames(VI)
VI.mat$gene <- r
VI.mat <- VI.mat[order(VI.mat$MeanDecreaseAccuracy, decreasing = T),]
VI.mat <- as.data.table(VI.mat)


# VI.mat
VI.mat %>% mutate_at(1:4, round, 3) %>% head() %>% 
  mutate_if(is.numeric, function(x) {
    cell_spec(x, bold = T, 
              color = spec_color(x, end = 0.9),
              font_size = spec_font_size(x))
  }) %>%
  kable(escape = F, align = "c", row.names = F, caption = "Table of selected genes used for classification of subgroups defined by circRNAs variation. For each class of response variable there is a OOB error rate of classification. In the 4th column there is the importance of the variable in the growing of the the random forest") %>%
  kable_styling(c("striped"), full_width = F)


## ----eval=FALSE, include=FALSE------------------------------------------------
#  library(randomForestExplainer)
#  min_depth_frame <- min_depth_distribution(gene_class_hipk3$RF)
#  importance_frame <- measure_importance(gene_class_hipk3$RF)
#  
#  plot_multi_way_importance(importance_frame, x_measure = "accuracy_decrease", y_measure = "gini_decrease",
#                            size_measure = "no_of_nodes", no_of_labels = 5)

## ----message=FALSE, warning=FALSE---------------------------------------------
library(enrichplot)

library(DOSE)

library(clusterProfiler)

organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)


## ----message=FALSE, warning=FALSE, include=FALSE------------------------------
#subset gene symbol deregulated using the interesting circRNA marker as stratificator
geneList <- gene_mark_hipk3$log2FoldChange
names(geneList) <- gene_mark_hipk3$gene_id

# order gene list by foldchange
geneList = sort(geneList, decreasing = TRUE)
# names(geneList) <- gene_mark$Gene[gene_mark$circIMPACT==markers.circrnas[5]]
geneList <- geneList[abs(geneList)>=1.89]
# library(gprofiler2)
# gostres2 <- gost(query = names(geneList)[names(geneList)!="."], 
#                  organism = "hsapiens", ordered_query = TRUE, 
#                  multi_query = FALSE, significant = FALSE, exclude_iea = FALSE, 
#                  measure_underrepresentation = FALSE, evcodes = TRUE, 
#                  user_threshold = 0.05, correction_method = "g_SCS", 
#                  domain_scope = "annotated", custom_bg = NULL, 
#                  numeric_ns = "", sources = NULL)
# 
# p <- gostplot(gostres2, capped = FALSE, interactive = TRUE)
# p

gse <- gseGO(geneList=geneList, 
             ont ="ALL", 
             keyType = "SYMBOL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")

df_gse <- gse@result
df_gse = as.data.table(df_gse)
dim(df_gse)

## -----------------------------------------------------------------------------
enrichplot::dotplot(gse, showCategory=10, split=".sign", title = "Enriched GO") + facet_grid(.~.sign)


## -----------------------------------------------------------------------------
upsetplot(gse, showCategory=10)


## -----------------------------------------------------------------------------
kable(head(df_gse[like(core_enrichment,"HPSE/"),]), escape = F, align = "c", row.names = F, caption = "Enrichment of GO including HPSE gene") %>%
  kable_styling(c("striped"), full_width = F)

## ----message=FALSE, warning=FALSE, include=FALSE------------------------------
original_gene_list <- geneList
ids<-bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)

de_ids = ids[!duplicated(ids[c("SYMBOL")]),]

Eids<-bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENSEMBL", OrgDb=organism)

# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
df <- gene_mark_hipk3
df$SYMBOL <- df$gene_id
df <- df[which(df$SYMBOL%in%de_ids$SYMBOL),]

# Create a new column with the corresponding ENTREZ IDs
df$Y = de_ids$ENTREZID[match(df$SYMBOL,de_ids$SYMBOL)]

# Create a vector of the gene unuiverse
kegg_gene_list <- df$log2FoldChange

# Name vector with ENTREZ ids
names(kegg_gene_list) <- df$Y

# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

kegg_organism = "hsa"
kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "kegg")
#> preparing geneSet collections...
#> GSEA analysis...
#> leading edge analysis...
#> done...
gene <- names(kegg_gene_list)
edox <- setReadable(kk2, 'org.Hs.eg.db', 'ENTREZID')
df_kegg = as.data.table(edox@result)

## -----------------------------------------------------------------------------
dotplot(kk2, showCategory=10, split=".sign", title = "Enriched Pathways") + facet_grid(.~.sign)


## -----------------------------------------------------------------------------
upsetplot(kk2, showCategory = 10)


## ----message=FALSE, warning=FALSE---------------------------------------------
library(enrichplot)

library(DOSE)

library(clusterProfiler)

organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)


## ----message=FALSE, warning=FALSE---------------------------------------------
#subset gene symbol deregulated using the interesting circRNA marker as stratificator
# geneList <- gene_mark$log2FoldChange[gene_mark$circIMPACT==markers.circrnas[5]]
geneList <- gene_mark_hipk3$log2FoldChange[gene_mark_hipk3$gene_id%in%rownames(gene_class_hipk3$VI)]
names(geneList) <- gene_mark_hipk3$gene_id[gene_mark_hipk3$gene_id%in%rownames(gene_class_hipk3$VI)]

# order gene list by foldchange
geneList = sort(geneList, decreasing = TRUE)
# names(geneList) <- gene_mark$Gene[gene_mark$circIMPACT==markers.circrnas[5]]


gse <- gseGO(geneList=geneList, 
             ont ="ALL", 
             keyType = "SYMBOL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")

df_gse <- gse@result
df_gse = as.data.table(df_gse)

## -----------------------------------------------------------------------------
enrichplot::dotplot(gse, showCategory=10, split=".sign", title = "Enriched GO") + facet_grid(.~.sign)


## -----------------------------------------------------------------------------
upsetplot(gse, showCategory=10)


## -----------------------------------------------------------------------------
kable(head(df_gse), escape = F, align = "c", row.names = F, caption = "Enrichment of GO including HPSE gene") %>%
  kable_styling(c("striped"), full_width = F)

## ----eval=FALSE, message=FALSE, warning=FALSE, include=FALSE------------------
#  original_gene_list <- geneList
#  ids<-bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)
#  
#  de_ids = ids[!duplicated(ids[c("SYMBOL")]),]
#  
#  Eids<-bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENSEMBL", OrgDb=organism)
#  
#  df <- gene_mark_hipk3[gene_mark_hipk3$gene_id%in%rownames(gene_class_hipk3$VI),]
#  df$SYMBOL <- df$gene_id
#  df <- df[which(df$SYMBOL%in%de_ids$SYMBOL),]
#  
#  # Create a new column with the corresponding ENTREZ IDs
#  df$Y = de_ids$ENTREZID[match(df$SYMBOL,de_ids$SYMBOL)]
#  
#  # Create a vector of the gene unuiverse
#  kegg_gene_list <- df$log2FoldChange
#  
#  # Name vector with ENTREZ ids
#  names(kegg_gene_list) <- df$Y
#  
#  # omit any NA values
#  kegg_gene_list<-na.omit(kegg_gene_list)
#  
#  # sort the list in decreasing order (required for clusterProfiler)
#  kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
#  
#  kegg_organism = "hsa"
#  kk2 <- gseKEGG(geneList     = kegg_gene_list,
#                 organism     = kegg_organism,
#                 nPerm        = 10000,
#                 minGSSize    = 3,
#                 maxGSSize    = 800,
#                 pvalueCutoff = 0.1,
#                 pAdjustMethod = "none",
#                 keyType       = "kegg")
#  
#  gene <- names(kegg_gene_list)
#  edox <- setReadable(kk2, 'org.Hs.eg.db', 'ENTREZID')
#  df_kegg = as.data.table(edox@result)

## ----eval=FALSE---------------------------------------------------------------
#  dotplot(kk2, showCategory=10, split=".sign", title = "Enriched Pathways") + facet_grid(.~.sign)
#  

## ----eval=FALSE---------------------------------------------------------------
#  upsetplot(kk2, showCategory = 10)
#  

