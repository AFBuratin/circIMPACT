## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval = FALSE-------------------------------------------------------------
#  devtools::install_github("ABuratin/CircTarget")

## ----setup, warning = FALSE, message = FALSE----------------------------------
library(CircTarget)
library(DESeq2)
library(data.table)
library(plyr)
library(Rtsne)

## -----------------------------------------------------------------------------
data("circularData")
data("count.matrix")
data("meta")
data("coldata.df")
## define sample/condition colors
intgroup.dt <- meta[, .(.N), by = .(sample_id, 
                                      condition)][order(sample_id), 
                                                  .(sample_id), by = condition]
samples.per.condition <- data.frame(intgroup.dt[, .N, by = .(condition)], row.names = "condition")
  
n.conditions <- nrow(samples.per.condition)
hues <- CircTarget::gg_color_hue(n.conditions)
for(i in 1:n.conditions){
      n.hues <- samples.per.condition[i, "N"] + 2
      col.hues <- colorRampPalette(colors = c(hues[i], 
                                              "white"))(n.hues)[1:(n.hues-2)]
      
      intgroup.dt[condition == rownames(samples.per.condition)[i], `:=`(color = col.hues,
                                                                        hue = hues[i])]
    }
## create a deseq dataset object 
dds.circular <- suppressMessages(DESeqDataSetFromMatrix(countData = ceiling(circularData[, coldata.df$sample_id[order(coldata.df$condition)]]),
                                   colData = coldata.df[order(coldata.df$condition),],
                                   design = ~ condition))
dds.circular <- suppressMessages(estimateSizeFactors(dds.circular))
sf <- sizeFactors(dds.circular)

## -----------------------------------------------------------------------------
data.filt <- circularData[rowSums(circularData >= 5) >= 3,]
dds.filt.expr <- suppressMessages(DESeqDataSetFromMatrix(countData = ceiling(data.filt[,coldata.df$sample_id[order(coldata.df$condition)]]),
                                   colData = coldata.df[order(coldata.df$condition),],
                                   design = ~ condition))
dds.filt.expr <- suppressMessages(estimateSizeFactors(dds.filt.expr))
sf.filt <- sizeFactors(dds.filt.expr)


## ----message=FALSE, warning=FALSE---------------------------------------------

circtarget <- marker.selection(dat = data.filt, dds = dds.filt.expr, sf = sf, p.cutoff = 0.1, lfc.cutoff = 1.5, 
                               method.d = "euclidean", method.c = "average", k = 2)


## -----------------------------------------------------------------------------
markers.circrnas <- circtarget$circ.targetIDS
mat.filt.mark <- circularData[markers.circrnas, ]

dds.filt.mark <- DESeqDataSetFromMatrix(countData = ceiling(mat.filt.mark[,coldata.df$sample_id[order(coldata.df$condition)]]),
                                   colData = coldata.df[order(coldata.df$condition),],
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
tsne_data <- Rtsne(mydist, pca = F, perplexity=5, max_iter=5000)

## getting the two dimension matrix
d_tsne_1 = as.data.frame(tsne_data$Y)
rownames(d_tsne_1) <- colnames(norm.counts.filt)


