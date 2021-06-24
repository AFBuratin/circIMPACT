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

