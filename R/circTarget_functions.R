#' Define colors of condition
#'
#' @param n number of condion
#'
#' @return data.frame of color hues per condition
#' 
#'
#' @export
gg_color_hue <- function(n) {
  hues = seq(20, 360, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
  # hcl_palettes("sequential (multi-hue)", n = n, palette = "Plasma")
}

#' Define colors of condition
#'
#' @param meta meta file
#'
#' @return data.frame of color hex and hue per condition
#' 
#'
#' @export
get.color.hues <- function(meta){
  intgroup.dt <- meta[, .(.N), by = .(sample_id, 
                                      condition)][order(sample_id), 
                                                  .(sample_id), by = condition]
  samples.per.condition <- intgroup.dt[, .N, by = .(condition)]
  rownames(samples.per.condition) <- samples.per.condition$condition
  n.conditions <- nrow(samples.per.condition)
  
  if(n.conditions <= 12){ 
    palette <- "Set3"
    if(n.conditions <= 8) palette <- "Set2"
    hues <- brewer_pal(palette = palette)(n.conditions)
    if(palette == "Set3"){
      hues[2] <- brewer_pal(palette = "Set1")(1)
      if(length(hues) == 12){
        hues[12] <- brewer_pal(palette = "Set1")(9)[7]
      }
    }
  }else{
    hues <- gg_color_hue(n.conditions)
  }
  
  if(nrow(intgroup.dt) == 1){
    intgroup.dt[, `:=`(color = hues[1],  hue = hues[1])]
  }else{
    for(i in 1:n.conditions){
      n.hues <- samples.per.condition[i, "N"] + 2
      col.hues <- colorRampPalette(colors = c(hues[i], 
                                              "white"))(n.hues$N)[1:(n.hues$N-2)]
      
      intgroup.dt[condition == rownames(samples.per.condition)[i], `:=`(color = col.hues,
                                                                        hue = hues[i])]
    }
  }
  intgroup.dt[]
}



#' Detect circRNA markers
#' To establish a possible impact of circRNA in gene expression
#' we define two groups of samples defined by similarity in circRNAs expression.
#' K-means algorithm is used to define these two groups and then DESeq test (adj. p-value<.01) 
#' have been used to detect circRNA markers.
#'
#' @param circ.m circRNA id
#' @param dds DeseqDataSet object
#' @param sf sizefactor  
#' @param k number of group in K-means algorithm 
#' @param method.d method for calculate matrix distance
#' @param method.c method for clustering
#' 
#' @return data.frame of padj
#' 
#'
#' @export
marker.detection <- function(circ.m, dds, sf, method.d, method.c, k){
  
  # if(median(circ.m)==0){
  #   group <- ifelse(circ.m<mean(circ.m), "DW", "UP")
  # } else {
  #   group <- ifelse(circ.m<median(circ.m), "DW", "UP")
  # }
  
  dist_mat <- dist(as.data.frame(scale(circ.m)), method = method.d)
  hclust_avg <- hclust(dist_mat, method = method.c)
  # kmeans <- kmeans(dist_mat, k, iter.max = 10, nstart = 1)
  cut_avg <- cutree(hclust_avg, k = k)
  group <- ifelse(cut_avg==1, "g1", "g2")
  # group <- ifelse(kmeans$cluster==1, "g1", "g2")
  colData(dds)$group <- as.factor(group)
  design(dds) <- ~group
  sizeFactors(dds) <- sf
  dds <- estimateDispersionsGeneEst(dds)
  dispersions(dds) <- mcols(dds)$dispGeneEst
  dds <- nbinomWaldTest(dds)
  padj =  as.data.frame(results(dds))
  group = t(as.data.frame(colData(dds)$group))
  
  return(res = cbind(padj, group))
}

#' Select circRNA markers
#' To establish a possible impact of circRNA in gene expression
#' we select a subset of circRNA defined as markers: 
#' only circRNA which expression can discretize two groups 
#' (DESeq adj. p-value<.01) have been selected.
#'
#' @param dat circRNA data set
#' @param dds DeseqDataSet object filterd
#' @param sf sizefactor  
#' @param p.cutoff select only circRNAs signficantly differentially expressed
#' @param lfc.cutoff log fold change cutoff to select circRNAs highly deregulated
#' @param k number of group in K-means algorithm 
#' @param method.d method for calculate matrix distance
#' @param method.c method for clustering
#' 
#' 
#' @return list of two elements: circRNAs estimate, circRNA-target ids
#' 
#'
#' @export
marker.selection <- function(dat, dds, sf, p.cutoff, lfc.cutoff, method.d, method.c, k){
  library(doParallel)
  library(dplyr)
  library(tidyr)
  no_cores <- detectCores() - 1  
  registerDoParallel(cores=no_cores)  
  
  circ_mark = foreach::foreach(i=1:nrow(dat), .combine=rbind) %dopar% {
  
  #make a for loop to estimate log2FC and p.adj to select marker circRNAs
  results.temp <- marker.detection(circ.m = dat[i,], dds = dds[i,], 
                                               sf = sf, method.d = method.d, 
                                               method.c = method.c, k = k)

  }

  circ_mark_selection <- circ_mark %>% dplyr::mutate(circ_id = rownames(dat)) %>% 
    tidyr::drop_na() %>% 
    dplyr::filter(padj<=p.cutoff) %>% dplyr::filter(log2FoldChange>lfc.cutoff) %>% 
    arrange(abs(log2FoldChange))

  markers.circrnas = circ_mark_selection$circ_id
  group = data.frame(circ_mark[,-c(1:6)], circ_id = rownames(circ_mark))
  group.df = gather(group, key="sample_id", value = "group", -circ_id) %>% arrange(circ_id)
 
  return(list(circ.mark = circ_mark, circ.targetIDS = markers.circrnas, group.df = group.df))
}
  
#' Detect dergulated genes using circRNA-markers as stratificator of samples
#'
#' @param circ_idofinterest circRNA-marker ID
#' @param circRNAs normalized read count of circRNAs
#' @param linearRNAs row read count of gene
#' @param colData coldata dataframe
#' @param covariates additinal covariates for batch effect
#' @param padj alpha thershold in deseq test 
#' @param expr A TRUE or FALSE variable controlling whether the median expression should be used to define the two group comparison in DEG analysis.
#' 
#' @return data.frame of gene deregulated per circRNA-markers
#' 
#'
#' @export
geneexpression <- function(circ_idofinterest, circRNAs, linearRNAs, group, colData, covariates, padj, expr=FALSE){
  
  if(expr){
    ## define covariates of interest: circlow vs circhigh
    circ_sample <- circRNAs %>% as.data.table() %>% dplyr::filter(circ_id==circ_idofinterest) %>% 
      dplyr::select_if(is.numeric) %>%
      gather(sample_id, sample_val) %>% dplyr::filter(sample_val>=median(sample_val))
    circ_sample <- merge(circ_sample, colData, by = "sample_id")
  
    colData <- colData %>% mutate(circKD = if_else(sample_id%in%circ_sample$sample_id, paste0(circ_idofinterest, "_high"), paste0(circ_idofinterest, "_low")))
    
    colData$circKD <- factor(colData$circKD)
    colData$condition <- factor(colData$condition)
    rownames(colData) <- colData$sample_id
    
    ## make deseqdataset for test
    dds <- DESeqDataSetFromMatrix(countData = ceiling(filt.mat[, rownames(colData)]),
                                  colData = colData,
                                  design = as.formula(paste0("~", covariates, "+ circKD")))
    
    ## performe DEGs
    dds <- DESeq(dds, fitType = "local")
    res <- results(dds)
    n.degs <- sum(res$padj <= padj, na.rm = T)
    degs <- res[which(res$padj <= padj), ]
    res.dt <- as.data.table(res[which(res$padj <= padj), ], keep.rownames = "gene_id")
    rownames(res.dt) <- res.dt$gene_id
    
    } else {
      colData = colData %>% 
        mutate(circM = if_else(sample_id%in%group$sample_id[group$group=="g1"], paste0(circ_idofinterest, "g1"), paste0(circ_idofinterest, "g2")))
  
      # colData$circKD <- meta$circKD[match(colData$sample_id, meta$sample_id)]
      colData$circM <- factor(colData$circM)
      colData$condition <- factor(colData$condition)
      rownames(colData) <- colData$sample_id
  
      ## make deseqdataset for test
      dds <- DESeqDataSetFromMatrix(countData = ceiling(filt.mat[, rownames(colData)]),
                                    colData = colData,
                                    design = as.formula(paste0("~", covariates, "+ circM")))
      
      ## performe DEGs
      dds <- DESeq(dds, fitType = "local")
      res <- results(dds)
      n.degs <- sum(res$padj <= padj, na.rm = T)
      degs <- res[which(res$padj <= padj), ]
      res.dt <- as.data.table(res[which(res$padj <= padj), ], keep.rownames = "gene_id")
      rownames(res.dt) <- res.dt$gene_id
      }
  return(as.data.frame(cbind(res.dt, n.degs)))
}
