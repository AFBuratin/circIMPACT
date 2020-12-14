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
marker.detection <- function(circ_id, circ.m, dds, sf, method.d, method.c, k){
  
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
  # group = t(as.data.frame(colData(dds)$group))
  group = plotCounts(dds, gene = circ_id, normalized = T,
                         intgroup = c("group"), 
                         returnData = TRUE)
  group$sample_id = rownames(group)
  group$circ_id = circ_id
  res = merge(group, padj, by.x="circ_id", by.y="row.names")
  return(res = res)
}

#' Define formatter function for plot table 
#' 
#' 
color_tile3 <- function(fun = "comma", digits = 3, palette = 'PiYG', n = 10) {
  fun=match.fun(fun)
  stopifnot(n >= 5)
  
  Thresh = 0
  nHalf = n/2
  
  ## Make vector of colors for values below threshold
  rc1 = colorRampPalette(colors = c("#4DAC26", "white"), space="Lab")(nHalf)    
  ## Make vector of colors for values above threshold
  rc2 = colorRampPalette(colors = c("white", "#D01C8B"), space="Lab")(nHalf)
  rampcols = c(rc1, rc2)
  ## In your example, this line sets the color for values between 49 and 51. 
  rampcols[c(nHalf, nHalf+1)] = rgb(t(col2rgb("#e1db86")), maxColorValue=256)
  
  # rb1 = seq(Min, Thresh, length.out=nHalf+1)
  # rb2 = seq(Thresh, Max, length.out=nHalf+1)[-1]
  # rampbreaks = c(rb1, rb2)
  
  return_cut <- function(y) 
    cut(y, breaks = c(seq(min(y), Thresh, length.out=nHalf+1), 
                      seq(Thresh, max(y), length.out=nHalf+1)[-1]), labels = 1:n+1, ordered_result = T)
  
  return_col <- function(y) 
    rampcols[as.integer(return_cut(y))]
  
  formattable::formatter("span", x ~ fun(x, digits = digits),
                           style = function(y) formattable::style(
                           display = "block",
                           padding = "0 4px",
                           "border-radius" = "4px",
                           "color" = csscolor("balck"),
                           "background-color" = return_col(y)
                           )
  )
}


#' Define formatter function for plot table 
#' 
#'
color_tile4 <- function(fun = "comma", digits = 0, palette = 'YlGnBu', n = 9) {
  fun=match.fun(fun)
  stopifnot(n >= 5)
  
  return_cut <- function(y) 
    cut(y, breaks = quantile(y, probs = 0:n/n, na.rm = T),
        labels = 1:n, ordered_result = T, include.lowest = T)
  
  return_col <- function(y) 
    RColorBrewer::brewer.pal(n, palette)[as.integer(return_cut(y))]
  
  formatter("span", x ~ fun(x, digits = digits),
            style = function(y) formattable::style(
              display = "block",
              padding = "0 4px",
              "border-radius" = "4px",
              "color" = ifelse( return_cut(y) %in% c(n-3, n-2, n-1, n),
                                csscolor("white"), csscolor("black")),
              "background-color" = return_col(y)
            )
  )
}


#' Define formatter function for plot table 
#' 
#'
stoplighttile <- function(cut1 = .01, cut2 = .05, cut3 = 0.1, fun = "comma", digits = 4) {
  fun=match.fun(fun)
  formattable::formatter("span", x ~ fun(x, digits = digits),
                         style = function(y) formattable::style(              
                           display = "block",              
                           padding = "0 4px",              
                           "border-radius" = "4px",              
                           "color" = ifelse( y <= cut3, csscolor("black"), csscolor("grey")),
                           "background-color" = ifelse( y <= cut1, csscolor("#d11b3a"),
                                                        ifelse( y <= cut2, csscolor("#d11b70"),
                                                                ifelse( y <= cut3, csscolor("#d11ba4"),
                                                                        csscolor("#731bd1"))))
                           )
  )
}

## Include this function somewhere
# sparkline_hist <- function(x, breaks, type="bar"){
#   as.character(
#     htmltools::as.tags(
#       sparkline(
#         hist(
#           x,
#           breaks=breaks,
#           plot=FALSE
#         )$density,
#         type = type
#       )
#     )
#   )
# }

#' Define formatter function for plot table 
#' 
#'
sparkline_hist <- function(x, g){
  data = data.frame(count=x, group=g)
    gg=ggplot(data, aes(x=count, color=group, fill=group)) +
  geom_density(alpha=0.3) + 
  # geom_vline(data=mu, aes(xintercept=grp.mean, color=group),
  #            linetype="dashed") +
  scale_fill_brewer(palette="Dark2") + 
  scale_color_brewer(palette="Dark2") + 
  # labs(title=paste0("circMarker (", circMark, ")", " counts density curve"), x = "Normalized read counts", y = "Density") + 
  theme_classic()
  plotly=plotly::ggplotly(gg)
as.character(htmltools::as.tags(as_widget(plotly)))
  }


#' Export a Formattable as PNG, PDF, or JPEG
#'
#' @param f A formattable.
#' @param file Export path with extension .png, .pdf, or .jpeg.
#' @param width Width specification of the html widget being exported.
#' @param height Height specification of the html widget being exported.
#' @param background Background color specification.
#' @param delay Time to wait before taking webshot, in seconds.
#'
#' @importFrom formattable as.htmlwidget
#' @importFrom htmltools html_print
#' @importFrom webshot webshot
#'
#' @export
export_formattable <- function(f, file, width = "100%", height = "100%", 
                               background = "white", delay = 10)
{
  w <- as.htmlwidget(f, width = width, height = height)
  path <- html_print(w, background = background, viewer = NULL)
  url <- paste0("file:///", gsub("\\\\", "/", normalizePath(path)))
  webshot(url,
          file = file,
          selector = ".formattable_widget",
          delay = delay)
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
marker.selection <- function(dat, dds, sf, p.cutoff, lfc.cutoff, method.d, method.c, k, plot=FALSE, n=9){
  library(doParallel)
  library(dplyr)
  library(tidyr)
  no_cores <- detectCores() - 1  
  registerDoParallel(cores=no_cores)  
  
  circ_mark = foreach::foreach(i=1:nrow(dat), .combine=rbind) %dopar% {

  circ_id <- rownames(dat)[i]
  #make a for loop to estimate log2FC and p.adj to select marker circRNAs
  results.temp <- marker.detection(circ_id = circ_id, circ.m = dat[circ_id,], dds = dds[circ_id,], 
                                   sf = sf, method.d = method.d, 
                                   method.c = method.c, k = k)
  
  }

  if(!is.null(p.cutoff)){
    if(!is.null(lfc.cutoff)){
      circ_mark_selection <- circ_mark %>%
        tidyr::drop_na() %>% 
        dplyr::filter(padj<=p.cutoff) %>% dplyr::filter(abs(log2FoldChange)>=lfc.cutoff) %>% 
        arrange(abs(log2FoldChange))
    } else {
      circ_mark_selection <- circ_mark %>% 
        tidyr::drop_na() %>% 
        dplyr::filter(padj<=p.cutoff) %>% 
        arrange(abs(log2FoldChange))
    } 
  } 
  
  if(!is.null(lfc.cutoff)) {
    circ_mark_selection <- circ_mark %>% 
      tidyr::drop_na() %>% 
      dplyr::filter(abs(log2FoldChange)>=lfc.cutoff) %>% 
      arrange(abs(log2FoldChange))
  } else {
    circ_mark_selection <- circ_mark %>% 
      tidyr::drop_na()
  } 
  
  markers.circrnas = unique(circ_mark_selection$circ_id)
  
    tab_merge <- circ_mark_selection %>%
      dplyr::select(circ_id, log2FoldChange, padj, group, count)

    tab_merge$Marker <- "FALSE"
    tab_merge$Marker[tab_merge$padj<=0.1] <- "TRUE"
    tab_merge$count <- as.numeric(tab_merge$count)
    
    tab_plot = tab_merge %>% 
      # tidyr::spread(group, count) %>% 
      dplyr::group_by(circ_id) %>% 
      dplyr::summarise(
        logFC=round(mean(log2FoldChange),4),
        p.adj=round(mean(padj),4),
        Marker=unique(Marker),
        mean.G1=round(mean(count[group=="g1"], na.rm=TRUE), 4),
        mean.G2=round(mean(count[group=="g2"], na.rm=TRUE), 4)
        # circ.counts=sparkline_hist(
        #   counts,
        #   hist(tab_merge$counts,plot=FALSE)$breaks
        # )
        # count.density = sparkline_hist(
        #     x = counts,
        #     g = group
        # )
        ) %>% formattable::formattable(., align = c("c","c","c","c","c"), list(
          circ_id = formattable::formatter("span", style = ~ formattable::style(color = "grey", font.weight = "bold")),
          logFC = color_tile3(digits = 3, n = 18, fun = "comma", palette = "PiYG"),
          p.adj = stoplighttile(cut1 = 0.01, cut2 = 0.05, cut3 = 0.1, fun = "comma", digits = 4),
          Marker = formattable::formatter("span", 
                             style = x ~ formattable::style(color = ifelse(x, "orange", "gray")), 
                             x ~ formattable::icontext(ifelse(x, "ok", "remove"), ifelse(x, "Yes", "No"))),
          mean.G1 = color_tile4(digits = 3, fun = "comma"),
          mean.G2 = color_tile4(digits = 3, fun = "comma")))
    
    # export_formattable(tab_plot, file = "formattble_table.png")
    print(tab_plot)
    # formattable::formattable(tab_reduce, align = c("c","c","c","c"), list(
    #   circ_id = formattable::formatter("span", style = ~ style(color = "grey", font.weight = "bold")),
    #   log2FoldChange = color_tile3(digits = 3, n = 18),
    #   padj = stoplighttile(cut1 = 0.01, cut2 = 0.05, cut3 = 0.1, fun = "comma", digits = 4),
    #   Marker = formattable::formatter("span", 
    #                      style = x ~ formastyle(color = ifelse(x, "orange", "gray")), 
    #                      x ~ icontext(ifelse(x, "ok", "remove"), ifelse(x, "Yes", "No"))))) 

  return(list(plot=tab_plot, circ.mark = circ_mark, circ.targetIDS = markers.circrnas, group.df = circ_mark[,c("circ_id", "sample_id", "group")]))
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
  
    colData <- colData %>% mutate(circKD = if_else(sample%in%circ_sample$sample, paste0(circ_idofinterest, "_high"), 
                                                   paste0(circ_idofinterest, "_low")))
    
    colData$circKD <- factor(colData$circKD)
    colData$condition <- factor(colData$condition)
    rownames(colData) <- colData$sample
    
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
        mutate(circM = if_else(sample%in%group$sample[group$group=="g1"], paste0(circ_idofinterest, "g1"), paste0(circ_idofinterest, "g2")))
  
      # colData$circKD <- meta$circKD[match(colData$sample_id, meta$sample_id)]
      colData$circM <- factor(colData$circM)
      colData$condition <- factor(colData$condition)
      rownames(colData) <- colData$sample
  
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
