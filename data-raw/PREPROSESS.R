# Load raw data from .csv file
exampleData.gtf <- read.csv("data-raw/bks.counts.union.intersected_Li2017.csv", sep = ",")
# Apply preprocessing...
library(stringr)
library(data.table)
exampleData.gtf$circ_id <- paste0(exampleData.gtf$chr, ":", exampleData.gtf$start, "-", exampleData.gtf$end)
exampleData.gtf <- as.data.table(exampleData.gtf)

circularData <- as.matrix(dcast(exampleData.gtf, circ_id ~ sample_id, value.var = "read.count", fill = 0, fun.aggregate = sum),
                          rownames = "circ_id")
# Save the cleaned data in the required R package location
usethis::use_data(circularData, overwrite = TRUE)

# meta data
meta <- unique(fread("data-raw/meta.csv")[,-1])
usethis::use_data(meta, overwrite = TRUE)

# create a coldata data.frame
coldata.df <- as.data.frame(meta)
rownames(coldata.df) <- coldata.df$sample
coldata.df$condition <- factor(coldata.df$condition)
usethis::use_data(coldata.df, overwrite = TRUE)

# make a dds object for readme.Rmd
library(DESeq2)
dds.circular <- suppressMessages(DESeqDataSetFromMatrix(countData = ceiling(circularData[, coldata.df$sample[order(coldata.df$condition)]]),
                                                        colData = coldata.df[order(coldata.df$condition),],
                                                        design = ~ condition))
dds.circular <- suppressMessages(estimateSizeFactors(dds.circular))
usethis::use_data(dds.circular, overwrite = TRUE)

# Load raw data from .csv file  
gene_expression_TPM <- read.csv("data-raw/gene_expression_TPM_table_Li2017.csv")
gene_expression <- read.csv("data-raw/transcript_expression_rawcounts_table_Li2017.csv")
gene_expression.dt <- gene_expression %>% dplyr::select('gene_id', starts_with("SRR")) %>% as.data.table() %>% reshape2::melt(id.vars = "gene_id",  variable.name = "sample", value.name = "reads")

gene_expression.dt$gene <- gene_expression_TPM$Gene.Name[match(gene_expression.dt$gene_id, gene_expression_TPM$Gene.ID)]
count.matrix.dt <- dcast(gene_expression.dt, 
                         formula = gene ~ sample, 
                         value.var = "reads", 
                         fun.aggregate = sum,
                         fill = 0)
count.matrix <- as.matrix(count.matrix.dt, 
                          rownames = "gene")

# Save the cleaned data in the required R package location
usethis::use_data(count.matrix, overwrite = TRUE)
