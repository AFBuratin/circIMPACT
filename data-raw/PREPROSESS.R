# Load raw data from .csv file
exampleData.gtf <- read.csv("data-raw/circRNA_expression_per_sample.csv", sep = ",")
# Apply preprocessing...
library(stringr)
library(data.table)
fix.name.dot <- function(x){
  gsub(pattern = " ", replacement = "", 
       sub("^([0-9].*)", "X\\1", 
           gsub(pattern = "\\.|/", replacement = "_", 
                gsub("\\+", "p", x))))
}
exampleData.gtf$circ_id <- as.character(exampleData.gtf$circ_id)
colnames(exampleData.gtf)[c(7:31)] <- str_replace(fix.name.dot(colnames(exampleData.gtf)[c(7:31)]), "C", "")
exampleData.gtf <- as.data.table(exampleData.gtf)

circularData <- as.matrix(exampleData.gtf[,c(7:31)])
rownames(circularData) <- exampleData.gtf$circ_id
# Save the cleaned data in the required R package location
usethis::use_data(circularData)

# Load meta data
meta <- unique(fread("data-raw/meta_tall.csv")[, .(sample_id = fix.name(Donor), 
                                                condition = Condition,
                                                tissue = Tissue,
                                                Age = gsub(",", ".", Age),
                                                Sex,
                                                RNAQuality)])
meta$condition <- gsub("  ", " ", meta$condition)
usethis::use_data(meta, overwrite = TRUE)


## create a coldata data.frame
coldata.df <- as.data.frame(meta)
rownames(coldata.df) <- coldata.df$sample_id
coldata.df$condition <- factor(coldata.df$condition)
usethis::use_data(coldata.df)

# Load raw data from .csv file  
gene_expression_TPM <- fread("data-raw/gene_expression_TPM_table.csv")
gene_expression <- read.csv("data-raw/transcript_expression_rawcounts_table.csv")
gene_expression.dt <- gene_expression %>% dplyr::select('gene_id', starts_with("X2015")) %>% as.data.table() %>% reshape2::melt(id.vars = "gene_id",  variable.name = "sample_id", value.name = "reads")

gene_expression.dt$sample_id <-sub("C$", "", gene_expression.dt$sample_id)
gene_expression.dt$sample_id <-gsub(pattern = "-|/", replacement = "_", gene_expression.dt$sample_id)
# gene_expression.dt$gene <- SYMBOL_ID$SYMBOL_ID[match(gene_expression.dt$gene_id, rownames(SYMBOL_ID))]
gene_expression.dt$gene <- gene_expression_TPM$`Gene Name`[match(gene_expression.dt$gene_id, gene_expression_TPM$`Gene ID`)]
count.matrix.dt <- dcast(gene_expression.dt, 
                         formula = gene ~ sample_id, 
                         value.var = "reads", 
                         fun.aggregate = sum,
                         fill = 0)

colnames(count.matrix.dt) <- fix.name.dot(colnames(count.matrix.dt))
count.matrix <- as.matrix(count.matrix.dt, 
                          rownames = "gene")

# Save the cleaned data in the required R package location
usethis::use_data(count.matrix)
