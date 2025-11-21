library(data.table)
library(Matrix)
library(openxlsx)
library(tidyverse)
select = dplyr::select
rename = dplyr::rename
filter = dplyr::filter

# generate the count matrix for run 241126-0144_Joshua_Bowen-8326 only -------

setwd('/cluster/tufts/patralab/rbator01/perlis_lab/cspr_apr25/analysis/')

# generate count matrix -----
matrix_dir <- "align/Solo.out/Gene/raw/"

# the non-dedup -----
f <- file(paste0(matrix_dir, "umiDedup-NoDedup.mtx"), "r")
mat <- as.data.frame(as.matrix(readMM(f)))
close(f)
feature.names = fread(paste0(matrix_dir, "features.tsv"), header = FALSE,stringsAsFactors = FALSE, data.table = F)
barcode.names = fread(paste0(matrix_dir, "barcodes.tsv"), header = FALSE, stringsAsFactors = FALSE, data.table = F)
colnames(mat) <- barcode.names$V1
rownames(mat) <- feature.names$V1

fwrite(mat, file = "processed/umi.trim.nondedup.counts.txt", sep = "\t", quote = F, row.names = T, col.names = T)

# the dedup -----
f <- file(paste0(matrix_dir, "umiDedup-1MM_All.mtx"), "r")
mat <- as.data.frame(as.matrix(readMM(f)))
close(f)
feature.names = fread(paste0(matrix_dir, "features.tsv"), header = FALSE,stringsAsFactors = FALSE, data.table = F)
barcode.names = fread(paste0(matrix_dir, "barcodes.tsv"), header = FALSE, stringsAsFactors = FALSE, data.table = F)
colnames(mat) <- barcode.names$V1
rownames(mat) <- feature.names$V1

fwrite(mat, file = "processed/umi.trim.1mm_all_dedup.counts.txt", sep = "\t", quote = F, row.names = T, col.names = T)


# combine with meta data -----
meta = read.xlsx("metadata_batch1/241210_Joy30_Metadata.xlsx") %>%
  mutate(gene_barcode = paste0(KO.Target,"|", Barcode))
write.xlsx(meta, "metadata_batch1/241210_Joy30_Metadata_format.xlsx")
table(meta$KO.Target)

# format mat -----
feature.names = fread(paste0(matrix_dir, "features.tsv"), header = FALSE,stringsAsFactors = FALSE, data.table = F)
colnames(feature.names) = c("ensembl_gene_id","hgnc_symbol","type")

mat = read.table("processed/umi.trim.1mm_all_dedup.counts.txt")
read_per_sample = colSums(mat)
mean(read_per_sample)

hist(read_per_sample)

mat_w_genes = mat %>%
  as.data.frame()

mat_w_genes$ensembl_gene_id = rownames(mat_w_genes)

mat_w_genes = mat_w_genes %>%
  left_join(feature.names , by="ensembl_gene_id")

mat_w_genes =  mat_w_genes %>%
  mutate(ens_gene = paste0(ensembl_gene_id,"|", hgnc_symbol))

mat_w_genes = mat_w_genes %>%
  column_to_rownames("ens_gene") %>%
  dplyr::select(-c(ensembl_gene_id, hgnc_symbol, type))

# change the column names from barcode to gene_barcode ----
convert_names = data.frame(Barcode = colnames(mat_w_genes))
head(convert_names)

convert_names = convert_names %>%
  left_join(meta, by="Barcode")

head(convert_names)

colnames(mat_w_genes) = convert_names$gene_barcode
head(mat_w_genes)

saveRDS(mat_w_genes, "processed/umi.trim.1mm_all_dedup.counts.rds")


# combine with meta data -----
mat = read.table("processed/umi.trim.nondedup.counts.txt")
mat_w_genes = mat
mat_w_genes$ensembl_gene_id = rownames(mat_w_genes)

mat_w_genes = mat_w_genes %>%
  left_join(feature.names, by="ensembl_gene_id")

mat_w_genes =  mat_w_genes %>%
  mutate(ens_gene = paste0(ensembl_gene_id,"|", hgnc_symbol))

mat_w_genes = mat_w_genes %>%
  column_to_rownames("ens_gene") %>%
  dplyr::select(-c(ensembl_gene_id, hgnc_symbol))

# change the column names from barcode to gene_barcode ----
convert_names = data.frame(Barcode = colnames(mat_w_genes))
head(convert_names)

convert_names = convert_names %>%
  left_join(meta, by="Barcode")

head(convert_names)

colnames(mat_w_genes) = convert_names$gene_barcode
head(mat_w_genes)

saveRDS(mat_w_genes, "processed/umi.trim.nondedup.counts.rds")
