LIB='/cluster/tufts/patralab/rbator01/R_libs/4.4.0'
.libPaths(c("",LIB))

library(DESeq2)
library(apeglm)
library(pheatmap)
library(tidyverse)
library(dplyr)
library(ComplexHeatmap)
library(openxlsx)
library(ggforce)
library(data.table)
library(clusterProfiler)
library('org.Hs.eg.db')
library(ggpubr)
library(ggrepel)

select = dplyr::select
rename = dplyr::rename
filter = dplyr::filter

# Paths & setup -----
setwd('/cluster/tufts/patralab/rbator01/perlis_lab/cspr_apr25/')

dds = readRDS("analysis/processed/umi.trim.1mm_all_dedup.counts.dds.rds")

samples_rm = read.csv("analysis/de/removed_samples_control_batch_q1q3_filtered_batch1_batch2.csv")

## Build count matrix & metadata -----
mat <- counts(dds, normalize=F)
genes_rm = c('RBFOX3','FCGR3B','PLD4','FCAR')

meta = colData(dds) %>%
  as.data.frame() %>%
  filter(!(KO.Target %in% genes_rm))

meta$name = rownames(meta)

ctrl = colnames(mat)[grep('NTC1',colnames(mat))]

##Define KO targets to test ----
all_res = NULL
genes_to_test = setdiff(unique(meta$KO.Target),'NTC1')

for (g in genes_to_test){
  print(g)

  select = meta %>% 
    filter(KO.Target %in% c(g,'NTC1')) %>% .$name
  
  select = setdiff(select, samples_rm$name)

  print(select)
  if (length(select) > 5) {
    mat_i = mat[, select]
    head(mat_i)
    
    meta_i = meta %>% filter(name %in% select)
    mat_i = mat_i[, rownames(meta_i)]

    meta_i$KO.Target = factor(meta_i$KO.Target, levels=c("NTC1", g))
    dds <- DESeqDataSetFromMatrix(countData = mat_i, colData = meta_i, design = ~ KO.Target)

    keep <- rowSums(counts(dds) >= 10) >= 2
    dds_filtered <- dds[keep, ]
    dds_filtered$KO.Target <- relevel(dds_filtered$KO.Target, "NTC1")
    dds_filtered <- DESeq(dds_filtered)

    saveRDS(dds_filtered, paste0("analysis/de/rds/", g, "_vs_ntc1_batch1_rmlow_redo7oct25.rds"))
    resultsNames(dds_filtered)
    name <- gsub("\\-","\\.", paste0("KO.Target_", g, "_vs_NTC1"))
    print(name)

    res_raw <- results(dds_filtered, name = name, independentFiltering = TRUE,
                       cooksCutoff = TRUE, alpha = 0.1)
    res_shrunk <- lfcShrink(dds_filtered, coef = name, type = "apeglm", res = res_raw)

    out <- as.data.frame(res_shrunk)
    out$log2FoldChange_shrunk <- out$log2FoldChange
    out$lfcSE_shrunk          <- out$lfcSE
    out$log2FoldChange_raw    <- res_raw$log2FoldChange
    out$lfcSE_raw             <- res_raw$lfcSE
    out$pvalue                <- res_raw$pvalue
    out$padj                  <- res_raw$padj
    out$gene                  <- rownames(res_shrunk)
    out$ko_gene               <- g
    rownames(out) <- NULL

    out <- out %>%
      dplyr::select(gene, ko_gene,
                    log2FoldChange_raw, lfcSE_raw, pvalue, padj,
                    log2FoldChange_shrunk, lfcSE_shrunk,
                    baseMean, starts_with("stat"), everything())

    fwrite(out, paste0("analysis/de/", g, "_vs_ntc1_batch1_rmlow_redo7oct25.tsv"), sep = "\t", quote = FALSE)
  }
}


## Combine per-KO results into a single table ----
ff = list.files(path = "analysis/de/", pattern = "*_vs_ntc1_batch1_rmlow_redo7oct25.tsv")
ff
all_res=NULL
for (f in ff){
  d = gsub("_vs_ntc1_batch1_rmlow.tsv","",f)
  res = fread(paste0("analysis/de/",f))
  all_res = rbind(all_res, res)
}

all_res  =   all_res %>%
  separate(gene, into=c("ens","gene_name"), sep="\\|", remove = F) %>%
  rename(feature=gene)

fwrite(all_res, "analysis/de/all_res_vs_ntc1_batch1_q1q3_filter_rmlow_redo7oct25.tsv",sep="\t",quote=F)


## Filter by padj and number of DEGs per KO ----
all_res = fread("analysis/de/all_res_vs_ntc1_batch1_q1q3_filter_rmlow_redo7oct25.tsv")

filter = all_res%>%
  filter(padj<0.1)

fwrite(filter, "analysis/de/all_res_vs_ntc1_padj_0.1_batch1_q1q3_filter_rmlow_redo7oct25.tsv",sep="\t",quote=F)

table(filter$ko_gene)

how_many_deg = filter %>%
  group_by(ko_gene) %>%
  summarise(count = n()) %>%
  filter(count>10)

filter_analyze = filter %>%
  filter(ko_gene %in% how_many_deg$ko_gene)
unique(filter_analyze$ko_gene)

fwrite(filter_analyze, "analysis/de/all_res_vs_ntc1_padj_0.1_10deg_batch1_q1q3_filter_rmlow_redo7oct25.tsv",sep="\t",quote=F)


###############################
## Example: Inspect a single KO target (ITGB2)
################################
#
# # Load the DESeq2 object for the KO target of interest
# dds <- readRDS("analysis/de/rds/ITGB2_vs_ntc1_batch1_rmlow.rds")
# 
# resultsNames(dds)
# # Get the raw (unshrunken) results
# res_raw <- results(dds, name = "KO.Target_ITGB2_vs_NTC1", independentFiltering = FALSE)
# 
# # Shrunken LFC (same type you used before)
# res_shrink <- lfcShrink(dds, coef = "KO.Target_ITGB2_vs_NTC1", type = "apeglm", res=res_raw)
# 
# # Normalized counts
# norm_counts <- counts(dds, normalized = TRUE)
# 
# # Extract values for ITGB2
# gene_id <- "ENSG00000160255|ITGB2"
# raw_lfc <- res_raw[gene_id, "log2FoldChange"]
# shrunk_lfc <- res_shrink[gene_id, "log2FoldChange"]
# 
# data.frame(
#   type = c("Raw LFC", "Shrunken LFC"),
#   log2FC = c(raw_lfc, shrunk_lfc)
# )
# 
# # Build a tidy counts table without accidental vector -> "." column
# counts_df <- tibble::tibble(
#   sample_name = colnames(norm_counts),
#   norm_count  = as.numeric(norm_counts[gene_id, , drop = TRUE])  # keep as numeric vector
# )
# 
# # Get sample metadata with sample IDs as a column
# col_df <- as.data.frame(colData(dds)) %>%
#   tibble::rownames_to_column("sample_name")
# 
# # Join successfully
# counts_df <- dplyr::left_join(counts_df, col_df, by = "sample_name")
# 
# counts_df
# 
# # Plot normalized counts
# ggplot(counts_df, aes(x = KO.Target, y = norm_count, color = KO.Target)) +
#   geom_jitter(width = 0.2, size = 3) +
#   stat_summary(fun = mean, geom = "point", shape = 4, size = 4, color = "black") +
#   scale_y_log10() +
#   labs(title = paste("Normalized counts for", gene_id),
#        y = "Normalized counts (log scale)", x = "Condition")
# 
# res_raw[gene_id, c("log2FoldChange", "lfcSE", "pvalue", "padj")]
# res_shrink[gene_id, c("log2FoldChange", "lfcSE", "pvalue", "padj")]

