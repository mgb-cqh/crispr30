LIB='/cluster/tufts/patralab/rbator01/R_libs/4.4.0'
.libPaths(c("",LIB))

library(DESeq2)
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
library(limma)

select = dplyr::select
rename = dplyr::rename
filter = dplyr::filter

# Paths and setup -----

setwd('/cluster/tufts/patralab/rbator01/perlis_lab/cspr_apr25/')
dedup_string="1mm_all_dedup"

dds1 = readRDS("analysis/processed/umi.trim.1mm_all_dedup.counts.dds.rds")

meta1 <- colData(dds1) %>%
  as.data.frame() 

target_genes <- unique(meta1$KO.Target) %>% 
  setdiff("NTC1")
target_genes

# prepare vst and normalized counts ------
vst1 <- vst(dds1, blind = TRUE)
vst1_mat <- assay(vst1) %>%
  as.data.frame() %>%
  rownames_to_column("feature") %>%
  tidyr::separate(feature, into = c("ens", "gene_name"), sep = "\\|")

norm_counts1 <- counts(dds1, normalized = TRUE) %>%
  as.data.frame() %>%
  rownames_to_column("feature") %>%
  tidyr::separate(feature, into = c("ens", "gene_name"), sep = "\\|")

# Helper functions ----

# long_vst_gene()
# ----------------
# For a single target gene g:
#   - Pull VST values for KO and NTC1 control samples in batch 1
#   - Return a long data frame ready for plotting

long_vst_gene <- function(g, keep_names) {
  # columns for KO and control within batch 1 (by naming convention)
  ko_cols <- grep(paste0("^", g, "\\|"), colnames(vst1), value = TRUE)
  c_cols  <- grep("^NTC1\\|", colnames(vst1), value = TRUE)
  
  # if no matching cols, return empty tibble
  if (length(ko_cols) + length(c_cols) == 0) return(tibble())
  
  df_g <- vst1_mat %>%
    dplyr::filter(gene_name == g)
  
  # KO long
  ko_long <- df_g %>%
    dplyr::select(all_of(ko_cols)) %>%
    pivot_longer(cols = everything(), names_to = "name", values_to = "value") %>%
    mutate(group = "ko", label = sub(paste0("^", g, "\\|"), "", name))
  
  # control long
  c_long <- df_g %>%
    dplyr::select(all_of(c_cols)) %>%
    pivot_longer(cols = everything(), names_to = "name", values_to = "value") %>%
    mutate(group = "c", label = sub("^NTC1\\|", "", name))
  
  bind_rows(ko_long, c_long) %>%
    mutate(color_jitter = ifelse(group == "ko" & !(name %in% keep_names), "remove", "keep"),
           gene = g)
}

# compute_kept_names()
# ---------------------
# Decide which KO samples to keep for a given gene g based on normalized counts:
#   1. Compute control (NTC1) quartiles
#   2. Check if KO mean is above or below control mean
#   3. Keep KO samples that are "clearly" shifted relative to the control distribution
compute_kept_names <- function(g) {
  ko_cols <- grep(paste0("^", g, "\\|"), colnames(norm_counts1 %>% dplyr::select(-ens, -gene_name)), value = TRUE)
  c_cols  <- grep("^NTC1\\|",  colnames(norm_counts1 %>% dplyr::select(-ens, -gene_name)), value = TRUE)
  
  if (length(ko_cols) + length(c_cols) == 0) return(character(0))
  
  nl <- norm_counts1 %>%
    dplyr::filter(gene_name == g) %>%
    dplyr::select(ens, gene_name, all_of(c(ko_cols, c_cols))) %>%
    pivot_longer(cols = -c(ens, gene_name), names_to = "name", values_to = "value") %>%
    mutate(group = ifelse(grepl(paste0("^", g, "\\|"), name), "ko", "c"))
  
  # control quartiles and KO vs control mean direction
  control_q1 <- quantile(nl$value[nl$group == "c"], 0.25, na.rm = TRUE)
  control_q3 <- quantile(nl$value[nl$group == "c"], 0.75, na.rm = TRUE)
  ko_mean    <- mean(nl$value[nl$group == "ko"], na.rm = TRUE)
  c_mean     <- mean(nl$value[nl$group == "c"],  na.rm = TRUE)
  ko_above   <- ko_mean > c_mean
  
  # keep rule 
  nl <- nl %>%
    mutate(keep = dplyr::case_when(
      group == "c" ~ TRUE,
      ko_above     ~ !(value <= control_q3 & value >= control_q1 | value <= control_q1),
      TRUE         ~ !(value <= control_q3 & value >= control_q1 | value >= control_q3)
    ))
  
  nl %>% dplyr::filter(group == "ko", keep) %>% pull(name) %>% unique()
}

## Version 1: Boxplots with KO outliers flagged (red) -----

plot_list <- vector("list", length(target_genes))
names(plot_list) <- target_genes

target_genes <- unique(meta1$KO.Target) %>% setdiff("NTC1")
target_genes <- sort(target_genes)   # <- add this

for (i in seq_along(target_genes)) {
  g <- target_genes[i]
  
  keep_names <- compute_kept_names(g)
  df_vst_g   <- long_vst_gene(g, keep_names)
  
  if (nrow(df_vst_g) == 0) {
    plot_list[[i]] <- ggplot() + theme_void() + ggtitle(paste0(g, " (no data)"))
    next
  }
  
  # Add display labels for x-axis only
  df_vst_g <- df_vst_g %>%
    mutate(group_disp = dplyr::recode(group,
                                      "c"  = "Control (NTC1)",
                                      "ko" = "Targeted (KO)"),
           group_disp = factor(group_disp, levels = c("Control (NTC1)", "Targeted (KO)")))
  
  p_g <- ggplot(df_vst_g, aes(x = group_disp, y = value)) +
    geom_boxplot(outlier.shape = NA, aes(color = group), fill = NA) +   # color still maps to 'group'
    geom_jitter(aes(color = color_jitter), width = 0.15, height = 0, size = 1.8) +
    ggrepel::geom_text_repel(data = df_vst_g %>% dplyr::filter(label != ""),
                             aes(label = label),
                             size = 2.5, max.overlaps = Inf) +
    scale_color_manual(values = c("keep" = "black", "remove" = "red", "ko" = "blue", "c" = "black")) +
    labs(title = g, x = NULL, y = "VST") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    theme_bw(base_size = 10) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  plot_list[[i]] <- p_g
}

# arrange all 31 into a single PDF (5 x 7 grid)
dir.create("analysis/de/plots/final_boxplots/", showWarnings = FALSE, recursive = TRUE)
combined <- ggpubr::ggarrange(plotlist = plot_list, ncol = 5, nrow = 7, align = "hv")

ggsave(filename = "analysis/de/plots/final_boxplots/batch1_only_VST_all31_genes_boxplots.pdf",
       plot = combined, width = 20, height = 28, units = "in")

#  Version 2: Boxplots without KO removal (no red) -----
plot_list <- vector("list", length(target_genes))
names(plot_list) <- target_genes

for (i in seq_along(target_genes)) {
  g <- target_genes[i]
  
  # no filtering logic needed now
  keep_names <- character(0)
  df_vst_g   <- long_vst_gene(g, keep_names)
  
  if (nrow(df_vst_g) == 0) {
    plot_list[[i]] <- ggplot() + theme_void() + ggtitle(paste0(g, " (no data)"))
    next
  }
  
  df_vst_g <- df_vst_g %>%
    mutate(group_disp = dplyr::recode(group,
                                      "c"  = "Control (NTC1)",
                                      "ko" = "Targeted (KO)"),
           group_disp = factor(group_disp, levels = c("Control (NTC1)", "Targeted (KO)")))
  
  # NOTE: color now maps ONLY to 'group' (ko vs c); no "remove" category
  p_g <- ggplot(df_vst_g, aes(x = group_disp, y = value)) +
    geom_boxplot(outlier.shape = NA, aes(color = group), fill = NA) +
    geom_jitter(aes(color = group), width = 0.15, height = 0, size = 1.8) +
    ggrepel::geom_text_repel(data = df_vst_g %>% dplyr::filter(label != ""),
                             aes(label = label),
                             size = 2.5, max.overlaps = Inf) +
    scale_color_manual(values = c("ko" = "blue", "c" = "black")) +
    labs(title = g, x = NULL, y = "VST") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    theme_bw(base_size = 10) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  plot_list[[i]] <- p_g
}

# arrange and save
dir.create("analysis/de/plots/final_boxplots/", showWarnings = FALSE, recursive = TRUE)
combined <- ggpubr::ggarrange(plotlist = plot_list, ncol = 5, nrow = 7, align = "hv")

ggsave(filename = "analysis/de/plots/final_boxplots/batch1_only_VST_all31_genes_boxplots_noRed.pdf",
       plot = combined, width = 20, height = 28, units = "in")
