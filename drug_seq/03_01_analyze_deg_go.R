LIB='/cluster/tufts/patralab/rbator01/R_libs/4.4.0'
.libPaths(c("",LIB))

library(pheatmap)
library(tidyverse)
library(ComplexHeatmap)
library(openxlsx)
library(ggforce)
library(data.table)
library(clusterProfiler)
library('org.Hs.eg.db')
library(data.table)
library(circlize)
library(grid)
library(forcats)
library(stringr)
library(janitor)
library(DESeq2)
library(ggbreak)

select = dplyr::select
rename = dplyr::rename
filter=dplyr::filter

# paths and setup ----
setwd('/cluster/tufts/patralab/rbator01/perlis_lab/cspr_apr25/')

# DEG thresholds
padj_thresh = 0.1 
praw_thresh = 0.1

# load metadata and counts ----
meta = read.xlsx("analysis/metadata_batch1/241210_Joy30_Metadata_format.xlsx")
ctrl_counts = read.xlsx("analysis/nct1_expression/nct1_cpm_ko_genes.xlsx")

rm_samples = fread("analysis/de/removed_samples_control_batch_q1q3_filtered_batch1_batch2.csv") %>%
  mutate(tmp = name) %>%
  separate(tmp, into=c("ko_gene","well","batch"), sep="\\|") 

n_rm = rm_samples %>%
  filter(batch == "b1") %>%
  group_by(ko_gene) %>%
  mutate(n_rm = n()) %>%
  ungroup() %>%
  as.data.frame()

# deg results b1 -----
b1 = fread("analysis/de/all_res_vs_ntc1_batch1_q1q3_filter_rmlow.tsv")
head(b1)

res_sig_b1 = b1 %>%
  filter(pvalue < praw_thresh)

# ko_gene in gene_name at praw < praw_thresh
ko_in_praw_df_b1 <- res_sig_b1 %>%
  group_by(ko_gene) %>%
  summarise(ko_in_praw_deg_b1 = as.integer(ko_gene %in% gene_name), .groups = "drop")%>%
  distinct()

# ko_gene in gene_name at padj < padj_thresh
ko_in_padj_df_b1 <- res_sig_b1 %>%
  filter(padj < padj_thresh) %>%
  group_by(ko_gene) %>%
  summarise(ko_in_padj_deg_b1 = as.integer(ko_gene %in% gene_name), .groups = "drop")%>%
  distinct()

# count DEGs
t_b1 = as.data.frame(table((b1 %>% filter(padj < padj_thresh))$ko_gene))
colnames(t_b1) = c("ko_gene","ndeg_padj_0.1_b1")

top_b1 <- b1 %>%
  group_by(ko_gene) %>%
  arrange(padj) %>%
  slice_min(order_by = padj, n = 10, with_ties = FALSE) %>%
  select(ko_gene, gene_name) %>%
  rename(top_b1 = gene_name)

top_b1_collapsed <- aggregate(top_b1 ~ ko_gene, data = top_b1, FUN = function(x) paste(x, collapse = ","))

# assemble summary workbook ---------

final = ctrl_counts %>%
  rename(ko_gene=gene_name) %>%
  mutate(coverage_pass = ifelse(mean_cpm > 5, 1,0)) %>%
  full_join(n_rm %>% dplyr::select(ko_gene, n_rm), by="ko_gene") %>%
  left_join(t_b1, by = "ko_gene")  %>%
  left_join(ko_in_padj_df_b1, by = "ko_gene") %>%
  left_join(ko_in_praw_df_b1, by = "ko_gene") %>%
  left_join(top_b1_collapsed, by = "ko_gene") 


# Replace NA with 0
final <- final %>%
  mutate(across(starts_with("ko_in_padj_deg_"), ~ if_else(is.na(.), 0L, .))) %>%
  mutate(across(starts_with("ko_in_praw_deg_"), ~ if_else(is.na(.), 0L, .))) %>%
  mutate(across(starts_with("ndeg_padj_0.1_"), ~ if_else(is.na(.), 0L, as.integer(.))))

final = final %>% 
  distinct()
write.xlsx(final, "analysis/de/summary_batch1_q1q3_filter_rmlow_22aug25.xlsx")

# Omics Fig 1A: DEG target overlap, batch 1 only =======
dds = readRDS("analysis/processed/umi.trim.1mm_all_dedup.counts.dds.rds")
genes_rm = c('RBFOX3','FCGR3B','PLD4','FCAR','NTC1')

meta = colData(dds) %>%
  as.data.frame() %>%
  filter(!(KO.Target %in% genes_rm))

all_targets = unique(meta$KO.Target)
all_targets = factor(all_targets, levels = all_targets)

## Helper: build KO–KO log2FC heatmap with asterisks for significant entries
make_heat <- function(deg_file, 
                      title_text, 
                      out_pdf, 
                      lfc_type = c("raw", "shrunken"),
                      sig_col = c("padj", "pvalue"),
                      sig_thresh = 0.1) {
  
  lfc_type <- match.arg(lfc_type)     # "raw" or "shrunken"
  sig_col  <- match.arg(sig_col)      # "padj" or "pvalue"
  
  # Append info to title for clarity
  #title_text <- paste0(title_text, " ", lfc_type, " (", sig_col, " < ", sig_thresh, ")")
  
  # read & standardize column names
  deg <- data.table::fread(deg_file) %>%
    dplyr::rename(gene = gene_name)
  
  # ensure chosen significance column is present
  if (!sig_col %in% colnames(deg)) {
    stop(sprintf("Column '%s' not found in '%s'.", sig_col, deg_file))
  }
  
  # choose LFC column
  lfc_col <- if (lfc_type == "raw") "log2FoldChange_raw" else "log2FoldChange"
  
  # keep targets, aggregate duplicates:
  # - LFC columns: mean
  # - star_flag: TRUE if any chosen sig_col < sig_thresh
  agg <- deg %>%
    dplyr::select(gene, ko_gene, log2FoldChange_raw, log2FoldChange, pvalue, padj) %>%
    dplyr::filter(gene %in% all_targets, ko_gene %in% all_targets) %>%
    dplyr::group_by(gene, ko_gene) %>%
    dplyr::summarise(
      log2FoldChange_raw = mean(log2FoldChange_raw, na.rm = TRUE),
      log2FoldChange     = mean(log2FoldChange,     na.rm = TRUE),
      star_flag          = any(.data[[sig_col]] < sig_thresh, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(dplyr::across(c(log2FoldChange_raw, log2FoldChange),
                                ~ ifelse(is.nan(.x), NA_real_, .x)))
  
  # complete to full target grid; impute chosen LFC to 0
  heat_long <- agg %>%
    tidyr::complete(
      gene    = all_targets,
      ko_gene = all_targets
    ) %>%
    dplyr::mutate(
      value = .data[[lfc_col]],
      value = tidyr::replace_na(value, 0),
      star_flag = tidyr::replace_na(star_flag, FALSE),
      gene    = factor(gene, levels = all_targets),
      ko_gene = factor(ko_gene, levels = all_targets)
    )
  
  # value matrix
  heat_mat <- heat_long %>%
    dplyr::select(gene, ko_gene, value) %>%
    tidyr::pivot_wider(names_from = ko_gene, values_from = value) %>%
    tibble::column_to_rownames("gene") %>%
    as.matrix()
  
  # asterisk (logical) matrix
  star_mat <- heat_long %>%
    dplyr::select(gene, ko_gene, star_flag) %>%
    tidyr::pivot_wider(names_from = ko_gene, values_from = star_flag) %>%
    tibble::column_to_rownames("gene") %>%
    as.matrix()
  
  ## cluster columns (correlation distance)
  d_col <- as.dist(1 - stats::cor(heat_mat, use = "pairwise.complete.obs"))
  col_clust <- stats::hclust(d_col, method = "average")
  col_order_names <- colnames(heat_mat)[col_clust$order]
  
  # row order aligned to clustered columns (append non-overlapping rows)
  common <- intersect(col_order_names, rownames(heat_mat))
  first_part <- match(common, rownames(heat_mat))
  row_order <- c(first_part, setdiff(seq_len(nrow(heat_mat)), first_part))
  
  row_gp    <- grid::gpar(fontfamily = "Arial", fontsize = 8)
  col_gp    <- grid::gpar(fontfamily = "Arial", fontsize = 8)
  title_gp  <- grid::gpar(fontfamily = "Arial", fontsize = 12, fontface = "bold")
  legend_gp <- grid::gpar(fontfamily = "Arial")
  
  ht <- ComplexHeatmap::Heatmap(
    heat_mat,
    name = "log2FC", #ifelse(lfc_type == "raw", "log2FC_raw", "log2FC_shrunken"),
    col = circlize::colorRamp2(
      c(min(heat_mat, na.rm = TRUE), 0, max(heat_mat, na.rm = TRUE)),
      c("blue", "white", "red")
    ),
    cluster_columns = as.dendrogram(col_clust),
    cluster_rows = FALSE,
    row_order = row_order,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = row_gp,
    column_names_gp = col_gp,
    column_title_gp = title_gp,
    heatmap_legend_param = list(
      labels_gp = legend_gp,
      title_gp  = grid::gpar(fontfamily = "Arial", fontface = "bold")
    ),
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (!is.na(star_mat[i, j]) && star_mat[i, j]) {
        grid::grid.text("*", x = x, y = y,
                        gp = grid::gpar(fontfamily = "Arial", fontsize = 10, fontface = "bold"))
      }
    }
  )
  
  # Use Cairo to embed Arial cleanly in the PDF
  Cairo::CairoPDF(out_pdf, height = 5, width = 6, family = "Arial")
  on.exit(dev.off(), add = TRUE)
  ComplexHeatmap::draw(ht)
}

# KO–KO overlap using raw LFC and p-value < 0.1

make_heat(
  deg_file  = "analysis/de/all_res_vs_ntc1_batch1_q1q3_filter_rmlow.tsv",
  title_text = "",
  out_pdf   = "analysis/de/plots/deg_target_overlap_b1_lfcraw_pval_0.1.pdf",
  lfc_type = "raw",
  sig_col = "pvalue",
  sig_thresh = 0.1 )


make_heat(
  deg_file  = "analysis/de/all_res_vs_ntc1_batch1_q1q3_filter_rmlow.tsv",
  title_text = "",
  out_pdf   = "analysis/de/plots/deg_target_overlap_b1_lfcraw_padj_0.1.pdf",
  lfc_type = "raw",
  sig_col = "padj",
  sig_thresh = 0.1 )


# Barplot NUMBER DEG =======

tab <- b1 %>%
  group_by(ko_gene) %>%
  summarise(n = sum(padj < 0.1, na.rm = TRUE), .groups = "drop")

# Color/category per KO from the *diagonal* row (ko_gene == gene_name)
diag_map <- b1 %>%
  filter(ko_gene == gene_name) %>%
  transmute(
    ko_gene,
    p_raw = pvalue,
    lfc   = dplyr::coalesce(log2FoldChange, log2FoldChange_raw),
    direction = case_when(
      is.finite(p_raw) & p_raw < 0.1 & is.finite(lfc) & lfc > 0 ~ "sig_up",
      is.finite(p_raw) & p_raw < 0.1 & is.finite(lfc) & lfc < 0 ~ "sig_down",
      TRUE ~ "other"
    )
  )

tab <- tab %>%
  left_join(diag_map[, c("ko_gene", "direction")], by = "ko_gene") %>%
  mutate(direction = tidyr::replace_na(direction, "other"),
         direction = factor(direction, levels = c("sig_up", "sig_down", "other")))


# Find outlier and compute a sensible break
top_row    <- tab %>% slice_max(n, n = 1, with_ties = FALSE)
top_val    <- top_row$n
max_others <- max(tab$n[tab$ko_gene != top_row$ko_gene], na.rm = TRUE)

gap       <- max(5, floor(0.05 * (top_val - max_others)))   # ~5% gap or ≥5
lower_brk <- max_others + gap
upper_brk <- top_val - gap

# Axis ticks: show normal ticks below break + a single tick at 1100 above break
lower_ticks <- pretty(c(0, max_others))
upper_tick  <- 1100
all_breaks  <- unique(c(lower_ticks, upper_tick))

p <- ggplot(tab, aes(x = fct_reorder(ko_gene, n), y = n, fill = direction)) +
  geom_col(color = "black") +
  coord_flip(clip = "off") +
  labs(
    x = "Target Gene",
    y = "n (DEGs with padj < 0.1)",
    title = "",
    fill = "DEG status of Target Gene"
  ) +
  scale_fill_manual(
    values = c(sig_up = "red", sig_down = "blue", other = "black"),
    labels = c(
      sig_up   = "p < 0.1 & log2FC > 0",
      sig_down = "p < 0.1 & log2FC < 0",
      other    = "otherwise"
    )
  ) +
  # Broken axis
  ggbreak::scale_y_break(c(lower_brk, upper_brk)) +
  # Only label lower ticks + one tick (1100) in the upper segment
  scale_y_continuous(
    breaks = all_breaks,
    labels = function(b) ifelse(b == upper_tick, "1100", b)
  ) +
  theme_minimal(base_size = 12, base_family = "Arial") +
  theme(
    text = element_text(family = "Arial"),
    plot.title = element_text(face = "bold"),
    legend.position = "right",
    plot.margin = margin(5.5, 20, 5.5, 5.5)
  ) +
  # ensure the 1100 tick is in range even if top_val < 1100
  expand_limits(y = max(top_val, upper_tick) * 1.02)

systemfonts::match_font("Arial")

show(p)
ggsave("analysis/de/plots/deg_per_ko_gene.pdf", p,width = 7.5, height = 5, device = cairo_pdf, family = "Arial",onefile = FALSE)

# Same barplot but all blue (no categories in legend)
p <- ggplot(tab, aes(x = fct_reorder(ko_gene, n), y = n)) +
  geom_col(color = "black", fill = "blue") +
  coord_flip(clip = "off") +
  labs(
    x = "Target Gene",
    y = "n (DEGs with padj < 0.1)",
    title = ""
  ) +
  ggbreak::scale_y_break(c(lower_brk, upper_brk)) +
  scale_y_continuous(
    breaks = all_breaks,
    labels = function(b) ifelse(b == upper_tick, "1100", b)
  ) +
  theme_minimal(base_size = 12, base_family = "Arial") +
  theme(
    text = element_text(family = "Arial"),
    plot.title = element_text(face = "bold"),
    legend.position = "none",          # hide legend
    plot.margin = margin(5.5, 20, 5.5, 5.5)
  ) +
  expand_limits(y = max(top_val, upper_tick) * 1.02)

print(p)
ggsave("analysis/de/plots/deg_per_ko_gene_blue.pdf", p,
       width = 5, height = 5, device = cairo_pdf, family = "Arial", onefile = FALSE)


# GO ENRICHMENT for targets that have >10 DEG =================

b1_filter = fread("analysis/de/all_res_vs_ntc1_padj_0.1_10deg_batch1_q1q3_filter_rmlow.tsv")


ck<- compareCluster(geneCluster = ens ~ ko_gene,
                    data = b1_filter,
                    OrgDb = org.Hs.eg.db,
                    keyType="ENSEMBL",
                    fun = "enrichGO",
                    ont="BP",
                    readable=T)

saveRDS(ck, "analysis/de/go_bp_batch1_padj_0.1_10deg.rds")
write.xlsx(ck, "analysis/de/go_bp_batch1_padj_0.1_10deg.xlsx")

ck<- compareCluster(geneCluster = ens ~ ko_gene,
                    data = b1_filter,
                    OrgDb = org.Hs.eg.db,
                    keyType="ENSEMBL",
                    fun = "enrichGO",
                    ont="CC",
                    readable=T)

saveRDS(ck, "analysis/de/go_cc_batch1_padj_0.1_10deg.rds")
write.xlsx(ck, "analysis/de/go_cc_batch1_padj_0.1_10deg.xlsx")

ck<- compareCluster(geneCluster = ens ~ ko_gene,
                    data = b1_filter,
                    OrgDb = org.Hs.eg.db,
                    keyType="ENSEMBL",
                    fun = "enrichGO",
                    ont="MF",
                    readable=T)

saveRDS(ck, "analysis/de/go_mf_batch1_padj_0.1_10deg.rds")
write.xlsx(ck, "analysis/de/go_mf_batch1_padj_0.1_10deg.xlsx")

## Quick sanity plots for GO enrichments (BP / MF / CC)

ck = readRDS("analysis/de/go_bp_batch1_padj_0.1.rds")
dotplot(ck, show=3)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_y_discrete(labels = ~ str_wrap(as.character(.x), 100))

ck = readRDS("analysis/de/go_mf_batch1_padj_0.1.rds")
dotplot(ck, show=3)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_y_discrete(labels = ~ str_wrap(as.character(.x), 100))

ck = readRDS("analysis/de/go_cc_batch1_padj_0.1.rds")
dotplot(ck, show=3)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_y_discrete(labels = ~ str_wrap(as.character(.x), 100))

## GO BP: filter by keyword families ----

keywords <- c(
  "phago", "autophag",
  "neuro",
  "macrophage", "activation",
  "adhesion", "attachment", "integrin", "actin",
  "chemotaxis", "chemokine", "cytokine",
  "inflamm", "stress",
  "lyso",
  "metabolic", "oxydat",
  "mitochondria"
)

ck = readRDS("analysis/de/go_bp_batch1_padj_0.1_10deg.rds")

# Collapse keywords into regex pattern
pattern <- paste(keywords, collapse = "|")

# Filter the dataframe based on keyword match in the specified column
ck_filter = ck  %>%
  filter(grepl(pattern, ck@compareClusterResult$Description, ignore.case = TRUE))

dotplot(ck_filter, show=3)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_y_discrete(labels = ~ str_wrap(as.character(.x), 100))


# Select terms from collaborators for initial combination -----
#BP
ck = readRDS("analysis/de/go_bp_batch1_padj_0.1_10deg.rds")

s = read.xlsx("analysis/de/go_bp_batch1_padj_0.1_10deg_sds.xlsx") %>%
  mutate(select = ifelse(select == 1, 's',select)) %>%
  filter(select == 's') %>%
  select(ID, select)

l = read.xlsx("analysis/de/go_bp_batch1_padj_0.1_10deg_LM.xlsx") %>%
  filter(!is.na(LM)) %>%
  rename(select = LM) %>%
  mutate(select = 'l') %>%
  select(ID, select)

c = read.xlsx("analysis/de/go_bp_batch1_padj_0.1_10deg.xlsx") %>%
  filter(!is.na(CAMILLA_list)) %>%
  rename(select = CAMILLA_list ) %>%
  select(ID, select)


all_select = rbind(s,l,c) %>%
  distinct() 

dedup_select <- all_select %>%
  mutate(
    is_num   = grepl("^\\d+$", select),
    # convert only the numeric entries; letters -> NA (no warning)
    num_val  = as.integer(replace(select, !is_num, NA_character_)),
    letter_score = case_when(
      !is_num & select == "s" ~ 0L,   # "s" beats other letters
      !is_num                  ~ 1L,
      TRUE                     ~ NA_integer_
    ),
    key1 = if_else(is_num, 0L, 1L),                # numbers (0) beat letters (1)
    key2 = if_else(is_num, num_val, letter_score)  # among numbers: smallest wins; among letters: "s" wins
  ) %>%
  arrange(ID, key1, key2) %>%
  distinct(ID, .keep_all = TRUE) %>%   # one row per GO ID
  select(ID, select)

ck_select = ck@compareClusterResult %>%
  left_join(dedup_select, by="ID") %>%
  filter(p.adjust < 0.05)


lab_map <- c(
  `1` = "Phagocytosis / Complement activation",
  `2` = "Cell adhesion",
  `3` = "Chemotaxis",
  `4` = "PH regulation",
  `5` = "Chemokines/Cytokines",
  `6` = "Protein production",
  l   = "liam",
  s   = "steve"
)

ck_select$select_label <- lab_map[as.character(ck_select$select)]

write.xlsx(ck_select, "analysis/de/go_bp_batch1_padj_0.1_10deg_allselect.xlsx")

# MF ------
ck = readRDS("analysis/de/go_mf_batch1_padj_0.1")

l = read.xlsx("analysis/de/go_mf_batch1_padj_0.1_10deg_LM.xlsx") %>%
  filter(!is.na(LM)) %>%
  rename(select = LM) %>%
  mutate(select = 'l') %>%
  select(ID, select)

c = read.xlsx("analysis/de/go_mf_batch1_padj_0.1_10deg.xlsx") %>%
  filter(!is.na(Camilla_list)) %>%
  rename(select = Camilla_list ) %>%
  select(ID,select)

all_select = rbind(l,c) %>%
  distinct() 

dedup_select <- all_select %>%
  mutate(
    is_num   = grepl("^\\d+$", select),
    # convert only the numeric entries; letters -> NA (no warning)
    num_val  = as.integer(replace(select, !is_num, NA_character_)),
    letter_score = case_when(
      !is_num & select == "s" ~ 0L,   # "s" beats other letters
      !is_num                  ~ 1L,
      TRUE                     ~ NA_integer_
    ),
    key1 = if_else(is_num, 0L, 1L),                # numbers (0) beat letters (1)
    key2 = if_else(is_num, num_val, letter_score)  # among numbers: smallest wins; among letters: "s" wins
  ) %>%
  arrange(ID, key1, key2) %>%
  distinct(ID, .keep_all = TRUE) %>%   # one row per GO ID
  select(ID, select)

ck_select = ck %>%
  filter(ID %in% dedup_select$ID)

ck_select_sim <- simplify(ck_select, cutoff = 0.7, by = "p.adjust", select_fun = min)

dotplot(ck_select_sim, show=200)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_y_discrete(labels = ~ str_wrap(as.character(.x), 100))

ck_select_sim_df = ck_select_sim@compareClusterResult %>%
  left_join(dedup_select, by="ID")

lab_map <- c(
  `1` = "Phagocytosis / Complement activation",
  `2` = "Cell adhesion",
  `3` = "Chemotaxis",
  `4` = "PH regulation",
  `5` = "Chemokines/Cytokines",
  `6` = "Protein production",
  l   = "liam",
  s   = "steve"
)

ck_select_sim_df$select_label <- lab_map[as.character(ck_select_sim_df$select)]

write.xlsx(ck_select_sim_df, "analysis/de/go_mf_batch1_padj_0.1_10deg_allselect.xlsx")

# CC -----
ck = readRDS("analysis/de/go_cc_batch1_padj_0.1_10deg.rds")

l = read.xlsx("analysis/de/go_cc_batch1_padj_0.1_10deg_LM.xlsx") %>%
  filter(!is.na(LM)) %>%
  rename(select = LM) %>%
  mutate(select = 'l') %>%
  select(ID, select)

c = read.xlsx("analysis/de/go_cc_batch1_padj_0.1_10deg.xlsx") %>%
  filter(!is.na(Camilla_list)) %>%
  rename(select = Camilla_list ) %>%
  select(ID,select)


all_select = rbind(l,c) %>%
  distinct() 

dedup_select <- all_select %>%
  mutate(
    is_num   = grepl("^\\d+$", select),
    # convert only the numeric entries; letters -> NA (no warning)
    num_val  = as.integer(replace(select, !is_num, NA_character_)),
    letter_score = case_when(
      !is_num & select == "s" ~ 0L,   # "s" beats other letters
      !is_num                  ~ 1L,
      TRUE                     ~ NA_integer_
    ),
    key1 = if_else(is_num, 0L, 1L),                # numbers (0) beat letters (1)
    key2 = if_else(is_num, num_val, letter_score)  # among numbers: smallest wins; among letters: "s" wins
  ) %>%
  arrange(ID, key1, key2) %>%
  distinct(ID, .keep_all = TRUE) %>%   # one row per GO ID
  select(ID, select)

ck_select = ck %>%
  filter(ID %in% dedup_select$ID)

ck_select_sim <- simplify(ck_select, cutoff = 0.7, by = "p.adjust", select_fun = min)

dotplot(ck_select_sim, show=200)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_y_discrete(labels = ~ str_wrap(as.character(.x), 100))

ck_select_sim_df = ck_select_sim@compareClusterResult %>%
  left_join(dedup_select, by="ID")

lab_map <- c(
  `1` = "Phagocytosis / Complement activation",
  `2` = "Cell adhesion",
  `3` = "Chemotaxis",
  `4` = "PH regulation",
  `5` = "Chemokines/Cytokines",
  `6` = "Protein production",
  l   = "liam",
  s   = "steve"
)

ck_select_sim_df$select_label <- lab_map[as.character(ck_select_sim_df$select)]

write.xlsx(ck_select_sim_df, "analysis/de/go_cc_batch1_padj_0.1_10deg_allselect.xlsx")


# make a plot of the CC individual results with collaborator selections

# ---------- helpers ----------
safe_min <- function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
strip_go_brackets <- function(x) sub("\\s*\\[[^\\]]+\\]$", "", x)

# choose one 'select' per y_id:
#   - if any numbers → pick the smallest number
#   - else letters only → prefer "s", otherwise first other letter
choose_select <- function(x) {
  xo <- as.character(x)
  nums <- suppressWarnings(as.integer(xo))
  if (any(!is.na(nums))) return(as.character(min(nums, na.rm = TRUE)))
  if (any(xo == "s")) return("s")
  setdiff(unique(xo), "s")[1]
}

# Start from ck_select_sim_df; map columns ----------

base_df <- ck_select_sim_df %>%
  transmute(
    comparison = Cluster,
    GO_ID      = ID,
    GO_Term    = Description,
    padj       = p.adjust,
    Significant_genes_count = Count,
    select
  )

# OPTIONAL: union of Top-10 terms per comparison by padj
top_ids <- base_df %>%
  group_by(comparison) %>%
  slice_min(order_by = padj, n = 100, with_ties = FALSE) %>%
  pull(GO_ID) %>% unique()
base_df <- base_df %>% filter(GO_ID %in% top_ids)

#  Build plot_df ----------
eps <- 1e-300
cap <- 10

# y_id is now JUST the term
plot_df <- base_df %>%
  mutate(
    mlog10 = pmin(-log10(pmax(padj, eps)), cap),
    y_id   = GO_Term
  ) %>%
  group_by(y_id, comparison) %>%
  summarise(
    mlog10 = if (all(is.na(mlog10))) NA_real_ else max(mlog10, na.rm = TRUE),
    Significant_genes_count = if (all(is.na(Significant_genes_count)))
      NA_real_ else max(Significant_genes_count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # annotation: collapse select per DESCRIPTION
  left_join(
    base_df %>%
      transmute(y_id = GO_Term, select) %>%
      group_by(y_id) %>%
      summarise(select = choose_select(select), .groups = "drop"),
    by = "y_id"
  )

# order y by best padj across comparisons
y_levels <- base_df %>%
  group_by(GO_Term) %>%
  summarise(best_padj = safe_min(padj), .groups = "drop") %>%
  arrange(best_padj) %>%
  pull(GO_Term) %>% unique()

plot_df <- plot_df %>% mutate(y_id = factor(y_id, levels = y_levels))

#  DOTPLOT of −log10(padj) ----------
gg_dot <- ggplot(plot_df, aes(x = comparison, y = y_id)) +
  geom_point(aes(size = Significant_genes_count, color = mlog10), na.rm = TRUE) +
  scale_size_continuous(name = "Significant genes", range = c(1.8, 8)) +
  scale_color_viridis_c(name = expression(-log[10](padj)), option = "C",
                        direction = 1, na.value = "grey85") +
  labs(title = "Selected GO Terms",
       subtitle = sprintf("Colored by −log10(padj), capped at %s;\n size = Significant_genes_count", cap),
       x = "Comparison", y = NULL) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 9),
        axis.text.x = element_text(angle = 45, hjust = 1))
print(gg_dot)

# FINAL plot with all annotations.  =====

# Target order from functional centered on NEG
well_summary=read.xlsx('20250415_Joshua_Bowen_8511/2025822_functionaldata_celllevel_exp31_only_to_show_ipa_go.xlsx')

exp31 = well_summary %>%
  filter(Experiment == "EXP31")

control <- "NEG"

average_by_experiment_treatment_well <- function(cell_data) {
  cell_data %>%
    group_by(Experiment, Treatment, Metadata_Well) %>%
    summarise(
      across(AreaShape_Area:Cell_Level_Phago_Index, mean, na.rm = TRUE),
      .groups = "drop"
    )
}

neg_mean <- well_summary %>%
  filter(Treatment == control) %>%
  summarise(m = mean(Cell_Level_Phago_Index, na.rm = TRUE), .groups = "drop") %>%
  pull(m)

centered_order_tbl <- well_summary %>%
  group_by(Treatment) %>%
  summarise(mean_clpi = mean(Cell_Level_Phago_Index, na.rm = TRUE), .groups = "drop") %>%
  mutate(delta_clpi = mean_clpi - neg_mean) %>%
  filter(Treatment != control) %>%
  arrange(desc(delta_clpi))

# read GO dataframe ------

df_ann <- openxlsx::read.xlsx("analysis/de/go_bp_batch1_padj_0.1_10deg_allselect_rebeccaedit.xlsx")

df <- df_ann %>% 
  janitor::clean_names() %>%
  mutate(
    select_label = as.character(select_label),
    select_label = stringr::str_squish(select_label),
    select_label = dplyr::na_if(select_label, "")   # treat "" as NA
  ) %>%
  filter(!is.na(select_label)) %>%
  mutate(
    padj_safe = pmax(p_adjust, 1e-300)  # guard tiny zeros
  )

df_plot0 <- df %>%
  transmute(
    Pathway = description,
    Target  = ko_gene,
    Count   = count,
    Group   = select_label,
    colval  = padj_safe       
  )

# Keep only targets present in well_summary order
go_targets <- unique(df_plot0$Target)

target_levels <- centered_order_tbl %>%
  filter(Treatment %in% go_targets) %>%
  pull(Treatment)

df_plot0 <- df_plot0 %>% 
  filter(Target %in% target_levels)

# count filter
df_plot0 <- df_plot0 %>% filter(Count > 2)

# ----------------------------
#   Collapse "X" with response/regulation variants
#    Keep, per base X: most treatments, then strongest mean -log10(padj),
#    then preference: exact < response < cellular response < positive < negative < regulation
# ----------------------------
strip_variants_to_base <- function(s) {
  s %>%
    # response variants
    stringr::str_replace(regex("^cellular\\s+response\\s+to\\s+", ignore_case = TRUE), "") %>%
    stringr::str_replace(regex("^response\\s+to\\s+",             ignore_case = TRUE), "") %>%
    # regulation variants
    stringr::str_replace(regex("^positive\\s+regulation\\s+of\\s+", ignore_case = TRUE), "") %>%
    stringr::str_replace(regex("^negative\\s+regulation\\s+of\\s+", ignore_case = TRUE), "") %>%
    stringr::str_replace(regex("^regulation\\s+of\\s+",            ignore_case = TRUE), "") %>%
    stringr::str_squish()
}

variant_type <- function(s) {
  if (stringr::str_detect(s, regex("^cellular\\s+response\\s+to\\s+", TRUE))) return("resp_cellular")
  if (stringr::str_detect(s, regex("^response\\s+to\\s+", TRUE)))             return("response")
  if (stringr::str_detect(s, regex("^positive\\s+regulation\\s+of\\s+", TRUE))) return("pos")
  if (stringr::str_detect(s, regex("^negative\\s+regulation\\s+of\\s+", TRUE))) return("neg")
  if (stringr::str_detect(s, regex("^regulation\\s+of\\s+", TRUE)))           return("reg")
  return("exact")
}

df_plot0 <- df_plot0 %>%
  mutate(
    base_name = strip_variants_to_base(Pathway),
    var_type  = vapply(Pathway, variant_type, character(1))
  )

# metrics per GO term across all targets
path_metrics <- df_plot0 %>%
  group_by(Pathway, base_name) %>%
  summarise(
    n_treat = n_distinct(Target),
    mean_sig = mean(colval, na.rm = TRUE),  # smaller is more significant
    .groups = "drop"
  )

# preference tier: exact (1) < response (2) < cellular response (3) < positive (4) < negative (5) < regulation (6)
path_metrics <- path_metrics %>%
  mutate(tier = dplyr::case_when(
    grepl("^cellular\\s+response\\s+to\\s+", Pathway, ignore.case = TRUE) ~ 3L,
    grepl("^response\\s+to\\s+",            Pathway, ignore.case = TRUE) ~ 2L,
    grepl("^positive\\s+regulation\\s+of\\s+", Pathway, ignore.case = TRUE) ~ 4L,
    grepl("^negative\\s+regulation\\s+of\\s+", Pathway, ignore.case = TRUE) ~ 5L,
    grepl("^regulation\\s+of\\s+",          Pathway, ignore.case = TRUE) ~ 6L,
    TRUE ~ 1L
  ))

winners <- path_metrics %>%
  arrange(base_name, desc(n_treat), mean_sig, tier, Pathway) %>%  
  group_by(base_name) %>%
  dplyr::slice(1) %>%
  ungroup()

# keep only the winning GO term for each base X
df_plot0 <- df_plot0 %>%
  semi_join(winners %>% dplyr::select(Pathway), by = "Pathway")  %>%
  mutate(Group = factor(Group, levels = unique(Group)))

# Keep only certain "metabolic process" categories

keep_processes <- c(
  "ATP metabolic process",
  "regulation of lipid metabolic process",
  "superoxide metabolic process"
)

df_plot0 <- df_plot0 %>%
  filter(
    # keep all terms that don't end in "metabolic process"
    !stringr::str_detect(Pathway, "metabolic process$") |
      # OR are one of the specific ones you care about
      Pathway %in% keep_processes
  )

# Keep union of top 10 per target
n_keep <- 12
top10_per_treat <- df_plot0 %>%
  group_by(Target) %>%
  arrange(colval, Pathway, .by_group = TRUE) %>%  # <- ASC padj
  slice_head(n = n_keep) %>%
  ungroup()

# Collect the unique GO terms that are top-10 in at least one target
keepers <- unique(top10_per_treat$Pathway)

# Filter the *full* table to those GO terms for *all* Treatments where they exist
df_plot0 <- df_plot0 %>%
  dplyr::filter(Pathway %in% keepers)

# Y ordering by max -log10(padj) within Group + dividers

ord <- df_plot0 %>%
  group_by(Group, Pathway) %>%
  summarise(ordkey = suppressWarnings(min(colval, na.rm = TRUE)), .groups = "drop") %>%  # <- MIN padj
  mutate(ordkey = ifelse(is.infinite(ordkey), NA_real_, ordkey)) %>%
  group_by(Group) %>%
  arrange(ordkey, Pathway, .by_group = TRUE) %>%  # <- ASC padj
  mutate(y_in_group = row_number()) %>%
  ungroup()

if (nrow(ord) == 0) stop("No rows left to plot after filters and collapsing.")

sizes <- ord %>% 
  dplyr::count(Group, name = "n") %>%
  dplyr::mutate(offset = dplyr::lag(cumsum(n), default = 0))

ord <- ord %>% left_join(sizes, by = "Group") %>%
  mutate(y = y_in_group + offset)

y_top    <- max(ord$y) + 0.6
y_bottom <- min(ord$y) - 0.4

# positions for right-side labels and horizontal dividers
label_pos <- ord %>%
  dplyr::group_by(Group) %>%
  dplyr::summarise(y = (min(y) + max(y))/2, .groups = "drop")

# compute one boundary after each group, then drop ONLY the very last boundary by y
dividers <- ord %>%
  dplyr::group_by(Group) %>%
  dplyr::summarise(y = max(y) + 0.5, .groups = "drop") %>%
  dplyr::filter(y < max(y))     # <- robust: keep all but the largest boundary

df_plot <- df_plot0 %>%
  left_join(ord %>% dplyr::select(Group, Pathway, y), by = c("Group","Pathway"))

write.xlsx(df_plot, "analysis/de/go_bp_batch1_padj_0.1_10deg_allselect_rebeccaedit_finalcatstoplot.xlsx")
# Plot (x order = ΔCLPI-descending; size=count; color=-log10(padj))
left_tag  <- ".LEFT"
right_tag <- ".RIGHT"
x_levels  <- target_levels

df_plot <- df_plot %>%
  mutate(Target2 = factor(Target, levels = c(left_tag, x_levels, right_tag)))

label_pos <- label_pos %>% mutate(x_lab = right_tag)

p <- ggplot(df_plot, aes(x = Target2, y = y, size = Count, colour = colval)) +
  geom_hline(data = dividers, aes(yintercept = y),
             inherit.aes = FALSE, linewidth = 0.4, colour = "grey70") +
  geom_point(alpha = 0.9) +
  geom_text(data = label_pos, aes(x = x_lab, y = y, label = Group),
            inherit.aes = FALSE, hjust = 0, vjust = 0.5, fontface = "bold") +
  scale_x_discrete(
    drop   = FALSE,
    limits = c(left_tag, x_levels, right_tag),
    breaks = x_levels,
    labels = x_levels,
    expand = expansion(add = c(0, 0.65))
  ) +
  scale_y_reverse(
    limits = c(y_top, y_bottom),
    breaks = ord$y,
    labels = ord$Pathway,
    expand = c(0, 0)
  ) +
  scale_size_area(max_size = 11, name = "Gene Count") +
  # clusterProfiler-like red↔blue (red = more significant)
  scale_colour_gradient(
    name = "padj",
    low  = "red",    # smaller padj = red
    high = "blue",
    oob  = scales::squish,
    guide = guide_colorbar(reverse = TRUE)   # <- lower numbers appear at the top
  ) +
  labs(x = "Target Gene", y = NULL) +
  coord_cartesian(clip = "off") +
  theme_bw(base_size = 12) +
  theme(
    plot.margin = margin(5.5, 70, 5.5, 5.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

p

# make a dotplot annotated with a boxplot --------
# Pick a global font family
base_family <- "Arial"

# Keep the SAME x limits as dotplot (including sentinel columns)
x_limits_full <- c(left_tag, x_levels, right_tag)

# Compute Treatment-level mean (across wells), subtract NEG mean
.delta_bar_data <- function(data, feature, control = "NEG") {
  # Only keep NEG + treatments that appear in the dotplot order
  present <- intersect(c(control, x_levels), unique(data$Treatment))
  
  means <- data %>%
    filter(Treatment %in% present) %>%
    group_by(Treatment) %>%
    summarise(mean_val = mean(.data[[feature]], na.rm = TRUE), .groups = "drop")
  
  neg_mean <- means %>% filter(Treatment == control) %>% pull(mean_val)
  if (length(neg_mean) == 0 || is.na(neg_mean)) {
    stop("NEG control not found in the provided data for feature: ", feature)
  }
  
  # Delta = gene mean - NEG mean, for non-NEG treatments in the dotplot order
  bars <- tibble(Treatment = x_levels) %>%
    left_join(means, by = "Treatment") %>%
    mutate(delta = mean_val - neg_mean)
  
  # Build a full x vector including sentinels (no bars for them)
  df <- tibble(Treatment2 = factor(x_limits_full, levels = x_limits_full)) %>%
    left_join(
      bars %>%
        transmute(Treatment2 = factor(Treatment, levels = x_limits_full),
                  delta),
      by = "Treatment2"
    ) %>%
    mutate(sign = case_when(
      is.na(delta)        ~ "zero",  # nothing drawn
      delta > 0           ~ "pos",
      delta < 0           ~ "neg",
      TRUE                ~ "zero"
    ))
  df
}

# Use Arial in the two delta bar panels
make_delta_bar_layer <- function(data, feature, panel_title, control = "NEG") {
  df <- .delta_bar_data(data, feature, control = control)
  
  ggplot(df, aes(x = Treatment2, y = delta, fill = sign)) +
    geom_col(width = 0.7, na.rm = TRUE) +
    geom_hline(yintercept = 0, linewidth = 0.4, color = "grey50") +
    scale_fill_manual(values = c(pos = "red", neg = "blue", zero = "grey80"), guide = "none") +
    scale_x_discrete(drop = FALSE, limits = x_limits_full, breaks = NULL, expand = expansion(add = c(0, 0.65))) +
    labs(title = panel_title, x = NULL, y = NULL) +
    theme_minimal(base_size = 12, base_family = base_family) +    # <- Arial
    theme(
      plot.title   = element_text(hjust = 0, face = "bold", size = 11, margin = margin(b = 2)),
      axis.text.x  = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none",
      plot.margin = margin(5.5, 70, 0, 5.5)
    )
}

p <- p +
  theme_bw(base_size = 12, base_family = base_family) +           # <- Arial
  theme(plot.margin = margin(5.5, 70, 5.5, 5.5))

p <- p +
  geom_text(
    data = label_pos,
    aes(x = x_lab, y = y, label = Group),
    inherit.aes = FALSE,
    hjust = 0, vjust = 0.5,
    fontface = "bold",
    family = base_family,
    size = 3    
  )

# Rebuild your two annotation bars with plotmath titles using Δ (delta)
p_bar_pi  <- make_delta_bar_layer(
  exp31, "Cell_Level_Phago_Index",
  panel_title = expression(Delta~" well averaged PI"),
  control = "NEG"
)
p_bar_ecc <- make_delta_bar_layer(
  exp31, "AreaShape_Eccentricity",
  panel_title = expression(Delta~" well-averaged Eccentricity"),
  control = "NEG"
)

# make sure label_pos has the x column for the right margin labels
label_pos <- label_pos %>% mutate(x_lab = right_tag)

# 1) Remove any existing text/label geoms from the dotplot
strip_text_layers <- function(g) {
  g$layers <- Filter(function(lyr) !inherits(lyr$geom, c("GeomText", "GeomLabel")), g$layers)
  g
}
p <- strip_text_layers(p)

# Add ONE right-side label layer (Arial + smaller size)
base_family <- "Arial"
p <- p +
  geom_text(
    data = label_pos,
    aes(x = x_lab, y = y, label = Group),
    inherit.aes = FALSE,
    hjust = 0, vjust = 0.5,
    fontface = "bold",
    family = base_family,
    size = 3.2
  )

p <- p +
  theme(
    axis.text.x = element_text(
      angle = 45, hjust = 1, vjust = 1,
      family = "Arial", size = 12
    ),
    # give a bit more bottom margin so angled labels don't clip
    plot.margin = margin(5.5, 70, 10, 5.5)
  )

# combine and save 
combined <- p_bar_pi / p_bar_ecc / p + patchwork::plot_layout(heights = c(0.7, 0.7, 5))

dir.create("analysis/de", recursive = TRUE, showWarnings = FALSE)
dev <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
ggplot2::ggsave(
  filename = "analysis/de/go_final_dotplot_barplot.pdf",
  plot     = combined,
  width    = 10, height = 14, units = "in",
  device   = dev
)

# SINGLE CATEGORY gene heatmap ======
plot_group_heatmaps <- function(df,
                                deg_path = "analysis/de/all_res_vs_ntc1_batch1_q1q3_filter_rmlow.tsv",
                                group    = "Lysosome",
                                out_dir  = "analysis/plots/single_cat/",
                                alpha = 0.1,
                                lfc_thresh = 0.2,
                                centered_order_tbl = NULL,  
                                min_genes_after_deg = 1,
                                min_treats_after_filter = 1,
                                cluster_rows = TRUE,
                                cluster_cols = TRUE,
                                plot_height = 17,
                                plot_width  = 5,
                                base_family = "Arial",
                                lfc_limits = c(-4, 2)
) {
  
  deg_path = "analysis/de/all_res_vs_ntc1_batch1_q1q3_filter_rmlow.tsv"
  out_dir  = "analysis/plots/single_cat/"
  alpha = 0.1
  lfc_thresh = 0.2
  min_genes_after_deg = 1
  min_treats_after_filter = 1
  cluster_rows = TRUE
  cluster_cols = TRUE
  base_family = "Arial"
  lfc_limits = c(-4, 2)
  df = df_phago
  group = "Phagocytosis"
  centered_order_tbl = centered_order_tbl
  plot_width  = sz$width_in
  plot_height = sz$height_in

  # --- Safety & inputs ---
  stopifnot(all(c("select_label","description","gene_id","ko_gene") %in% names(df)))
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Read DEG (only once)
  deg <- data.table::fread(deg_path) %>% 
    dplyr::mutate(gene = gene_name)
  
  # Terms in this group
  terms <- df %>%
    dplyr::filter(select_label == group) %>%
    dplyr::pull(description) %>%
    unique()
  if (length(terms) == 0) stop(sprintf("No terms found with select_label == '%s'.", group))
  
  # Helper to sanitize filenames
  sanitize <- function(x) gsub("[^A-Za-z0-9_\\-]+", "_", x)
  
  # PDF device that embeds fonts
  .pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  
  # Plot each term
  out_files <- list()
  for (term_desc in terms) {
    message("Plotting [", group, "]: ", term_desc)
    
    # genes for this term
    genes <- df %>%
      dplyr::filter(select_label == group, description == term_desc) %>%
      dplyr::pull(gene_id) %>%
      stringr::str_split("/", simplify = FALSE) %>%
      unlist(use.names = FALSE) %>%
      stringr::str_trim() %>%
      unique()
    genes <- genes[genes != ""]
    if (!length(genes)) { warning("No genes parsed for: ", term_desc); next }
    
    # treatments that actually list this term in df
    treats_in_term_df <- df %>%
      dplyr::filter(select_label == group, description == term_desc) %>%
      dplyr::pull(ko_gene) %>%
      unique()
    if (!length(treats_in_term_df)) { warning("No treatments listed for this term in df: ", term_desc); next }
    
    # subset DEG
    deg_sel <- deg %>%
      dplyr::filter(gene %in% genes) %>%
      dplyr::select(gene, log2FoldChange, padj, ko_gene)
    if (!nrow(deg_sel)) { warning("No DEG rows matched for: ", term_desc); next }
    
    # wide -> matrices
    lfc_wide <- deg_sel %>%
      dplyr::select(gene, ko_gene, log2FoldChange) %>%
      tidyr::pivot_wider(names_from = ko_gene, values_from = log2FoldChange)
    p_wide <- deg_sel %>%
      dplyr::select(gene, ko_gene, padj) %>%
      tidyr::pivot_wider(names_from = ko_gene, values_from = padj)
    
    lfc_raw <- lfc_wide %>% tibble::column_to_rownames("gene") %>% as.matrix()
    p_raw   <- p_wide   %>% tibble::column_to_rownames("gene") %>% as.matrix()
    
    # drop cols not present in df for this term
    lfc_raw <- lfc_raw[, colnames(lfc_raw) %in% treats_in_term_df, drop = FALSE]
    p_raw   <- p_raw[,   colnames(p_raw)   %in% treats_in_term_df, drop = FALSE]
    
    # drop all-NA columns, then order by ΔCLPI (only remaining)
    nonempty_cols <- colSums(!is.na(lfc_raw)) > 0
    lfc_raw <- lfc_raw[, nonempty_cols, drop = FALSE]
    p_raw   <- p_raw[,   nonempty_cols, drop = FALSE]
    
    desired_order <- intersect(centered_order_tbl$Treatment, colnames(lfc_raw))
    lfc_raw <- lfc_raw[, desired_order, drop = FALSE]
    p_raw   <- p_raw[,   desired_order, drop = FALSE]
    if (ncol(lfc_raw) < min_treats_after_filter) { warning("Too few treatments after filtering for: ", term_desc); next }
    
    # fill NA (after pruning)
    lfc_mat <- lfc_raw; lfc_mat[is.na(lfc_mat)] <- 0
    p_mat   <- p_raw;   p_mat[is.na(p_mat)]     <- 1
    
    # drop genes with no data across remaining treatments
    nonempty_rows <- rowSums(!is.na(lfc_raw)) > 0
    lfc_mat <- lfc_mat[nonempty_rows, , drop = FALSE]
    p_mat   <- p_mat[nonempty_rows,   , drop = FALSE]
    if (nrow(lfc_mat) < min_genes_after_deg) { warning("Too few genes with data left for: ", term_desc); next }
    
    heatmap_max <- suppressWarnings(max(lfc_mat, na.rm = TRUE)); if (!is.finite(heatmap_max)) heatmap_max <- 1
    heatmap_min <- suppressWarnings(min(lfc_mat, na.rm = TRUE)); if (!is.finite(heatmap_min)) heatmap_min <- -1
    
    # star overlay (use Arial)
    cell_fun_star <- function(j, i, x, y, w, h, fill) {
      if (!is.na(p_mat[i, j]) && p_mat[i, j] < alpha &&
          !is.na(lfc_mat[i, j]) && abs(lfc_mat[i, j]) > lfc_thresh) {
        grid::grid.text("*", x, y, gp = grid::gpar(fontfamily = base_family))
      }
    }
    heatmap_min <- lfc_limits[1]
    heatmap_max <- lfc_limits[2]
    
    # keep white at 0 (diverging)
    col_fun <- circlize::colorRamp2(
      c(heatmap_min, 0, heatmap_max),
      c("darkblue", "white", "red")
    )
    
    # legend ticks (top→bottom shows 2, 0, -2, -4). If yours renders bottom→top,
    # swap to at = c(-4, -2, 0, 2) to match your preference.
    legend_at     <- c(2, 0, -2, -4)
    legend_labels <- c("2", "0", "-2", "-4")
    
    pdf_file <- file.path(out_dir, paste0(sanitize(group), "__", sanitize(term_desc), ".pdf"))
    .pdf_device(pdf_file, height = plot_height, width = plot_width, family = base_family)
    
    view(lfc_mat)
    ht <- ComplexHeatmap::Heatmap(
      lfc_mat,
      name = "Log2FC",
      col  = col_fun,
      cluster_columns = cluster_cols,
      cluster_rows    = cluster_rows,
      column_names_gp = grid::gpar(fontsize = 10, fontfamily = base_family),
      row_names_gp    = grid::gpar(fontsize = 10, fontfamily = base_family),
      heatmap_legend_param = list(
        title     = "Log2FoldChange\nKO vs. Control",
        title_gp  = grid::gpar(fontfamily = base_family, fontface = "bold"),
        labels_gp = grid::gpar(fontfamily = base_family),
        at        = legend_at,
        labels    = legend_labels
      ),
      column_title    = term_desc,
      column_title_gp = grid::gpar(fontfamily = base_family, fontface = "bold"),
      cell_fun        = cell_fun_star
    )

    ComplexHeatmap::draw(ht)
    grDevices::dev.off()
    message("Saved: ", pdf_file)
    out_files[[term_desc]] <- pdf_file
  }
  invisible(out_files)
}

estimate_heatmap_dims <- function(df_plot, deg_path, group, centered_order_tbl) {
  # genes for this group/term set (split "gene_id" on '/')
  genes <- df_plot %>%
    filter(select_label == group) %>%
    pull(gene_id) %>%
    str_split("/", simplify = FALSE) %>%
    unlist(use.names = FALSE) %>%
    str_trim() %>%
    unique()
  genes <- genes[genes != ""]
  
  # treatments mentioned for these terms
  treats_in_term_df <- df_plot %>%
    filter(select_label == group) %>%
    pull(ko_gene) %>%
    unique()
  
  # read DEG and subset to these genes/treatments
  deg <- data.table::fread(deg_path) %>% mutate(gene = gene_name)
  deg_sel <- deg %>%
    filter(gene %in% genes, ko_gene %in% treats_in_term_df) %>%
    select(gene, ko_gene, log2FoldChange)
  
  # wide, then drop all-NA cols and keep centered order
  lfc_wide <- deg_sel %>%
    pivot_wider(names_from = ko_gene, values_from = log2FoldChange)
  if (nrow(lfc_wide) == 0) return(list(nr = 0L, nc = 0L))
  
  lfc_raw <- lfc_wide %>% column_to_rownames("gene") %>% as.matrix()
  nonempty_cols <- colSums(!is.na(lfc_raw)) > 0
  lfc_raw <- lfc_raw[, nonempty_cols, drop = FALSE]
  
  desired_order <- intersect(centered_order_tbl$Treatment, colnames(lfc_raw))
  lfc_raw <- lfc_raw[, desired_order, drop = FALSE]
  
  # drop genes with all NA across remaining treatments
  nonempty_rows <- rowSums(!is.na(lfc_raw)) > 0
  lfc_raw <- lfc_raw[nonempty_rows, , drop = FALSE]
  
  list(nr = nrow(lfc_raw), nc = ncol(lfc_raw))
}

heatmap_device_size <- function(nr, nc,
                                cell_w_mm = 6, 
                                cell_h_mm = 6.5,     
                                row_names_mm = 35, 
                                col_names_mm = 12,
                                row_dend_mm = 4,  
                                col_dend_mm = 4,
                                legend_mm = 18,   
                                gap_mm = 6) {
  mm_to_in <- function(mm) mm / 25.4
  heatmap_w_in <- mm_to_in(cell_w_mm * nc)
  heatmap_h_in <- mm_to_in(cell_h_mm * nr)
  width_in  <- heatmap_w_in + mm_to_in(row_names_mm + row_dend_mm + gap_mm + legend_mm)
  height_in <- heatmap_h_in + mm_to_in(col_names_mm + col_dend_mm + gap_mm)
  list(width_in = width_in, height_in = height_in)
}


# Example: Phagocytosis panel
df_phago <- df %>% filter(description == "phagocytosis")
head(df_phago)

dims <- estimate_heatmap_dims(df_phago, deg_path = "analysis/de/all_res_vs_ntc1_batch1_q1q3_filter_rmlow.tsv",
                              group = "Phagocytosis", centered_order_tbl = centered_order_tbl)

sz <- heatmap_device_size(dims$nr, dims$nc)

files_phago <- plot_group_heatmaps(
  df_phago,
  group = "Phagocytosis",
  centered_order_tbl = centered_order_tbl,
  plot_width  = sz$width_in,
  plot_height = sz$height_in
)

# Lysosome terms
df_plot = df %>%
  filter(description == "lysosome organization")

dims <- estimate_heatmap_dims(df_plot, 
                              deg_path = "analysis/de/all_res_vs_ntc1_batch1_q1q3_filter_rmlow.tsv",
                              group = "Lysosome",
                              centered_order_tbl = centered_order_tbl)

sz <- heatmap_device_size(dims$nr, dims$nc)
files_lyso <- plot_group_heatmaps(
  df_plot,
  group = "Lysosome",
  centered_order_tbl = centered_order_tbl,
  plot_width  = sz$width_in,
  plot_height = sz$height_in
)

# Activation terms
df_plot = df %>%
  filter(description == "macrophage activation")


head(df_plot)
dims <- estimate_heatmap_dims(df_plot, 
                              deg_path = "analysis/de/all_res_vs_ntc1_batch1_q1q3_filter_rmlow.tsv",
                              group = "Inflammation/Activation",
                              centered_order_tbl = centered_order_tbl)

sz <- heatmap_device_size(dims$nr, dims$nc)
files_act <- plot_group_heatmaps(
  df_plot,
  group = "Inflammation/Activation",
  centered_order_tbl = centered_order_tbl,
  plot_width  = sz$width_in,
  plot_height = sz$height_in
)

# Actin terms
df_plot = df %>%
  filter(description == "regulation of cell shape")
dims <- estimate_heatmap_dims(df_plot, 
                              deg_path = "analysis/de/all_res_vs_ntc1_batch1_q1q3_filter_rmlow.tsv",
                              group = "Actin/Cell Shape",
                              centered_order_tbl = centered_order_tbl)

sz <- heatmap_device_size(dims$nr, dims$nc)
files_act <- plot_group_heatmaps(
  df_plot,
  group = "Actin/Cell Shape",
  centered_order_tbl = centered_order_tbl,
  plot_width  = sz$width_in,
  plot_height = sz$height_in
)

# Actin term
df_plot = df %>%
  filter(description == "actin crosslink formation")
dims <- estimate_heatmap_dims(df_plot, 
                              deg_path = "analysis/de/all_res_vs_ntc1_batch1_q1q3_filter_rmlow.tsv",
                              group = "Actin/Cell Shape",
                              centered_order_tbl = centered_order_tbl)

sz <- heatmap_device_size(dims$nr, dims$nc)
files_act <- plot_group_heatmaps(
  df_plot,
  group = "Actin/Cell Shape",
  centered_order_tbl = centered_order_tbl,
  plot_width  = sz$width_in,
  plot_height = sz$height_in
)

df_plot = df %>%
  filter(description == "actin filament bundle assembly")
dims <- estimate_heatmap_dims(df_plot, 
                              deg_path = "analysis/de/all_res_vs_ntc1_batch1_q1q3_filter_rmlow.tsv",
                              group = "Actin/Cell Shape",
                              centered_order_tbl = centered_order_tbl)

sz <- heatmap_device_size(dims$nr, dims$nc)
files_act <- plot_group_heatmaps(
  df_plot,
  group = "Actin/Cell Shape",
  centered_order_tbl = centered_order_tbl,
  plot_width  = sz$width_in,
  plot_height = sz$height_in
)


df_plot = df %>%
  filter(description == "regulation of actin cytoskeleton organization")
dims <- estimate_heatmap_dims(df_plot, 
                              deg_path = "analysis/de/all_res_vs_ntc1_batch1_q1q3_filter_rmlow.tsv",
                              group = "Actin/Cell Shape",
                              centered_order_tbl = centered_order_tbl)

sz <- heatmap_device_size(dims$nr, dims$nc)
files_act <- plot_group_heatmaps(
  df_plot,
  group = "Actin/Cell Shape",
  centered_order_tbl = centered_order_tbl,
  plot_width  = sz$width_in,
  plot_height = sz$height_in
)

