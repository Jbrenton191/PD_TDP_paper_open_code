#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Figure 5 - sig genes and proteins in Rio data overlap

# Setup

library(magrittr)
library(patchwork)
library(tidyverse)
library(ggpubr)
library(rstatix)

# 0. Setup ----

main_path<-here::here() # main project level i.e TDP43_paper_open_code level

data_path<-file.path(main_path, "Data")

figure_out_path<-file.path(main_path, "Plots/Main_Figures/Figure_5")

# create figure output directory if not present
if (!dir.exists(figure_out_path)) {
  dir.create(figure_out_path, recursive = T)
}

extrafont::font_import(paths = data_path, 
                       pattern = 'Roboto', prompt = F)

theme_jb <- function(text_size=12) {
  theme_bw(base_family = "Roboto", base_size = text_size) +
    theme(
      legend.position = "top",
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.title.y = element_text(vjust = 0.6),
      panel.spacing = unit(0.1, "lines")
    )
}

analysis_colours<-c("#0c5394", "#e06666")

lrrk2_colours<-c("#C3B3A1", "#125C5D")

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Load Proteomic data: Gene level counts with library normalisation added - t.tests across all genes.
lib_normalised_gene_level_summed_counts <- read_csv(file.path(data_path, "LRRK2_HS_PSC/mDN1/lib_normalised_gene_level_summed_counts.csv"))

# Load RNAseq Differential gene expression data  
DESeq_Rio_hPSCs_LRRK2_vs_Ctrl_all_genes <- read_csv(file.path(data_path, "LRRK2_HS_PSC/mDN1/DESeq_Rio_hPSCs_LRRK2_vs_Ctrl_all_genes.csv"))

# combine the two
Log2FC_tibble<-inner_join(
  lib_normalised_gene_level_summed_counts %>% dplyr::select(Genes, Log2FC) %>% dplyr::rename(prot_Log2FC=Log2FC),
  DESeq_Rio_IPSCs_LRRK2_vs_Ctrl_all_genes %>% dplyr::select(gene_name, log2FoldChange) %>% 
    dplyr::rename(DE_Log2FC=log2FoldChange, Genes=gene_name),
  by="Genes"
)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Load significant splicing data:
full_results<-readRDS(file.path(data_path, 'LRRK2_HS_PSC/mDN1/lrrk2_hpsc_results_1.RDS'))
lrrk2_sig_leaf<-full_results$all_results$significant_clusters_0.05_filter

# to get exon - exon junction need to subtract 1 from start of leaf results as are in sj.out format with one added to end when leafcutter run
lrrk2_sig_leaf$start<-as.character(as.numeric(lrrk2_sig_leaf$start)-1)
# separate stand and cluster 
lrrk2_sig_leaf$strand<-gsub(x = lrrk2_sig_leaf$cluster, pattern = "clu.+_(.)", replacement = "\\1")
lrrk2_sig_leaf$cluster<-gsub(x = lrrk2_sig_leaf$cluster, pattern = "(clu.+_).", replacement = "\\1")

# load tdp 43 targets - already in exon - exon format
tdp_junc_list<-readRDS(file.path(data_path, "TDP-43_HS_KD/Significant_TDP_KD_Database_ex_ex_annotated_df.RDS"))

# Distinct junctions across experiments
tdp_junc_list %<>% filter(p.adjust<0.05) %>% 
  mutate(start=as.character(start), end=as.character(end)) %>%
  group_by(start, end, chr, strand) %>%
  distinct(start, end,chr, strand, .keep_all = T) 


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2. Finding signficant TDP-43 junctions in significant leafcutter analysis ----
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sig_leaf_tdp<-lrrk2_sig_leaf %>% distinct(genes,.keep_all = T) %>%# not distinct as multiple brain areas 
  inner_join(., tdp_junc_list %>% distinct(genes),
             by=c("genes")) %>% distinct()


# find novel junctions that are upregulated in LRRK2
novel_sig_leaf_tdp <- sig_leaf_tdp

# testing in sig genes
Log2FC_tibble<-inner_join(
  lib_normalised_gene_level_summed_counts %>% 
    mutate(t.test_fdr=p.adjust(t.test_p_val, "fdr")) %>%
    filter(t.test_fdr<0.05) %>%
    dplyr::select(Genes, Log2FC) %>% dplyr::rename(prot_Log2FC=Log2FC),
  DESeq_Rio_hPSCs_LRRK2_vs_Ctrl_all_genes %>% 
    filter(padj <0.05) %>%
    dplyr::select(gene_name, log2FoldChange) %>% 
    dplyr::rename(DE_Log2FC=log2FoldChange, Genes=gene_name),
  by="Genes"
)

ranked_df<-Log2FC_tibble %>% arrange(desc(prot_Log2FC)) %>% rownames_to_column(var = "Order") %>% 
  
  mutate(is_tdpkd=case_when(Genes %in% unique(novel_sig_leaf_tdp$genes) ~ T,
                            .default = F)) %>% 
  mutate(Order=as.numeric(Order))


# Graph: ----

x <- ranked_df %>%
  filter(is_tdpkd == TRUE) %>%
  mutate(tdp_kd_regulation = case_when(
    DE_Log2FC > 0 & prot_Log2FC > 0 ~ "DE_up / Prot_up",
    DE_Log2FC > 0 & prot_Log2FC < 0 ~ "DE_up / Prot_down",
    DE_Log2FC < 0 & prot_Log2FC > 0 ~ "DE_down / Prot_up",
    DE_Log2FC < 0 & prot_Log2FC < 0 ~ "DE_down / Prot_down"
  )) %>%
  count(tdp_kd_regulation) %>%
  mutate(
    percent = n / sum(n) * 100,
    tdp_kd_regulation = fct_reorder(tdp_kd_regulation, percent)
  )  %>% arrange(desc(percent)) %>%
  # mutate(id = row_number())
  mutate(id = c(2,3,4,1))



ring_labels <- data.frame(
  x = rep(0.51, 4),  # just pick one category to anchor the labels
  y = c(25, 50, 75, 100),
  label = c("25%", "50%", "75%", "100%")
)

p1<-ggplot(x %>% mutate(back=rep(100,4)), aes(x = id, y = percent, 
                                          fill = tdp_kd_regulation,
)) +
  geom_col(width = 1, show.legend = TRUE, alpha=0.8, color="black") +
  geom_col(aes(x = id, y = back, 
               fill = tdp_kd_regulation), 
           width = 1, show.legend = TRUE, alpha=0.15) +
  geom_hline(
    aes(yintercept = y),
    data.frame(y = seq(0, 100, 25)),
    color = "lightgrey"
  ) +
  geom_text(
    aes(label = paste0(round(percent), "%"), y = percent + 7),
    color = "black",
    size = 4.5
  ) +
  coord_polar(start = pi / 4, clip = "off") +  # Rotate to center first bar
  # scale_x_continuous(breaks = x$id, labels = x$tdp_kd_regulation) +
  scale_x_continuous(breaks = x$id, labels = c(
    "Downregulated mRNA,\nDownregulated Protein",
    "Downregulated mRNA,\nUpregulated Protein",
    "Upregulated mRNA,\nUpregulated Protein",  
    "Upregulated mRNA,\nDownregulated Protein"
  )) +
  theme_void()+
  theme(
    # Remove axis ticks and text
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(
      color = "gray12",
      size = 12),
    legend.position = "bottom",
    )+
  scale_fill_manual(values=
                      c(  "#9B5EA3","#CA0020", "#836FA8", "#0571B0")
  )+ 
    geom_text(
    data = ring_labels,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    angle = 0,
    color = "gray12",
    size = 3
  )+
  guides(fill="none")


ggsave(
  plot=p1,
  filename = file.path(figure_out_path, "Figure_5d.png"),
  width = 3900,
  height = 3250,
  dpi = 500,
  units = 'px'
)
