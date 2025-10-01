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

lrrk2_group_order<-c("Control", "LRRK2\nG2019S")


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### ++++++ ++++++ ++++++ ++++++ ++++++ ++++++ ++++++ ++++++ ++++++ ++++++ ++++++ 
# Graphing ----
### ++++++ ++++++ ++++++ ++++++ ++++++ ++++++ ++++++ ++++++ ++++++ ++++++ ++++++ 


# get matching results from novel (non-uniprot database) MS ions aligned to translated novel transcript sequences
matching_results <-
  read_csv(
    file.path(
      data_path,
      "LRRK2_HS_PSC/mDN1/novel_peps_detection_rate_min_added_LIB_size_norm_Upreg_tdp_tsl_nmd_included_NOVEL_MS_seq_matching_to_novel_junction_glm_and_t.tests_23.6.25.csv"
    )
  )

# filter for amino acid sequences/MS ions that aren't too small or too large to generate artificial matches. 
# Also that have a delta alignment score (blast score to novel transcript peptide sequence - blast score to annotated transcript peptide sequence)
stat_results <- matching_results %>% 
  filter(
    number_matches>=7,
         number_matches<=30
  ) %>% 
  filter(delta_alignment_score_by_full_query_length>=0.5)

  
# Make graphing dataframe with only Significant Upregulated novel peptides
# If multiple possible significant novel peptides per gene - average of these is taken for Log2FC
top_proteins <- stat_results %>%
  filter(
    t.test_fdr<0.05,
    Log2FC>0
  ) %>%
  group_by(gene, transcript_status) %>%
  mutate(avg_Log2FC=mean(Log2FC)) %>%
  arrange(desc(avg_Log2FC)) %>%
  distinct(gene,.keep_all = T) %>%
  mutate(gene = factor(gene, levels = gene)) %>%
  mutate(junction_type=case_when(junction_type=="novel_donor" ~ "Novel donor",
                                 junction_type=="novel_acceptor" ~ "Novel acceptor",
  )) %>%
  mutate(junction_type=factor(junction_type, levels=c("Novel donor", "Novel acceptor"))) 


# Get ordered df by log2fc for plot
gene_order_df <- top_proteins %>%
  ungroup() %>%
  select(gene, avg_Log2FC) %>%
  arrange(avg_Log2FC)


# Create Log2FC heatmap of genes with significant novel peptides
p1<-top_proteins %>%
  mutate(blank="") %>% # Make a blank column title to use for x axis
  mutate(Gene_ordered = factor(gene, levels = gene_order_df$gene))  %>%
  ggplot(aes(x = Gene_ordered, y = junction_type, fill=Log2FC))+
  geom_tile( width = 1, show.legend = T)+
  geom_hline(yintercept = seq(0.5, length(unique(top_proteins$junction_type))+1), color="black", linewidth = 0.2) +
  geom_vline(xintercept = seq(0.5, length(top_proteins[, 'Log2FC', T])+1), color="black",
             linewidth = 0.1) +
  scale_y_discrete(limits=rev, expand = expansion(add = 0.5)) +
  scale_x_discrete(drop = F, expand = expansion(add = 0.5))+
  scale_fill_viridis_c(name = "log"[2]~"(Fold Change)",
                       na.value = "white",
                                             option = "viridis",
                       begin = 0.3,
                       end = 1
                       ) +
  coord_cartesian()+
  theme_void(base_family = "Roboto")+
  labs(x="", y="")+
  theme(legend.position = "top",
        axis.text.y = element_text(
          size = 12, margin = margin(r = 2)
        ),
        legend.key.size = unit(0.7, 'cm'),
        legend.text = element_text(size=9),
        legend.title = element_text(size=9),
        legend.spacing.x =  unit(0.7, 'cm'),
        legend.spacing.y = unit(0.7, 'cm')
  )+
  theme(plot.margin = margin(l = 8, r = 8, t=3, b=3, unit = "pt"),
        axis.text.x = element_text(angle = 45, vjust = 0.97, 
                                   hjust=0.97, margin = margin(t = 1)),
        )

# save graph
ggsave(
  plot=p1,
  filename = file.path(figure_out_path, "Figure_5b.png"),
  width = 8000,
  height = 1750,
  dpi = 600,
  units = 'px'
)

