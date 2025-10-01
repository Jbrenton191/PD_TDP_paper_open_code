library(magrittr)
library(patchwork)
library(tidyverse)
library(foreach)
library(doParallel)

main_path<-here::here() # main project level i.e TDP43_paper_open_code level

data_path<-file.path(main_path, "Data")

figure_out_path<-file.path(main_path, "Plots/Main_Figures/Figure_4")

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


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# 1. Load splicing data and TDP KD database ----

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# Load PD mouse leafcutter results
PD_mice_full_res <- readRDS(file.path(data_path, "LRRK2_MM/LRRK2_vs_WT_leafcutter_results.RDS"))

# load all successfully tested clusters and introns, i.e those that pass leafcutter QC filters.
all_results <- PD_mice_full_res$all_successful_clusters

all_results$strand<-gsub(x = all_results$cluster, pattern = "clu.+_(.)", replacement = "\\1")
all_results$cluster<-gsub(x = all_results$cluster, pattern = "(clu.+_).", replacement = "\\1")

# adjust to get exon start.
all_results$start<-as.character(as.numeric(all_results$start)-1)
# get fdr significant clusters
sig_leaf <- all_results %>% filter(p.adjust<0.05)

# get background cluster/gene list
full_junc_universe_together<-all_results
full_junc_universe_together %<>% distinct(start, end, chr, .keep_all = T)
PD_mouse_bg_uni<-full_junc_universe_together

## TDP targets download / prepare TDP KD list
# Load tdp43 kd junctions in mouse neuronal cell lines
tdp_junc_list<-readRDS(file.path(data_path, "TDP-43_KD_MM/TDP_KD_DB_Mouse.Rds"))

tdp_junc_list$start<-as.character(as.numeric(tdp_junc_list$start)-1)
tdp_junc_list$strand<-gsub(x = tdp_junc_list$cluster, pattern = "clu.+_(.)", replacement = "\\1")
tdp_junc_list$cluster<-gsub(x = tdp_junc_list$cluster, pattern = "(clu.+_).", replacement = "\\1")

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# 2. Run Enrichment tests ----

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## 2.1 Junction Enrichment ---- 

  
  # number of TDP43 matches in sig leafcutter list
  sig_in<-sig_leaf %>% distinct(start, end, chr, .keep_all = T) %>% inner_join(tdp_junc_list %>% 
                                                                                 distinct(start, end, chr, .keep_all = T), by=c('start', 'end', "chr"))
  
  # number of TDP43 matches NOT in sig leafcutter list
  sig_out<-sig_leaf %>% distinct(start, end, chr, .keep_all = T) %>% anti_join(tdp_junc_list %>% 
                                                                                 distinct(start, end, chr, .keep_all = T), by=c('start', 'end', "chr")) %>% arrange(genes)
  
  # number of TDP43 matches in leafcutter universe list
  back_in<-PD_mouse_bg_uni %>% 
    distinct(start, end, chr, .keep_all = T) %>%
    # removed filter to keep list larger - rather than brain area specific list - can revert back - but this is same as GO
    inner_join(tdp_junc_list %>%  
                 distinct(start, end, chr, .keep_all = T), by=c('start', 'end', "chr"))
  
  # number of TDP43 matches NOT in leafcutter universe list
  back_out<-PD_mouse_bg_uni %>%
    distinct(start, end, chr, .keep_all = T) %>%
    # removed filter to keep list larger - rather than brain area specific list - can revert back - but this is same as GO
    anti_join(tdp_junc_list %>%  
                distinct(start, end, chr,.keep_all = T), by=c('start', 'end', "chr"))
  
  dat <- data.frame(
    "In TDP Target List" = c(dim(sig_in)[1], dim(back_in)[1]),
    "Not In TDP Target List" = c(dim(sig_out)[1], dim(back_out)[1]),
    row.names = c( "Significantly Spliced Junction List", "Background Junctions"),
    stringsAsFactors = FALSE
  )
  
  print(knitr::kable(dat))
  
  test <- fisher.test(dat, alternative = 'greater')
  
  # test$p.value
  
  fisher_res<-tibble(fisher_statistic=test$estimate,
                     fisher_conf_int=list(test$conf.int),
                     fisher_p_value=test$p.value,
                     tdp_junction_ratio=length(sig_in)/length(sig_out))
  
  fisher_res_df_junctions<-fisher_res
  
  fisher_res_df_junctions$g_or_j<-"Junction"
  

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  

## 2.2 Gene Enrichment ----

  # get universe of genes
universe<-all_results %>% pull(genes)
# full_junc_universe_together<-all_results
# full_junc_universe_together %<>% distinct(start, end, chr, .keep_all = T)

# number of TDP43 matches in sig leafcutter list
sig_in<- unique(sig_leaf$genes)[unique(sig_leaf$genes) %in% unique(tdp_junc_list$genes)]

# number of TDP43 matches NOT in sig leafcutter list
sig_out<-unique(sig_leaf$genes)[!unique(sig_leaf$genes) %in% unique(tdp_junc_list$genes)]

# number of TDP43 matches in leafcutter universe list
back_in<-unique(universe)[unique(universe) %in% unique(tdp_junc_list$genes)]
# number of TDP43 matches NOT in leafcutter universe list
back_out<-unique(universe)[!unique(universe) %in% unique(tdp_junc_list$genes)]

dat <- data.frame("In TDP target list" = c(length(sig_in), length(back_in)),
                  "Not in TDP target list" = c(length(sig_out), length(back_out)),
                  row.names = c("Genes with Differentially Spliced Clusters", "Background Genes (Any cluster found across brain tissues)"),
                  stringsAsFactors = FALSE
)

test <- fisher.test(dat, alternative = 'greater')

fisher_res<-tibble(fisher_statistic=test$estimate,
                   fisher_conf_int=list(test$conf.int),
                   fisher_p_value=test$p.value,
                   tdp_junction_ratio=length(sig_in)/length(sig_out))


fisher_res_df_genes<-fisher_res

fisher_res_df_genes$g_or_j<-"Gene"


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# 3. Graph Results ----

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# merge gene and junction results
fisher_res_merge<-rbind(fisher_res_df_junctions, fisher_res_df_genes)

# add colour and significance columns
fisher_res_merge <-
  fisher_res_merge %>% mutate(
    comparison = "LRRK2 G2019S\nKI vs WT"
    ) %>% mutate(stage_colour=case_when(
  fisher_p_value < 0.05 ~ "#125C5D",
  fisher_p_value > 0.05 ~ "grey"
),
is_sig=factor(case_when(fisher_p_value < 0.05 ~ TRUE,
                        fisher_p_value > 0.05 ~ FALSE), levels=c(FALSE, TRUE))
)

# set order to appear in
fisher_res_merge$g_or_j<-factor(fisher_res_merge$g_or_j, levels=c("Junction", "Gene"))

# graph
p2<-
  ggplot(data = fisher_res_merge,
           aes(x = comparison , y = -log10(fisher_p_value),
               colour=stage_colour, 
               alpha=is_sig
           ))+
  facet_wrap(~g_or_j)+
  geom_point(size=8)+
  geom_hline(yintercept = -log10(0.05), linetype='dashed', colour='black')+
  labs(x="", y=expression("-log"[10]*"(p)"), size="TDP KD\nEnrichment\nRatio")+
  scale_color_identity()+
  scale_alpha_manual(values = c(0.8))+
  guides(
    size="none",
    color = "none",   # Remove legend for color
    shape = "none",
    alpha="none"
  )+
  theme_jb()+
  theme_bw(base_family = "Roboto")+
  theme(
    text = element_text(size = 16),
    legend.position = "right",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.spacing = unit(3, "lines")
  )


ggsave(
  plot=p2,
  filename = file.path(figure_out_path, "Figure_4b.png"),
  width = 2000,
  height = 1500,
  dpi = 350,
  units = 'px'
)
