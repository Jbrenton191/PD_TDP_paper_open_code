library(magrittr)
library(patchwork)
library(tidyverse)
library(foreach)
library(doParallel)
library(ggh4x)

# 0 Setup ----

main_path<-here::here() # main project level i.e TDP43_paper_open_code level

data_path<-file.path(main_path, "Data")

figure_out_path<-file.path(main_path, "Plots/Main_Figures/Figure_3")

# create figure output directory if not present
if (!dir.exists(figure_out_path)) {
  dir.create(figure_out_path, recursive = T)
}

analysis_colours<-c("#0c5394", "#e06666")

stage_colours<-c("#91D1C2", "#4B8D86")

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

## TDP target download ----

# tdp 43 targets - exon to exon junction format
tdp_junc_list<-as_tibble(readRDS(file = file.path(data_path, "Significant_TDP_KD_Database_ex_ex_annotated_df.RDS")))
# Keep only distinct significant junctions across experiments
tdp_junc_list %<>% filter(p.adjust<0.05) %>%
  mutate(start=as.character(start), end=as.character(end)) %>%
  distinct(start, end,chr, strand,.keep_all = T) 


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# LRRK2 p.G2019S post-mortem brain analysis: Junctions ----

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# 
# Load LRRK2 G2019S post-mortem brain splicing data
# Using 4 covariates - see QC steps - select that for the cluster and intron usage dataframes
lrrk2_full_results<-readRDS(file = file.path(data_path, "LRRK2_PM_HS/lrrk2_pm_hs_results.RDS"))

# Extract splicing function: ----

fisher_test_junctions<-function(results, dataset, 
                                tdp_junc_list){
  # generate dataframe  
  lc_all<-results$cluster_significance
  int_all<-results$intron_usage
  
  all_results <- lc_all %>% 
    # filter for successful tests and remove clusters overlapping multiple genes
    dplyr::filter(!str_detect(genes, ",")) %>%
    tidyr::separate(cluster, into = c("chr", "cluster"), sep = ":") %>% 
    dplyr::inner_join(int_all %>% 
                        tidyr::separate(intron, into = c("chr", "start", "end", "cluster"), sep = ":")) %>% 
    dplyr::mutate(direction_of_effect = case_when(deltapsi > 0 ~ "upregulated",
                                                  deltapsi < 0 ~ "downregulated",
                                                  deltapsi == 0 ~ "no_change"))
  
  all_results$strand<-gsub(x = all_results$cluster, pattern = "clu.+_(.)", replacement = "\\1")
  all_results$cluster<-gsub(x = all_results$cluster, pattern = "(clu.+_).", replacement = "\\1")
  
  # adjust to get exon start
  all_results$start<-as.character(as.numeric(all_results$start)-1)
  
  # get universe of any junction tested in either brain area
  universe<-all_results
  full_junc_universe_together<-all_results
  full_junc_universe_together %<>% distinct(start, end, chr, .keep_all = T)
  
  sig_leaf <- all_results %>% filter(p.adjust<0.05)
  
  
  # fisher tests ----
  
  fisher_res_df_merge<-c()
  
  fisher_res_df<-c()
  
  # run across each brain area/comparison
  for (comp in unique(sig_leaf$comparison)) {
    
    # number of TDP43 matches in sig leafcutter list
    sig_in<-sig_leaf %>% filter(comparison==comp) %>% 
      distinct(start, end, chr, strand, .keep_all = T) %>% 
      inner_join(tdp_junc_list %>% distinct(start, end, chr, strand,.keep_all = T), by=c('start', 'end', "chr", "strand"))
    
    # number of TDP43 matches NOT in sig leafcutter list
    sig_out<-sig_leaf %>% filter(comparison==comp) %>% distinct(start, end, chr,strand, .keep_all = T) %>% 
      anti_join(tdp_junc_list %>% distinct(start, end, chr, strand,.keep_all = T), by=c('start', 'end', "chr", "strand")) 
    
    # number of TDP43 matches in leafcutter universe list
    back_in<-  full_junc_universe_together %>% 
      # filter(comparison==comp) %>% 
      distinct(start, end, chr,strand, .keep_all = T) %>%
      inner_join(tdp_junc_list %>%  distinct(start, end, chr, strand,.keep_all = T), by=c('start', 'end', "chr", "strand"))
    
    # number of TDP43 matches NOT in leafcutter universe list
    back_out<-full_junc_universe_together %>% 
      # filter(comparison==comp) %>%
      distinct(start, end, chr,strand, .keep_all = T) %>%
      # removed filter to keep list larger - rather than brain area specific list - can revert back - but this is same as GO
      anti_join(tdp_junc_list %>%  distinct(start, end, chr, strand,.keep_all = T), by=c('start', 'end', "chr", "strand"))
    
    dat <- data.frame(
      "In TDP Target List" = c(dim(sig_in)[1], dim(back_in)[1]),
      "Not In TDP Target List" = c(dim(sig_out)[1], dim(back_out)[1]),
      row.names = c( "Significantly Spliced Junction List", "Background Junctions"),
      stringsAsFactors = FALSE
    )
    
    test <- fisher.test(dat, alternative = 'greater')
    
    
    fisher_res<-tibble(fisher_statistic=test$estimate,
                       fisher_conf_int=list(test$conf.int),
                       fisher_p_value=test$p.value, comparison=comp,
                       tdp_junction_ratio=dim(sig_in)[1]/dim(sig_out)[1])
    
    fisher_res_df<-rbind(fisher_res_df, fisher_res)
    
  }
  
  
  fisher_res_df_merge<-rbind(fisher_res_df,fisher_res_df_merge)
  
  # FDR correcting for two brain areas tested
  fisher_res_df_merge$fdr<-NA
  for (i in 1:length(fisher_res_df_merge$fdr)) {
    fisher_res_df_merge$fdr[i]<-p.adjust(p = fisher_res_df_merge$fisher_p_value[i], 
                                         n=length(fisher_res_df$fisher_p_value),
                                         method = "fdr")
  }
  
  
  fisher_res_df_merge$leafcutter_dataset<-dataset
  
  return(fisher_res_df_merge)
  
}

lrrk2_fisher_res<-fisher_test_junctions(results = lrrk2_full_results, 
                                        dataset = "LRRK2 p.G2019S", 
                                        tdp_junc_list =  tdp_junc_list)

lrrk2_fisher_res <-
  lrrk2_fisher_res %>% mutate(
    comparison = case_when(
        comparison == "CC" ~ "Ant.\nCingulate",
            comparison == "FC" ~ "Frontal",
      TRUE ~ "DOUBLE CHECK"
    )
  )

lrrk2_fisher_res$comparison <-
  factor(
    lrrk2_fisher_res$comparison,
    levels = c("Ant.\nCingulate",
               "Frontal"
               # , "Combined"
               )
  )

lrrk2_fisher_res %<>% mutate(stage_colour=case_when(
fisher_p_value < 0.05 ~ "#125C5D",
 fisher_p_value > 0.05 ~ "grey"
),
is_sig=case_when(fisher_p_value < 0.05 ~ TRUE,
                 fisher_p_value > 0.05 ~ FALSE)
)


# dot graphs
lrrk2_fisher_res$experiment<-"Post-mortem brain"

p1<-ggplot(data = lrrk2_fisher_res,
           aes(x = comparison , y = -log10(fdr),
               colour=stage_colour, 
               # size=round(tdp_junction_ratio*100,3), 
               alpha=is_sig
           ))+
  geom_point( size=8)+
  ggh4x::facet_grid2(~experiment,
                     scales = "free_x", 
                     strip = ggh4x::strip_themed(background_y = elem_list_rect(fill = "white"),
                                          background_x =
                                            elem_list_rect(fill = scales::alpha("grey", 0.2))))+
  geom_hline(yintercept = -log10(0.05), linetype='dashed', colour='black')+
  labs(x="", y=expression("-log"[10]*"(FDR)"), size="TDP KD\nEnrichment\nRatio")+
  scale_alpha_manual(values = c(0.8))+
  scale_color_identity()+
  scale_alpha_manual(values = c(0.8))+ # only one alpha as all sig
  guides(
    size="none",
    color = "none",   # Remove legend for color
    shape = "none",
    alpha="none"
  )+
  theme_jb()+

theme_bw(base_family = "Roboto")+
  theme(axis.ticks.x = element_blank(),
        plot.title = element_text(size=18, hjust = 0.5),
        legend.position = "right",
        axis.title.y = element_text(size = 16),
        text = element_text(size = 18),
        axis.text.x = element_blank(),
        strip.text = element_text(size=16),
        strip.background = element_rect(fill = scales::alpha("grey", 0.2)),
        panel.spacing = unit(2, "lines")
  )


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# LRRK2 p.G2019S post-mortem brain analysis: Genes ----

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


fisher_test_genes<-function(results, dataset, 
                            tdp_junc_list){
  # generate dataframe  
  lc_all<-results$cluster_significance
  int_all<-results$intron_usage
  
  all_results <- lc_all %>% 
    dplyr::filter(!str_detect(genes, ",")) %>%
    tidyr::separate(cluster, into = c("chr", "cluster"), sep = ":") %>% 
    dplyr::inner_join(int_all %>% 
                        tidyr::separate(intron, into = c("chr", "start", "end", "cluster"), sep = ":")) %>% 
    dplyr::mutate(direction_of_effect = case_when(deltapsi > 0 ~ "upregulated",
                                                  deltapsi < 0 ~ "downregulated",
                                                  deltapsi == 0 ~ "no_change"))
  
  all_results$strand<-gsub(x = all_results$cluster, pattern = "clu.+_(.)", replacement = "\\1")
  all_results$cluster<-gsub(x = all_results$cluster, pattern = "(clu.+_).", replacement = "\\1")
  
  # adjust to get exon start
  all_results$start<-as.character(as.numeric(all_results$start)-1)
  
  # find any gene with junction tested
  universe<-all_results %>% pull(genes)
  
  sig_leaf <- all_results %>% filter(p.adjust<0.05)
  
  # fisher tests ----
  fisher_res_df_merge<-c()
 
  # Running over each comparison
  
  fisher_res_df<-c()
  
  for (comp in unique(sig_leaf$comparison)) {
    
    sig_leaf <- all_results %>% filter(p.adjust<0.05) %>% filter(comparison==comp)
    
    
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
                      stringsAsFactors = FALSE)
    
    test <- fisher.test(dat, alternative = 'greater')
    
    fisher_res<-tibble(fisher_statistic=test$estimate,
                       fisher_conf_int=list(test$conf.int),
                       fisher_p_value=test$p.value, 
                       comparison=comp,
                       tdp_junction_ratio=length(sig_in)/length(sig_out))
    
    fisher_res_df<-rbind(fisher_res_df, fisher_res)
    
  }
  
  
  fisher_res_df_merge<-rbind(fisher_res_df,fisher_res_df_merge)
  
  # only FDR correcting the mass brain area p value
  fisher_res_df_merge$fdr<-NA
  for (i in 1:length(fisher_res_df_merge$fdr)) {
    fisher_res_df_merge$fdr[i]<-p.adjust(p = fisher_res_df_merge$fisher_p_value[i], 
                                         n=length(fisher_res_df$fisher_p_value),
                                         method = "fdr")
  }
  
  
  fisher_res_df_merge$leafcutter_dataset<-dataset
  
  return(fisher_res_df_merge)
  
}

lrrk2_fisher_res_genes<-fisher_test_genes(results = lrrk2_full_results, 
                                        dataset = "LRRK2 p.G2019S", 
                                        tdp_junc_list =  tdp_junc_list)

lrrk2_fisher_res_genes <-
  lrrk2_fisher_res_genes %>% mutate(
    comparison = case_when(
      comparison == "CC" ~ "Ant.\nCingulate",
      comparison == "FC" ~ "Frontal",
      TRUE ~ "DOUBLE CHECK"
    )
  )

lrrk2_fisher_res_genes$comparison <-
  factor(
    lrrk2_fisher_res_genes$comparison,
    levels = c("Ant.\nCingulate",
               "Frontal"
               )
  )

lrrk2_fisher_res_genes %<>% mutate(stage_colour=case_when(
  fisher_p_value < 0.05 ~ "#125C5D",
  fisher_p_value > 0.05 ~ "grey"
),
is_sig=factor(case_when(fisher_p_value < 0.05 ~ TRUE,
                 fisher_p_value > 0.05 ~ FALSE), levels=c(FALSE, TRUE))
)

lrrk2_fisher_res_genes$experiment<-"Post-mortem brain"

p2<-ggplot(data = lrrk2_fisher_res_genes,
           aes(x = comparison , y = -log10(fdr),
               colour=stage_colour, 
               alpha=is_sig
           ))+
  geom_point(size=8)+
  ggh4x::facet_grid2(~experiment,
                     scales = "free_x", 
                     strip = ggh4x::strip_themed(background_y = elem_list_rect(fill = "white"),
                                          background_x =
                                            elem_list_rect(fill = scales::alpha("grey", 0.2))))+
   geom_hline(yintercept = -log10(0.05), linetype='dashed', colour='black')+
  labs(x="", y=expression("-log"[10]*"(FDR)"), size="TDP KD\nEnrichment\nRatio")+
  scale_color_identity()+
  scale_alpha_manual(values = c(0.8))+
  # scale_size_binned(n.breaks = 3)+
  guides(
    size="none",
    color = "none",   # Remove legend for color
    shape = "none",
    alpha="none"
  )+
  theme_jb()+
  # ggtitle("Enrichment of TDP-43 KD Genes")+
  theme_bw(base_family = "Roboto")+
  theme(
    text = element_text(size = 18),
    # plot.title = element_text(size=18, hjust = 0.5),
    legend.position = "right",
    # strip.background = element_rect(fill = stage_colours[c(1,2)]),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(
      angle = 60,
      hjust = 1.05,
      vjust = 1.05,
      size = 18
    ),
      panel.spacing = unit(2, "lines"),
    strip.background = element_blank(),
    strip.text.y = element_blank(),
    strip.text =  element_blank()
  )

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# mDN1: Team Rio human embryonic stem cell LRRK2 G2019S mDN model Results - Junctions: ####

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


lrrk2_hpsc_results_1<-readRDS(file = file.path(data_path, 'LRRK2_HS_PSC/mDN1/lrrk2_hpsc_results_1.RDS'))

# selected from QC data: 3 covariates

# Fisher test code run again here - bit slimmed down and not run as a function as only need to run once as no multiple brain areas

# generate dataframe  
lc_all<-lrrk2_hpsc_results_1cluster_significance
int_all<-lrrk2_hpsc_results_1$intron_usage

all_results <- lc_all %>% 
  # filter for successful tests and remove clusters overlapping multiple genes
  dplyr::filter(!str_detect(genes, ",")) %>%
  tidyr::separate(cluster, into = c("chr", "cluster"), sep = ":") %>% 
  dplyr::inner_join(int_all %>% 
                      tidyr::separate(intron, into = c("chr", "start", "end", "cluster"), sep = ":")) %>% 
  dplyr::mutate(direction_of_effect = case_when(deltapsi > 0 ~ "upregulated",
                                                deltapsi < 0 ~ "downregulated",
                                                deltapsi == 0 ~ "no_change"))

all_results$strand<-gsub(x = all_results$cluster, pattern = "clu.+_(.)", replacement = "\\1")
all_results$cluster<-gsub(x = all_results$cluster, pattern = "(clu.+_).", replacement = "\\1")

# adjust to get exon start
all_results$start<-as.character(as.numeric(all_results$start)-1)

universe<-all_results
full_junc_universe_together<-all_results
full_junc_universe_together %<>% distinct(start, end, chr, .keep_all = T)

sig_leaf <- all_results %>% filter(p.adjust<0.05)

# fisher tests ----
fisher_res_df_merge<-c()
fisher_res_df<-c()

# number of TDP43 matches in sig leafcutter list
sig_in<-sig_leaf %>% distinct(start, end, chr,strand, .keep_all = T) %>%
  inner_join(tdp_junc_list %>% distinct(start, end, chr,strand), by=c('start', 'end', "chr", "strand"))

# number of TDP43 matches NOT in sig leafcutter list
sig_out<-sig_leaf %>% distinct(start, end, chr,strand, .keep_all = T) %>%
  anti_join(tdp_junc_list %>% distinct(start, end, chr,strand), by=c('start', 'end', "chr", "strand"))

# number of TDP43 matches in leafcutter universe list
back_in<-full_junc_universe_together %>%
  inner_join(tdp_junc_list %>%
               distinct(start, end, chr,strand), by=c('start', 'end', "chr", "strand"))

# number of TDP43 matches NOT in leafcutter universe list
back_out<-full_junc_universe_together %>%
  anti_join(tdp_junc_list %>%
              distinct(start, end, chr,strand), by=c('start', 'end', "chr", "strand")) 

dat <- data.frame("In TDP target list" = c(dim(sig_in)[1], dim(back_in)[1]),
                  "Not in TDP target list" = c(dim(sig_out)[1], dim(back_out)[1]),
                  row.names = c("Genes with Differentially Spliced Clusters", "Background Genes (Any cluster found across brain tissues)"),
                  stringsAsFactors = FALSE
)

test <- fisher.test(dat, alternative = 'greater')

fisher_res<-tibble(fisher_statistic=test$estimate,
                   fisher_conf_int=list(test$conf.int),
                   fisher_p_value=test$p.value,  comparison="mDN1", 
                   tdp_junction_ratio=dim(sig_in)[1]/dim(sig_out)[1])


junc_res_lrrk2_hpsc_1<-fisher_res

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# mDN1: Team Rio human embryonic stem cell LRRK2 G2019S mDN model Results - Genes: ####

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

universe<-all_results %>% pull(genes)

sig_leaf <- all_results %>% filter(p.adjust<0.05)

# fisher tests ----
fisher_res_df_merge<-c()
fisher_res_df<-c()

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
                   fisher_p_value=test$p.value,  comparison="mDN1", 
                   tdp_junction_ratio=length(sig_in)/length(sig_out))

gene_res_lrrk2_hpsc_1<-fisher_res

# ****************************************************************************************************************************************

# mDN2: Dataset from Gandhi lab - hIPSC derived LRRK2 pG2019S midbrain neurons - Junctions: ----

# ****************************************************************************************************************************************

lrrk2_hpsc_results_2<-readRDS(file = file.path(data_path, 'LRRK2_HS_PSC/mDN2/combGrp_clust_junct_anno_singleGene.tsv'))


lrrk2_hpsc_results_2 %<>% separate(col = intron, into = c("chr", "start", "end", "cluster"), sep=":")

lrrk2_hpsc_results_2$strand<-gsub(x = lrrk2_hpsc_results_2$cluster, pattern = "clu.+_(.)", replacement = "\\1")
lrrk2_hpsc_results_2$cluster<-gsub(x = lrrk2_hpsc_results_2$cluster, pattern = "(clu.+)_.", replacement = "\\1")

# adjust to get exon start
lrrk2_hpsc_results_2$start<-as.character(as.numeric(lrrk2_hpsc_results_2$start)-1)

universe<-lrrk2_hpsc_results_2
full_junc_universe_together<-lrrk2_hpsc_results_2
full_junc_universe_together %<>% distinct(start, end, chr, .keep_all = T)

sig_leaf <- lrrk2_hpsc_results_2 %>% filter(p.adjust<0.05)

# fisher tests ----
## Collapsed result
fisher_res_df_merge<-c()
fisher_res_df<-c()

# number of TDP43 matches in sig leafcutter list
sig_in<-sig_leaf %>% distinct(start, end, chr,strand, .keep_all = T) %>%
  inner_join(tdp_junc_list %>% distinct(start, end, chr,strand), by=c('start', 'end', "chr", "strand"))

# number of TDP43 matches NOT in sig leafcutter list
sig_out<-sig_leaf %>% distinct(start, end, chr,strand, .keep_all = T) %>%
  anti_join(tdp_junc_list %>% distinct(start, end, chr,strand), by=c('start', 'end', "chr", "strand"))

# number of TDP43 matches in leafcutter universe list
back_in<-full_junc_universe_together %>%
  inner_join(tdp_junc_list %>%
               distinct(start, end, chr,strand), by=c('start', 'end', "chr", "strand"))

# number of TDP43 matches NOT in leafcutter universe list
back_out<-full_junc_universe_together %>%
  anti_join(tdp_junc_list %>%
              distinct(start, end, chr,strand), by=c('start', 'end', "chr", "strand")) 

dat <- data.frame("In TDP target list" = c(dim(sig_in)[1], dim(back_in)[1]),
                  "Not in TDP target list" = c(dim(sig_out)[1], dim(back_out)[1]),
                  row.names = c("Genes with Differentially Spliced Clusters", "Background Genes (Any cluster found across brain tissues)"),
                  stringsAsFactors = FALSE
)

test <- fisher.test(dat, alternative = 'greater')

fisher_res<-tibble(fisher_statistic=test$estimate,
                   fisher_conf_int=list(test$conf.int),
                   fisher_p_value=test$p.value, comparison="mDN2", 
                   tdp_junction_ratio=dim(sig_in)[1]/dim(sig_out)[1])

junc_res_lrrk2_hpsc_2<-fisher_res


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# mDN2: Dataset from Gandhi lab - hIPSC derived LRRK2 pG2019S midbrain neurons - Genes: ----

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

universe<-lrrk2_hpsc_results_2 %>% pull(genes)


sig_leaf <- lrrk2_hpsc_results_2 %>% filter(p.adjust<0.05)

# fisher tests ----
fisher_res_df_merge<-c()
fisher_res_df<-c()

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
                   fisher_p_value=test$p.value,  comparison="mDN2", 
                   tdp_junction_ratio=length(sig_in)/length(sig_out))

gene_res_lrrk2_hpsc_2<-fisher_res


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

### Merge results hPSCs ----

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


junc_res_lrrk2_hpsc_both<-rbind(junc_res_lrrk2_hpsc_1, junc_res_lrrk2_hpsc_2) %>% mutate(stage_colour=case_when(
  fisher_p_value < 0.05 ~ "#125C5D",
  fisher_p_value > 0.05 ~ "grey"
),
is_sig=case_when(fisher_p_value < 0.05 ~ TRUE,
                 fisher_p_value > 0.05 ~ FALSE)
)

junc_res_lrrk2_hpsc_both$experiment<-"mDN"


# merged hPSC junction graph
hpsc_pj<-ggplot(data = junc_res_lrrk2_hpsc_both %>% mutate(Analysis="Junction"),
       aes(x = comparison , y = -log10(fisher_p_value),
           colour=stage_colour, 
           alpha=is_sig
       ))+
  ggh4x::facet_grid2(Analysis~experiment,
                     scales = "free_x", 
                     strip = strip_themed(background_y = elem_list_rect(fill = "white"),
                                          background_x =
                                            elem_list_rect(fill = scales::alpha("grey", 0.2))))+
  # facet_wrap(~experiment, scales = "free_x")+
  geom_point( size=8)+
  geom_hline(yintercept = -log10(0.05), linetype='dashed', colour='black')+
  labs(x="", y=expression("-log"[10]*"(p)"), size="TDP KD\nEnrichment\nRatio")+
  scale_alpha_manual(values = c(0.8))+
  scale_color_identity()+
  scale_alpha_manual(values = c(0.8))+ # only one alpha as all sig
  # scale_size_binned(n.breaks = 3)+
  guides(
    size="none",
    color = "none",   # Remove legend for color
    shape = "none",
    alpha="none"
  )+
  theme_jb()+
  # ggtitle("Enrichment of TDP-43 KD Genes")+
  theme_bw(base_family = "Roboto")+
   theme(axis.ticks.x = element_blank(),
        plot.title = element_text(size=18, hjust = 0.5),
        legend.position = "right",
        text = element_text(size = 18),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_blank(),
        strip.text.x =  element_text(size = 16),
        strip.text.y =  element_text(size = 14),
        # strip.background = element_rect(fill = scales::alpha("grey", 0.2)),
        panel.spacing = unit(2, "lines")
  )

gene_res_lrrk2_hpsc_both<-rbind(gene_res_lrrk2_hpsc_1, gene_res_lrrk2_hpsc_2) %>% mutate(stage_colour=case_when(
  fisher_p_value < 0.05 ~ "#125C5D",
  fisher_p_value > 0.05 ~ "grey"
),
is_sig=case_when(fisher_p_value < 0.05 ~ TRUE,
                 fisher_p_value > 0.05 ~ FALSE)
)

gene_res_lrrk2_hpsc_both$experiment<-"mDN"

# merged hPSC gene graph

hpsc_pg<-ggplot(data = gene_res_lrrk2_hpsc_both %>% mutate(Analysis="Gene"),
       aes(x = comparison , y = -log10(fisher_p_value),
           colour=stage_colour, 
           # size=round(tdp_junction_ratio*100,3), 
           alpha=is_sig
       ))+
  ggh4x::facet_grid2(Analysis~experiment,
                     scales = "free_x", 
                     strip = strip_themed(background_y = elem_list_rect(fill = "white"),
                                          background_x =
                                            elem_list_rect(fill = scales::alpha("grey", 0.2))))+
  geom_point( size=8)+
  geom_hline(yintercept = -log10(0.05), linetype='dashed', colour='black')+
  labs(x="", y=expression("-log"[10]*"(p)"), size="TDP KD\nEnrichment\nRatio")+
  scale_alpha_manual(values = c(0.8))+
  scale_color_identity()+
  scale_alpha_manual(values = c(0.8))+ # only one alpha as all sig
  # scale_size_binned(n.breaks = 3)+
  guides(
    size="none",
    color = "none",   # Remove legend for color
    shape = "none",
    alpha="none"
  )+
  theme_jb()+
  # ggtitle("Enrichment of TDP-43 KD Genes")+
  theme_bw(base_family = "Roboto")+
  theme( axis.text.x = element_text(
    angle = 60,
    hjust = 1.05,
    vjust = 1.05,
    size = 18
  ),
        plot.title = element_text(size=18, hjust = 0.5),
        legend.position = "right",
        text = element_text(size = 18),
  axis.title.y = element_text(size = 16),
        # axis.text.x = element_blank(),
        strip.text.x =  element_blank(),
        strip.text.y =  element_text(size = 14),
        # strip.background = element_rect(fill = scales::alpha("grey", 0.2)),
        panel.spacing = unit(2, "lines")
  )


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Putting plots together with different y axes labels ----

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

p3<-((p1+theme(plot.margin = unit(c(0,0,50,0), "pt")))/p2+plot_layout(axes = 'collect'))

p4<-((hpsc_pj+theme(plot.margin = unit(c(0,0,50,0), "pt")))/hpsc_pg+plot_layout(axes = 'collect'))



ggsave(
  plot=p3+theme(plot.margin = unit(c(0,50,0,0), "pt"))|p4,
  filename = file.path(figure_out_path, "Figure_3b.png"),
  width = 7000,
  height = 6000,
  dpi = 550,
  units = 'px'
)


