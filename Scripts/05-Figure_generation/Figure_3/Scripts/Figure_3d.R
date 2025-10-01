
library(magrittr)
library(patchwork)
library(tidyverse)
library(foreach)
library(doParallel)
library(ggh4x)

main_path<-here::here() # main project level i.e TDP43_paper_open_code level

data_path<-file.path(main_path, "Data")

figure_out_path<-file.path(main_path, "Plots/Main_Figures/Figure_3")

# create figure output directory if not present
if (!dir.exists(figure_out_path)) {
  dir.create(figure_out_path, recursive = T)
}


analysis_colours<-c("#0c5394", "#e06666")

stage_colours<-c("#125C5D")

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

# 1. Load splicing data, merge brain areas and hPSC datasets, and find overlapping genes in TDP KD database ----

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# 1.1. Get LRRK2 HS DS Genes - combined/ unique genes across both areas

# Load LRRK2 G2019S post-mortem brain splicing data
# Using 4 covariates - see QC steps - select that for the cluster and intron usage dataframes
lrrk2_full_results<-readRDS(file = file.path(data_path, "LRRK2_PM_HS/lrrk2_pm_hs_results.RDS"))

results<-lrrk2_full_results

lc_all<-results$cluster_significance
int_all<-results$intron_usage

lrrk2_results <- lc_all %>% 
  # filter for successful tests and remove clusters overlapping multiple genes
  dplyr::filter(!str_detect(genes, ",")) %>%
  tidyr::separate(cluster, into = c("chr", "cluster"), sep = ":") %>% 
  dplyr::inner_join(int_all %>% 
                      tidyr::separate(intron, into = c("chr", "start", "end", "cluster"), sep = ":")) %>% 
  dplyr::mutate(direction_of_effect = case_when(deltapsi > 0 ~ "upregulated",
                                                deltapsi < 0 ~ "downregulated",
                                                deltapsi == 0 ~ "no_change"))

universe_hs_pm<-lrrk2_results %>% distinct(genes) %>% pull(genes)

# Taking significant genes across both brain areas - combining them
sig_leaf_hs_pm <- lrrk2_results %>% filter(p.adjust<0.05) %>% distinct(genes) %>% mutate(Dataset="LRRK2 p.G2019S\nPost-mortem\nbrain")

# 1.2.Load Rio team hPSC results

# Rio LRRK2 Results - Junctions:

lrrk2_hpsc_results_1<-readRDS(file = file.path(data_path, 'LRRK2_HS_PSC/mDN1/lrrk2_hpsc_results_1.RDS'))

lc_all<-lrrk2_hpsc_results_1$cluster_significance
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

universe_hpsc_1<-all_results %>% distinct(genes) %>% pull(genes)

sig_leaf_hpsc_1 <- all_results %>% filter(p.adjust<0.05) %>% distinct(genes) %>% mutate(Dataset="LRRK2 p.G2019S\nmDNs")

## 1.3. Load Gandhi lab IPSC results

lrrk2_ipsc_results_2<-read_tsv(file = "/home/jbrenton/TDP43_PD_Paper_Figures/Drafting_figures_code/Figure_3/Data/combGrp_clust_junct_anno_singleGene.tsv")


universe_ipsc_2<-lrrk2_ipsc_results_2 %>% distinct(genes) %>% pull(genes)


sig_leaf_ipsc_2 <- lrrk2_ipsc_results_2 %>% filter(p.adjust<0.05) %>% distinct(genes) %>% mutate(Dataset='LRRK2 p.G2019S\nmDNs')

## 1.4. Merge significantly spliced genes and universes together

sig_leaf_hpsc<-rbind(sig_leaf_hpsc_1, sig_leaf_ipsc_2) %>% distinct(genes,.keep_all = T)

lrrk2_sig_leaf<-rbind(sig_leaf_hs_pm, sig_leaf_hpsc)

universe_hs_pm<-tibble(genes=unique(c(universe_hs_pm)), Dataset="LRRK2 p.G2019S\nPost-mortem\nbrain")
universe_hpsc<-tibble(genes=unique(c(universe_hpsc_1, universe_ipsc_2)), Dataset='LRRK2 p.G2019S\nmDNs')

universe<-rbind(universe_hs_pm, universe_hpsc)

# 1.5. find overlap with genes changed upon TDP KD or absence from nucleus ----

tdp_junc_list<-readRDS(file = "/home/jbrenton/TDP_KD_HS/Merged_junction_database/Data/Significant_TDP_KD_Database_ex_ex_annotated_df.RDS")
# Distinct junctions across experiments
tdp_junc_list %<>% filter(p.adjust<0.05) %>%
  mutate(start=as.character(start), end=as.character(end)) %>%
  # mutate(location=str_c(chr, start, end, strand, sep = ":")) %>% # make one value for chr, start, end, strand - sanity check
  distinct(start, end,chr, strand,.keep_all = T) 

sig_overlap<- lrrk2_sig_leaf[lrrk2_sig_leaf$genes %in% unique(tdp_junc_list$genes),] %>% mutate(SYMBOL=genes)
## Checking PDPRS pathways in GO data overlap with sporadic splicing

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
# Load gencode gene map
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
gene_map_path <- file.path(data_path, "References/gencode_txid_to_geneid.txt")

gencode_txid_to_geneid <- vroom::vroom(gene_map_path) %>%
  `colnames<-`(c("tx_id", "gene_id","gene_name", "description")) %>%
  dplyr::mutate(tx_id = sub("\\..+", "", tx_id),
                gene_id = sub("\\..+", "", gene_id)) %>% 
  dplyr::select(gene_id, gene_name) %>% distinct(gene_id, .keep_all=T)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Load PDPRS pathways downloaded from https://pdgenetics.shinyapps.io/pathwaysbrowser/ , and adjusted into df format
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

PDPRSpathways_ensg <- readRDS(data_path, "Gene_lists/PDPRSpathways_ensg.rds")
# Convert to long tibble format
PDPRSpathways_ensg_tibble<-enframe(PDPRSpathways_ensg, 
                                   name = "Pathway", value = "Gene") %>% 
  unnest(Gene) %>% mutate(Gene_no_version=gsub(x=Gene, pattern="\\..*", replacement=""))

# Missing genes between the two:
PDPRSpathways_ensg_tibble[!PDPRSpathways_ensg_tibble$Gene_no_version %in% gencode_txid_to_geneid$gene_id,]
# works out as 15 distinct genes across 15 pathways - probably because using MSigDB database v7.2?

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# map gene names to Gene in PDPRS pathways
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

PDPRSpathways_ensg_tibble %<>% inner_join(., gencode_txid_to_geneid, by=c("Gene_no_version"= "gene_id"))


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Load sig data splicing data and get backgrounds ----
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

asap_sig_leaf<-sig_overlap

universe_df <- universe

fisher_res_df_merge<-c()

  
for (dataset in unique(asap_sig_leaf$Dataset)) {
  
  asap_sig_leaf2 <- asap_sig_leaf %>% filter(Dataset==dataset)
  
  asap_sig_leaf_filt <- asap_sig_leaf2 %>% distinct(genes) %>% pull(genes)
  
  universe_df_filt<-universe_df %>% filter(Dataset==dataset) %>% distinct(genes) %>% pull(genes)
  
  fisher_res_df<-c()
  
  for (i in unique(PDPRSpathways_ensg_tibble$Pathway)) {
    
    pathway<-PDPRSpathways_ensg_tibble %>% filter(Pathway==i) %>% pull(gene_name)
    
    # number of pathway matches in sig leafcutter list
    sig_in<-asap_sig_leaf_filt[asap_sig_leaf_filt %in% pathway]
    
    # number of pathway matches NOT in sig leafcutter list
    sig_out<-asap_sig_leaf_filt[!asap_sig_leaf_filt %in% pathway]
    
    # number of matches in leafcutter universe list
    back_in<-universe_df_filt[universe_df_filt %in% pathway]
    
    # number of matches NOT in leafcutter universe list
    back_out<-universe_df_filt[!universe_df_filt %in% pathway]
    
    dat <- data.frame("In Pathway" = c(length(sig_in), length(back_in)),
                      "Not in Pathway" = c(length(sig_out), length(back_out)),
                      row.names = c("Significant Genes", "Background Genes (Any successfully tested junction found across brain tissues)"),
                      stringsAsFactors = FALSE
    )
    
    test_f <- fisher.test(dat, alternative = "greater")
    test_c <- chisq.test(dat)
    
    
    fisher_res<-tibble(fisher_statistic=test_f$estimate, fisher_p_value=test_f$p.value, fisher_fdr=test_f$fdr,
                       chisq_statistic=test_c$statistic, chisq_p_value=test_c$p.value, fisher_chisq=test_c$fdr,
                       sig_in_list=length(sig_in), sig_not_in_list=length(sig_out), background_in_list=length(back_in), 
                       background_not_in_list=length(back_out),
                       pathway=i)
    
    fisher_res_df<-rbind(fisher_res_df, fisher_res)
    
    # print(paste(i, test$p.value))
    
  } # end of tdp sig dataset loop
  
  fisher_res_df$fisher_fdr<-p.adjust(p = fisher_res_df$fisher_p_value,
                                     method = "fdr")
  fisher_res_df$chisq_fdr<-p.adjust(p = fisher_res_df$chisq_p_value, 
                                    method = "fdr")
  
  fisher_res_df$Dataset<-dataset

  fisher_res_df_merge<-rbind(fisher_res_df_merge, fisher_res_df)
  
}

# arrange factor order
fisher_res_df_merge %<>% 
  mutate(Dataset=case_when(Dataset=='LRRK2 p.G2019S\nPost-mortem\nbrain' ~ 'LRRK2 p.G2019S\nPost-mortem brain',
         .default=Dataset)) %>%
  mutate(Dataset=factor(Dataset,
                      levels=c("LRRK2 p.G2019S\nPost-mortem brain",
                               'LRRK2 p.G2019S\nmDNs'))) %>% 
  arrange(desc(fisher_fdr))
  
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Graphing of  pathways ----
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Graph of nominal significant pathways - with fdr sig mark


PRS_plot<-
  fisher_res_df_merge %>% filter(fisher_p_value<0.05) %>%
  mutate(is_fdr_sig=case_when(
    fisher_fdr<0.05 ~ TRUE,
    fisher_fdr>0.05 ~ FALSE
  )) %>% 
  ungroup() %>%
  # mutate(pathway=str_to_sentence(str_trim(pathway))) %>% # change off all caps and remove blank space at end
  mutate(pathway=str_trim(pathway)) %>%
  mutate(pathway=factor(str_wrap(pathway, width=20))) %>%
  # mutate(pathway=fct_reorder(pathway, -log10(fisher_p_value))) %>%
mutate(pathway = fct_reorder2(pathway, Dataset, -log10(fisher_p_value), .desc = T)) %>%
  # mutate(pathway=reorder_within(pathway,by=log10(fisher_p_value), within=Dataset)) %>%
  # mutate(pathway = fct_reorder(pathway, -log10(fisher_p_value),.desc = TRUE)) %>%
  ggplot(aes(-log10(fisher_p_value), pathway, fill=Dataset, group=Dataset
  ))+
  # facet_wrap(~Dataset, drop = T, scales="free_x", nrow=2)+
  facet_wrap(~Dataset, ncol = 1, strip.position = "top", drop = TRUE) +
  theme(strip.placement = "outside")+
  geom_col(position = position_dodge2(), alpha=0.9)+
  geom_text(aes(label = ifelse(is_fdr_sig, "*", "")),
            size = 7, 
            # nudge_y = c(rep(0,3),0.2, rep(0,4)),
            nudge_x = 0.1
  ) +
  scale_fill_manual(values = rep(stage_colours,2))+
  # scale_alpha_manual(values = c(0.6, 0.9))+
  labs(y="", x="-log"[1][0] ~ "(p)", fill="")+
  guides(fill="none")+
  theme_bw(base_size = 14, base_family = "Roboto")+
  # scale_y_reordered()+
  theme(text = element_text(size=15),
        legend.position = "top",
        strip.text = element_text(size=18),
        axis.text.x = element_text(size=13, angle=60, vjust = 1, hjust=1)
  )+
  coord_flip()+
  theme(
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0)
  )



ggsave(
  plot =  PRS_plot,
  file.path(
    figure_out_path,
    "Figure_3d.png"
  ),
  # width = 8500,
  # height = 4000,
  height = 4000,
  width = 3500,
  dpi = 300,
  units = 'px'
)




