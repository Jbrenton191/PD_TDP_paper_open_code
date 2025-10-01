
library(magrittr)
library(patchwork)
library(tidyverse)
library(foreach)
library(doParallel)
library(ggh4x)

# 0 Setup ----

main_path<-here::here() # main project level i.e TDP43_paper_open_code level

data_path<-file.path(main_path, "Data")

figure_out_path<-file.path(main_path, "Plots/Main_Figures/Figure_2")

# create figure output directory if not present
if (!dir.exists(figure_out_path)) {
  dir.create(figure_out_path, recursive = T)
}

analysis_colours<-c("#0c5394", "#e06666")

stage_colours<-c("#91D1C2", "#4B8D86")

Braak_stage_order<-c("Mid-stage", "Late-stage")

mid_brain_area_order<-c("Substantia nigra", "Caudate", "Putamen", "Ant.\nCingulate",
                        "Parahippocampal", 'Temporal',
                        "Frontal", "Parietal", "Combined")

late_brain_area_order<-c("Ant.\nCingulate", 'Temporal ',
                         "Frontal", "Parietal", "Combined")

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


## Checking PDPRS pathways in GO data overlap with sporadic splicing

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Load pathways downloaded from https://pdgenetics.shinyapps.io/pathwaysbrowser/ , and adjusted into df format
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

PDPRSpathways_ensg <- readRDS(file.path(data_path, "Gene_lists/PDPRSpathways_ensg.rds"))
# Convert to long tibble format
PDPRSpathways_ensg_tibble<-enframe(PDPRSpathways_ensg, name = "Pathway", value = "Gene") %>% unnest(Gene) %>% mutate(Gene_no_version=gsub(x=Gene, pattern="\\..*", replacement=""))

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Load sig data splicng data and get backgrounds ----
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1. Get Mid-stage (Team Wood splicing results) DS genes and their overlap with TPD KD DS genes
# Load Braak stage 3/4, Team Wood splicing results
wood_full_results<-readRDS(file = file.path(data_path, "Mid_stage/full_results_list_for_default_params.RDS"))


results<-wood_full_results$all_results
# cluster df
lc_all<-results$cluster_significance
# intron df
int_all<-results$intron_usage

# conbine cluster and intron df
mid_results <- lc_all %>% 
  # filter for successful tests and remove clusters overlapping multiple genes
  dplyr::filter(!str_detect(genes, ",")) %>%
  tidyr::separate(cluster, into = c("chr", "cluster"), sep = ":") %>% 
  dplyr::inner_join(int_all %>% 
                      tidyr::separate(intron, into = c("chr", "start", "end", "cluster"), sep = ":")) %>% 
  dplyr::mutate(direction_of_effect = case_when(deltapsi > 0 ~ "upregulated",
                                                deltapsi < 0 ~ "downregulated",
                                                deltapsi == 0 ~ "no_change"))

# create universe of all genes tested across brain areas
mid_universe<-mid_results %>% distinct(genes) %>% pull(genes)
# keep only significant clusters/genes 
asap_sig_leaf <- mid_results %>% filter(p.adjust<0.05)

# tdp 43 targets - exon to exon junction format
tdp_junc_list<-as_tibble(readRDS(file = file.path(data_path, "TDP-43_HS_KD/Significant_TDP_KD_Database_ex_ex_annotated_df.RDS")))
# Keep only distinct singificant junctions across experiments
tdp_junc_list %<>% filter(p.adjust<0.05) %>%
  mutate(start=as.character(start), end=as.character(end)) %>%
  distinct(start, end,chr, strand,.keep_all = T)

# keep genes with significant splicing alterations in mid-stage PD (across any brain area) and those that show significant splicing alterations upon TDP-43 knockdown
mid_sig_overlap<- unique(asap_sig_leaf$genes[asap_sig_leaf$genes %in% unique(tdp_junc_list$genes)])

# Late-stage Overlap ----

# Load Braak stage 6, Team Hardy splicing results
hardy_full_results<-readRDS(file = file.path(data_path, "Late_stage/full_results_list_for_default_params.RDS"))

results<-hardy_full_results$all_results

lc_all<-results$cluster_significance
int_all<-results$intron_usage

late_results <- lc_all %>% 
  # filter for successful tests and remove clusters overlapping multiple genes
  dplyr::filter(!str_detect(genes, ",")) %>%
  tidyr::separate(cluster, into = c("chr", "cluster"), sep = ":") %>% 
  dplyr::inner_join(int_all %>% 
                      tidyr::separate(intron, into = c("chr", "start", "end", "cluster"), sep = ":")) %>% 
  dplyr::mutate(direction_of_effect = case_when(deltapsi > 0 ~ "upregulated",
                                                deltapsi < 0 ~ "downregulated",
                                                deltapsi == 0 ~ "no_change"))
# create universe of all genes tested across brain areas
late_universe<-late_results %>% distinct(genes) %>% pull(genes)
# keep only significant clusters/genes
asap_sig_leaf <- late_results %>% filter(p.adjust<0.05)

# keep genes with significant splicing alterations in late-stage PD (across any brain area) and those that show significant splicing alterations upon TDP-43 knockdown
late_sig_overlap<- unique(asap_sig_leaf$genes[asap_sig_leaf$genes %in% unique(tdp_junc_list$genes)])

# Merge the two datasets
mid_sig<-tibble(SYMBOL=mid_sig_overlap, Dataset="Mid-stage")
late_sig<-tibble(SYMBOL=late_sig_overlap, Dataset="Late-stage")
merge_df_genes<-rbind(mid_sig, late_sig)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
# Load gencode gene map
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
gene_map_path <- file.path(data_path, "References/gencode_txid_to_geneid.txt")

gencode_txid_to_geneid <- vroom::vroom(gene_map_path) %>%
  `colnames<-`(c("tx_id", "gene_id","gene_name", "description")) %>%
  dplyr::mutate(tx_id = sub("\\..+", "", tx_id),
                gene_id = sub("\\..+", "", gene_id)) %>% 
  dplyr::select(gene_id, gene_name) %>% distinct(gene_id, .keep_all=T)

# Missing genes between the two:
PDPRSpathways_ensg_tibble[!PDPRSpathways_ensg_tibble$Gene_no_version %in% gencode_txid_to_geneid$gene_id,]
# works out as 15 distinct genes across 15 pathways - probably because using MSigDB database v7.2?

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# map gene names to Gene in PDPRS pathways
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

PDPRSpathways_ensg_tibble %<>% inner_join(., gencode_txid_to_geneid, by=c("Gene_no_version"= "gene_id"))

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# run enrichment analysis
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Mid-stage overlap ----
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

asap_sig_leaf <- merge_df_genes %>% filter(Dataset=="Mid-stage")

fisher_res_df<-c()

fisher_res_df_merge<-c()

for (i in unique(PDPRSpathways_ensg_tibble$Pathway)) {
  
  pathway<-PDPRSpathways_ensg_tibble %>% filter(Pathway==i) %>% pull(gene_name)
  
  asap_sig_leaf_filt <- asap_sig_leaf %>% distinct(SYMBOL,.keep_all = T) %>% pull(SYMBOL)
  
  # number of pathway matches in sig leafcutter list
  sig_in<-asap_sig_leaf_filt[asap_sig_leaf_filt %in% pathway]
  
  # number of pathway matches NOT in sig leafcutter list
  sig_out<-asap_sig_leaf_filt[!asap_sig_leaf_filt %in% pathway]
  
  # number of matches in leafcutter universe list
  back_in<-mid_universe[mid_universe %in% pathway]
  
  # number of matches NOT in leafcutter universe list
  back_out<-mid_universe[!mid_universe %in% pathway]
  
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

fisher_res_df$fisher_fdr<-p.adjust(p = fisher_res_df$fisher_p_value, method = "fdr")
fisher_res_df$chisq_fdr<-p.adjust(p = fisher_res_df$chisq_p_value, method = "fdr")

fisher_res_df$Dataset<-"Mid-stage"

fisher_res_df_merge<-rbind(fisher_res_df_merge,fisher_res_df)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Late-stage enrichment ----
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

asap_sig_leaf <- merge_df_genes %>% filter(Dataset=="Late-stage")

fisher_res_df<-c()

for (i in unique(PDPRSpathways_ensg_tibble$Pathway)) {
  
  pathway<-PDPRSpathways_ensg_tibble %>% filter(Pathway==i) %>% pull(gene_name)
  
  asap_sig_leaf_filt <- asap_sig_leaf %>% distinct(SYMBOL,.keep_all = T) %>% pull(SYMBOL)
  
  # number of pathway matches in sig leafcutter list
  sig_in<-asap_sig_leaf_filt[asap_sig_leaf_filt %in% pathway]
  
  # number of pathway matches NOT in sig leafcutter list
  sig_out<-asap_sig_leaf_filt[!asap_sig_leaf_filt %in% pathway]
  
  # number of matches in leafcutter universe list
  back_in<-late_universe[late_universe %in% pathway]
  
  # number of matches NOT in leafcutter universe list
  back_out<-late_universe[!late_universe %in% pathway]
  
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
  
} # end of tdp sig dataset loop

fisher_res_df$fisher_fdr<-p.adjust(p = fisher_res_df$fisher_p_value, method = "fdr")
fisher_res_df$chisq_fdr<-p.adjust(p = fisher_res_df$chisq_p_value, method = "fdr")

fisher_res_df$Dataset<-"Late-stage"

fisher_res_df_merge<-rbind(fisher_res_df_merge,fisher_res_df) %>% mutate(Dataset=factor(Dataset, levels=Braak_stage_order))


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Graph of nominal significant pathways with FDR significance highlighted by a *

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


stage_colours<-c("#91D1C2", "#4B8D86")

PRS_plot<-fisher_res_df_merge %>% filter(fisher_p_value<0.05) %>%
  mutate(is_fdr_sig=case_when(
   fisher_fdr<0.05 ~ TRUE,
   fisher_fdr>0.05 ~ FALSE
  )) %>% 
  
  # mutate(pathway=str_to_sentence(str_trim(pathway))) %>% # change off all caps and remove blank space at end
  mutate(pathway=(str_trim(pathway))) %>% # remove blank space at end
  mutate(pathway=factor(str_wrap(pathway, width=20))) %>%
  mutate(pathway=fct_reorder(pathway, fisher_p_value)) %>%
  ggplot(aes(-log10(fisher_p_value), pathway, fill=Dataset, group=Dataset
             ))+
  geom_col(position = position_dodge2(), alpha=0.9)+
  # geom_col(position =  position_dodge2(preserve = "single"), alpha=0.9)+
  geom_text(aes(label = ifelse(is_fdr_sig, "*", "")),
            size = 7, 
            nudge_y = c(rep(0,3),-0.22, rep(0,4)),
            nudge_x = 0.1
            ) +
  scale_fill_manual(values = stage_colours)+
  # scale_alpha_manual(values = c(0.6, 0.9))+
  labs(y="", x="-log"[1][0] ~ "(p)", fill="")+
  theme_bw(base_size = 12, base_family = "Roboto")+
  theme(legend.position = "top",
        axis.text.x = element_text(size=10, angle=60, vjust = 1.08, hjust=1.1)
        )+
  coord_flip()


ggsave(
  plot =  PRS_plot,
  file.path(
    figure_out_path,
    "Figure_2d.png"
  ),
  width = 2000,
  height = 3000,
  dpi = 350,
  units = 'px'
)




