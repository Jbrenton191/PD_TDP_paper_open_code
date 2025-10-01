# 0 Setup ----
library(magrittr)
library(patchwork)
library(tidyverse)
library(foreach)
library(doParallel)
library(ggh4x)

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

#=====================================================================================================================================================================

# 1. TDP KD junctions and PD splicing analsyes loading ----

#=====================================================================================================================================================================

# tdp 43 targets - exon to exon junction format
tdp_junc_list<-as_tibble(readRDS(file = file.path(data_path, "TDP-43_HS_KD/Significant_TDP_KD_Database_ex_ex_annotated_df.RDS")))
# Keep only distinct singificant junctions across experiments
tdp_junc_list %<>% filter(p.adjust<0.05) %>%
  mutate(start=as.character(start), end=as.character(end)) %>%
  distinct(start, end,chr, strand,.keep_all = T) 


# Load Braak stage 6, Team Hardy splicing results
hardy_full_results<-readRDS(file = file.path(data_path, "Late_stage/full_results_list_for_default_params.RDS"))

# Load Braak stage 3/4, Team Wood splicing results
wood_full_results<-readRDS(file = file.path(data_path, "Mid_stage/full_results_list_for_default_params.RDS"))


#=====================================================================================================================================================================

# 2. Run junction enrichment analyses: ----

#=====================================================================================================================================================================

## Extract splicing junctions and run enrichment function: ----

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
  
  # adjust to get exon start, leafcutter adds 1 to the end but not doesn't adjust start as expecting regtools input (which has done this) and not SJ.outs to be used.
  all_results$start<-as.character(as.numeric(all_results$start)-1)
  
  # get all the junctions tested across brain areas
  universe<-all_results
  full_junc_universe_together<-all_results
  full_junc_universe_together %<>% distinct(start, end, chr, .keep_all = T)
  
  asap_sig_leaf <- all_results %>% filter(p.adjust<0.05)
  
  
  # fisher tests ----
  ## Collapsed result
  fisher_res_df_merge<-c()
  fisher_res_df<-c()
  # Running over each comparison

  for (comp in unique(asap_sig_leaf$comparison)) {
    
    # number of TDP43 matches in sig leafcutter list
    sig_in<-asap_sig_leaf %>% 
      filter(comparison==comp) %>%
      distinct(start, end, chr, strand, .keep_all = T) %>% 
      inner_join(tdp_junc_list %>% distinct(start, end, chr, strand,.keep_all = T), by=c('start', 'end', "chr", "strand"))
    
    # number of TDP43 matches NOT in sig leafcutter list
    sig_out<-asap_sig_leaf %>% 
      filter(comparison==comp) %>%
      distinct(start, end, chr,strand, .keep_all = T) %>% 
      anti_join(tdp_junc_list %>% distinct(start, end, chr, strand,.keep_all = T), by=c('start', 'end', "chr", "strand")) %>% arrange(genes)
    
    # number of TDP43 matches in leafcutter universe list
    back_in<-full_junc_universe_together %>% 
      distinct(start, end, chr,strand, .keep_all = T) %>%
      # removed filter to keep list larger - rather than brain area specific list - can revert back - but this is same as GO
      inner_join(tdp_junc_list %>%  distinct(start, end, chr, strand,.keep_all = T), by=c('start', 'end', "chr", "strand"))
    
    # number of TDP43 matches NOT in leafcutter universe list
    back_out<-full_junc_universe_together %>%
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
    
    
    fisher_res<-tibble(
      fisher_statistic=test$estimate,
                       fisher_conf_int=list(test$conf.int),
                       fisher_p_value=test$p.value, 
                       comparison=comp,
                       sig_in=dim(sig_in)[1], 
                       sig_out=dim(sig_out)[1],
                       back_in=dim(back_in)[1],
                       back_out=dim(back_out)[1]
                       )
    
    
    tdp_junction_ratio=(dim(sig_in)[1]/(dim(sig_in)[1]+dim(sig_out)[1]))
    
    fisher_res$tdp_junction_ratio<-tdp_junction_ratio
    
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

# Adjust comparison to read brain areas only in cluster df
hardy_full_results$all_results$cluster_significance %<>% 
  mutate(comparison=str_remove(comparison,"\\..*"))

# Adjust comparison to read brain areas only in intron df
hardy_full_results$all_results$intron_usage %<>% 
  mutate(comparison=str_remove(comparison,"\\..*"))

# Run enrichment function
hardy_fisher_res<-fisher_test_junctions(results = hardy_full_results$all_results, 
                                        dataset = "Late-stage", 
                                        tdp_junc_list =  tdp_junc_list)

# Adjust comparison to read brain areas only in cluster df
wood_full_results$all_results$cluster_significance %<>% 
  mutate(comparison=str_remove(comparison, "_Control.*"))

# Adjust comparison to read brain areas only in intron df
wood_full_results$all_results$intron_usage %<>% 
  mutate(comparison=str_remove(comparison, "_Control.*"))

# Run enrichment function
wood_fisher_res<-fisher_test_junctions(results = wood_full_results$all_results, 
                                       dataset = "Mid-stage", 
                                       tdp_junc_list =  tdp_junc_list)

# Merge results together
merge_df<-rbind(hardy_fisher_res, wood_fisher_res)

# Add brain area column and remove other text in comparisons columns 
merge_df %<>% mutate(dataset_ba=str_c(leafcutter_dataset, comparison, sep = ": ")) 


# add significance T/F column
merge_df %<>% mutate(
  is_sig = case_when(fdr <= 0.05 ~ TRUE,
                     fdr > 0.05 ~ FALSE)
)

# rename comparisons/brain areas to full names
merge_df <-
  merge_df %>% mutate(
    comparison = case_when(
      leafcutter_dataset == "Mid-stage" &
        comparison == "C_CTX" ~ "Ant.\nCingulate",
      leafcutter_dataset == "Mid-stage" &
        comparison == "F_CTX" ~ "Frontal",
      leafcutter_dataset == "Mid-stage" &
        comparison == "P_CTX" ~ "Parietal",
      leafcutter_dataset == "Mid-stage" &
        comparison == "T_CTX" ~ "Temporal",
      leafcutter_dataset == "Mid-stage" &
        comparison == "CAU" ~ "Caudate",
      leafcutter_dataset == "Mid-stage" &
        comparison == "PUT" ~ "Putamen",
      leafcutter_dataset == "Mid-stage" &
        comparison == "SN" ~ "Substantia nigra",
      leafcutter_dataset == "Mid-stage" &
        comparison == "PARA" ~ "Parahippocampal",
      leafcutter_dataset == "Late-stage" &
        comparison == "ACG" ~ "Ant.\nCingulate",
      leafcutter_dataset == "Late-stage" & comparison ==
        "MTG" ~ "Temporal",
      leafcutter_dataset == "Late-stage" & comparison ==
        "MFG" ~ "Frontal",
      leafcutter_dataset == "Late-stage" & comparison ==
        "IPL" ~ "Parietal",
      TRUE ~ "DOUBLE CHECK"
    )
  )

# Factor based off approx braak stage order
merge_df$comparison <-
  factor(
    merge_df$comparison,
    levels = c("Substantia nigra", "Caudate", "Putamen", "Ant.\nCingulate",
               "Parahippocampal", 'Temporal',
               "Frontal", "Parietal"
               )
    )
  

# factor order of braak stages
merge_df$leafcutter_dataset<-factor(merge_df$leafcutter_dataset, 
                                    levels = c("Mid-stage", "Late-stage"))

# add colours to each significance x braak stage
merge_df %<>% mutate(stage_colour=case_when(
  leafcutter_dataset == "Mid-stage" & fisher_p_value < 0.05 ~ "#91D1C2",
  leafcutter_dataset == "Mid-stage" & fisher_p_value > 0.05 ~ "grey",
  leafcutter_dataset == "Late-stage" & fisher_p_value < 0.05 ~ "#4B8D86",
  leafcutter_dataset == "Late-stage" & fisher_p_value > 0.05 ~ "grey"
)
                    )



#=====================================================================================================================================================================

# 3. Gene overlap ----

#=====================================================================================================================================================================



# extract genes from results and run enrichment function

fisher_test_genes<-function(results, dataset, 
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
  
  # filter for all genes tested and extract vector of gene names
  universe<-all_results %>% pull(genes)
  
  # filter for all significant genes
  asap_sig_leaf <- all_results %>% filter(p.adjust<0.05)
  
  
  # fisher tests ----
  fisher_res_df_merge<-c()
  fisher_res_df<-c()
  
 
  fisher_res_df<-c()
  
  for (comp in unique(asap_sig_leaf$comparison)) {
    
    asap_sig_leaf <- all_results %>% filter(p.adjust<0.05) %>% filter(comparison==comp)
    
    # number of TDP43 matches in sig leafcutter list
    sig_in<- unique(asap_sig_leaf$genes)[unique(asap_sig_leaf$genes) %in% unique(tdp_junc_list$genes)]
    
    # number of TDP43 matches NOT in sig leafcutter list
    sig_out<-unique(asap_sig_leaf$genes)[!unique(asap_sig_leaf$genes) %in% unique(tdp_junc_list$genes)]
    
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
                       fisher_p_value=test$p.value, comparison=comp,
                       sig_in=length(sig_in), 
                       sig_out=length(sig_out),
                       back_in=length(back_in),
                       back_out=length(back_out)
                       )
    
    tdp_gene_ratio=length(sig_in)/(length(sig_in)+length(sig_out))
    
    fisher_res$tdp_gene_ratio<-tdp_gene_ratio 
    
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

hardy_fisher_genes<-fisher_test_genes(results = hardy_full_results$all_results, 
                                      dataset = "Late-stage",
                                      tdp_junc_list =  tdp_junc_list)

wood_fisher_genes<-fisher_test_genes(results = wood_full_results$all_results, 
                                     dataset = "Mid-stage", 
                                     tdp_junc_list =  tdp_junc_list)

# Merge results together
merge_df_genes<-rbind(hardy_fisher_genes, wood_fisher_genes)

# Add brain area column and remove other text in comparisons columns 
merge_df_genes %<>% mutate(dataset_ba=str_c(leafcutter_dataset, comparison, sep = ": ")) 


# add significance T/F column
merge_df_genes %<>% mutate(
  is_sig = case_when(fdr <= 0.05 ~ TRUE,
                     fdr > 0.05 ~ FALSE)
)

# rename comparisons/brain areas to full names
merge_df_genes <-
  merge_df_genes %>% mutate(
    comparison = case_when(
      leafcutter_dataset == "Mid-stage" &
        comparison == "C_CTX" ~ "Ant.\nCingulate",
      leafcutter_dataset == "Mid-stage" &
        comparison == "F_CTX" ~ "Frontal",
      leafcutter_dataset == "Mid-stage" &
        comparison == "P_CTX" ~ "Parietal",
      leafcutter_dataset == "Mid-stage" &
        comparison == "T_CTX" ~ "Temporal",
      leafcutter_dataset == "Mid-stage" &
        comparison == "CAU" ~ "Caudate",
      leafcutter_dataset == "Mid-stage" &
        comparison == "PUT" ~ "Putamen",
      leafcutter_dataset == "Mid-stage" &
        comparison == "SN" ~ "Substantia nigra",
      leafcutter_dataset == "Mid-stage" &
        comparison == "PARA" ~ "Parahippocampal",
      leafcutter_dataset == "Late-stage" &
        comparison == "ACG" ~ "Ant.\nCingulate",
      leafcutter_dataset == "Late-stage" & comparison ==
        "MTG" ~ "Temporal",
      leafcutter_dataset == "Late-stage" & comparison ==
        "MFG" ~ "Frontal",
      leafcutter_dataset == "Late-stage" & comparison ==
        "IPL" ~ "Parietal",
      TRUE ~ "DOUBLE CHECK"
    )
  )

merge_df_genes$comparison <-
  factor(
    merge_df_genes$comparison,
    levels = c("Substantia nigra", "Caudate", "Putamen", "Ant.\nCingulate",
               "Parahippocampal", 'Temporal',
               "Frontal", "Parietal"
               )
  )
merge_df_genes$leafcutter_dataset<-factor(merge_df_genes$leafcutter_dataset, 
                                    levels = c("Mid-stage", "Late-stage"))

merge_df_genes %<>% mutate(stage_colour=case_when(
  leafcutter_dataset == "Mid-stage" & fisher_p_value < 0.05 ~ "#91D1C2",
  leafcutter_dataset == "Mid-stage" & fisher_p_value > 0.05 ~ "grey",
  leafcutter_dataset == "Late-stage" & fisher_p_value < 0.05 ~ "#4B8D86",
  leafcutter_dataset == "Late-stage" & fisher_p_value > 0.05 ~ "grey"
)
)


  
#=====================================================================================================================================================================

# 3. Graphing junctions and genes together: ----

#=====================================================================================================================================================================

# Junction graph

p3<-ggplot(data = merge_df %>% mutate(Analysis="Junction"),
       aes(x = comparison , y = -log10(fdr),
           colour=stage_colour, 
           alpha=is_sig,
       ))+
  geom_point( 
    size=7
    )+
  ggh4x::facet_grid2(Analysis~leafcutter_dataset,
                     scales = "free_x", 
                     strip = strip_themed(background_y = elem_list_rect(fill = "white"),
                                          background_x =
                                            elem_list_rect(fill = stage_colours)))+
  geom_hline(yintercept = -log10(0.05), linetype='dashed', colour='black')+
  labs(x="", y=expression("-log"[10]*"(FDR)"), size="Enrichment\nRatio (%)")+
  scale_alpha_manual(values =  merge_df$sig_alpha)+
  scale_color_identity()+
  scale_alpha_manual(values = c(0.8))+ # one size needed as all comparisons are significant
  scale_size_continuous(limits = c(0,100), range=c(1,10),breaks = c(25,50,75,100), transform = "identity")+
  guides(
    size="none",
    color = "none",   # Remove legend for color
    shape = "none",
    alpha="none"
  )+
  theme_jb()+
  theme_bw(base_family = "Roboto")+
  theme(
    text = element_text(size = 14),
    legend.text = element_text(size=8),
    legend.title = element_text(size=10),
    legend.key.size = unit(0.1, 'cm'),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 14),
    legend.position = "right",
    axis.text.x = element_text(
      angle = 60,
      hjust = 1.05,
      vjust = 1.05,
      size = 14
    ),
    plot.margin=unit(c(0,0,1.5,0), 'cm'),
    strip.text.x =  element_text(size = 16),
    strip.text.y =  element_text(size = 14),
    panel.spacing = unit(2, "lines")
  )


# gene graph
p4<-ggplot(data = merge_df_genes %>% mutate(Analysis="Gene"),
       aes(x = comparison , y = -log10(fdr),
           colour=stage_colour, 
           alpha=is_sig
       ))+
  geom_point(
    size=7
    )+
  ggh4x::facet_grid2(Analysis~leafcutter_dataset,
                     scales = "free_x", 
                     strip = strip_themed(background_y = elem_list_rect(fill = "white"),
                                          background_x =
                                            elem_list_rect(fill = stage_colours)))+
  geom_hline(yintercept = -log10(0.05), linetype='dashed', colour='black')+
  labs(x="", y=expression("-log"[10]*"(FDR)"), size="Enrichment\nRatio (%)")+
  scale_color_identity()+
  scale_alpha_manual(values = c(0.8))+ # one size needed as all comparisons are significant
  scale_size_continuous(limits = c(0,100), range=c(1,10),
                        breaks = c(25,50,75,100), transform = "identity")+
  guides(
    size="none",
    color = "none",   # Remove legend for color
    shape = "none",
    alpha="none"
  )+
  theme_jb()+
  theme_bw(base_family = "Roboto")+
  theme(
    text = element_text(size = 14),
    legend.text = element_text(size=8),
    legend.title = element_text(size=10),
    legend.key.size = unit(0.1, 'cm'),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 14),
    legend.position = "right",
    axis.text.x = element_text(
      angle = 60,
      hjust = 1.05,
      vjust = 1.05,
      size = 14
    ),
    strip.background = element_blank(),
    strip.text.y =  element_text(size = 16),
    strip.text.x =  element_blank(),
    # ,
    panel.spacing = unit(2, "lines")
  )


ggsave(
  plot =  
    (p3/p4)+plot_layout(axes = 'collect'),
  filename = file.path(figure_out_path, "Figure_2b.png"),
  width = 6500,
  height = 5250,
  dpi = 535,
  units = 'px'
)



