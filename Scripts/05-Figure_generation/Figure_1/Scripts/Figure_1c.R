# 0. Running enrichment for ALS and PD open target risk genes ----
# 1. Setup ----
library(tidyverse)
library(magrittr)
library(doParallel)
library(UpSetR)
library(patchwork)
library(ggnewscale)

main_path<-here::here() # main project level i.e TDP43_paper_open_code level


data_path<-file.path(main_path, "Data")

figure_out_path<-file.path(main_path, "Plots/Main_Figures/Figure_1")

# create figure output directory if not present
if (!dir.exists(figure_out_path)) {
  dir.create(figure_out_path, recursive = T)
}


analysis_colours<-c("#0c5394", "#e06666")

stage_colours<-c("#91D1C2", "#4B8D86")

Braak_stage_order<-c("Mid-stage", "Late-stage")

mid_brain_area_order<-c("Substantia\nnigra", "Caudate", "Putamen", "Ant.\nCingulate",
                        "Parahippocampal", 'Temporal',
                        "Frontal", "Parietal", "Combined")

late_brain_area_order<-c("Ant.\nCingulate", 'Temporal',
"Frontal", "Parietal", "Combined")

extrafont::font_import(paths = data_path, 
                       pattern = 'Roboto', prompt = F)

# New theme based off Regina theme:

theme_jb <- function(text_size=12) {
  theme_bw(base_family = "Roboto", base_size = text_size) +
    theme(
      legend.position = "top",
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.title.y = element_text(vjust = 0.6),
      panel.spacing = unit(0.1, "lines")
    )
}


# set for all comparisons - looking for an association across Opentargets risk and mendelian genes
genetic_association_filter<-0.1

# Load Wood DEG data

dge_all <- readRDS(file.path(data_path, "Mid_stage/dge_all.rds"))

dge_all<-dge_all[-grep(x = dge_all$tissue, pattern = "Collapse", ignore.case = 'T'),]
  
# Load gene_map *************

# preprocessed GTF for specific columns: transcript id, gene id, gene name and transcript description
gene_map_path <- file.path(data_path, "References/gencode_txid_to_geneid.txt")

gencode_txid_to_geneid <- vroom::vroom(gene_map_path) %>%
  `colnames<-`(c("tx_id", "gene_id","gene_name", "description")) %>%
  dplyr::mutate(tx_id = sub("\\..+", "", tx_id),
                gene_id = sub("\\..+", "", gene_id)) %>% 
  dplyr::select(gene_id, gene_name) %>% distinct(gene_id, .keep_all=T)

# **********************
# 2. Run Braak 3/4 dataset/Team Wood data
dge_all %<>% inner_join(gencode_txid_to_geneid) 

dge_sig<-dge_all %>% filter(adj.P.Val<0.05)

dge_background<-dge_all %>% distinct(gene_name) %>% pull()

# Load Braak stage 3/4,Wood splicing data
full_res<-readRDS(file = file.path(data_path, "Mid_stage/full_results_list_for_default_params.RDS"))

merged_summaries_gene_clust_lists<-full_res$all_results$significant_clusters_0.05_filter %>% 
  mutate(comparison=str_remove(comparison, "_Control.*"))

sig_splice<-merged_summaries_gene_clust_lists
# use all distinct genes found in any Brain area as background 
background<-unique(full_res$universe_list_of_genes_tested$genes)

## 2.1 Run ALS risk gene lst ----

# Get ALS merged list

ALS_gene_table <-
  vroom::vroom( 
    file.path(data_path, "Gene_lists/OT-MONDO_0004976-associated-targets-30_04_2025-v25_03.tsv")) %>% 
  dplyr::select(1, 3:11) %>% 
    pivot_longer(cols = -c(1),names_to = "Target_analysis", values_to = "Score") %>% 
  filter(Score!="No data", Score>genetic_association_filter) %>%
distinct(symbol) %>% pull()


# Run Enrichment

# all must be distinct vectors of gene symbols or IDs rather dataframes - as long as sig list, background and test_list in same 

gene_enrichment<-function(sig_genes, background, test_gene_list){
  
  sig_in<-length(sig_genes[sig_genes %in% test_gene_list])
  sig_out<-length(sig_genes[!sig_genes %in% test_gene_list])
  back_in<-length(background[background %in% test_gene_list])
  back_out<- length(background[!background %in% test_gene_list])
  
  dat <- data.frame("In Gene Set" = c(sig_in, sig_out),
                    "Not In Gene Set" = c(back_in, back_out),
                    row.names = c("Significantly Spliced Genes", "Backround Genes"),
                    stringsAsFactors = FALSE
  )
  
  test <- fisher.test(dat, alternative = "greater")
  
  
  res<-tibble(p_value=test$p.value, odds_ratio=test$estimate, sig_ratio=sig_in/(sig_in+sig_out), bg_ratio=back_in/(back_in+back_out))
  
  return(res)
  
}

# Run on all DEG data - only brain areas in Wood data

res_table_degs<-c()

for (i in unique(dge_all$tissue)) {
  
  degs<-dge_sig %>% filter(tissue==i) %>% distinct(gene_name) %>% pull()
  
  enrich_out<-gene_enrichment(degs, dge_background, ALS_gene_table)
  
  enrich_out$comparison<-i
  
  res_table_degs<-rbind(res_table_degs, enrich_out)
  
}


# Merge
ALS_risk_degs<-res_table_degs

ALS_risk_degs$Analysis<-"Differential expression"

ALS_risk_degs$Braak_stage<-"Mid-stage"

ALS_risk_degs$target_list<-"ALS risk genes"


# Run on splicing data

res_table_splice<-c()

for (i in unique(sig_splice$comparison)) {
  
  sig_splice_filt<-sig_splice %>% filter(comparison==i) %>% distinct(genes) %>% pull()
  
  enrich_out<-gene_enrichment(sig_splice_filt, background, ALS_gene_table)
  
  enrich_out$comparison<-i
  
  res_table_splice<-rbind(res_table_splice, enrich_out)
  
}


ALS_risk_splice<-res_table_splice

ALS_risk_splice$Analysis<-"Differential splicing"

ALS_risk_splice$Braak_stage<-"Mid-stage"

ALS_risk_splice$target_list<-"ALS risk genes"


## 2.2. PD Risk gene list ----

PD_gene_table <-
  vroom::vroom(
    file.path(data_path, "Gene_lists/OT-MONDO_0005180-associated-targets-30_04_2025-v25_03.tsv")) %>%
    dplyr::select(1, 3:11) %>% 
      pivot_longer(cols = -c(1),names_to = "Target_analysis", values_to = "Score") %>% 
      filter(Score!="No data", Score>genetic_association_filter) %>%
      distinct(symbol) %>% pull()

res_table_degs<-c()

for (i in unique(dge_all$tissue)) {
  
  degs<-dge_sig %>% filter(tissue==i) %>% distinct(gene_name) %>% pull()
  
  enrich_out<-gene_enrichment(degs, dge_background, PD_gene_table)
  
  enrich_out$comparison<-i
  
  res_table_degs<-rbind(res_table_degs, enrich_out)
  
}

# merge

PD_risk_degs<-res_table_degs

PD_risk_degs$Analysis<-"Differential expression"

PD_risk_degs$Braak_stage<-"Mid-stage"

PD_risk_degs$target_list<-"Parkinson's risk genes"

# Run on splicing data

res_table_splice<-c()

for (i in unique(sig_splice$comparison)) {
  
  sig_splice_filt<-sig_splice %>% filter(comparison==i) %>% distinct(genes) %>% pull()
  
  enrich_out<-gene_enrichment(sig_splice_filt, background, PD_gene_table)
  
  enrich_out$comparison<-i
  
  res_table_splice<-rbind(res_table_splice, enrich_out)
  
}



PD_risk_splice<-res_table_splice

PD_risk_splice$Analysis<-"Differential splicing"

PD_risk_splice$Braak_stage<-"Mid-stage"

PD_risk_splice$target_list<-"Parkinson's risk genes"

# run collapse again of all regions
 
b34_merged_risk_lists<-rbind(ALS_risk_degs,ALS_risk_splice,
  PD_risk_degs,PD_risk_splice)


# 3. Run Braak stage 6 analysis -------
# *******************************************************************************************************************

dge_all<- read_csv(file.path(data_path, 
                             "Late_stage/dge_limma.csv"))

# Load gene_map *************

# **********************

dge_all %<>% inner_join(gencode_txid_to_geneid) 

dge_all<-dge_all[-grep(x = dge_all$tissue, pattern = "Collapse", ignore.case = 'T'),]

dge_sig<-dge_all %>% filter(adj.P.Val<0.05)

dge_background<-dge_all %>% distinct(gene_name) %>% pull()



# Load Hardy splicing data

full_res<-readRDS(file = file.path(data_path, "Late_stage/full_results_list_for_default_params.RDS"))

sig_splice<-full_res$all_results$significant_clusters_0.05_filter %>% mutate(comparison=str_remove(comparison,"\\..*"))
# use all distinct genes found in any Brain area as background - already unique but there if copied code over to new data
background<-unique(full_res$universe_list_of_genes_tested$genes)

## 3.1 Run ALS risk gene lst ----

# Get ALS merged list

# Run Enrichment

# all must be distinct vectors of gene symbols or IDs rather dataframes - as long as sig list, background and test_list in same 

gene_enrichment<-function(sig_genes, background, test_gene_list){
  
  sig_in<-length(sig_genes[sig_genes %in% test_gene_list])
  sig_out<-length(sig_genes[!sig_genes %in% test_gene_list])
  back_in<-length(background[background %in% test_gene_list])
  back_out<- length(background[!background %in% test_gene_list])
  
  dat <- data.frame("Significantly Spliced Genes" = c(sig_in, sig_out),
                    "Background Genes" = c(back_in, back_out),
                    row.names = c("In Gene Set","Not In Gene Set"),
                    stringsAsFactors = FALSE
  )
  
  test <- fisher.test(dat, alternative = "greater")
  
  
  res<-tibble(p_value=test$p.value, odds_ratio=test$estimate, sig_ratio=sig_in/sig_out, bg_ratio=back_in/back_out)
  
  return(res)
  
}

# Run on all DEG data - only brain areas in Wood data

res_table_degs<-c()

for (i in unique(dge_sig$tissue)) {
  
  degs<-dge_sig %>% filter(tissue==i) %>% distinct(gene_name) %>% pull()
  
  enrich_out<-gene_enrichment(degs, dge_background, ALS_gene_table)
  
  enrich_out$comparison<-i
  
  res_table_degs<-rbind(res_table_degs, enrich_out)
  
}

# Merge
ALS_risk_degs<-res_table_degs

ALS_risk_degs$Analysis<-"Differential expression"

ALS_risk_degs$Braak_stage<-"Late-stage"

ALS_risk_degs$target_list<-"ALS risk genes"


# Run on splicing data

res_table_splice<-c()

for (i in unique(sig_splice$comparison)) {
  
  sig_splice_filt<-sig_splice %>% filter(comparison==i) %>% distinct(genes) %>% pull()
  
  enrich_out<-gene_enrichment(sig_genes = sig_splice_filt, 
                              background = background, 
                              test_gene_list = ALS_gene_table)
  
  enrich_out$comparison<-i
  
  res_table_splice<-rbind(res_table_splice, enrich_out)
  
}


# correct p value

ALS_risk_splice<-res_table_splice

ALS_risk_splice$Analysis<-"Differential splicing"

ALS_risk_splice$Braak_stage<-"Late-stage"

ALS_risk_splice$target_list<-"ALS risk genes"


# run collapse of all distinct genes

## 3.2 PD Risk gene list ----

res_table_degs<-c()

for (i in unique(dge_all$tissue)) {
  
  degs<-dge_sig %>% filter(tissue==i) %>% distinct(gene_name) %>% pull()
  
  enrich_out<-gene_enrichment(degs, dge_background, PD_gene_table)
  
  enrich_out$comparison<-i
  
  res_table_degs<-rbind(res_table_degs, enrich_out)
  
}


PD_risk_degs<-res_table_degs

PD_risk_degs$Analysis<-"Differential expression"

PD_risk_degs$Braak_stage<-"Late-stage"

PD_risk_degs$target_list<-"Parkinson's risk genes"

# Run on splicing data

res_table_splice<-c()

for (i in unique(sig_splice$comparison)) {
  
  sig_splice_filt<-sig_splice %>% filter(comparison==i) %>% distinct(genes) %>% pull()
  
  enrich_out<-gene_enrichment(sig_splice_filt, background, PD_gene_table)
  
  enrich_out$comparison<-i
  
  res_table_splice<-rbind(res_table_splice, enrich_out)
  
}



PD_risk_splice<-res_table_splice

PD_risk_splice$Analysis<-"Differential splicing"

PD_risk_splice$Braak_stage<-"Late-stage"

PD_risk_splice$target_list<-"Parkinson's risk genes"

# run collapse again of all regions

b6_merged_risk_lists<-rbind(ALS_risk_degs,ALS_risk_splice,
                             PD_risk_degs,PD_risk_splice)

# Change Brain areas in each stage df

# Braak 3 + 4

b34_merged_risk_lists <-
  b34_merged_risk_lists %>% mutate(
    comparison = case_when(
      comparison == "C_CTX" ~ "Ant.\nCingulate",
      comparison == "F_CTX" ~ "Frontal",
      comparison == "P_CTX" ~ "Parietal",
      comparison == "T_CTX" ~ 'Temporal',
      comparison == "CAU" ~ "Caudate",
      comparison == "PUT" ~ "Putamen",
      comparison == "SN" ~ "Substantia\nnigra",
      comparison == "PARA" ~ "Parahippocampal",
      # comparison == "Collapse" ~ "Combined",
      TRUE ~ "DOUBLE CHECK"
    )
  )



b34_merged_risk_lists$comparison <-
  factor(
    b34_merged_risk_lists$comparison,
    levels = c(
      "Substantia\nnigra",
      "Caudate",
      "Putamen",
      "Ant.\nCingulate",
      "Parahippocampal",
      'Temporal',
      "Frontal",
      "Parietal"
         )
  )


# Braak 6

b6_merged_risk_lists <-
  b6_merged_risk_lists %>% mutate(
    comparison = case_when(
      comparison == "ACG" ~ "Ant.\nCingulate",
      comparison ==
        "MTG" ~ "Temporal",
      comparison ==
        "MFG" ~ "Frontal",
      comparison ==
        "IPL" ~ "Parietal"
     
    )
  )



b6_merged_risk_lists$comparison <-
  factor(
    b6_merged_risk_lists$comparison,
    levels = c(
      "Ant.\nCingulate",
      "Temporal",
      "Frontal",
      "Parietal"
    
    )
  )

# ****************************************************************************************************************************************************


# merge the two path stages

path_merged_risk_lists <- rbind(b34_merged_risk_lists,
                                b6_merged_risk_lists)

path_merged_risk_lists %<>% mutate(Braak_stage = factor(Braak_stage, levels =
                                                          Braak_stage_order))

path_merged_risk_lists %<>% mutate(
  target_list = case_when(
    target_list == "ALS risk genes" ~ "Amyotrophic lateral sclerosis",
    target_list == "Parkinson's risk genes" ~ "Parkinson's disease"
  )
)

# order analysis
path_merged_risk_lists$Analysis <-
  factor(
    path_merged_risk_lists$Analysis,
    levels = c("Differential expression",
               "Differential splicing")
  )

# add significance bar
path_merged_risk_lists %<>% mutate(is_sig = case_when(p_value <= 0.05 ~ "*",
                                                      p_value > 0.05 ~ ""))

# fdr correct each sig gene list by 2 - i.e the two pathways tested for each sig gene set.
path_merged_risk_lists %<>% group_by(Braak_stage, comparison, Analysis) %>% mutate(fdr=p.adjust(p_value, method="BH")) 

# sanity check for dplyr grouping code:
# path_merged_risk_lists %>% filter(Braak_stage=="Mid-stage", comparison=="Ant.\nCingulate", Analysis=="Differential splicing") %>% mutate(fdr=p.adjust(p_value, method="BH")) %>% select(fdr, everything())
# same outcome

# graph of combined mid-stage DEG and Splicing
p1 <- ggplot() +
  theme_bw(base_family = "Roboto") +
  facet_grid(Analysis ~ Braak_stage, drop = TRUE, scales = 'free') +
  # Plot DEGs with scale 1
  geom_tile(
    data = path_merged_risk_lists %>%
      filter(Braak_stage == "Mid-stage", Analysis == "Differential expression") %>%
      mutate(
        Braak_stage = droplevels(Braak_stage),
        Analysis = droplevels(Analysis)
      ),
    aes(
      x = comparison,
      y = target_list,
      fill = -log10(fdr)
    ),
    color = "black",
    width = 1,
    show.legend = TRUE
  ) +
  scale_fill_gradient(
    limits = c(0, 4),
    n.breaks = 3,
    low = "white",
    high = "#0c5394",
    name = "-log"[1][0] ~ "(FDR)"

  ) +
  
  new_scale_fill() +  # Add a new fill scale
  
  # Plot second subset with scale 2 - mid stage splicing
  geom_tile(
    data = path_merged_risk_lists %>%
      filter(Braak_stage == "Mid-stage", Analysis == "Differential splicing") %>%
      mutate(
        Braak_stage = droplevels(Braak_stage),
        Analysis = droplevels(Analysis)
      ),
    aes(
      x = comparison,
      y = target_list,
      fill = -log10(fdr)
    ),
    color = "black",
    width = 1,
    show.legend = TRUE
  ) +
  scale_fill_gradient(
    limits = c(0, 4),
    n.breaks = 3,
    low = "white",
    high = "#e06666",
    name = "-log"[1][0] ~ "(FDR)"
  ) +
  geom_text(
    data = path_merged_risk_lists %>%
      filter(Braak_stage == "Mid-stage") %>%
      mutate(
        Braak_stage = droplevels(Braak_stage),
        Analysis = droplevels(Analysis)
      ), size=5,

    aes(x = comparison, y = target_list, label = is_sig)
  ) +
  # Add labels and theme elements
  labs(x = "", y = "", fill = "hola") +
  theme_jb()+
  guides(size="none")+
  theme(
    legend.position = "right",
    strip.background = element_rect(fill = stage_colours[1]),
    axis.text.x = element_text(
      angle = 60,
      hjust = 1,
      vjust = 1,
      size = 10
    ),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(), # removed ticks and y row names
    text = element_text(size = 10),
    strip.text.y = element_blank(),
    strip.text = element_text(size = 12),
    panel.background = element_rect(fill = "white"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid = element_blank()
  )+
  theme(legend.key.size = unit(0.4, 'cm'))



# graph of combined late-stage DEG and Splicing
p2 <- ggplot() +
  theme_bw(base_family = "Roboto") +
  facet_grid(Analysis ~ Braak_stage, drop = TRUE, scales = 'free') +
  # Plot first DEG subset with scale 1
  geom_tile(
    data = path_merged_risk_lists %>%
      filter(
        Braak_stage == "Late-stage",
        Analysis == "Differential expression"
      ) %>%
      mutate(comparison = factor(
        comparison,
        levels = c(
          "Ant.\nCingulate",
          "Temporal",
          "Frontal",
          "Parietal"
        )
      )) %>%
      mutate(
        Braak_stage = droplevels(Braak_stage),
        Analysis = droplevels(Analysis)
      ),
    aes(
      x = comparison,
      y = target_list,
      fill = -log10(fdr)
    ),
    color = "black",
    width = 1,
    show.legend = TRUE
  ) +
  
  scale_fill_gradient(
    limits = c(0, 4),
    n.breaks = 3,
    low = "white",
    high = "#0c5394",
    name = "-log"[1][0] ~ "(FDR)"
    # ,
  ) +
  
  new_scale_fill() +  # Add a new fill scale
  
  # Plot second DS subset with scale 2
  geom_tile(
    data = path_merged_risk_lists %>%
      filter(Braak_stage == "Late-stage", Analysis == "Differential splicing") %>%
      mutate(comparison = factor(
        comparison,
        levels = c(
          "Ant.\nCingulate",
          "Temporal",
          "Frontal",
          "Parietal"
        )
      )) %>%
      mutate(
        Braak_stage = droplevels(Braak_stage),
        Analysis = droplevels(Analysis)
      ),
    aes(
      x = comparison,
      y = target_list,
      fill = -log10(fdr)
    ),
    color = "black",
    width = 1,
    show.legend = TRUE
  ) +
  scale_fill_gradient(
    limits = c(0, 4),
    n.breaks = 3,
    low = "white",
    high = "#e06666",
    name = "-log"[1][0] ~ "(FDR)"
    # ,
    # guide = guide_legend(order = 2)  # Set order of this legend
  ) +
  geom_text(
    data = path_merged_risk_lists %>%
      filter(Braak_stage == "Late-stage") %>%
      mutate(
        Braak_stage = droplevels(Braak_stage),
        Analysis = droplevels(Analysis)
      ), size=5,
    aes(x = comparison, y = target_list, label = is_sig)
  ) +
  # Add labels and theme elements
  labs(x = "", y = "", fill = "hola", size="") +
  scale_y_discrete(
    position = "right" # Move y-axis text to the right
  )+
  guides(size="none")+
  theme_jb()+
  theme(
    legend.position = "right",
    strip.background = element_rect(fill = stage_colours[2]),
    axis.text.x = element_text(
      angle = 60,
      hjust = 1.05,
      vjust = 1.05,
      size = 10
    ),
    axis.text.y = element_text(
      size = 8
    ),
    text = element_text(size = 10),
    strip.text.y = element_blank(),
    strip.text = element_text(size = 12),
    panel.background = element_rect(fill = "white"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid = element_blank()
    )+
  theme(legend.key.size = unit(0.4, 'cm'))



combined_plots<-p1+theme(legend.position = "none")+p2+plot_layout(ncol = 2, widths = c(9, 5))
# set width to number of x axis elements in each plot to make equal tile size

# save plot
ggsave(
  plot =  combined_plots,
  file.path(
    figure_out_path,
    "Figure_1c.png"
  ),
  width = 8000,
  height = 2500,
  dpi = 600,
  units = 'px'
)

