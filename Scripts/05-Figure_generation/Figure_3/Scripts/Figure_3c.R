
library(magrittr)
library(patchwork)
library(tidyverse)
library(foreach)
library(doParallel)
library(ggh4x)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggnewscale)
library(ggtree)
library(DOSE)
library(enrichplot)

# Functions:

## 1. GO Reduction function
updated_go_reduce<-function (pathway_df, orgdb = "org.Hs.eg.db", threshold = 0.7, 
                             scores = NULL, measure = "Wang") 
{
  if (!measure %in% c("Resnik", "Lin", "Rel", "Jiang", "Wang")) {
    stop("Chosen measure is not one of the recognised measures, c(\"Resnik\", \"Lin\", \"Rel\", \"Jiang\", \"Wang\").")
  }
  if (measure == "Wang") {
    computeIC <- FALSE
  }
  else {
    computeIC <- TRUE
  }
  ont <- pathway_df %>% .[["go_type"]] %>% unique()
  if (any(!ont %in% c("BP", "CC", "MF"))) {
    stop("Column go_type does not contain the recognised sub-ontologies, c(\"BP\", \"CC\", \"MF\")")
  }
  go_similarity <- setNames(object = vector(mode = "list", 
                                            length = length(ont)), nm = ont)
  for (i in 1:length(ont)) {
    print(stringr::str_c("Reducing sub-ontology: ", ont[i]))
    hsGO <- GOSemSim::godata(OrgDb = orgdb, ont = ont[i], 
                             computeIC = computeIC)
    terms <- pathway_df %>% dplyr::filter(.data$go_type == 
                                            ont[i]) %>% .[["go_id"]] %>% unique()
    sim <- GOSemSim::mgoSim(GO1 = terms, GO2 = terms, semData = hsGO, 
                            measure = measure, combine = NULL)
    go_similarity[[i]] <- rrvgo::reduceSimMatrix(simMatrix = sim, 
                                                 threshold = threshold, orgdb = orgdb) %>% 
      tibble::as_tibble() %>% dplyr::rename(parent_id = .data$parent, 
                                            parent_term = .data$parentTerm, parent_sim_score = .data$termUniqueness)
  }
  go_sim_df <- go_similarity %>% qdapTools::list_df2df(col1 = "go_type")
  pathway_go_sim_df <- pathway_df %>% dplyr::inner_join(go_sim_df %>% 
                                                          dplyr::select(.data$go_type, go_id = .data$go, contains("parent")), 
                                                        by = c("go_type", "go_id")) %>% dplyr::arrange(.data$go_type, 
                                                                                                       .data$parent_id, -.data$parent_sim_score)
  return(pathway_go_sim_df)
}

## 2. Generate dataframe from GO results for plotting function:

generatePlotDataframe <- function(formula_res){
  if("parent_term" %in% colnames(formula_res)){
    description_names <- formula_res %>% 
      dplyr::distinct(Description, parent_term) %>% 
      dplyr::group_by(parent_term) %>% 
      dplyr::count() %>%
      dplyr::mutate(mod_description = stringr::str_to_title(paste0(parent_term, " (", n, ")"))) %>%
      dplyr::mutate(mod_description = paste0(strwrap(mod_description, width = 50, prefix = "\n", initial = ""),
                                             collapse = "")) %>%
      dplyr::select(parent_term, n, mod_description)
    
    plot_df <- formula_res %>% 
      dplyr::group_by(Dataset, parent_term) %>%
      dplyr::filter(p.adjust == min(p.adjust)) %>%
      dplyr::filter(parent_sim_score == max(parent_sim_score)) %>%
      dplyr::ungroup() %>%
      dplyr::left_join(description_names, by = "parent_term") %>%
      dplyr::mutate(Description = mod_description) %>%
      # tidyr::separate(GeneRatio, into = c("gr1", "gr2"), convert = T) %>%
      # dplyr::mutate(geneRatio = gr1/gr2) %>%
      dplyr::group_by(Description) %>% 
      dplyr::mutate(min.padj = min(p.adjust)) %>% 
      dplyr::ungroup() %>%
      # dplyr::arrange(desc(n), tissue, p.adjust) %>%
      dplyr::arrange(desc(n), p.adjust) %>%
      dplyr::mutate(Description = factor(Description, levels = unique(Description)))
    # ,tissue = factor(tissue, levels =c("SN", "CAU", "PUT", "T_CTX", "PARA", "C_CTX", "F_CTX", "P_CTX")))
    
  }else{
    plot_df <- formula_res %>% 
      dplyr::rowwise() %>% 
      dplyr::mutate(Description = stringr::str_to_title(Description),
                    Description = paste0(strwrap(Description, width = 50, prefix = "\n", initial = ""), 
                                         collapse = "")) %>% 
      dplyr::group_by(Description) %>% 
      dplyr::mutate(n = n(),
                    min.padj = min(p.adjust)) %>% 
      dplyr::ungroup() %>%
      # dplyr::arrange(desc(n), min.padj, tissue, p.adjust) %>%
      dplyr::arrange(desc(n), p.adjust) %>%
      dplyr::mutate(Description = factor(Description, levels = unique(Description)))
  }
}


## 3. Tile GO plot funciton

drawPathwaysGO_tile <- function(formula_res_df, x='Description'){
  plot<-formula_res_df %>% 
    
    mutate(blank="") %>% # Make a blank column title to use for x axis
    mutate(GeneRatio=DOSE::parse_ratio(GeneRatio)) %>%
    ggplot(aes(
      # y = str_wrap(Description, width=20), 
      y = fct_reorder(str_wrap(Description, width=40), n,.desc = TRUE),
      x=blank, alpha = -log10(p.adjust))) +
    geom_tile(color = "black", width = 1, show.legend = T, fill="#125C5D") +
    geom_hline(yintercept = seq(0.5, length(plot_df$Description)), color="black", linewidth = 0.2) +
    geom_vline(xintercept = seq(0.5, length(plot_df[, x, T])+1), color="black", linewidth = 0.2) +
    scale_y_discrete(limits=rev, expand = expansion(add = x_expansion)) +
    scale_x_discrete(drop = drop_x, expand = expansion(add = x_expansion)) +
    coord_cartesian()+
    labs( x = "GO Term", y = "", fill = "hola") +
    theme(legend.position = "bottom",
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          strip.text.x = element_text(size = 24),
          text=element_text(size=22),
          panel.background = element_rect(fill = "white"),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid = element_blank())
  
  return(plot)
}


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

Braak_stage_order<-c("Mid Stage", "Late Stage")

mid_brain_area_order<-c("Substantia nigra", "Caudate", "Putamen", "Ant.\nCingulate",
                        "Parahippocampal", 'Temporal',
                        "Frontal", "Parietal", "Combined")

late_brain_area_order<-c("Ant.\nCingulate", 'Temporal ',
                         "Frontal", "Parietal", "Combined")

extrafont::font_import(paths = data_path, 
                       pattern = 'Roboto', prompt = F)

# GO variables
count_threshold <- 2
go_reduce_threshold <- 0.7

# GO plotting variables
x_expansion = c(0.5, 0.5)
x='Description'
drop_facet = T
drop_x = F

# Get df for convering between gene symbols and ENTREZ IDs
uniKeys <- AnnotationDbi::keys(org.Hs.eg.db, keytype = "ENSEMBL")
cols <- c("SYMBOL", "ENTREZID", "ENSEMBL")
cross_comp_list <- AnnotationDbi::select(org.Hs.eg.db, keys = uniKeys, columns = cols, keytype = "ENSEMBL")

# Remove duplicate symbols
non_ducplicates_idx <- which(duplicated(cross_comp_list$SYMBOL) == F)
cross_comp_list <- cross_comp_list[non_ducplicates_idx, ]




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

sig_overlap <- sig_overlap %>% 
  dplyr::left_join(cross_comp_list, by = c("SYMBOL"), relationship = "many-to-many") %>%
  dplyr::select(ENTREZID, SYMBOL, Dataset) %>%
  dplyr::filter(!is.na(ENTREZID))


universe_df<-universe %>% left_join(.,cross_comp_list, by=c("genes" = "SYMBOL")) %>%
  dplyr::filter(!is.na(ENTREZID))

#============================================================================================================================================

# 2. Run GO Enrichment ----

#============================================================================================================================================

formula_res_merge_path <- file.path(data_path, "GO_results/ASAP_LRRK2_HS_PM_hPSC_tdpKO_match_formula_GO.rds")
formula_res_merge_red_path <- file.path(data_path, "GO_results/ASAP_LRRK2_HS_PM_hPSC_tdpKO_match_formula_GO_Reduced.rds")

# set variable path so don't need to rerun everytime for new graphs - if need to rerun GO analysis: delete this files first
if(!file_test("-f", formula_res_merge_path, formula_res_merge_red_path)){
# Running each dataset with own universe
  
# human post-mortem
  formula_res_pm <- enrichGO(
    gene = sig_overlap %>% filter(Dataset=='LRRK2 p.G2019S\nPost-mortem\nbrain') %>% pull(ENTREZID),
    OrgDb = org.Hs.eg.db,
    ont = "ALL",
    universe = universe_df %>% filter(Dataset=='LRRK2 p.G2019S\nPost-mortem\nbrain') %>% pull(ENTREZID),
    readable = T
  )
  
# hpscs
  
  formula_res_ips <- enrichGO(
    gene = sig_overlap %>% filter(Dataset=='LRRK2 p.G2019S\nmDNs') %>% pull(ENTREZID),
    OrgDb = org.Hs.eg.db,
    ont = "ALL",
    universe = universe_df %>% filter(Dataset=='LRRK2 p.G2019S\nmDNs')%>% pull(ENTREZID),
    readable = T
  )
  
# Merge ontologies
    
  formula_res_pm@result %<>% dplyr::filter(Count >= count_threshold)
    formula_res_ips@result %<>% dplyr::filter(Count >= count_threshold)

    formula_res_merge <- clusterProfiler::merge_result(list('LRRK2 p.G2019S\nPost-mortem\nbrain' = formula_res_pm,
                                                            'LRRK2 p.G2019S\nmDNs' = formula_res_ips))
    
    names(formula_res_merge@compareClusterResult)[1] <- "Dataset"
    # change names to match previous code
    names(formula_res_merge@compareClusterResult)[which(names(formula_res_merge@compareClusterResult)=="ONTOLOGY")] <- "go_type"
    names(formula_res_merge@compareClusterResult)[which(names(formula_res_merge@compareClusterResult)=="ID")] <- "go_id"

   
  
  # Reduce terms  
  formula_res_merge_red<-formula_res_merge@compareClusterResult %<>% updated_go_reduce(threshold = go_reduce_threshold)

  # Save results
  formula_res_merge %>% saveRDS(formula_res_merge_path)
  formula_res_merge_red %>% saveRDS(formula_res_merge_red_path)
  } else {
  formula_res_merge <- readRDS(formula_res_merge_path)
  formula_res_merge_red <- readRDS(formula_res_merge_red_path)
  }


plot_df<-generatePlotDataframe(formula_res_merge_red) %>% arrange(qvalue)
plot_df$go_type<-factor(plot_df$go_type, levels=c("MF", "CC", "BP"))


plot_df %<>% 
  mutate(Dataset=factor(Dataset,
                        levels=c("LRRK2 p.G2019S\nPost-mortem\nbrain",'LRRK2 p.G2019S\nmDNs'))) %>% arrange(desc(n))

top_10_terms<-plot_df %>% arrange(desc(n)) %>% distinct(Description) %>% slice_head(n=25) %>% pull(Description)

p_tile<-drawPathwaysGO_tile(plot_df %>% 
                              mutate(blank="") %>% 
                              dplyr::filter(Description %in% top_10_terms) %>% arrange(desc(n))
                            )+ 
  facet_wrap(
  ~Dataset,
  drop = T,
)+
  labs(alpha="-log"[1][0]~"(FDR)")+
  theme_bw(base_family = "Roboto")+
  theme(legend.position = "right",
                 axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                 axis.title.x = element_blank(),
                 axis.title.y = element_blank(),
                 text=element_text(size=8),
        strip.text = element_text(size = 12),
                 panel.background = element_rect(fill = "white"),
                 panel.grid.minor = element_blank(),
                 panel.grid.major.y = element_blank(),
                 panel.grid.major.x = element_blank(),
                 panel.grid = element_blank(),
        axis.text.y = element_text(
          size = 14
        ),
        legend.key.size = unit(0.9, 'cm'),
        legend.text = element_text(size=10),
        legend.title = element_text(size=11),
        legend.spacing.x =  unit(0.9, 'cm'),
        legend.spacing.y = unit(0.9, 'cm')
        )+
  theme(
    strip.placement = "inside",
    panel.spacing = unit(0.01, "lines") # Remove space between facets
  )+
  theme(
    strip.text = element_text(size = 15, margin = margin(t = 7, b = 7)),
    strip.background = element_rect(fill = "#f0f0f0", color = "black")
  )



p_tile



ggsave(
  plot =  p_tile,
  file.path(figure_out_path, "Figure_3c.png"),
  width = 3500,
  height = 4500,
  dpi = 375,
  units = 'px'
)


