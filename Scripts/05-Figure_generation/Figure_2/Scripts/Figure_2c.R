
library(magrittr)
library(patchwork)
library(tidyverse)
library(foreach)
library(doParallel)
library(ggh4x)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggnewscale)

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
      dplyr::group_by(Description) %>% 
      dplyr::mutate(min.padj = min(p.adjust)) %>% 
      dplyr::ungroup() %>%
      dplyr::arrange(desc(n), p.adjust) %>%
      dplyr::mutate(Description = factor(Description, levels = unique(Description)))
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



# GO variables
count_threshold <- 2
go_reduce_threshold <- 0.7

# GO plotting variables
x_expansion = c(0.5, 0.5)
x='Description'
drop_facet = T
drop_x = F

#=====================================================================================================================================================================

# Loading splicing results and thresholds for functions

#=====================================================================================================================================================================

# make sure not to have small pathways where only one gene is found, at least 2 or more required
count_threshold <- 2

# similarity threshold for collapsing pathways
go_reduce_threshold <- 0.7

# Get df for convering between gene symbols and ENTREZ IDs
uniKeys <- AnnotationDbi::keys(org.Hs.eg.db, keytype = "ENSEMBL")
cols <- c("SYMBOL", "ENTREZID", "ENSEMBL")
cross_comp_list <- AnnotationDbi::select(org.Hs.eg.db, keys = uniKeys, columns = cols, keytype = "ENSEMBL")

# Remove duplicate symbols
non_ducplicates_idx <- which(duplicated(cross_comp_list$SYMBOL) == F)
cross_comp_list <- cross_comp_list[non_ducplicates_idx, ]

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

merge_df_genes %<>%
  mutate(Dataset = 
           factor(Dataset, 
                  levels = c("Mid-stage", "Late-stage")))

merge_df_genes %<>% mutate(stage_colour=case_when(
  Dataset == "Mid-stage" ~ "#91D1C2",
  Dataset == "Late-stage" ~ "#4B8D86"
)
)

# Combine with Entrez gene reference list for GO analysis and remove genes without an entrez id
merge_df_genes <- merge_df_genes %>% 
  dplyr::left_join(cross_comp_list, by = c("SYMBOL"), relationship = "many-to-many") %>%
  dplyr::select(ENTREZID, SYMBOL, Dataset) %>%
  dplyr::filter(!is.na(ENTREZID))

# change mid universe to valid entrez ids

mid_universe_df<-cross_comp_list[cross_comp_list$SYMBOL %in% mid_universe,] %>%
  dplyr::filter(!is.na(ENTREZID))

mid_universe<-mid_universe_df$ENTREZID

# change late universe to valid entrez ids

late_universe_df<-cross_comp_list[cross_comp_list$SYMBOL %in% late_universe,] %>%
  dplyr::filter(!is.na(ENTREZID))

late_universe<-late_universe_df$ENTREZID

## Run GO Enrichment ----

# set variable path so don't need to rerun everytime for new graphs - if need to rerun GO analysis: delete this files first
formula_res_merge_path <- file.path(data_path, "GO_results/ASAP_idiopathic_PD_tdpKO_match_formula_GO.rds")
formula_res_merge_red_path <- file.path(data_path, "GO_results/ASAP_idiopathic_PD_tdpKO_match_formula_GO_Reduced.rds")

# check if previous output still present
if(!file_test("-f", formula_res_merge_path, formula_res_merge_red_path)){

  # run GO enrichment on mid-stage overlap across all ontologies: BP, MF, CC
  formula_res_mid <- enrichGO(
    gene = merge_df_genes %>% filter(Dataset=='Mid-stage') %>% pull(ENTREZID),
    OrgDb = org.Hs.eg.db,
    ont = "ALL",
    universe = mid_universe,
    readable = T
  )
  
  
  # run GO enrichment on late-stage overlap genes cross all ontologies: BP, MF, CC
  
  formula_res_late <- enrichGO(
    gene = merge_df_genes %>% filter(Dataset=='Late-stage') %>% pull(ENTREZID),
    OrgDb = org.Hs.eg.db,
    ont = "ALL",
    universe = late_universe,
    readable = T
  )
  

  # filter out pathways with only one gene present
  formula_res_mid@result %<>% dplyr::filter(Count >= count_threshold)
  formula_res_late@result %<>% dplyr::filter(Count >= count_threshold)
  
  # Merge analyses
  formula_res_merge <- clusterProfiler::merge_result(list('Mid-stage' = formula_res_mid,
                                                          'Late-stage' = formula_res_late))
  
  # rename first column
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
  # if previous results present read those in
  formula_res_merge <- readRDS(formula_res_merge_path)
  formula_res_merge_red <- readRDS(formula_res_merge_red_path)
}


# run dataframe function
plot_df<-generatePlotDataframe(formula_res_merge_red) %>% arrange(qvalue)
# set order of ontologies
plot_df$go_type<-factor(plot_df$go_type, levels=c("MF", "CC", "BP"))
# set order of Datasets: mid then late stage
plot_df$Dataset<-factor(plot_df$Dataset, levels=c("Mid-stage", "Late-stage"))

#....................................................................................................

# Graphing ontologies separately to then combine with patchwork ----

#....................................................................................................


p_final_mf<-ggplot(plot_df %>% filter(go_type=="MF"),
                   aes(str_wrap(Description, width=15),
                       -log10(p.adjust), colour=Dataset, size=DOSE::parse_ratio(GeneRatio)*100))+
  geom_point()+
  facet_wrap(
    ~go_type, 
    drop = T, 
    scales = 'free_x', 
    labeller = labeller(go_type = c(
      "MF" = "Molecular function"
    ))
  )+
  scale_color_manual(values = stage_colours)+
  scale_size(range = c(4, 8))+
  # theme_bw()+
  theme_jb()+
  guides(color="none")+
  labs(x="", y="-log"[1][0]~"(FDR)", size="GeneRatio (%)", color="")+
  theme(
    strip.text = element_text(size=15),
    strip.background = element_rect(fill = scales::alpha("grey", 0.2)),
    legend.position = "top",
    text=element_text(size=14),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 12
    ),

    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.key.size = unit(0.5, "cm")
  )

p_final_cc<-ggplot(plot_df %>% filter(go_type=="CC"),
                   aes(fct_reorder(str_wrap(Description, width=15), n,.desc = TRUE),
                       -log10(p.adjust), colour=Dataset, size=DOSE::parse_ratio(GeneRatio)*100))+
  geom_point()+
  facet_wrap(
    ~go_type, 
    drop = T, 
    scales = 'free_x', 
    labeller = labeller(go_type = c(
      "CC" = "Cellular component"
    ))
  )+
  scale_color_manual(values = stage_colours)+
  scale_size(range = c(4, 8))+
  theme_jb()+
  guides(
    colour = guide_legend(override.aes = list(size = 6), order = 1),
    size = guide_legend(order = 2)
  )+ # increase point size in color legend
  labs(x="", y="-log"[1][0]~"(FDR)", size="GeneRatio (%)", color="")+
  theme(
    strip.text = element_text(size=15),
    strip.background = element_rect(fill = scales::alpha("grey", 0.2)),
    legend.position = "top",
    text=element_text(size=14),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 12
    ),
 
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.key.size = unit(0.7, "cm")
  )

p_final_bp<-ggplot(plot_df %>% filter(go_type=="BP"),
                   aes(fct_reorder(str_wrap(Description, width=22), n,.desc = TRUE),
                       # aes(Description, 
                       -log10(p.adjust), colour=Dataset, size=DOSE::parse_ratio(GeneRatio)*100))+
  geom_point()+
  facet_wrap(
    ~go_type, 
    drop = T, 
    scales = 'free_x', 
    labeller = labeller(go_type = c(
      "BP" = "Biological process"
    ))
  )+
  scale_color_manual(values = stage_colours)+
  scale_size(range = c(4, 8))+
  theme_jb()+
  guides(color="none")+
  labs(x="", y="-log"[1][0]~"(FDR)", size="GeneRatio (%)", color="")+
  theme(
    strip.text = element_text(size=15),
    strip.background = element_rect(fill = scales::alpha("grey", 0.2)),
    legend.position = "bottom",
    text=element_text(size=14),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 12
    ),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.key.size = unit(0.5, "cm")
  )


ggsave(
  plot =
    ((p_final_mf + p_final_cc + plot_layout(guides = 'collect')) +
       plot_annotation(theme = theme(
         legend.position = "top",
         legend.box.margin = margin(t = 0, r = 0, b = 0, l = 50, unit = "pt")
       ))) /
    plot_spacer() +
    p_final_bp +
    plot_layout(heights = c(5, -0.5, 5)),
  file.path(figure_out_path, "Figure_2c.png"),
  width = 10500,
  height = 9000,
  dpi = 600,
  units = 'px'
)




