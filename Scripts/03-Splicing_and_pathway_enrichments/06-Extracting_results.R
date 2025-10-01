library(tidyverse)
library(magrittr)



# specify folder with runs within - can be levels above and will search through them - check doesn't break gsub code!
# only final settings params
leafcutter_dir<-"/home/drihome/MRJonathanBrenton/Hardy_ASAP_bulk/Analysis/Hardy_Leafcutter_R_analysis/07-Rerun_of_final_PD_groups_collapsed/Leafcutter_run"
clusters <- list.files(file.path(leafcutter_dir), pattern = "cluster_", recursive = T, full.names = T)
intron <- list.files(file.path(leafcutter_dir), pattern = "effect_", recursive = T, full.names = T)

extract_leafcutter_results<-function(clusters, intron){
  
  
  for(i in 1:length(clusters)){
    
    comparison <- clusters[i] %>% 
      str_replace("/.*/", "") %>% 
      str_replace("_cluster_significance.txt", "")
    
    comparison<-paste(gsub(x = clusters[i], pattern = "^.*/(.*)/.*$", replacement = '\\1'), comparison, sep = ".")
    
    lc <- read_delim(clusters[i], delim = "\t") %>% 
      
      dplyr::mutate(comparison = comparison) %>% 
      dplyr::select(comparison, everything())
    
    if(i == 1){
      
      lc_all <- lc
      
    } else{
      
      lc_all <- lc_all %>% 
        bind_rows(lc)
      
    }
    
  }
  
  
  for(i in 1:length(intron)){
    
    comparison <- intron[i] %>% 
      str_replace("/.*/", "") %>% 
      str_replace("_effect_sizes.txt", "")
    
    comparison<-paste(gsub(x = intron[i], pattern = "^.*/(.*)/.*$", replacement = '\\1'), comparison, sep = ".")
    
    int <- read_delim(intron[i], delim = "\t") %>% 
      dplyr::mutate(comparison = comparison) %>% 
      dplyr::select(comparison, everything())
    
    if(i == 1){
      
      int_all <- int
      
    } else{
      
      int_all <- int_all %>% 
        bind_rows(int)
      
    }
    
  }
  
  print(unique(lc_all$comparison))
  
  cluster_list<- lc_all %>% 
    # filter for successful tests and remove clusters overlapping multiple genes
    dplyr::filter(status == "Success", p.adjust < 0.05, !str_detect(cluster, ","))
  
  
  gene_list<- lc_all %>% 
    # filter for successful tests and remove clusters overlapping multiple genes
    dplyr::filter(status == "Success", p.adjust < 0.05, !str_detect(genes, ",")) %>% group_by(comparison) %>% 
    distinct(genes, .keep_all = T)
  
  universe<-lc_all %>% 
    # filter for successful tests and remove clusters overlapping multiple genes
    dplyr::filter(!str_detect(genes, ",")) %>% 
    select(genes) %>% distinct(genes)
  
  # Filter for significant clusters
  # Add direction of effect
  lc_filtered <- lc_all %>% 
    # filter for successful tests and remove clusters overlapping multiple genes
    dplyr::filter(status == "Success", p.adjust < 0.05, !str_detect(genes, ",")) %>%
    tidyr::separate(cluster, into = c("chr", "cluster"), sep = ":") %>% 
    dplyr::inner_join(int_all %>% 
                        tidyr::separate(intron, into = c("chr", "start", "end", "cluster"), sep = ":")) %>% 
    dplyr::mutate(direction_of_effect = case_when(deltapsi > 0 ~ "upregulated",
                                                  deltapsi < 0 ~ "downregulated",
                                                  deltapsi == 0 ~ "no_change"))
  
  int_all$intron_cluster<-gsub(x = int_all$intron, pattern = ".*:.*:.*:(clu_[0-9]*)_.*", replacement = "\\1")
  
  x<-int_all %>% group_by(comparison) %>% count(intron_cluster)
  
  med_mean<-x %>% group_by(comparison) %>% mutate(Median=median(n), Mean=mean(n), Minimum=min(n), Maximum=max(n)) %>% 
    distinct(comparison, .keep_all = T) %>% select(-intron_cluster, -n)
  
  
  all_res<-setNames(list(lc_all, int_all, lc_filtered),
                    c("cluster_significance", "intron_usage", "significant_clusters_0.05_filter"))
  
  
  # rm(int_all)
  # rm(lc_all)
  
  # SUMMARISE ----
  
  summary <- all_res$cluster_significance %>% 
    dplyr::group_by(comparison) %>% 
    dplyr::summarise(n = n()) %>%
    tidyr::spread(key = comparison, value = n) %>% 
    dplyr::mutate(status = "Total clusters across all samples (read >= 30)") %>% 
    dplyr::select(status, everything()) %>% 
    dplyr::bind_rows(all_res$cluster_significance %>% 
                       dplyr::group_by(comparison, status) %>% 
                       dplyr::summarise(n = n()) %>% 
                       tidyr::spread(key = comparison, value = n)) %>% 
    dplyr::bind_rows(all_res$significant_clusters_0.05_filter %>% 
                       dplyr::distinct(comparison, cluster, genes) %>% 
                       dplyr::group_by(comparison) %>% 
                       summarise(n = n()) %>% 
                       tidyr::spread(key = comparison, value = n) %>% 
                       dplyr::mutate(status = "Differentially spliced clusters, p.adjust < 0.05") %>% 
                       dplyr::select(status, everything())) %>% 
    dplyr::bind_rows(all_res$cluster_significance %>% 
                       # filter for successful tests and remove clusters that overlap multiple genes
                       dplyr::filter(status == "Success", p.adjust < 0.05, !str_detect(genes, ",")) %>%
                       tidyr::separate(cluster, into = c("chr", "cluster_id"), sep = ":") %>% 
                       dplyr::inner_join(all_res$intron_usage %>% 
                                           tidyr::separate(intron, into = c("chr", "start", "end", "cluster_id"), sep = ":")) %>% 
                       dplyr::filter(abs(deltapsi) >= 0.1) %>% 
                       dplyr::distinct(comparison, cluster_id, genes) %>% 
                       dplyr::group_by(comparison) %>% 
                       summarise(n = n()) %>% 
                       tidyr::spread(key = comparison, value = n) %>% 
                       dplyr::mutate(status = "Differentially spliced clusters, p.adjust < 0.05, |dPSI| >= 0.1") %>% 
                       dplyr::select(status, everything())) %>% 
    dplyr::mutate(status = str_replace(status, "Success", "Successfully tested"))
  
  # Add propotions
  summary <- summary %>% 
    dplyr::bind_rows((((summary %>% 
                          dplyr::filter(status == "Successfully tested") %>% 
                          dplyr::select(-status))/(summary %>% 
                                                     dplyr::filter(status == "Total clusters across all samples (read >= 30)") %>% 
                                                     dplyr::select(-status))) *100) %>% 
                       dplyr::mutate(status = "Successfully tested/Total clusters (%)")) %>% 
    dplyr::bind_rows((((summary %>% 
                          dplyr::filter(status == "Differentially spliced clusters, p.adjust < 0.05") %>% 
                          dplyr::select(-status))/(summary %>% 
                                                     dplyr::filter(status == "Successfully tested") %>% 
                                                     dplyr::select(-status))) *100) %>% 
                       dplyr::mutate(status = "Differentially spliced clusters/Successfully tested (%)"))
  
  
  
  return(setNames(list(summary,med_mean, all_res, cluster_list, gene_list, universe),
                  c("Summary", "intron_stats", "all_results", "list_of_significant_clusters", "list_of_significant_genes", "universe_list_of_genes_tested")))
}

results_list<-extract_leafcutter_results(clusters, intron)

results_list$Summary

saveRDS(results_list, 
        file = "/home/drihome/MRJonathanBrenton/Hardy_ASAP_bulk/Analysis/Hardy_Leafcutter_R_analysis/07-Rerun_of_final_PD_groups_collapsed/Analysis_and_Results/full_results_list_for_default_params.RDS")
