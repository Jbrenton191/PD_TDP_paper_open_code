# Gene level differential expression testing of proteomic data
# used for comparisons of gene regulation at transcript and proteomic level

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1. Setup ----
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Load libraries
suppressPackageStartupMessages({
  library(magrittr)
  library(tidyverse)
  library(car)
  library(lme4)
  
})


conflicted::conflict_prefer_all(winner = "dplyr", quiet = T)


# set paths
main_path<-here::here() # main project level i.e TDP43_paper_open_code level

data_path<-file.path(main_path, "Data/LRRK2_HS_PSC/mDN1")

results_path<-file.path(main_path, "Data/LRRK2_HS_PSC/mDN1")


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2. Finding Gene targets in mass spec data ----
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if (!file.exists(file.path(data_path, "lib_normalised_gene_level_summed_counts.csv"))) {
  
# read in isoform level counts for each gene and summate across genes to get per gene depth
gene_ms_counts<-read_tsv(file.path(data_path, "DRio_LRRK2-iPSC_TP-Raw.txt"))

Gene_level_sum_counts<-gene_ms_counts %>%
  separate_rows(Genes, sep = ";") %>%
  select(starts_with("Ewt"), starts_with("LRRK2-G"), Genes, Protein.Names, Protein.Ids) %>%
  pivot_longer(cols = c(starts_with("Ewt"), starts_with("LRRK2-G")), 
               names_to = 'Sample', values_to = 'Counts') %>%
  group_by(Sample) %>% 
  mutate(Library_size=sum(Counts)) %>%
  mutate(med_lib_size=median(Library_size)) %>% 
  mutate(sample_scaling_factor=Library_size/med_lib_size) %>%
  mutate(normalised_counts=Counts*sample_scaling_factor) %>%
  group_by(Sample, Genes) %>% 
  mutate(sum_norm_counts=sum(normalised_counts)) %>%
  ungroup() %>%
  mutate(Line=gsub(x = Sample, pattern="\\.0[0-9]$", "")) %>% 
  mutate(Group=case_when(
  grepl("Ewt.*", Sample)~"Control",
  grepl("LRRK2.*", Sample)~"LRRK2",
           )
           ) %>% 
  group_by(Sample, Genes) %>% 
  distinct(Sample, Genes, .keep_all = T) %>%
  select(Genes, Sample, Line, Group, sum_norm_counts) %>% 
  ungroup() %>%
  drop_na(Genes)


results_df<-c()



for (i in unique(Gene_level_sum_counts$Genes)) {
  
  filt_Gene_level_sum_counts<-Gene_level_sum_counts %>% 
    filter(Genes == i) 
  
 
  t.test_res<-t.test(sum_norm_counts ~ 
                       Group, filt_Gene_level_sum_counts)
  
  # Get fold change of counts
  log2FC<-filt_Gene_level_sum_counts %>% group_by(Group) %>% summarise(mean=mean(sum_norm_counts)) %>% ungroup() %>% 
    mutate(Log2FC=log2(.$mean[2]/.$mean[1])) %>% distinct(Log2FC) %>% pull()
  
  results<-tibble(Genes=i, 
                  Log2FC=log2FC,
                  t.test_t_val=t.test_res$statistic,  
                  t.test_p_val=t.test_res$p.value, 
                  t.test_control_mean=t.test_res$estimate[[1]],
                  t.test_lrrk2_mean=t.test_res$estimate[[2]]
  )
  
  results_df<-rbind(results_df, results)  
  
}

results_df_all<-results_df

write_csv(results_df_all, file.path(data_path, "lib_normalised_gene_level_summed_counts.csv"))

} else {

results_df_all<-read_csv(file = file.path(data_path, "lib_normalised_gene_level_summed_counts.csv"))


}

