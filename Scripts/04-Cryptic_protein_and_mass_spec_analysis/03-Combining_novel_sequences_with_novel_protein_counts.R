# Differential testing of protein sequences

# Novel 

library(car)
library(lme4)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1. Setup ----
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(purrr)
  library(readxl)
  library(Biostrings)
  library(magrittr)
})


conflicted::conflict_prefer_all(winner = "dplyr", quiet = T)


# set paths
main_path<-here::here() # main project level i.e TDP43_paper_open_code level

data_path<-file.path(main_path, "Data/LRRK2_HS_PSC/mDN1")

results_path<-file.path(main_path, "Data/LRRK2_HS_PSC/mDN1")


# load peptide sequences
novel_sequences_donors <-
  readxl::read_xlsx(
    path = file.path(
      data_path,
      "DRio_NovelDonor_DIA-NN_Updated.xlsx"
    ),
    sheet = 2
  ) %>% arrange(Protein.Group)

novel_sequences_acceptor <-
  readxl::read_xlsx(
    path = file.path(
      data_path,
      "DRio_NovelAcceptor_DIA-NN_Updated.xlsx"
    ),
    sheet = 2
  ) %>% arrange(Protein.Group)


# read in isoform level counts for each gene and summate across genes to get per gene depth
gene_ms_counts<-read_tsv(file.path(data_path, "DRio_LRRK2-iPSC_TP-Raw.txt"))

Gene_level_sum_counts<-gene_ms_counts %>%
  separate_rows(Genes, sep = ";") %>%
  select(starts_with("Ewt"), starts_with("LRRK2-G"), Genes, Protein.Names, Protein.Ids) %>%
  pivot_longer(cols = c(starts_with("Ewt"), starts_with("LRRK2-G")), 
               names_to = 'Sample', values_to = 'Counts') %>%
  group_by(Sample) %>% 
  summarise(Library_size=sum(Counts)) %>%
  mutate(med_lib_size=median(Library_size)) %>% 
  mutate(sample_scaling_factor=Library_size/med_lib_size)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2. Change results for model testing structure ----
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


novel_sequences_donor_long<-novel_sequences_donors %>% select(`Gene names`, Protein.Group, Precursor.Charge, Stripped.Sequence, starts_with("Ewt"), starts_with("LRRK2")) %>% 
  pivot_longer(cols = -c(`Gene names`, Protein.Group,Precursor.Charge, Stripped.Sequence,), names_to = 'Sample', values_to = 'Counts') %>% mutate(Group=case_when(
    grepl("Ewt.*", Sample)~"Control",
    grepl("LRRK2.*", Sample)~"LRRK2",
  )
  ) %>% mutate(Line=gsub(x = Sample, pattern="\\.0[0-9]$", "")
  ) %>% rename(gene=`Gene names`)
  


novel_sequences_acceptor_long<-novel_sequences_acceptor %>% select(`Gene names`, Protein.Group, Precursor.Charge, Stripped.Sequence, starts_with("Ewt"), starts_with("LRRK2")) %>% 
  pivot_longer(cols = -c(`Gene names`, Protein.Group,Precursor.Charge, Stripped.Sequence,), names_to = 'Sample', values_to = 'Counts') %>% mutate(Group=case_when(
    grepl("Ewt.*", Sample)~"Control",
    grepl("LRRK2.*", Sample)~"LRRK2",
    )
    ) %>% mutate(Line=gsub(x = Sample, pattern="\\.0[0-9]$", "")
    ) %>% rename(gene=`Gene names`)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2.5 Filter novel sequences for 7 amino acid matches ----
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# get matching results
matching_results <- read_csv(file.path(data_path,"Upreg_tdp_tsl_nmd_included_ALL_MS_seq_matching_to_novel_junction_aa_differences_to_annotated_7.5.25.csv"))
# filter this for sequences that have more than 6 amino acid matches between peptide ion sequence and junction derived sequence that does not matching annotated transcript
filt_num_matches<-matching_results %>% filter(number_matches>6)

# get novel peptide half of the detection rate - lowest value you can account for  - 
# lowest normalised peptide value /2 - add to every value
# add a scaling factor - to avoid infinite log2 fold changes - no infinite values, add constant  to all sequences

detection_rate_min<-round(min(c(novel_sequences_acceptor_long %>% filter(Counts>0) %>% pull(Counts),
                                
                                novel_sequences_donor_long %>% filter(Counts>0) %>% pull(Counts)))/2)

# filter novel sequence counts for these ions and protein only - add in 0s for not found ions / NAs
novel_sequences_donor_long<-novel_sequences_donor_long %>% 
  inner_join(., filt_num_matches %>% filter(transcript_status=="Novel_donor"), 
             by=c("gene", 'Stripped.Sequence'='query_peptide_ion', 
                  'Precursor.Charge'='precursor_charge')) %>%
  mutate(Counts = replace(Counts, is.na(Counts), 0)) %>% 
  mutate(Counts = Counts+detection_rate_min)

# merge with gene level counts - exclude ions that don't have gene level counts
novel_sequences_donor_long %<>%
  inner_join(.,
             Gene_level_sum_counts %>%
               select(Sample, sample_scaling_factor), by=c("Sample")) %>%

  mutate(normalised_counts=Counts*sample_scaling_factor)


# filter novel sequence counts for these ions and protein only - add in 0s for not found ions / NAs
novel_sequences_acceptor_long<-novel_sequences_acceptor_long %>% 
  inner_join(., filt_num_matches %>% filter(transcript_status=="Novel_acceptor"), 
  by=c("gene", 'Stripped.Sequence'='query_peptide_ion', 'Precursor.Charge'='precursor_charge')) %>%
  mutate(Counts = replace(Counts, is.na(Counts), 0)) %>% 
  mutate(Counts = Counts+detection_rate_min)

# merge with gene level counts - exclude ions that don't have gene level counts
novel_sequences_acceptor_long %<>%
  inner_join(.,
             Gene_level_sum_counts %>%
               select(Sample, sample_scaling_factor), by=c("Sample")) %>%

  mutate(normalised_counts=Counts*sample_scaling_factor)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 3. Run loop to run glm with random effect of line and t tests ----
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



# 3.1 donors ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


results_df<-c()

novel_sequences_donor_long<-novel_sequences_donor_long %>% 
  mutate(search_col=str_c(gene, transcript, modified_sequence, Precursor.Charge, sep = "-"))

for (i in unique(novel_sequences_donor_long$search_col)) {
  
  filt_novel_sequences_donor_long<-novel_sequences_donor_long %>% 
    filter(search_col == i) 
  
  filt_novel_sequences_donor_long$Group <- relevel(factor(filt_novel_sequences_donor_long$Group), ref = "Control")
  
  # run t test
  
  t.test_res<-t.test(normalised_counts ~ 
                       Group, filt_novel_sequences_donor_long
                     )
  
  # Get fold change of counts
  log2FC<-filt_novel_sequences_donor_long %>% group_by(Group) %>% summarise(mean=mean(normalised_counts)) %>% ungroup() %>% 
    mutate(Log2FC=log2(.$mean[2]/.$mean[1])) %>% distinct(Log2FC) %>% pull()
  
  # Handy calculation of Log2FC error:
  # se_fold_change <- abs(mean_A / mean_B) * sqrt((se_A / mean_A)^2 + (se_B / mean_B)^2)
  # print(paste("Standard Error of Fold Change:", round(se_fold_change, 2)))
  
  results<-tibble(Protein=i,
                  transcript=unique(filt_novel_sequences_donor_long$transcript),
                  gene=unique(filt_novel_sequences_donor_long$gene),
                  Precursor.Charge=unique(filt_novel_sequences_donor_long$Precursor.Charge),
                  Stripped.Sequence=unique(filt_novel_sequences_donor_long$Stripped.Sequence),
                  modified_sequence=unique(filt_novel_sequences_donor_long$modified_sequence),
                  AA_seq=unique(filt_novel_sequences_donor_long$AA_seq),
                  Log2FC=log2FC,
                  t.test_t_val=t.test_res$statistic,  
                  t.test_p_val=t.test_res$p.value, 
                  t.test_control_mean=t.test_res$estimate[[1]],
                  t.test_lrrk2_mean=t.test_res$estimate[[2]]
  )
  
  results_df<-rbind(results_df, results)  
  
}

results_df_donor<-results_df


## 3.2 acceptors +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

results_df<-c()

novel_sequences_acceptor_long<-novel_sequences_acceptor_long %>% mutate(search_col=str_c(gene, transcript, modified_sequence, Precursor.Charge, sep = "-"))

for (i in unique(novel_sequences_acceptor_long$search_col)) {
  
  filt_novel_sequences_acceptor_long<-novel_sequences_acceptor_long %>% 
    filter(search_col == i) 
  
  filt_novel_sequences_acceptor_long$Group <- relevel(factor(filt_novel_sequences_acceptor_long$Group), ref = "Control")
  
  
  t.test_res<-t.test(normalised_counts ~ 
                       Group, filt_novel_sequences_acceptor_long)
  
  # Get fold change of counts
  log2FC<-filt_novel_sequences_acceptor_long %>% group_by(Group) %>% summarise(mean=mean(normalised_counts)) %>% ungroup() %>% 
    mutate(Log2FC=log2(.$mean[2]/.$mean[1])) %>% distinct(Log2FC) %>% pull()
  
  results<-tibble(Protein=i, 
                  transcript=unique(filt_novel_sequences_acceptor_long$transcript),
                  gene=unique(filt_novel_sequences_acceptor_long$gene),
                  Precursor.Charge=unique(filt_novel_sequences_acceptor_long$Precursor.Charge),
                  Stripped.Sequence=unique(filt_novel_sequences_acceptor_long$Stripped.Sequence),
                  modified_sequence=unique(filt_novel_sequences_acceptor_long$modified_sequence),
                  AA_seq=unique(filt_novel_sequences_acceptor_long$AA_seq),
                  Log2FC=log2FC,
                  t.test_t_val=t.test_res$statistic,  
                  t.test_p_val=t.test_res$p.value, 
                  t.test_control_mean=t.test_res$estimate[[1]],
                  t.test_lrrk2_mean=t.test_res$estimate[[2]]
  )
  
  results_df<-rbind(results_df, results)  
  
}

results_df_acceptor<-results_df


## 3.3 combine the two and arrange ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


combined_novel_ms_junc_stats <-
  rbind(
    results_df_donor %>% mutate(junction_type = "novel_donor"),
    results_df_acceptor %>% mutate(junction_type = "novel_acceptor")
  ) %>% arrange(t.test_p_val) %>% mutate(
    t.test_fdr = p.adjust(t.test_p_val, method = "fdr")
  ) %>% left_join(
    .,
    filt_num_matches %>% distinct(),
    by = c(
      "gene",
      'Stripped.Sequence' = 'query_peptide_ion',
      'Precursor.Charge' = 'precursor_charge',
      'modified_sequence',
      "transcript",
      "AA_seq"
    )
  )


write_csv(combined_novel_ms_junc_stats, file = file.path(results_path, "novel_peps_detection_rate_min_added_LIB_size_norm_Upreg_tdp_tsl_nmd_included_NOVEL_MS_seq_matching_to_novel_junction_glm_and_t.tests_23.6.25.csv"))


