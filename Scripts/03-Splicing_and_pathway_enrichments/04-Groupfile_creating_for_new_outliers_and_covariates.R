library(optparse)
library(tidyverse)
library(magrittr)

# need to change to new intron clustering file - 
count_file<-"~/Hardy_ASAP_bulk/Analysis/Hardy_Leafcutter_R_analysis/03-Checking_percent_cluster_reads_leafcutter_params/Clustering_introns/0.1_pcnt/leafcutter_perind_numers.counts.gz"
metadata_cols_path<-"~/Hardy_ASAP_bulk/Analysis/Hardy_Leafcutter_R_analysis/03-Checking_percent_cluster_reads_leafcutter_params//Metadata_and_selection/metadata_cols_selected.txt"
out_dir<-"~/Hardy_ASAP_bulk/Analysis/Hardy_Leafcutter_R_analysis/07-Rerun_of_final_PD_groups_collapsed/Creating_groupfiles/"

cf_sample_names <- data.table::fread(
  file.path(count_file) ) %>% dplyr::select(-V1) %>%
  colnames()

cf_sample_names

metadata_cols<-read.table(file = metadata_cols_path, header = T, sep = " ")

Sample_col_header<-names(metadata_cols)[1]
Group_col_header<-names(metadata_cols)[2]

md_samples<-metadata_cols[,1]

searched_md_samps<-sapply(md_samples, cf_sample_names, FUN = grep)

samp_names_key<-unique(names(unlist(searched_md_samps)))
order<-unique(unlist(searched_md_samps))

metadata<-unique(metadata_cols[which(metadata_cols[,1] %in% samp_names_key),])

meta_df <-
  readRDS(
    "~/Hardy_ASAP_bulk/Analysis/Hardy_Leafcutter_R_analysis/02-Run_with_updated_covariates_and_outliers/Metadata_and_selection/metadata_without_removal_of_correlated_tersm.rds"
  )

# Regina tibble way - check ordering doesn't matter but does appear to as Regina reordered hers
# https://rhreynolds.github.io/LBD-seq-bulk-analyses/quantification_splicing/Leafcutter.html#
# master <- 
#   tibble(lc_sample_name = cf_sample_names,
#          bxp_id_full = cf_sample_names %>% 
#            str_replace("_S.*", "")) %>% 
#   inner_join(metadata_cols, by = c('bxp_id_full')) %>% 
#   dplyr::select(-bxp_id_full) %>% 
#   arrange(Group)

# metadata$order<-order
metadata$lc_names<-cf_sample_names[order]

dir.create(str_c(out_dir, "ACG"))
dir.create(str_c(out_dir, "IPL"))
dir.create(str_c(out_dir, "MFG"))
dir.create(str_c(out_dir, "MTG"))

BAs<-c("ACG", "IPL", "MFG","MTG")
# BAs<-c("ACG", "IPL")

metadata$Group<-gsub(x = metadata$Group, pattern = " ", replacement = "_")
metadata$Group<-gsub(x = metadata$Group, pattern = ",", replacement = "")

metadata %<>% mutate(Group2 = case_when(Group == "short_no_dementia" ~ 'PD',
                                        Group == "long_no_dementia" ~ 'PD',
                                        Group == "short_with_dementia" ~ 'PD',
                                        Group == "long_with_dementia" ~ 'PD',
                                        TRUE ~ 'Control'))
metadata$Group<-metadata$Group2
metadata %<>% select(-Group2)

# Removing outliers - this was done prior to clustering so not needed here anymore
# sample_outliers16 <-
#   readRDS(
#     "~/Hardy_ASAP_bulk/Analysis/Hardy_Leafcutter_R_analysis/02-Run_with_updated_covariates_and_outliers/Creating_groupfiles/16_outliers/sample_outliers16.rds"
#   )
# 
# sample_outliers16 %<>% filter(outlier==T)
# 
# metadata<-metadata[!metadata$bxp_id_full %in% sample_outliers16$bxp_id_full,]


for (i in BAs) {
  x<-meta_df %>% filter(brain_region==i)
  metadata2<-metadata[metadata$bxp_id_full %in% x$bxp_id_full,]
  metadata_refined<-metadata2 %>% relocate(lc_names, Group_col_header)
  metadata_refined<-select(metadata_refined, -c(Sample_col_header))
  
  RNAseqProcessing::create_group_files_multi_pairwisecomp(df = metadata_refined, 
                                                          group_column_name = Group_col_header, 
                                                          output_path = str_c(out_dir, i))
}

