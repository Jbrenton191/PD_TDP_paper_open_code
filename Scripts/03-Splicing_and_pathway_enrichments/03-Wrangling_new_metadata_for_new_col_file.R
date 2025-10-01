library(tidyverse)
library(optparse)

# Adapted from select_metadata_cols_man_run_2.6.23.R in nextflow folder 
# moved metadata file from Ryten cluster that contained all metadata for all samples 
# used this version as Leafcutter_ds.R scales continuous cofounders - so don't want changed beforehand as did for DE analysis

meta_df <-
  readRDS(
    "~/Hardy_ASAP_bulk/Analysis/Hardy_Leafcutter_R_analysis/08-No_blacklist_removal_final_PD_groups_collapsed/Metadata/metadata_without_removal_of_correlated_terms.rds"
  )


# New covariates found through Guillermo analysis - taken from Onedrive - new covariates folder: 04-Hardy_DGE_report.html
# And added bxp_id_full and group
name_df<-c("bxp_id_full",
           "Group",
           'GC_NC_40_59',
           'PCT_UTR_BASES',
           'gender',
           'Total_Sequences',
           'INTRONIC_BASES',
           'tss_up_1kb_tag_pct')


# "GC_NC_40_59", "gender",
# "tss_up_1kb_tag_pct", "PCT_UTR_BASES",
# "Total_Sequences", "INTRONIC_BASES"

colnames(meta_df)<-gsub(pattern = "%", replacement = "percent_", colnames(meta_df))

colnames(meta_df)<-gsub(pattern = " ", replacement = "_", colnames(meta_df))


name_vec<-as.vector(unlist(name_df))

# selecting the metadata files
key<-c()

for(i in 1:length(name_vec)){
  key[i]<-which(names(meta_df) %in% name_vec[i])
}
selected_metadata<-meta_df[,key]

write.table(file = "~/Hardy_ASAP_bulk/Analysis/Hardy_Leafcutter_R_analysis/03-Checking_percent_cluster_reads_leafcutter_params/Metadata_and_selection/metadata_cols_selected.txt",
            x = selected_metadata,
            row.names = F)




